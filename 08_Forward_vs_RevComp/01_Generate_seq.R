# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(gtools)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(purrr)))
suppressPackageStartupMessages(suppressWarnings(library(furrr)))
suppressPackageStartupMessages(suppressWarnings(library(progressr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))
suppressPackageStartupMessages(suppressWarnings(library(arrow)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))

# source functions
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

pbapply::pboptions(char = "=", type = "txt")
options(future.seed = TRUE)
source("../lib/Breakpoints.R")
source("../lib/Features.R")
source("../lib/Preprocess.R")
Rcpp::sourceCpp("../lib/edlibFunction.cpp")
Rcpp::sourceCpp("../lib/eigenMatrixMultiply.cpp")

# randomly select sequences from human genome and get its reverse complement
seq_len <- 10000
genome <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0

# get the start and end positions of each chromosome
refseq.table <- as.data.frame(genome@seqinfo)
refseq.table <- refseq.table[grepl(
    pattern = "^([1-9]|1[0-9]|2[0-2])$", 
    x = rownames(refseq.table)
),]
refseq.table <- data.table(
    seqnames = rownames(refseq.table),
    len = refseq.table$seqlengths
)

# sample sequences
set.seed(1234)
N_sample <- 100
sample_dt <- refseq.table[sample(x = .N, size = N_sample, replace = TRUE)]
sample_dt[, `:=`(
    start = sapply(len, function(l) sample(x = 1:(l-1), size = 1)),
    width = seq_len,
    len = NULL
)]
setorder(sample_dt, seqnames, start)
sample_dt[, seqnames := paste0("chr", seqnames)]
seq_names <- paste0("Forward_", sample_dt$seqnames, "_", sample_dt$start)

# extract sequences
sample_dt <- plyranges::as_granges(sample_dt)

if(!any(grepl(pattern = "^chr", x = seqnames(genome)))){
    chr_names <- paste0("chr", seqnames(genome))
    seqnames(genome@seqinfo) <- seqnames(genome) <- chr_names
}
sequence <- Biostrings::getSeq(genome, sample_dt)

# sequence <- Biostrings::DNAStringSet(paste(c(rep("A", 1000), rep("C", 1000)), collapse = ""))
names(sequence) <- seq_names
sequence <- unique(sequence)

# get its reverse complement
rc_sequence <- Biostrings::reverseComplement(sequence)
names(rc_sequence) <- gsub("Forward", "RevComp", names(sequence))
sequence <- c(sequence, rc_sequence)

# save results
dir.create(path = "./data/", showWarnings = FALSE)
Biostrings::writeXStringSet(
    sequence,
    filepath = "./data/random_sequence.fasta.gz",
    format = "fasta",
    compress = TRUE
)

# extract features
seed <- 1234
k <- 8
break_score <- "zscore"
break_type <- "Biological"
ranges <- NULL
kmer_window <- 5
cores <- 1

future::plan(future::multicore)
data.table::setDTthreads(threads = cores)
set.seed(seed)
rng <- sample(100000, size = 20, replace = FALSE)
seeds <- rng[1]

for(x in 1:length(sequence)){
    cat(paste0("Processing sequence ", x, "/", length(sequence), ". \n"))
    
    preprocess <- Preprocess$new(
        fasta_sequence = sequence[[x]],
        label = names(sequence[x]),
        k = k,
        break_type = break_type,
        ranges = ranges
    )
    preprocess$preprocess_sequences()

    features <- Features$new(
        fasta_sequence = preprocess$fasta_sequence,
        label = preprocess$label,
        df_bp = preprocess$df_bp,
        k = k, 
        seed = seeds, 
        break_score = break_score,
        regression = FALSE,
        get_controls = FALSE,
        any_predictors = FALSE,
        only_breaks = TRUE
    )

    features$get_features(
        break_type = break_type, 
        ranges = NULL,
        FEAT_LONGRANGE_PREPARAMS = TRUE, 
        FEAT_G4_REGEX = TRUE, g4_type = "GPQS", 
        FEAT_GC_COUNT = TRUE,
        FEAT_KMER_COUNTS = TRUE, kmer_window = kmer_window,
        FEAT_VIENNA_RNA = FALSE, sliding_window = NULL, 
        nuc_type = "DNA", 
        RNAfold.CALL = "/home/imm/hert6114/anaconda3/bin/RNAfold",
        maxloopsize = 12,
        FEAT_DNA_SHAPE = FALSE
    )

    if(features$out == -1){
        cat("\n")
        cat("All short-range k-mers are Ns. Starting next chunk.")
        cat("\n")
        next
    }

    file_label <- paste0("../data/Forward_vs_RevComp/", preprocess$label)
    dir.create(
        path = file_label,
        showWarnings = FALSE,
        recursive = TRUE
    )

    write_parquet(
        features$feature_matrix, 
        paste0(file_label, "/feature_matrix.parquet")
    )
    write_parquet(
        features$true_bp_table,
         paste0(file_label, "/positions.parquet")
    )
}