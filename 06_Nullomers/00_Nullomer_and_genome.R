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

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
repeats <- as.numeric(args[2])
setwd(my.path)

pbapply::pboptions(char = "=", type = "txt")
options(future.seed = TRUE)
source("../lib/Breakpoints.R")
source("../lib/Features.R")
source("../lib/Preprocess.R")
Rcpp::sourceCpp("../lib/edlibFunction.cpp")
Rcpp::sourceCpp("../lib/eigenMatrixMultiply.cpp")

# import nullomer sequences from human genome
# Downloaded from https://www.nullomers.org/ on 1st Feb 2024.
# From paper: https://doi.org/10.1093/nar/gkab139.
df <- fread("./data/Genomic_MAWs.tsv")
ind <- grepl("Homo sapiens", df$`Species name`, ignore.case = TRUE)
df <- df[ind,]
nullomers <- df$MAW

genome <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0

# randomly generate sequences of fixed size of only nullomer sequences
sample_size <- 1000
seed <- 1234
k <- 8
break_score <- "zscore"
break_type <- "Biological"
ranges <- NULL
kmer_window <- 5
cores <- 1

round_to_nearest_even <- function(x) round(x/2)*2

set.seed(seed)
all_sequences <- pbapply::pblapply(1:repeats, function(x){
    nullomer_subset <- sample(
        nullomers, 
        size = sample_size, 
        replace = TRUE
    )
    
    nullomer_sequence <- paste0(
        nullomer_subset, 
        collapse = ""
    )

    letter_counts <- Biostrings::letterFrequency(
        Biostrings::DNAStringSet(nullomer_sequence),
        letters = "ACGT", OR = 0
    )
    letter_counts <- letter_counts * 100 / width(nullomer_sequence)

    return(list(
        "Nullomer" = nullomer_sequence,
        "base_content" = letter_counts
    ))
})
all_sequences <- do.call(c, all_sequences)

nullomer_ind <- which("Nullomer" == names(all_sequences))
base_ind <- which("base_content" == names(all_sequences))

nullomer_base_content <- do.call(rbind, all_sequences[base_ind])
nullomer_base_range <- round_to_nearest_even(nullomer_base_content)

# symmetrical ones, sum of 3 mutation rates, average heptameric to octameric
nullomer_sequence <- unlist(as(
    all_sequences[nullomer_ind], "DNAStringSetList"
), use.names = FALSE)
names(nullomer_sequence) <- paste0("Nullomer_", 1:length(nullomer_sequence))

# subsample lots of sequences from the human genome, retain those with the same GC content
genome_attr <- attr(genome, "seqinfo")
genome_lens <- attr(genome_attr, "seqlengths")[1:22]
N_subsample <- 1000000
seq_len <- round_to_nearest_even(mean(width(nullomer_sequence)))

which_chr <- sample(x = 1:22, size = N_subsample, replace = TRUE)
sample_dt <- data.table(seqnames = 1:22, len = genome_lens)
sample_dt <- sample_dt[match(which_chr, sample_dt$seqnames)]
start_pos <- sapply(1:nrow(sample_dt), function(x){
    sample(x = 1:(sample_dt$len[x]-seq_len), size = 1)
})
sample_dt[, `:=`(
    start = start_pos, 
    width = seq_len,
    len = NULL
)]
setorder(sample_dt, seqnames)
sample_dt[, seqnames := seqnames]
seq_names <- paste0("Genome_", sample_dt$seqnames, "_", sample_dt$start)

# extract sequences
sample_dt <- plyranges::as_granges(sample_dt)
ref_sequences <- Biostrings::getSeq(genome, sample_dt)
names(ref_sequences) <- seq_names
ref_sequences <- unique(ref_sequences)

# filter for the sequences within the tolerance of base composition values
ref_sequence_base <- Biostrings::letterFrequency(
    Biostrings::DNAStringSet(ref_sequences),
    letters = "ACGT", OR = 0
)
ref_sequence_base <- ref_sequence_base * 100 / width(ref_sequences)
ref_sequence_base <- round_to_nearest_even(ref_sequence_base)
unique_nullomer_base_range <- unique(nullomer_base_range)

to_keep <- which(
    (ref_sequence_base[,"A"] %in% unique_nullomer_base_range[,"A"]) &
    (ref_sequence_base[,"T"] %in% unique_nullomer_base_range[,"T"]) &
    (ref_sequence_base[,"G"] %in% unique_nullomer_base_range[,"G"]) &
    (ref_sequence_base[,"C"] %in% unique_nullomer_base_range[,"C"]) 
)

ref_sequences <- ref_sequences[to_keep]
seq_names <- seq_names[to_keep]

# retain the same number of samples as nullomer sequences
ref_sequences <- sample(ref_sequences, size = repeats, replace = FALSE)

# combine both sequences
sequence <- c(nullomer_sequence, ref_sequences)

# save results
Biostrings::writeXStringSet(
    sequence,
    filepath = "./data/nullomer_and_genome.fasta.gz",
    format = "fasta",
    compress = TRUE
)

# extract features
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
        RNAfold.CALL = "",
        maxloopsize = 12,
        FEAT_DNA_SHAPE = FALSE
    )

    if(features$out == -1){
        cat("\n")
        cat("All short-range k-mers are Ns. Starting next chunk.")
        cat("\n")
        next
    }

    file_label <- paste0("../data/nullomer_and_genome/", preprocess$label)
    unlink(file_label, recursive = TRUE)
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