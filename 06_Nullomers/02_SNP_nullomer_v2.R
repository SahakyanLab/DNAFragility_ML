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

# import long-range sequence effect
range_cutoff <- fread(
    "../data/Ranges_cutoffs.csv", 
    select = c("Cluster", "kmer_8")
)
LR <- range_cutoff[Cluster == "long.range", "kmer_8"][[1]]
round_to_nearest_even <- function(x) round(x/2)*2
LR <- round_to_nearest_even(LR / 2)

# randomly select sequences from human genome
seq_len <- 10000
N_sample <- 10000
seed <- 1234
genome <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0

# get the start and end positions of each chromosome
refseq.table <- as.data.frame(genome@seqinfo)
refseq.table <- refseq.table[grepl(
    pattern = "^([1-9]|1[0-9]|2[0-2])$", 
    x = rownames(refseq.table)
),]
refseq.table <- data.table(
    seqnames = rownames(refseq.table),
    len = refseq.table$seqlengths - LR
)

# import nullomer sequences from human genome
# Downloaded from https://www.nullomers.org/ on 1st Feb 2024.
# From paper: https://doi.org/10.1093/nar/gkab139.
df <- fread("./data/Genomic_MAWs.tsv")
ind <- grepl("Homo sapiens", df$`Species name`, ignore.case = TRUE)
df <- df[ind,]
nullomers <- df$MAW

# randomly subsample regions from the genome
set.seed(seed)
sample_dt <- refseq.table[sample(x = .N, size = N_sample, replace = TRUE)]
sample_dt[, `:=`(
    start = sapply(len, function(l) sample(x = LR:(l-1), size = 1)),
    nullomer = sample(nullomers, size = N_sample, replace = TRUE),
    len = NULL
)]
setorder(sample_dt, seqnames, start)
sample_dt <- as_tibble(sample_dt) %>% 
    dplyr::mutate(
        seqnames = seqnames, 
        nullomer_len = nchar(nullomer),
        end = start + nullomer_len - 1
    )

# extract left and right genome sequences
LR_grange_left <- sample_dt %>% 
    dplyr::mutate(end = start, start = start - LR) %>% 
    plyranges::as_granges(.)
dnastring_left <- getSeq(genome, LR_grange_left)

LR_grange_right <- sample_dt %>% 
    dplyr::mutate(start = end, end = start + LR) %>% 
    plyranges::as_granges(.)
dnastring_right <- getSeq(genome, LR_grange_right)

# extract middle sequence
middle_seq <- sample_dt %>% 
    plyranges::as_granges(.)
dnastring_middle <- getSeq(genome, middle_seq)

# replace mismatched sequence in the middle with nullomer
sequence_with_nullomer <- Biostrings::DNAStringSet(paste0(
    dnastring_left, 
    sample_dt$nullomer,
    dnastring_right
))
names(sequence_with_nullomer) <- paste0(
    "Nullomer_sample_", 
    1:length(sequence_with_nullomer)
)

# extract real sequences
sequence_without_nullomer <- Biostrings::DNAStringSet(paste0(
    dnastring_left, 
    dnastring_middle,
    dnastring_right
))
names(sequence_without_nullomer) <- paste0(
    "Present_sample_", 
    1:length(sequence_without_nullomer)
)

sequences <- list(
    Nullomer = sequence_with_nullomer,
    Present = sequence_without_nullomer
)
start_pos <- LR + ceiling(sample_dt$nullomer_len / 2)

# save results
Biostrings::writeXStringSet(
    c(sequence_with_nullomer, sequence_without_nullomer),
    filepath = "./data/SNP_nullomer_sequence.fasta.gz",
    format = "fasta",
    compress = TRUE
)

# extract features
seed <- 1234
k <- 8
break_score <- "zscore"
break_type <- "Biological"
SV_type <- "SNP"
ranges <- NULL
kmer_window <- 5
cores <- 1

future::plan(future::multicore)
data.table::setDTthreads(threads = cores)
set.seed(seed)
rng <- sample(100000, size = 20, replace = FALSE)
seeds <- rng[1]

file_label <- paste0("../data/SNP_nullomer_sequence_SNPs/SNPs")
unlink(file_label, recursive = TRUE)
dir.create(
    path = file_label,
    showWarnings = FALSE,
    recursive = TRUE
)

for(sequence_type in c("Nullomer", "Present")){
    cat(paste0("Processing ", sequence_type, " sequences.\n"))

    preprocess <- Preprocess$new(
        fasta_sequence = sequences[[sequence_type]],
        label = names(sequences[[sequence_type]]),
        k = k,
        break_type = break_type,
        SV_type = SV_type,
        ranges = ranges
    )
    preprocess$preprocess_sequences()

    preprocess$df_bp <- data.table(start.pos = start_pos)
    
    features <- Features$new(
        fasta_sequence = preprocess$fasta_sequence,
        label = preprocess$label,
        k = k, 
        seed = seeds, 
        break_score = break_score,
        df_bp = preprocess$df_bp,
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

    print(dim(features$feature_matrix))

    write_parquet(
        features$feature_matrix, 
        paste0(file_label, "/", sequence_type, "_feature_matrix.parquet")
    )
    write_parquet(
        features$true_bp_table, 
        paste0(file_label, "/", sequence_type, "_positions.parquet")
    )
}