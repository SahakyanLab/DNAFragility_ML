args <- commandArgs(trailingOnly = TRUE)
which_chunk <- as.numeric(args[1])

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

# source functions
setwd("/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/05_DeltaFragility")
pbapply::pboptions(char = "=", type = "txt")
options(future.seed = TRUE)
source("../lib/Breakpoints.R")
source("../lib/Features.R")
source("../lib/Preprocess.R")
Rcpp::sourceCpp("../lib/edlibFunction.cpp")
Rcpp::sourceCpp("../lib/eigenMatrixMultiply.cpp")

# import data set
clinvar_processed <- arrow::read_parquet("./data/all_labels.parquet")

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

# save individual sequences
before_sequences <- Biostrings::DNAStringSet(clinvar_processed$str_before)
names(before_sequences) <- labels <- clinvar_processed$labels

after_sequences <- Biostrings::DNAStringSet(clinvar_processed$str_after)
names(after_sequences) <- labels <- clinvar_processed$labels

sequences <- list(
    before = before_sequences,
    after = after_sequences
)

labels_dt <- clinvar_processed[, .(
    labels, 1:nrow(clinvar_processed), GeneSymbol, 
    Type, ClinicalSignificance, seqnames
)]
fwrite(labels_dt, "./data/labels.csv")

start <- seq(from = 1, to = length(labels), by = 10000)
end <- seq(from = 10000, to = length(labels), by = 10000)
end <- c(end, length(labels))
df_chunks <- data.frame(
    start = as.integer(start), 
    end = as.integer(end)
)
df_chunks_filter <- df_chunks[which_chunk,]

for(x in df_chunks_filter$start:df_chunks_filter$end){
    for(before_after in c("before", "after")){
        cat(paste0("Processing ", before_after, " sequence ", x, "/", df_chunks_filter$end, ".\n"))
        
        preprocess <- Preprocess$new(
            fasta_sequence = sequences[[before_after]][[x]],
            label = labels[x],
            k = k,
            break_type = break_type,
            ranges = ranges
        )
        preprocess$preprocess_sequences()

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

        file_label <- paste0(
            "/media/hert6114/Paddy_6TB/ProjectBoard_Patrick/",
            "04_DNAFragility/data/deltafragility/", 
            preprocess$label
        )
        dir.create(
            path = file_label,
            showWarnings = FALSE,
            recursive = TRUE
        )

        write_parquet(
            features$feature_matrix, 
            paste0(file_label, "/", before_after, "_feature_matrix.parquet")
        )
        write_parquet(
            features$true_bp_table, 
            paste0(file_label, "/", before_after, "_positions.parquet")
        )
    }
}