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
repeats <- 10
seed <- 1234
k <- 8
break_score <- "zscore"
break_type <- "Biological"
ranges <- NULL
kmer_window <- 5
cores <- 1

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

    random_shuffle <- strsplit(
        nullomer_sequence, 
        split = ""
    )[[1]]
    random_shuffles <- paste0(
        sample(random_shuffle), 
        collapse = ""
    )

    letter.counts <- Biostrings::letterFrequency(
        Biostrings::DNAStringSet(random_shuffles),
        letters = "ACGT", OR = 0
    )
    letter.counts.norm <- letter.counts / width(random_shuffles)
    gc.content <- letter.counts.norm[, "G"]+letter.counts.norm[, "C"]

    return(list(
        "Nullomer" = nullomer_sequence,
        "Shuffled" = random_shuffles,
        "GC_content" = gc.content
    ))
})
all_sequences <- do.call(c, all_sequences)

nullomer_ind <- which("Nullomer" == names(all_sequences))
shuffle_ind <- which("Shuffled" == names(all_sequences))
gc_ind <- which("GC_content" == names(all_sequences))

# unlist(all_sequences[gc_ind], use.names = FALSE) %>% range()
# symmetrical ones, sum of 3 mutation rates, average heptameric to octameric

nullomer_sequence <- unlist(as(
    all_sequences[nullomer_ind], "DNAStringSetList"
), use.names = FALSE)
names(nullomer_sequence) <- paste0("Nullomer_", 1:length(nullomer_sequence))

shuffle_sequence <- unlist(as(
    all_sequences[shuffle_ind], "DNAStringSetList"
), use.names = FALSE)
names(shuffle_sequence) <- paste0("Shuffle_", 1:length(shuffle_sequence))

sequence <- c(nullomer_sequence, shuffle_sequence)

# save results
Biostrings::writeXStringSet(
    sequence,
    filepath = "./data/nullomer_and_shuffle.fasta.gz",
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

    file_label <- paste0("../data/nullomer_and_shuffle/", preprocess$label)
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