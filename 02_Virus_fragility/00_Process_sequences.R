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

#' Download history:
#' Reference Viral DataBase (RVDB).
#' Website: https://rvdb.dbi.udel.edu/
#' Download Current Release (v27.0, Sept 20, 2023). 
#' Downloaded ClusteredDB on Wednesday December 6, 2023. 
sequences <- Biostrings::readDNAStringSet("./data/C-RVDBvCurrent.fasta.gz")

# process names
labels <- names(sequences)
split_labels <- strsplit(labels, "\\|")
labels <- data.table(
    accession = sapply(split_labels, `[`, 3),
    title = sapply(split_labels, `[`, 4),
    organism = sapply(split_labels, `[`, 5),
    virus = sapply(split_labels, `[`, 6)
)

# only keep complete genome and complete sequences
labels[, complete_genome := ifelse(
    grepl("complete genome|complete sequence", title), 
    TRUE, FALSE)
]
to_keep <- which((labels$complete_genome) == TRUE | (!is.na(labels$organism)))
labels[, complete_genome := NULL]

labels <- labels[to_keep, ]
sequences <- sequences[to_keep]

# remove potential presence of complete cds, genes, enzymes
labels[, genic := ifelse(
    grepl(paste0(
        "complete cds|gene|enzyme|merase|retro|repeat|site|associate|", 
        "terminal|unknown|untrans|common|intergenic|gene|region|defect|",
        "segment|isolate|unveri|sequencing|nearly|LTR|UTR|PCR|like|promoter|fragment|",
        "relate|sense|anti|clone|similar|structur|sequence|partial|end|cds"
    ), title, ignore.case = TRUE),
    TRUE, FALSE
)]
to_keep <- which(labels$genic == FALSE)
labels[, genic := NULL]

labels <- labels[to_keep, ]
sequences <- sequences[to_keep]

# remove complete genome or compelte sequence from the 
labels[, title := gsub(", complete genome|, complete sequence", "", title)]
labels[, title := gsub(" |\\/", "_", title)]
labels[, title := gsub("[[:punct:]]", "_", title)]
labels[, title := gsub("___+", "", title)]
labels[, title := gsub("__+", "_", title)]
labels[, length := width(sequences)]

# remove genome_assembly if in title
labels[, title := sub("_genome.*", "", title)]

# remove last underscore if exists
labels[, title := sub("_([^_]*)$", "\\1", title)]

# if no organism exists, then extract everything up to and incl. virus
labels[, organism := ifelse(
    is.na(organism), gsub("^(.*virus).*", "\\1", title),
    organism
)]

# final filtering process
labels[, to_discard := ifelse(
    grepl(paste0(
        "ORF|tissue|Ori|origin|hunan|protein|^_|esterase|domain|",
        "terminus|FMDV|junction|primer|primary|neuraminidase"
    ), title, ignore.case = TRUE),
    TRUE, FALSE
)]
to_keep <- which(labels$to_discard == FALSE)
labels[, to_discard := NULL]

labels[, `:=`(virus = NULL, length = NULL)]

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
big_to_small <- order(width(sequences), decreasing = TRUE)
sequences <- sequences[big_to_small]
labels <- labels[big_to_small]
labels[, length := width(sequences)]

rm_dup <- match(unique(labels$title), labels$title)
labels <- labels[rm_dup]
sequences <- sequences[rm_dup]

Biostrings::writeXStringSet(
    sequences,
    filepath = "./data/virus_sequences.fasta.gz", 
    format = "fasta",
    compress = TRUE
)
fwrite(labels, "./data/virus_labels.csv")

max_ind <- length(labels$length)
for(x in 1:max_ind){
    cat(paste0("Processing sequence ", x, "/", max_ind, ".\n"))

    preprocess <- Preprocess$new(
        fasta_sequence = sequences[[x]],
        label = labels$title[x],
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

    file_label <- paste0("../data/human_viruses/", preprocess$label)
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