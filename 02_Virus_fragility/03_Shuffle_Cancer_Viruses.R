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
RNAfold_path <- as.character(args[2])
setwd(my.path)

pbapply::pboptions(char = "=", type = "txt")
options(future.seed = TRUE)
source("../lib/Breakpoints.R")
source("../lib/Features.R")
source("../lib/Preprocess.R")
Rcpp::sourceCpp("../lib/edlibFunction.cpp")
Rcpp::sourceCpp("../lib/eigenMatrixMultiply.cpp")

# import virus sequences and labels
sequences <- Biostrings::readDNAStringSet(filepath = "./data/virus_sequences.fasta.gz")
labels <- fread("./data/virus_labels.csv", showProgress = FALSE)

# select the cancer-associated virus sequences
labels <- as_tibble(labels) %>% 
    dplyr::filter(!grepl(paste0(
        "ORF|tissue|Ori|origin|hunan|protein|^_|esterase|domain|",
        "terminus|FMDV|junction|primer|primary|neuraminidase"
    ), title), ignore.case = TRUE) %>% 
    dplyr::mutate(
        species = gsub(" ", "_", organism),
        species = gsub("^(.*virus).*", "\\1", species),
        species = stringr::str_to_title(species)
    ) %>%  
    dplyr::mutate(
        species = dplyr::case_when(
            species == "Hepatitis_B_virus" ~ "Hepatitis_B_virus",
            species == "Hepatitis_c_virus" ~ "Hepatitis_C_virus",
            species == "Hepacivirus" ~ "Hepatitis_C_virus",
            species == "Epstein_Barr_virus" ~ "Epstein_Barr_virus",
            grepl("Epstein", title, ignore.case = TRUE) ~ "Epstein_Barr_virus",
            TRUE ~ species
        )
    ) %>% 
    dplyr::select(-c(accession, organism))

ind_of_int <- which(grepl(paste0(
    "^Hepatitis_B_viruscomplete$|^Hepatitis_C_virus_HCVcomplete$|",
    "^Human_papillomavirus1$|^Epstein_Barr_virusEBV$"
), labels$title, ignore.case = TRUE))[1:4]

sequences <- sequences[ind_of_int]
labels <- labels[ind_of_int,]

# for each virus, take the sequence and shuffle the bases randomly 100 times. 
repeats <- 100
all_sequences <- pbapply::pblapply(1:nrow(labels), function(row){
    this_seq <- sequences[row]
    this_label <- labels[row,]
    bases <- unlist(strsplit(paste0(this_seq), ""))

    set.seed(1234)
    shuffled_seq <- sapply(1:repeats, function(x){
        shuffled_bases <- sample(bases)
        return(paste(shuffled_bases, collapse = ""))
    })
    all_sequences <- c(this_seq, shuffled_seq)
    all_sequences <- DNAStringSet(all_sequences)
    names(all_sequences) <- paste0(this_label$species, "_", 0:repeats)

    return(all_sequences)
})
all_sequences <- do.call(c, all_sequences)

Biostrings::writeXStringSet(
    all_sequences,
    filepath = "./data/virus_sequences_shuffle.fasta.gz", 
    format = "fasta",
    compress = TRUE
)
labels <- data.table(
    label = names(all_sequences),
    length = width(all_sequences)
)
fwrite(
    labels, 
    "./data/virus_labels_shuffle.csv",
    showProgress = FALSE
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

for(x in 1:length(all_sequences)){
    cat(paste0("Processing sequence ", x, "/", length(all_sequences), ".\n"))
    
    preprocess <- Preprocess$new(
        fasta_sequence = all_sequences[[x]],
        label = names(all_sequences[x]),
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
        RNAfold.CALL = RNAfold_path,
        maxloopsize = 12,
        FEAT_DNA_SHAPE = FALSE
    )

    if(features$out == -1){
        cat("\n")
        cat("All short-range k-mers are Ns. Starting next chunk.")
        cat("\n")
        next
    }

    file_label <- paste0("../data/human_viruses_shuffle/", preprocess$label)
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