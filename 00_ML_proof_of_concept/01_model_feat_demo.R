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
suppressPackageStartupMessages(suppressWarnings(library(DNAshapeR)))

# global variables
args <- commandArgs(trailingOnly = TRUE)
exp = as.character(args[1])
kmer_window = as.numeric(args[2])
crash_test = as.logical(args[3])
only_breaks = as.logical(args[4])
my.path <- as.character(args[5])
setwd(my.path)
RNAfold_path = as.character(args[6])
fast_matrix = as.logical(args[7])

# exp="K562_DMSO_DSBs"
# kmer_window=3
# crash_test=TRUE
# only_breaks=TRUE
# my.path="/Users/paddy/Documents/DPhil/DNAFragility/00_ML_proof_of_concept"
# RNAfold_path="RNAfold"
# fast_matrix=FALSE

# source functions
pbapply::pboptions(char = "=", type = "timer")
options(future.seed = TRUE)
source("../01_LGBM_FullGenome/lib/Breakpoints.R")
source("../01_LGBM_FullGenome/lib/Features.R")
Rcpp::sourceCpp("../lib/edlibFunction.cpp")

# if relevant files don't exist, don't run below.
if(fast_matrix) Rcpp::sourceCpp("../lib/eigenMatrixMultiply.cpp")

k = 8
statistic = "mean"
maxloopsize = 12
seed = 1234
break_score = "zscore"
true_prop = 1
num_cores = 1

if(crash_test){
    FEAT_G4_REGEX = FALSE
    FEAT_GC_COUNT = FALSE
    FEAT_VIENNA_RNA = FALSE
    FEAT_DNA_SHAPE = FALSE
} else {
    FEAT_G4_REGEX = TRUE
    FEAT_GC_COUNT = TRUE
    FEAT_VIENNA_RNA = TRUE
    FEAT_DNA_SHAPE = TRUE
}

future::plan(future::multicore)
data.table::setDTthreads(threads = num_cores)

# experiment specific
assembly = "hg19"
exp.to.remove.in.table = ""

# # get short, medium, and long-range for this experiment
df_ranges <- fread(paste0(
    "../data/ranges/",
    "kmer_8_Ranges_cutoffs_from_clustering_all-exp.csv"
))

ranges <- as_tibble(df_ranges) %>% 
    dplyr::rename(full_exp_name = exp) %>%  
    dplyr::filter(stringr::str_detect(
        string = full_exp_name,
        pattern = paste0(
            "K562_Top2_mediated_DSBs/",
            ifelse(grepl("DMSO", exp), "DMSO", "ETO")
        )
    )) %>% 
    dplyr::filter(rowid == "ranges") %>%
    dplyr::select(dplyr::contains("range")) %>% 
    tidyr::pivot_longer(cols = everything()) %>% 
    dplyr::filter(!is.na(value)) %>% 
    dplyr::arrange(value) %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarise(value = mean(value)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(value) %>% 
    dplyr::pull(value)

# in case I need to run more than one experiment.
set.seed(seed)
rng <- sample(100000, size = 20, replace = FALSE)
seeds <- rng[1]

features <- Features$new(
    fast_matrix = fast_matrix,
    k = k, 
    exp = exp, 
    seed = seeds, 
    break_score = break_score,
    assembly = assembly,
    scores_with_kmers = FALSE,
    exp_to_remove_in_table = exp.to.remove.in.table,
    regression = FALSE,
    get_controls = TRUE,
    true_prop = true_prop,
    which_chr = 1:22,
    only_breaks = only_breaks
)

#' Logic flow:
#' 1. Get all features and the full data set.
#' 2. Split the full data set into 70% training and 30% testing.
#' 3. Keep the columns of interest.
#' 4. Fit model to 70% training with the subset of columns.

# 1. get all features and the full data set
features$get_features(
    break_type = "High_frequency", ranges = ranges,
    FEAT_LONGRANGE_PREPARAMS = TRUE, 
    FEAT_G4_REGEX = FEAT_G4_REGEX, g4_type = "GPQS", 
    FEAT_GC_COUNT = FEAT_GC_COUNT,
    FEAT_KMER_COUNTS = TRUE, kmer_window = kmer_window,
    FEAT_VIENNA_RNA = FEAT_VIENNA_RNA, sliding_window = NULL, 
    nuc_type = "DNA", 
    RNAfold.CALL = RNAfold_path,
    maxloopsize = maxloopsize,
    FEAT_DNA_SHAPE = FEAT_DNA_SHAPE,
    SAVE_OUTPUT = FALSE
)

dir.create("../data/ML_demo/", showWarnings = FALSE)
file_name <- paste0(
    "../data/ML_demo/", 
    exp, "_",
    "kmerwindow-", kmer_window,
    "crash_test-", crash_test,
    "only_breaks-", only_breaks,
    "kmer-", features$k, "_",
    "seed-", features$seed, "_", 
    break_score, "_"
)

# save feature_matrix
arrow::write_parquet(
    features$feature_matrix,
    paste0(file_name, "features.parquet")
)

# save column names
df_col_names <- lapply(1:length(features$column_names), function(x){
    return(tibble(
        group = names(features$column_names)[x],
        columns = features$column_names[[x]]
    ))
})
df_col_names <- do.call(dplyr::bind_rows, df_col_names)

df_col_names <- as_tibble(df_col_names) %>% 
    dplyr::mutate(
        group = dplyr::case_when(
            grepl("dEhof", columns) ~ "QM_PARAMETERS",
            grepl("G4", columns) ~ "G4MAP",
            TRUE ~ group
        )
    ) %>% 
    as.data.table()

fwrite(df_col_names, paste0(file_name, "colnames.csv"))

# 2. Split the full data set into 70% training and 30% testing.
train_ind <- caret::createDataPartition(
    features$feature_matrix$predictor, 
    p = 0.7, 
    list = FALSE,
    times = 1
)
train_matrix <- features$feature_matrix[train_ind, ]
test_matrix <- features$feature_matrix[-train_ind, ]
arrow::write_parquet(train_matrix, paste0(file_name, "train.parquet"))
arrow::write_parquet(test_matrix, paste0(file_name, "test.parquet"))