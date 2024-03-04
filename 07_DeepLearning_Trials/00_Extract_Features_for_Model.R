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
regression <- as.logical(args[2])
bw <- as.numeric(args[3])
setwd(my.path)

pbapply::pboptions(char = "=", type = "timer")
options(future.seed = TRUE)
source("../01_LGBM_FullGenome/lib/Breakpoints.R")
source("../01_LGBM_FullGenome/lib/Features.R")
Rcpp::sourceCpp("../01_LGBM_FullGenome/lib/edlibFunction.cpp")
Rcpp::sourceCpp("../01_LGBM_FullGenome/lib/eigenMatrixMultiply.cpp")

extract_features <- function(
    getTable, k, seed, num_cores,
    break_score, statistic, maxloopsize,
    exp, assembly, cols, regression, 
    break_type, long_range, bw, true_prop,
    get_controls, kmer_window, which_chr,
    ranges, file_chunk=1
    ){  
    future::plan(future::multicore)
    data.table::setDTthreads(threads = num_cores)

    # in case I need to run more than one experiment.
    set.seed(seed)
    rng <- sample(100000, size = 20, replace = FALSE)
    seeds <- rng[1]

    features <- Features$new(
        k = k, 
        exp = exp, 
        seed = seeds, 
        break_score = break_score,
        assembly = assembly,
        scores_with_kmers = FALSE,
        regression = regression,
        get_controls = get_controls,
        true_prop = true_prop,
        num_cores = num_cores,
        which_chr = which_chr,
        only_breaks = TRUE
    )

    features$get_features(
        break_type = break_type, ranges = ranges,
        FEAT_LONGRANGE_PREPARAMS = TRUE, 
        FEAT_G4_REGEX = TRUE, g4_type = "GPQS", 
        FEAT_GC_COUNT = TRUE,
        FEAT_KMER_COUNTS = TRUE, kmer_window = kmer_window,
        FEAT_VIENNA_RNA = FALSE, sliding_window = NULL, 
        nuc_type = "DNA", 
        RNAfold.CALL = "/home/imm/hert6114/anaconda3/bin/RNAfold",
        maxloopsize = maxloopsize,
        FEAT_DNA_SHAPE = FALSE,
        SAVE_OUTPUT = FALSE
    )
    if(features$out == -1){
        cat("\n")
        cat("All short-range k-mers are Ns. Starting next chunk.")
        cat("\n")
        return(-1)
    }

    # full data set
    path_to_dir <- "../07_DeepLearning_Trials/data/"

    df_positions <- plyranges::as_granges(features$true_bp_table)
    df_positions <- plyranges::stretch(
        plyranges::anchor_center(df_positions),
        long_range-1
    )
    df_positions <- as.data.table(df_positions)
    df_positions[, `:=`(
        seqnames = as.character(seqnames), 
        width = NULL, 
        strand = NULL
    )]

    # train/test split for LightGBM model training process
    index_2 <- caret::createDataPartition(
        y = features$feature_matrix$predictor,
        p = 0.7,
        list=FALSE
    )
    arrow::write_parquet(
        features$feature_matrix[index_2, ], 
        paste0(
            path_to_dir, "/train_", 
            "regression-", regression, ".parquet"
        )
    )
    arrow::write_parquet(
        features$feature_matrix[-index_2, ], 
        paste0(
            path_to_dir, "/test_", 
            "regression-", regression, ".parquet"
        )
    )

    # for hyena-dna and CNN models
    if(!regression){
        df_positions[, Breaks := ifelse(Breaks == "YES", 1, 0)]
    }
    
    fwrite(
        df_positions[index_2, ], 
        paste0(
            path_to_dir, "/train_", 
            "regression-", regression, ".csv"
        )
    )
    fwrite(
        df_positions[-index_2, ], 
        paste0(
            path_to_dir, "/test_",
            "regression-", regression, ".csv"
        )
    )
}

# global vars
get_controls <- ifelse(regression, FALSE, TRUE)
num_cores <- 1
break_type <- "Biological"
chr <- 1:22

# get biological long range effect
ranges <- fread(paste0(
    "../data/range_effects/MaxValuesFromClustersByType.csv"
), select = c("short", "medium", "long"))
ranges <- apply(ranges, 2, max)

if(regression){
    # Bin width | Overlapping bins | Overlap factor | Remove Zero breaks
    ranges[3] <- long_range <- bw <- as.integer(bw)
    system(paste(
        "Rscript ../05_Cosmic/scripts/02_Bin_breaks_full_genome.R", 
        bw, "FALSE 10 FALSE", my.path, "../05_Cosmic/scripts/"
    ))
    exp <- "FullGenome_DeepLearning"
} else {
    bw <- as.integer(bw)
    long_range <- max(ranges)
    exp <- "COSMIC"
}

extract_features(
    getTable=FALSE,
    k=8,
    seed=1234,
    num_cores=num_cores,
    break_score="zscore",
    statistic="mean",
    maxloopsize=12,
    exp=exp,
    assembly="telo_to_telo",
    cols="COSMIC_columns",
    regression=regression,
    break_type=break_type,
    long_range=long_range,
    ranges=ranges,
    true_prop=1,
    bw=bw,
    get_controls=get_controls,
    kmer_window=3,
    which_chr=chr
)