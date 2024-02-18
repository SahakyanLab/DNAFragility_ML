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
my.path <- as.character(args[3])
setwd(my.path)
fast_matrix <- as.logical(args[5])

pbapply::pboptions(char = "=", type = "timer")
options(future.seed = TRUE)
source("./lib/Breakpoints.R")
source("./lib/Features.R")
Rcpp::sourceCpp("../lib/edlibFunction.cpp")
if(fast_matrix) Rcpp::sourceCpp("../lib/eigenMatrixMultiply.cpp")

# Breaks <- Breakpoints$new(
#     exp="", 
#     seed=1234, 
#     scores_with_kmers=FALSE, 
#     assembly="hg38", 
#     break_score="zscore", 
#     break_score_all=FALSE
# )
# Breaks$get_extended_tables()

extract_features <- function(
    fast_matrix, getTable, k, seed, num_cores,
    break_score, statistic, maxloopsize,
    exp, assembly, cols, regression, 
    break_type, long_range, bw, true_prop,
    get_controls, kmer_window, which_chr,
    training_process, RNAfold_path, file_chunk=1
    ){  
    future::plan(future::multicore)
    data.table::setDTthreads(threads = num_cores)

    # in case I need to run more than one experiment.
    set.seed(seed)
    rng <- sample(100000, size = 20, replace = FALSE)
    seeds <- rng[1]

    # if(getTable){
    #     source("../lib/KmerTable.R")
    #     for(k in seq(2,8,2)){
    #         tables <- KmerTable$new(
    #             k = k, 
    #             exp = "", 
    #             break_score = break_score,
    #             statistic = statistic,
    #             all_exp = TRUE,
    #             group_exp = TRUE
    #         )

    #         tables$generate_querytable()
    #     }
    # }

    features <- Features$new(
        fast_matrix = fast_matrix,
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
        break_type = break_type, 
        ranges = NULL,
        FEAT_LONGRANGE_PREPARAMS = TRUE, 
        FEAT_G4_REGEX = TRUE, g4_type = "GPQS", 
        FEAT_GC_COUNT = TRUE,
        FEAT_KMER_COUNTS = TRUE, kmer_window = kmer_window,
        FEAT_VIENNA_RNA = FALSE, sliding_window = NULL, 
        nuc_type = "DNA", 
        RNAfold.CALL = RNAfold_path,
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

    # # only keep columns of interest
    # features$select_columns(cols = cols)

    # full data set
    path_to_dir <- "../data/models/python/lightgbm"

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

    if(training_process){
        # train/test split for LightGBM model training process
        index_2 <- caret::createDataPartition(
          y = features$feature_matrix$predictor,
          p = 0.7,
          list=FALSE
        )
        arrow::write_parquet(
            features$feature_matrix[index_2, ], 
            paste0(path_to_dir, "/train.parquet")
        )
        arrow::write_parquet(
            features$feature_matrix[-index_2, ], 
            paste0(path_to_dir, "/test.parquet")
        )

        fwrite(
            df_positions, 
            paste0(path_to_dir, "/FullGenome_Ranges_", bw, ".csv"),
            showProgress = FALSE
        )
    } else {
        dir.create(
            path = paste0(path_to_dir, "/chr", which_chr, "/"),
            showWarnings = FALSE,
            recursive = TRUE
        )

        arrow::write_parquet(
            features$feature_matrix, 
            paste0(
                path_to_dir, "/chr", which_chr, "/FullGenome_FeatureMatrix_", 
                bw, "_file_chunk_", file_chunk, ".parquet"
            )
        )

        arrow::write_parquet(
            df_positions, 
            paste0(
                path_to_dir, "/chr", which_chr, "/FullGenome_Ranges_", 
                bw, "_file_chunk_", file_chunk, ".parquet"
            )
        )
    }
}

# global vars
training_process <- as.logical(args[1])
RNAfold_path <- as.character(args[4])
bw <- as.integer(1)
num_cores <- 1
break_type <- "Biological"

# get biological long range effect
df_long_range <- fread(paste0(
    "../data/range_effects/", 
    "MaxValuesFromClustersByType.csv"
))
which_row <- which(df_long_range$break_type == break_type)
long_range <- df_long_range$long[which_row]

if(training_process){
    # use below 2 arguments for train/testing LGBM model
    chr <- 1:22
    exp <- "COSMIC"

    extract_features(
        fast_matrix=fast_matrix,
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
        regression=FALSE,
        break_type=break_type,
        long_range=long_range,
        true_prop=1.5,
        bw=bw,
        get_controls=TRUE,
        kmer_window=5,
        which_chr=chr,
        RNAfold_path=RNAfold_path,
        training_process=training_process
    )
} else {
    # use below 3 arguments for testing on full genome..
    chr <- as.integer(args[2])
    exp <- "FullGenomeChunks"

    if(chr == 1){
        #' params:
        #'  Bin width = 1
        #'  Overlapping bins = TRUE
        #'  Overlap factor = 1
        #'  Remove Zero breaks = FALSE
        # system(paste(
        #     "Rscript ../../05_Cosmic/scripts/02_Bin_breaks_full_genome.R", 
        #     bw, "TRUE 1 FALSE"
        # ))

        message("Already binned genome.")
    }

    # get testing data sets and move each chunk over
    all.files <- list.files(
        path = paste0(
            "../data/experiments/", exp, 
            "/breakpoint_positions",
            "/chr", chr
        ),
        pattern = "chunk_.*\\.csv$",
        full.names = TRUE
    )
    all.files <- stringr::str_sort(all.files, numeric = TRUE)
    chunk_num <- stringr::str_extract(
        string = all.files, 
        pattern = "(?<=chunk_)[0-9]+(?=\\.csv)"
    ) %>% as.integer()

    for(f in 1:length(all.files)){
        print(paste(
            "Running experiment:", all.files[f], 
            "of", max(chunk_num)
        ))

        file.mv <- file.rename(
            all.files[f], 
            paste0(dirname(all.files[f]), ".csv")
        )

        extract_features(
            fast_matrix=fast_matrix,
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
            regression=FALSE,
            break_type=break_type,
            long_range=long_range,
            bw=bw,
            get_controls=FALSE,
            kmer_window=5,
            which_chr=chr,
            training_process=training_process,
            RNAfold_path=RNAfold_path,
            file_chunk=chunk_num[f]
        )

        gc(verbose = FALSE)
    }
}