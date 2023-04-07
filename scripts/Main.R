setwd("/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/scripts/")

# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(gtools)))
suppressPackageStartupMessages(suppressWarnings(library(matrixStats)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(DNAshapeR)))
data.table::setDTthreads(threads = 1) # prevents segmentation faults

# source functions
source("../lib/Breakpoints.R")
source("../lib/Features.R")

# global variables
k = 8
statistic = "mean"
exp = "K562_Top2_mediated_DSBs"
exp.to.remove.in.table = "K562_cells_Top2_ETO"
maxloopsize = 12
seed = 1234
break_score = "zscore"

# extract features
set.seed(seed)
rng <- sample(100000, size = 1, replace = FALSE)

for(seeds in rng){
    features <- Features$new(
        k = k, 
        exp = exp, 
        seed = seeds, 
        break_score = break_score,
        assembly = "hg19",
        scores_with_kmers = FALSE,
        exp_to_remove_in_table = exp.to.remove.in.table
    )

    features$get_features(
        break_type = "biological", ranges = NULL,
        FEAT_G4_REGEX = TRUE, g4_type = "GPQS", 
        FEAT_GC_COUNT = TRUE,
        FEAT_KMER_COUNTS = TRUE, kmer_window = 3,
        FEAT_VIENNA_RNA = TRUE, sliding_window = NULL, 
        nuc_type = "DNA", 
        RNAfold.CALL = "/home/imm/hert6114/anaconda3/bin/RNAfold",
        maxloopsize = maxloopsize,
        FEAT_DNA_SHAPE = TRUE,
        SAVE_OUTPUT = TRUE
    )    
}

###########################################################
#' uncomment the below if going line-by-line in code
# self=private=NULL
# self$k <- k
# private$bp_dir <- paste0("../data/experiments/", exp)
# self$seed <- seed
# assembly <- "hg19"
# private$exp_to_remove_in_table <- exp.to.remove.in.table
# scores_with_kmers = FALSE
# break_type = "biological"
# ranges = NULL
# g4_type = "GPQS"
# kmer_window = 3
# nuc_type = "DNA"
# crash_test = FALSE
# self=features
# private$break_score <- break_score
# private$break_score_all <- FALSE
# private$exp=exp
###########################################################