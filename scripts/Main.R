setwd("/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/scripts/")

# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(gtools)))
suppressPackageStartupMessages(suppressWarnings(library(matrixStats)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(DNAshapeR)))

# source functions
source("../lib/KmerTable.R")
source("../lib/Breakpoints.R")
source("../lib/Features.R")

# global variables
k = 8
statistic = "mean"
exp = "DSBCapture"
# exp = "K562_Top2_mediated_DSBs"
maxloopsize = 12
seed = 1234
break_score = "all"

# main script
# generate query table of k-mers
dir.create("../data/kmertone/QueryTable/", showWarnings = FALSE)
for(table.kmer in c(2,4,6,8)){
    for(scores in c("zscore", "ratio")){
        table <- KmerTable$new(
            k = table.kmer, 
            break_score = scores, 
            statistic = statistic, 
            exp = exp,
            all_exp = FALSE,
            group_exp = TRUE
        )
        table$generate_querytable()
    }
}

# extract features
features <- Features$new(k = k, exp = exp, 
                         seed = seed, 
                         break_score = break_score,
                         assembly = "hg19",
                         scores_with_kmers = FALSE)

features$get_features(
    FEAT_G4_REGEX = TRUE, g4_type = "GPQS", 
    FEAT_GC_COUNT = TRUE,
    FEAT_KMER_COUNTS = TRUE, kmer_window = 3,
    FEAT_VIENNA_RNA = TRUE, sliding_window = NULL, 
    nuc_type = "DNA", 
    RNAfold.CALL = "/home/imm/hert6114/anaconda3/bin/RNAfold",
    FEAT_DNA_SHAPE = TRUE,
    FEAT_TFBS_EUCLIDEAN_DISTANCE = FALSE,
    FEAT_OCCUPANCY_SCORES = TRUE,
    FEAT_QM_PARAMETERS = TRUE,  
    SAVE_OUTPUT = TRUE
)