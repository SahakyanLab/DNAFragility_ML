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

# import long-range sequence effect
range_cutoff <- fread(
    "../data/Ranges_cutoffs.csv", 
    select = c("Cluster", "kmer_8")
)
LR <- range_cutoff[Cluster == "long.range", "kmer_8"][[1]]
round_to_nearest_even <- function(x) round(x/2)*2
LR <- round_to_nearest_even(LR / 2)

# randomly select sequences from human genome and get its reverse complement
genome <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0

# import nullomer sequences from human genome
# Downloaded from https://www.nullomers.org/ on 1st Feb 2024.
# From paper: https://doi.org/10.1093/nar/gkab139.
df <- fread("./data/Genomic_MAWs.tsv")
ind <- grepl("Homo sapiens", df$`Species name`, ignore.case = TRUE)
df <- df[ind,]
nullomers <- df$MAW

# find sequence in human genome that match by 1 mismatch.
all_sequences <- pbapply::pblapply(1:length(nullomers), function(x){
    query_nullomer_string <- nullomers[x]
    query_nullomer <- Biostrings::DNAString(query_nullomer_string)

    matches <- Biostrings::vmatchPattern(
        pattern = query_nullomer, 
        subject = genome, 
        max.mismatch = 1,
        min.mismatch = 1
        # with.indels = TRUE
    )

    matches <- matches %>%
        dplyr::filter(seqnames %in% 1:22)
        
    subset_sequence <- getSeq(genome, matches)

    # find position of mismatch
    matched_sequences <- paste0(subset_sequence)
    
    find_diff_ind <- function(x, y){
        len1 <- nchar(x)
        len2 <- nchar(y)
        min_len <- min(len1, len2)

        len1 <- strsplit(x, "")[[1]][1:min_len]
        len2 <- strsplit(y, "")[[1]][1:min_len]
        
        diff_ind <- which(len1 != len2)[1]
        if(is.na(diff_ind)) diff_ind <- min_len + 1
            
        return(diff_ind)
    }

    ind_match <- mapply(
        find_diff_ind, 
        query_nullomer_string, matched_sequences,
        USE.NAMES = FALSE
    )

    # extract the long-range sequence effects
    LR_grange_left <- matches %>% 
        dplyr::mutate(end = start, start = start - LR)
    dnastring_left <- getSeq(genome, LR_grange_left)

    LR_grange_right <- matches %>% 
        dplyr::mutate(start = end, end = start + LR)
    dnastring_right <- getSeq(genome, LR_grange_right)

    # replace mismatched sequence in the middle with nullomer
    dnastring_with_nullomer <- Biostrings::DNAStringSet(paste0(
        dnastring_left, 
        query_nullomer,
        dnastring_right
    ))
    names(dnastring_with_nullomer) <- paste0(
        "QuerySequence_", query_nullomer_string,
        "_Nullomer_", 
        "_sample_", 1:length(dnastring_with_nullomer)
    )

    # extract real sequences
    dnastring_with_mismatch <- Biostrings::DNAStringSet(paste0(
        dnastring_left, 
        subset_sequence,
        dnastring_right
    ))
    names(dnastring_with_mismatch) <- paste0(
        "QuerySequence_", query_nullomer_string,
        "_Present_", 
        "_sample_", 1:length(dnastring_with_mismatch)
    )

    return(list(
        "Nullomer" = dnastring_with_nullomer, 
        "Present" = dnastring_with_mismatch,
        "start" = LR + ind_match
    ))
})
all_sequences <- do.call(c, all_sequences)

nullomer_ind <- which("Nullomer" == names(all_sequences))
present_ind <- which("Present" == names(all_sequences))
start_ind <- which("start" == names(all_sequences))

sequence_with_nullomer <- unlist(as(
    all_sequences[nullomer_ind], "DNAStringSetList"
), use.names = FALSE)
sequence_without_nullomer <- unlist(as(
    all_sequences[present_ind], "DNAStringSetList"
), use.names = FALSE)

sequences <- list(
    Nullomer = sequence_with_nullomer,
    Present = sequence_without_nullomer
)
start_pos <- unlist(all_sequences[start_ind], use.names = FALSE)

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