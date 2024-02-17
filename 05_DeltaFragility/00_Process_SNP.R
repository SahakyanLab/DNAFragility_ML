# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.UCSC.hg38)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.UCSC.hg19)))
suppressPackageStartupMessages(suppressWarnings(library(arrow)))
pbapply::pboptions(char = "=", type = "txt")

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

hg19_chain <- import.chain("../data/liftover/hg19ToHg38.over.chain")
hg38 <- BSgenome.Hsapiens.UCSC.hg38
hg19 <- BSgenome.Hsapiens.UCSC.hg19

get_ref_genome_lens <- function(version){
    refseq <- switch(version,
        "hg19" = BSgenome.Hsapiens.UCSC.hg19,
        "hg38" = BSgenome.Hsapiens.UCSC.hg38
    )

    # get the start and end positions of each chromosome
    refseq.table <- as.data.frame(refseq@seqinfo)
    refseq.table <- refseq.table[grepl(
        pattern = "^chr([1-9]|1[0-9]|2[0-2])$", 
        x = rownames(refseq.table)
    ), ]
    refseq.table <- data.table(
        chromName = rownames(refseq.table),
        chromSize = refseq.table$seqlengths
    )

    return(refseq.table)
}
hg19_table <- get_ref_genome_lens(version = "hg19")
hg38_table <- get_ref_genome_lens(version = "hg38")

#' Download ClinVar dataset.
#' Download: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz.
#' Accessed Wednesday December 6, 2023.
t1 <- Sys.time()
cur.msg <- "Processing and cleaning the ClinVar variant summary table"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

clinvar_table <- fread('../04_ClinVar/variant_summary.txt', showProgress = FALSE)

# only keep single nucleotide mutations
clinvar <- clinvar_table["single nucleotide variant" == Type,]
clinvar <- clinvar[grepl(
    paste0("pathogenic|benign"), ClinicalSignificance, ignore.case = TRUE
)]

# columns to keep
clinvar <- clinvar[, .(
    Type, GeneSymbol, ClinicalSignificance, OriginSimple, 
    Assembly, Chromosome, Start, Stop,
    ReferenceAlleleVCF, AlternateAlleleVCF
)]
setnames(clinvar, c(
    "Type", "GeneSymbol", "ClinicalSignificance", 
    "Origin", "Assembly", "seqnames", "start", "end",
    "str_before", "str_after"
))
clinvar[, seqnames := paste0("chr", seqnames)]

# keep GRCh37 and GRCh38 genome versions
clinvar <- clinvar[grepl("GRCh37|GRCh38", Assembly)]

# separate out start and end positions
clinvar <- dplyr::bind_rows(
        as_tibble(clinvar) %>% 
            dplyr::select(-end),
        as_tibble(clinvar) %>% 
            dplyr::select(-start) %>%
            dplyr::rename(start = end)
    ) %>% 
    dplyr::filter(
        !grepl(" ", GeneSymbol) & 
        !grepl(";", GeneSymbol) &
        GeneSymbol != "-" & 
        grepl("^chr([1-9]|1[0-9]|2[0-2])$", seqnames) &
        nchar(str_before) == 1 &
        nchar(str_after) == 1
    ) %>%
    dplyr::mutate(
        width = 1,
        ClinicalSignificance = dplyr::case_when(
            grepl("pathogenic", ClinicalSignificance, ignore.case = TRUE) ~ "Pathogenic",
            grepl("benign", ClinicalSignificance, ignore.case = TRUE) ~ "Benign"
        )
    )

# only get unique breakpoint positions by genes
clinvar <- clinvar %>%
    dplyr::group_by(GeneSymbol, ClinicalSignificance, seqnames, start) %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::filter(str_before != "na" | str_after != "na") %>%
    dplyr::select(-Type) %>%
    as.data.table()

# get short range of influence
range_cutoff <- fread("../data/range_effects/MaxValuesFromClustersByType.csv")
ind <- which.max(range_cutoff$long)
LR <- range_cutoff$long[ind]
break_type <- range_cutoff$break_type[ind]

round_to_nearest_even <- function(x) round(x/2)*2
LR <- round_to_nearest_even(LR / 2)

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")

####################################################################################################
# Processing the sequences of all SNVs
####################################################################################################
t1 <- Sys.time()
cur.msg <- "Processing single nucleotide variants"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

clinvar_left <- copy(clinvar); clinvar_right <- copy(clinvar)

clinvar_left[, `:=`(
    end = start - 1,
    start = start - LR - 1,
    width = NULL
)]
clinvar_left[, str_left := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_left)]

clinvar_right[, `:=`(
    start = start + 1,
    width = LR
)]
clinvar_right[, str_right := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_right)]

clinvar[, str_before := paste0(
    clinvar_left$str_left, 
    clinvar$str_before, 
    clinvar_right$str_right
)]
clinvar[, str_after := paste0(
    clinvar_left$str_left, 
    clinvar$str_after, 
    clinvar_right$str_right
)]
clinvar[, `:=`(Assembly = NULL, width = NULL)]

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")

####################################################################################################
# Bind all breakage types together
####################################################################################################
t1 <- Sys.time()
cur.msg <- "Saving all results"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

clinvar[, labels := paste0(
    1:nrow(clinvar), "_",
    clinvar$GeneSymbol, "_",
    clinvar$ClinicalSignificance, "_",
    clinvar$seqnames
)]

# only keep sequences with ATGC bases
str_before_bases <- Biostrings::alphabetFrequency(
    Biostrings::DNAStringSet(
        clinvar$str_before
), as.prob = FALSE)
no_atgc_left <- !(colnames(str_before_bases) %in% c("A","T","G","C","N"))
cols_to_rm_left <- which(rowSums(str_before_bases[, no_atgc_left]) != 0)

str_after_bases <- Biostrings::alphabetFrequency(
    Biostrings::DNAStringSet(
        clinvar$str_after
), as.prob = FALSE)
no_atgc_right <- !(colnames(str_after_bases) %in% c("A","T","G","C","N"))
cols_to_rm_right <- which(rowSums(str_after_bases[, no_atgc_right]) != 0)

cols_to_rm <- union(cols_to_rm_left, cols_to_rm_right)
clinvar <- clinvar[-cols_to_rm,]
clinvar_processed <- copy(clinvar)
clinvar_processed[, `:=`(str_before = NULL, str_after = NULL)]

arrow::write_parquet(clinvar_processed, "./data/all_SNP.parquet")

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")

####################################################################################################
# Extract features only looking at the middle base to evaluate the breakage
####################################################################################################
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(gtools)))
suppressPackageStartupMessages(suppressWarnings(library(purrr)))
suppressPackageStartupMessages(suppressWarnings(library(furrr)))
suppressPackageStartupMessages(suppressWarnings(library(progressr)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))

pbapply::pboptions(char = "=", type = "txt")
options(future.seed = TRUE)
source("../lib/Breakpoints.R")
source("../lib/Features.R")
source("../lib/Preprocess.R")
Rcpp::sourceCpp("../lib/edlibFunction.cpp")
Rcpp::sourceCpp("../lib/eigenMatrixMultiply.cpp")

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
data.table::setDTthreads(threads = 1)
set.seed(seed)
rng <- sample(100000, size = 20, replace = FALSE)
seeds <- rng[1]

# save individual sequences
before_sequences <- Biostrings::DNAStringSet(clinvar$str_before)
names(before_sequences) <- labels <- clinvar$labels

after_sequences <- Biostrings::DNAStringSet(clinvar$str_after)
names(after_sequences) <- labels <- clinvar$labels

sequences <- list(
    before = before_sequences,
    after = after_sequences
)

labels_dt <- clinvar[, .(
    labels, 1:nrow(clinvar), GeneSymbol, 
    ClinicalSignificance, seqnames
)]

#' Logic flow is the following within Preprocess function.
#' 1. Make a new column called start.pos which is the same for all sequences. 
#'      It's just the middle where the SNP happens.
#' 2. Process all sequences at the same time.
#' 3. Extract features.
#' 4. Predict. 

file_label <- paste0(
    "../data/deltafragility_SNPs/SNPs"
)
dir.create(
    path = file_label,
    showWarnings = FALSE,
    recursive = TRUE
)

for(before_after in c("before", "after")){
    cat(paste0("Processing ", before_after, " sequences.\n"))

    preprocess <- Preprocess$new(
        fasta_sequence = sequences[[before_after]],
        label = names(sequences[[before_after]]),
        k = k,
        break_type = break_type,
        SV_type = SV_type,
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

    print(dim(features$feature_matrix))

    write_parquet(
        features$feature_matrix, 
        paste0(file_label, "/", before_after, "_feature_matrix.parquet")
    )
    write_parquet(
        features$true_bp_table, 
        paste0(file_label, "/", before_after, "_positions.parquet")
    )
}