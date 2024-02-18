# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
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

# only keep insertion / deletion / indels
clinvar <- clinvar_table[grepl(paste0(
    "^Indel$|^Deletion$|^Duplication$|^Insertion$|",
    "^Inversion$|^Translocation$|^Tandem duplication"
    ), Type, ignore.case = TRUE
)]

clinvar <- clinvar[grepl(
    paste0("^Pathogenic|^Likely pathogenic|^Benign|^Likely benign"), 
    ClinicalSignificance, 
    ignore.case = TRUE
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
        grepl("^chr([1-9]|1[0-9]|2[0-2])$", seqnames)
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
    as.data.table()

# get short range of influence
range_cutoff <- fread("../data/Ranges_cutoffs.csv", select = c("Cluster", "kmer_8"))
LR <- range_cutoff[Cluster == "long.range", "kmer_8"][[1]]
round_to_nearest_even <- function(x) round(x/2)*2
LR <- round_to_nearest_even(LR / 2)

#' Thought process:
#' Before mutation: count the number of breaks by gene. No model needs to be applied as the 
#' before mutation sequence is the ground truth experimental detail. In other words, these 
#' sites were *definitely* broken in the experiment. Doesn't make sense to either:
#'      1) use the model to predict the number of breaks along the sequence.
#'      2) use combination of model prediction and true breakpoint locations. 
#'
#' After mutation: modify the sequence based on the type of modification, then expand by the 
#' maximum long range sequence effect, and extract features and predict per sequence. The theory is 
#' that a change in the sequence may cause the region around the modification to become more/less fragile.
clinvar_split <- split(clinvar, by = "Type")

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")

over_LR_length <- rbind(
    as.data.table(table(nchar(clinvar$str_before))), 
    as.data.table(table(nchar(clinvar$str_after)))
    ) %>% 
    dplyr::group_by(V1) %>% 
    dplyr::summarise(N = sum(N)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
        frac = N / sum(N) * 100,
        V1 = as.numeric(V1)
    ) %>% 
    dplyr::filter(V1 > (LR * 2)) %>% 
    dplyr::pull(frac) %>% 
    sum(.)
over_LR_length <- signif(over_LR_length, 3)
cat(paste0("Percentage of SVs longer than the LR seq-context: ", over_LR_length, "%\n"))

####################################################################################################
# Insertion
####################################################################################################
t1 <- Sys.time()
cur.msg <- "Processing sequence variant: insertions"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

clinvar_ins <- clinvar_split[[match("Insertion", names(clinvar_split))]]
clinvar_ins_left <- copy(clinvar_ins); clinvar_ins_right <- copy(clinvar_ins)

clinvar_ins_left[, `:=`(
    end = start,
    start = start - LR + 1,
    width = NULL
)]
clinvar_ins_left[, str_after := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_ins_left)]

clinvar_ins_right[, `:=`(
    # because str_before contains the first letter of the str_after.
    start = start + 1,
    width = LR - nchar(str_after)
)]
clinvar_ins_right[, str_after := {
    if(unique(Assembly) == "GRCh37"){
        paste0(str_after, paste0(getSeq(hg19, plyranges::as_granges(.SD))))
    } else {
        paste0(str_after, paste0(getSeq(hg38, plyranges::as_granges(.SD))))
    }
}, by = .(Assembly), .SDcols = names(clinvar_ins_right)]
clinvar_ins[, str_after := paste0(clinvar_ins_left$str_after, clinvar_ins_right$str_after)]
clinvar_ins[, str_after := substring(text = str_after, first = 1, last = LR * 2)]

clinvar_ins[, `:=`(start = start - LR, width = LR * 2)]
clinvar_ins[, str_before := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_ins)]

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")
####################################################################################################
# Deletion
####################################################################################################
t1 <- Sys.time()
cur.msg <- "Processing sequence variant: deletions"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

clinvar_del <- clinvar_split[[match("Deletion", names(clinvar_split))]]
clinvar_del_left <- copy(clinvar_del); clinvar_del_right <- copy(clinvar_del)

clinvar_del_left[, `:=`(
    # start-1 is the first base of str_before.
    end = start - 2, 
    start = start - LR,
    width = NULL
)]
clinvar_del_left[, str_after := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_del_left)]

clinvar_del_right[, `:=`(
    # because str_before contains the first letter of the str_after.
    start = start + nchar(str_before) - 1,
    width = LR
)]
clinvar_del_right[, str_after := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_del_right)]
clinvar_del[, str_after := paste0(clinvar_del_left$str_after, clinvar_del_right$str_after)]
clinvar_del[, str_after := substring(text = str_after, first = 1, last = LR * 2)]

clinvar_del[, `:=`(start = start - LR, width = LR * 2)]
clinvar_del[, str_before := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_del)]

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")
####################################################################################################
# Indel
# str_before: sequence deleted and str_after sequence inserts. 
####################################################################################################
t1 <- Sys.time()
cur.msg <- "Processing sequence variant: insertion-deletions"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

clinvar_indel <- clinvar_split[[match("Indel", names(clinvar_split))]]
clinvar_indel_left <- copy(clinvar_indel); clinvar_indel_right <- copy(clinvar_indel)

clinvar_indel_left[, `:=`(
    end = start - 2, # start-1 is the first base of str_before.
    start = start - LR,
    width = NULL
)]
clinvar_indel_left[, str_after := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_indel_left)]

clinvar_indel_right[, `:=`(
    # because str_before contains the first letter of the str_after.
    start = start + nchar(str_before),
    width = LR
)]
clinvar_indel_right[, str_after := {
    if(unique(Assembly) == "GRCh37"){
        paste0(str_after, paste0(getSeq(hg19, plyranges::as_granges(.SD))))
    } else {
        paste0(str_after, paste0(getSeq(hg38, plyranges::as_granges(.SD))))
    }
}, by = .(Assembly), .SDcols = names(clinvar_indel_right)]
clinvar_indel[, str_after := paste0(clinvar_indel_left$str_after, clinvar_indel_right$str_after)]
clinvar_indel[, str_after := substring(text = str_after, first = 1, last = LR * 2)]

clinvar_indel[, `:=`(start = start - LR, width = LR * 2)]
clinvar_indel[, str_before := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_indel)]

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")
####################################################################################################
# Duplication
####################################################################################################
t1 <- Sys.time()
cur.msg <- "Processing sequence variant: duplications"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

clinvar_dup <- clinvar_split[[match("Duplication", names(clinvar_split))]]
clinvar_dup_left <- copy(clinvar_dup); clinvar_dup_right <- copy(clinvar_dup)

clinvar_dup_left[, `:=`(
    end = start + nchar(str_after) - 1,
    start = start - LR,
    width = NULL
)]
clinvar_dup_left[, str_after := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_dup_left)]

clinvar_dup_right[, `:=`(
    start = start + nchar(str_after),
    width = LR
)]
clinvar_dup_right[, str_after := {
    if(unique(Assembly) == "GRCh37"){
        paste0(str_after, paste0(getSeq(hg19, plyranges::as_granges(.SD))))
    } else {
        paste0(str_after, paste0(getSeq(hg38, plyranges::as_granges(.SD))))
    }
}, by = .(Assembly), .SDcols = names(clinvar_dup_right)]
clinvar_dup[, str_after := paste0(clinvar_dup_left$str_after, clinvar_dup_right$str_after)]
clinvar_dup[, str_after := substring(text = str_after, first = 1, last = LR * 2)]

clinvar_dup[, `:=`(start = start - LR, width = LR * 2)]
clinvar_dup[, str_before := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_dup)]

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")
####################################################################################################
# Inversion
####################################################################################################
t1 <- Sys.time()
cur.msg <- "Processing sequence variant: inversions"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

clinvar_inv <- clinvar_split[[match("Inversion", names(clinvar_split))]]
clinvar_inv_left <- copy(clinvar_inv); clinvar_inv_right <- copy(clinvar_inv)

clinvar_inv_left[, `:=`(
    end = start - 1,
    start = start - LR,
    width = NULL
)]
clinvar_inv_left[, str_after := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_inv_left)]

clinvar_inv_right[, `:=`(
    start = start + 2,
    width = LR
)]
clinvar_inv_right[, str_after := {
    if(unique(Assembly) == "GRCh37"){
        paste0(str_after, paste0(getSeq(hg19, plyranges::as_granges(.SD))))
    } else {
        paste0(str_after, paste0(getSeq(hg38, plyranges::as_granges(.SD))))
    }
}, by = .(Assembly), .SDcols = names(clinvar_inv_right)]
clinvar_inv[, str_after := paste0(clinvar_inv_left$str_after, clinvar_inv_right$str_after)]
clinvar_inv[, str_after := substring(text = str_after, first = 1, last = LR * 2)]

clinvar_inv[, `:=`(start = start - LR, width = LR * 2)]
clinvar_inv[, str_before := {
    if(unique(Assembly) == "GRCh37"){
        paste0(getSeq(hg19, plyranges::as_granges(.SD)))
    } else {
        paste0(getSeq(hg38, plyranges::as_granges(.SD)))
    }
}, by = .(Assembly), .SDcols = names(clinvar_inv)]

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")
####################################################################################################
# Bind all breakage types together
####################################################################################################
t1 <- Sys.time()
cur.msg <- "Binding all tables together and saving results"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

clinvar_processed <- rbindlist(list(
    clinvar_ins, clinvar_del, 
    clinvar_indel, clinvar_dup, 
    clinvar_inv
))
clinvar_processed[, width := NULL]
clinvar_processed[, labels := paste0(
    1:nrow(clinvar_processed), "_",
    clinvar_processed$GeneSymbol, "_",
    clinvar_processed$Type, "_",
    clinvar_processed$ClinicalSignificance, "_",
    clinvar_processed$seqnames
)]
arrow::write_parquet(clinvar_processed, "./data/all_labels.parquet")

clinvar_labels <- copy(clinvar_processed)
clinvar_labels[, `:=`(str_before = NULL, str_after = NULL)]
arrow::write_parquet(clinvar_labels, "./data/all_SV_labels.parquet")

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")