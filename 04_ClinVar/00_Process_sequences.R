# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))
pbapply::pboptions(char = "=")

# source functions
setwd("/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/04_ClinVar")
hg38_chain <- import.chain("../../05_Cosmic/data/liftover/hg38-chm13v2.over.chain")
hg19_chain <- import.chain("../../05_Cosmic/data/liftover/hg19ToHg38.over.chain")

# get gene lengths
file_to_load <- "../../05_Cosmic/data/annotations/group_genes.csv"
group_names <- stringr::str_extract(
    string = file_to_load,
    pattern = "(?<=group_)[^.]+(?=\\.csv)"
)
genes <- fread(file_to_load, showProgress = FALSE)
genes <- genes[type == "Genes"]

get_ref_genome <- function(version){
    refseq <- switch(version,
        "hg38" = BSgenome.Hsapiens.UCSC.hg38,
        "t2t" = BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
    )

    if(!any(grepl(pattern = "^chr", x = seqnames(refseq)))){
        chr_names <- paste0("chr", seqnames(refseq))
        seqnames(refseq@seqinfo) <- seqnames(refseq) <- chr_names
    }
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
    human_genome_len <- sum(refseq.table$chromSize)

    return(list(refseq, refseq.table, human_genome_len))
}

# get thresholds and column names
thresholds <- fread("../data/models/python/lightgbm/best_LGBM_model_thresholds.csv")
col_names <- thresholds$FPR
col_names <- gsub("\\.", "_", paste0("Pred_", col_names))

# t2t genome version
ref_out <- get_ref_genome(version = "t2t")
refseq <- ref_out[[1]]
refseq.table <- ref_out[[2]]
human_genome_len <- ref_out[[3]]

# import predicted breakpoints of the full human genome
path_to_pred <- "../data/Whole_genome_pred"
df_pred_all <- pbapply::pblapply(1:length(col_names), function(x){
    df_pred <- arrow::read_parquet(
        paste0(path_to_pred, "/", col_names[x], ".parquet")
    )
    df_pred[, width := 1]
    df_pred <- plyranges::as_granges(df_pred) 
    df_pred <- unique(df_pred)

    return(df_pred)
})
names(df_pred_all) <- col_names

hex_codes <- c("#990000", "#fd952c", "#669e85", "#394d7e")
names(hex_codes) <- col_names

# calculate the average human genome fragility
avg_human_fragility <- lengths(df_pred_all) / human_genome_len * 100
avg_human_fragility <- as.data.frame(as.matrix(avg_human_fragility))
avg_human_fragility <- avg_human_fragility %>% 
    tibble::rownames_to_column() %>% 
    as_tibble() %>% 
    dplyr::rename_with(~c("Key", "Value")) %>% 
    dplyr::mutate(hex = hex_codes[Key])

# Download CNV data from the TCGA.
# Data Portal: Data Release 39.0 - December 04, 2023. Accessed Wednesday December 6, 2023. 
# download: https://portal.gdc.cancer.gov/repository
# Data category: copy number variation.
# Data type: Copy number segment. 
# Platform: affymetrix snp 6.0.
# Workflow type: DNAcopy.
# Access: Open

#' Download ClinVar dataset.
#' Download: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz.
#' Accessed Wednesday December 6, 2023.
clinvar <- fread('./variant_summary.txt', showProgress = TRUE)

# only keep insertion / deletion / indels
clinvar <- clinvar[grepl(paste0(
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
    Assembly, Chromosome, Start, Stop
)]
setnames(clinvar, c(
    "Type", "GeneSymbol", "ClinicalSignificance", 
    "Origin", "Assembly", "seqnames", "start", "end"
))
clinvar[, seqnames := paste0("chr", seqnames)]

# keep GRCh37 and GRCh38 genome versions
clinvar <- clinvar[grepl("GRCh37|GRCh38", Assembly)]

# liftover from hg38 to T2T genome version.
clinvar_hg38 <- liftOver(
    as_tibble(clinvar) %>% 
        dplyr::filter(Assembly == "GRCh38") %>% 
        dplyr::filter(grepl("^chr([1-9]|1[0-9]|2[0-2])$", seqnames)) %>% 
        dplyr::filter(end >= start) %>% 
        plyranges::as_granges(.), 
    hg38_chain
)
clinvar_hg38 <- unlist(as(clinvar_hg38, "GRangesList"))

# liftover from hg37 to T2T genome version.
clinvar_hg37 <- liftOver(
    as_tibble(clinvar) %>% 
        dplyr::filter(Assembly == "GRCh37") %>% 
        dplyr::filter(grepl("^chr([1-9]|1[0-9]|2[0-2])$", seqnames)) %>% 
        dplyr::filter(end >= start) %>% 
        plyranges::as_granges(.), 
    hg19_chain
)
clinvar_hg37 <- unlist(as(clinvar_hg37, "GRangesList"))
clinvar_hg37 <- liftOver(clinvar_hg37, hg38_chain)
clinvar_hg37 <- unlist(as(clinvar_hg37, "GRangesList"))
clinvar_granges <- plyranges::bind_ranges(clinvar_hg38, clinvar_hg37)
clinvar_granges <- clinvar_granges %>% dplyr::select(-Assembly)

# separate out start and end positions
clinvar_granges <- dplyr::bind_rows(
        as_tibble(clinvar_granges) %>% 
            dplyr::select(-end, -width, -strand),
        as_tibble(clinvar_granges) %>% 
            dplyr::select(-start, -width, -strand) %>%
            dplyr::rename(start = end)
    ) %>% 
    dplyr::mutate(width = 1) %>% 
    dplyr::filter(!grepl(" ", GeneSymbol)) %>%
    dplyr::filter(!grepl(";", GeneSymbol)) %>%
    dplyr::filter(GeneSymbol != "-") %>%
    dplyr::mutate(
        ClinicalSig_Simple = dplyr::case_when(
            grepl("pathogenic", ClinicalSignificance, ignore.case = TRUE) ~ "Pathogenic",
            grepl("benign", ClinicalSignificance, ignore.case = TRUE) ~ "Benign"
        )
    ) %>%
    plyranges::as_granges()

# only get unique breakpoint positions by genes
clinvar_granges <- as_tibble(clinvar_granges) %>%
    dplyr::group_by(GeneSymbol, ClinicalSig_Simple, seqnames, start) %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    # dplyr::mutate(
    #     GeneLength = genes$width[match(GeneSymbol, genes$gene_id)],
    #     GeneSeqnames = genes$seqnames[match(GeneSymbol, genes$gene_id)],
    #     GeneStart = genes$start[match(GeneSymbol, genes$gene_id)],
    #     GeneEnd = genes$end[match(GeneSymbol, genes$gene_id)],
    # ) %>%
    # dplyr::filter(!is.na(GeneLength)) %>%
    # dplyr::mutate(
    #     rmGene = dplyr::case_when(
    #         seqnames == GeneSeqnames & start >= GeneStart & start <= GeneEnd ~ TRUE,
    #         .default = FALSE
    #     )
    # ) %>%
    # dplyr::filter(rmGene == TRUE) %>%
    # dplyr::select(-rmGene) %>%
    plyranges::as_granges()

# number of benign and pathogenic breaks by genes
min_const <- 1e-10
clinvar_true <- as_tibble(clinvar_granges) %>%
    dplyr::select(-end, -strand) %>%
    dplyr::group_by(GeneSymbol, ClinicalSig_Simple) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    # dplyr::mutate(GeneLength = genes$width[match(GeneSymbol, genes$gene_id)]) %>%
    # dplyr::mutate(count = (count + min_const) / (GeneLength + min_const)) %>%
    tidyr::pivot_wider(
        names_from = "ClinicalSig_Simple",
        values_from = "count"
    ) %>%
    dplyr::mutate(Key = "True") %>%
    dplyr::select(GeneSymbol, Key, Benign, Pathogenic) %>%
    dplyr::mutate(
        Benign = ifelse(is.na(Benign), 0, Benign),
        Pathogenic = ifelse(is.na(Pathogenic), 0, Pathogenic)
    ) %>%
    dplyr::mutate(PB_ratio = (Pathogenic + min_const) / sum(Pathogenic, Benign, min_const)) %>%
    dplyr::mutate(PB_ratio = PB_ratio / max(PB_ratio)) %>%
    # dplyr::mutate(Pathogenic_norm = Pathogenic / max(Pathogenic)) %>%
    # dplyr::mutate(Pathogenic = Pathogenic / max(Pathogenic)) %>%
    dplyr::arrange(desc(PB_ratio)) %>%
    suppressMessages()

# get short range of influence
range_cutoff <- fread("../data/Ranges_cutoffs.csv", select = c("Cluster", "kmer_8"))
range_cutoff <- range_cutoff[Cluster == "short.range", "kmer_8"][[1]]
round_to_nearest_even <- function(x) round(x/2)*2
range_cutoff <- round_to_nearest_even(range_cutoff)

unique_clinical_sig <- unique(clinvar_granges$ClinicalSig_Simple)
all_perc_overlaps <- lapply(1:length(unique_clinical_sig), function(x){
    label <- unique_clinical_sig[x]
    df_bps <- clinvar_granges %>% 
        dplyr::filter(ClinicalSignificance == unique_clinical_sig[x])
    # df_bps <- plyranges::stretch(plyranges::anchor_center(df_bps), range_cutoff)

    perc_overlaps <- pbapply::pblapply(1:length(df_pred_all), function(i){
        #' instead of predicting whether the exact position is 
        #' predicted to be broken, evaluate whether they overlap 
        #' within the medium range of influence
        overlaps <- GenomicRanges::findOverlaps(df_bps, df_pred_all[[i]])

        # count by gene
        gene_counts <- mcols(df_bps)$GeneSymbol[queryHits(overlaps)]
        gene_counts <- as.data.table(table(gene_counts))
        setnames(gene_counts, c("GeneSymbol", "Count"))
        gene_counts[, Count := (Count + min_const)]

        # add genes missing
        missing_genes <- data.table(
            GeneSymbol = setdiff(unique(clinvar_true$GeneSymbol), gene_counts$GeneSymbol),
            Count = 0
        )
        gene_counts <- rbind(gene_counts, missing_genes)

        # normalise each gene by its length
        # gene_counts[, GeneLength := genes$width[match(gene_counts$GeneSymbol, genes$gene_id)]]
        # gene_counts[, Count := (Count + min_const) / (GeneLength + min_const)]
        # gene_counts[, GeneLength := NULL]
        setnames(gene_counts, c("GeneSymbol", col_names[i]))

        return(gene_counts)
    })
    perc_overlaps <- do.call(cbind, perc_overlaps)
    perc_overlaps <- cbind(perc_overlaps[, 1], perc_overlaps[, ..col_names])
    perc_overlaps[, ClinicalSig_Simple := label]
})

gene_symbols <- all_perc_overlaps[[1]][["GeneSymbol"]]
dt_all_perc_overlaps <- all_perc_overlaps[[1]][, ..col_names] / all_perc_overlaps[[2]][, ..col_names]
dt_all_perc_overlaps[, `:=`(GeneSymbol = gene_symbols, ClinicalSig_Simple = "Pathogenic_Benign_ratio")]
setcolorder(dt_all_perc_overlaps, c("GeneSymbol", col_names, "ClinicalSig_Simple"))

# dt_all_perc_overlaps <- rbindlist(all_perc_overlaps)
# setorder(dt_all_perc_overlaps, ClinicalSig_Simple)

dt_all_perc_overlaps <- as_tibble(dt_all_perc_overlaps) %>%
    tidyr::pivot_longer(
        cols = -c(GeneSymbol, ClinicalSig_Simple),
        names_to = "Key",
        values_to = "count"
    ) %>%
    dplyr::mutate(count = dplyr::case_when(
        is.na(count) ~ 0,
        is.infinite(count) ~ 0,
        .default = count
    )) %>%
    tidyr::pivot_wider(
        names_from = "ClinicalSig_Simple",
        values_from = "count"
    ) %>%
    dplyr::group_by(Key) %>%
    # dplyr::mutate(Pathogenic_norm = Pathogenic / max(Pathogenic)) %>%
    # dplyr::arrange(desc(Pathogenic_norm)) %>%
    dplyr::mutate(PB_ratio = Pathogenic_Benign_ratio / max(Pathogenic_Benign_ratio)) %>%
    dplyr::arrange(desc(PB_ratio)) %>%
    dplyr::ungroup() %>%
    suppressMessages()

# combine with true data set
hex_codes <- c("#990000", "#fd952c", "#669e85", "#394d7e")
names(hex_codes) <- col_names
hex_codes <- c(True = "#000000", hex_codes)

dt_all <- dplyr::bind_rows(
    clinvar_true %>%
        dplyr::select(GeneSymbol, Key, PB_ratio), 
    dt_all_perc_overlaps %>%
        dplyr::select(GeneSymbol, Key, PB_ratio)
)
dt_all <- dt_all %>% 
    dplyr::mutate(
        # GeneLength = genes$width[match(GeneSymbol, genes$gene_id)],
        hex = hex_codes[Key],
        Key = factor(Key, levels = names(hex_codes))
    )

# get the top 20 per group
dt_all_split <- dplyr::group_split(dt_all, Key)
plot_all_groups <- lapply(1:length(dt_all_split), function(x){
    temp <- dt_all_split[[x]] %>%
        dplyr::slice_max(n = 20, order_by = PB_ratio) %>%
        dplyr::arrange(PB_ratio) %>%
        dplyr::mutate(GeneSymbol = forcats::fct_inorder(GeneSymbol))

    plot_title <- stringr::str_replace(
        string = names(hex_codes)[x], 
        pattern = "^(Pred)_0_(.*)$", 
        replacement = "\\1: 0.\\2 (FPR)"
    )

    p1 <- temp %>% 
        ggplot2::ggplot(aes(x = PB_ratio, y = GeneSymbol)) +
        ggplot2::geom_segment(aes(
            x = 0, xend = PB_ratio, 
            y = GeneSymbol, yend = GeneSymbol),
            linewidth = 1.1,
            col = 'grey80'
        ) +
        ggplot2::geom_point(size = 2, col = 'black') +
        ggplot2::theme_bw() + 
        ggplot2::theme_classic() +
        suppressWarnings(ggplot2::theme(
            text = element_text(size = 17),
            axis.text.y = element_text(colour = hex_codes[[x]])
        )) +
        ggplot2::coord_cartesian(xlim = c(0, NA)) + 
        ggplot2::labs(
            subtitle = plot_title,
            x = "", y = ""
        )

    return(p1)
})

dir.create(path = "../figures/ClinVar/", showWarnings = FALSE)
pdf(
    file = paste0(
        "../figures/ClinVar/",
        "Top_genes_by_pathogen_benign_ratio.pdf"
    ),
    height = 7, width = 20
)
do.call(gridExtra::grid.arrange, c(plot_all_groups, nrow = 1))
plots.saved <- dev.off()

# # plot as side-by-side barplot
# dt_all_perc_overlaps <- as_tibble(dt_all_perc_overlaps) %>% 
#     tidyr::pivot_longer(
#         cols = -clinical_sig,
#         names_to = "Key",
#         values_to = "Value"
#     ) %>% 
#     dplyr::mutate(
#         hex = hex_codes[Key],
#         Key = factor(Key, levels = col_names),
#         clinical_sig = as.factor(clinical_sig)
#     ) 

# p1 <- dt_all_perc_overlaps %>% 
#     ggplot2::ggplot(aes(x = clinical_sig, y = Value, fill = Key)) + 
#     ggplot2::geom_bar(stat = "identity", position = position_dodge(), col = "black") +
#     ggplot2::theme_bw() + 
#     ggplot2::theme_classic() +
#     ggplot2::scale_fill_manual(values = hex_codes) + 
#     ggplot2::theme(text = element_text(size = 17)) +
#     ggplot2::coord_cartesian(ylim = c(0, 100)) + 
#     ggplot2::labs(
#         subtitle = paste0(
#             "Accuracy of ClinVar predictions within the medium range of influence: ",
#             range_cutoff, " bases."
#         ),
#         x = "", y = "Accuracy, %"
#     )

# dir.create(path = "../figures/ClinVar/", showWarnings = FALSE)
# ggsave(
#     filename = "../figures/ClinVar/Clinical_Significance_Mid_range_pred.pdf",
#     plot = p1,
#     height = 7, width = 13
# )

# unique_type <- unique(clinvar$type)
# all_perc_overlaps <- lapply(1:length(unique_type), function(x){
#     label <- unique_type[x]
#     df_bps <- clinvar_liftover %>% 
#         dplyr::filter(type == unique_type[x])
#     df_bps <- plyranges::stretch(plyranges::anchor_center(df_bps), range_cutoff)

#     perc_overlaps <- lapply(df_pred_all, function(pred){
#         #' instead of predicting whether the exact position is 
#         #' predicted to be broken, evaluate whether they overlap 
#         #' within the medium range of influence
#         overlaps <- GenomicRanges::findOverlaps(df_bps, pred)
#         perc <- length(unique(queryHits(overlaps))) / length(df_bps) * 100
#         return(perc)
#     })
#     perc_overlaps <- as.data.table(perc_overlaps)
#     perc_overlaps[, type := label]
# })
# dt_all_perc_overlaps <- rbindlist(all_perc_overlaps)

# # plot as side-by-side barplot
# dt_all_perc_overlaps <- as_tibble(dt_all_perc_overlaps) %>% 
#     tidyr::pivot_longer(
#         cols = -type,
#         names_to = "Key",
#         values_to = "Value"
#     ) %>% 
#     dplyr::mutate(
#         hex = hex_codes[Key],
#         Key = factor(Key, levels = col_names),
#         type = as.factor(type)
#     ) 

# p1 <- dt_all_perc_overlaps %>% 
#     ggplot2::ggplot(aes(x = type, y = Value, fill = Key)) + 
#     ggplot2::geom_bar(stat = "identity", position = position_dodge(), col = "black") +
#     ggplot2::theme_bw() + 
#     ggplot2::theme_classic() +
#     ggplot2::scale_fill_manual(values = hex_codes) + 
#     ggplot2::theme(text = element_text(size = 17)) +
#     ggplot2::coord_cartesian(ylim = c(0, 100)) + 
#     ggplot2::labs(
#         subtitle = paste0(
#             "Accuracy of ClinVar predictions within the medium range of influence: ",
#             range_cutoff, " bases."
#         ),
#         x = "", y = "Accuracy, %"
#     )

# ggsave(
#     filename = "../figures/ClinVar/Break_Type_Mid_range_pred.pdf",
#     plot = p1,
#     height = 7, width = 13
# )