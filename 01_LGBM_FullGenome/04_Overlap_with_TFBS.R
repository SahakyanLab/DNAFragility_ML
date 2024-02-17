# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))
suppressPackageStartupMessages(suppressWarnings(library(arrow)))
suppressPackageStartupMessages(suppressWarnings(library(rtracklayer)))
pbapply::pboptions(char = "=", type = "txt")

# source functions
my.path = "/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/01_LGBM_FullGenome"
setwd(my.path)

cat(paste0("Processing overlaps with TFBS...\n"))

########################################################################
# load in all the motifs and positions from JASPAR
path_to_dat <- "../data/TFBS"

load_motif_mapping <- function(path){
    mapping <- readr::read_file(path) %>%
        str_split(., ">") %>%
        sapply(., function(x) {
            str_extract_all(x, "[.:A-Za-z0-9]*\\t[.:A-Za-z0-9]*(?=\\n)")
        }) %>%
        unlist() %>%
        sapply(., function(x) {
            str_split(x, "\\t")
        })
    names(mapping) <- sapply(mapping, function(x) {
        x[1]
    })
    mapping <- sapply(mapping, function(x) {
        x[2]
    })
    mapping
}

# Load mapping between matrix ID and motif name
jaspar_motifs <- load_motif_mapping(
    "../../03_Breakpoints_v2/data/02_JASPAR/JASPAR2024_CORE_vertebrates_redundant_pfms_jaspar.txt"
)
hg38_chain <- import.chain("../../05_Cosmic/data/liftover/hg38-chm13v2.over.chain")
hg19_chain <- import.chain("../../05_Cosmic/data/liftover/hg19ToHg38.over.chain")

# Load motifs bams into single object, and make method for anottations
bed_path <- "../data/TFBS/jaspar_beds"
bed_files <- list.files(bed_path, ".bed", full.names = TRUE)
all_bed_names <- paste0(bed_path, "/", names(jaspar_motifs), ".bed")
valid_files <- intersect(all_bed_names, bed_files)
bed_ranges <- pbapply::pbsapply(valid_files, plyranges::read_bed)

# only keep hg19 and hg38 versions
assemblies <- lapply(1:length(bed_ranges), function(x){
  jaspar_label <- stringr::str_extract(
      names(bed_ranges[x]), 
      "jaspar_beds/([A-Za-z0-9:.]*).bed", group = 1
  )
  tfbs_label <- as.character(jaspar_motifs[jaspar_label])

  # if any punctuations
  tfbs_label <- gsub("[[:punct:]]", "_", tfbs_label)
  tfbs_label <- gsub("___+", "", tfbs_label)
  tfbs_label <- gsub("__+", "_", tfbs_label)

  # capitalise
  tfbs_label <- toupper(tfbs_label)

  data.table(
    ind = x,
    assembly = stringr::str_extract(mcols(bed_ranges[[x]])$name[[1]], ".+?(?=_)"),
    jaspar_name = jaspar_label,
    tfbs_label = tfbs_label
  )
})
assemblies <- as.data.table(do.call(rbind, assemblies))

# if multiple TFBS exist, keep the first occurring one.
assemblies <- as_tibble(assemblies) %>% 
  dplyr::filter(grepl("hg19|hg38", assembly, ignore.case = TRUE)) %>% 
  dplyr::arrange(tfbs_label, jaspar_name) %>% 
  dplyr::group_by(tfbs_label) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()

# liftover from hg37 to T2T genome version.
liftover_granges <- pbapply::pblapply(1:nrow(assemblies), function(x){
    genome_assembly <- stringr::str_extract(mcols(bed_ranges[[assemblies$ind[x]]])$name[[1]], ".+?(?=_)")
    bed_grange <- bed_ranges[[assemblies$ind[x]]]
    bed_grange_label <- assemblies$tfbs_label[x]
    
    # liftover to t2t genome
    if(grepl("hg19", genome_assembly)){
        liftover_hg19 <- suppressMessages(liftOver(bed_grange, hg19_chain))
        bed_grange <- unlist(as(liftover_hg19, "GRangesList"))
    } 
    liftover_hg38 <- suppressMessages(liftOver(bed_grange, hg38_chain))
    liftover_hg38 <- unlist(as(liftover_hg38, "GRangesList"))
    liftover_hg38 <- liftover_hg38 %>% 
        dplyr::select(-name, -score) %>% 
        plyranges::filter(seqnames %in% paste0("chr", 1:22)) %>%
        reduce(.) %>%
        dplyr::mutate(label = bed_grange_label)

    return(liftover_hg38)
})
liftover_granges_all <- do.call(plyranges::bind_ranges, liftover_granges)
saveRDS(liftover_granges_all, "../data/TFBS/t2t_jaspar_beds.RData")

TFBS_labels <- unique(mcols(liftover_granges_all)$label)
TFBS_sums <- pbapply::pblapply(TFBS_labels, function(l){
    TFBS_sums <- liftover_granges_all %>% 
        dplyr::filter(label == l) %>%
        width() %>%
        sum()
    return(TFBS_sums)
})
TFBS_table <- data.table(
    TF = TFBS_labels,
    len = unlist(TFBS_sums)
)
setorder(TFBS_table, -len)

refseq <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
if(!any(grepl(pattern = "^chr", x = seqnames(refseq)))){
    chr_names <- paste0("chr", seqnames(refseq))
    seqnames(refseq@seqinfo) <- seqnames(refseq) <- chr_names
}

# get the start and end positions of each chromosome
refseq.table <- as.data.frame(refseq@seqinfo)
refseq.table <- refseq.table[grepl(
    pattern = "^chr([1-9]|1[0-9]|2[0-2])$", 
    x = rownames(refseq.table)
),]
refseq.table <- data.table(
    chr = rownames(refseq.table),
    end = refseq.table$seqlengths
)

# Import true cosmic breakpoint data
get_cosmic_breaks <- function(kmer){
    df.parsed <- fread(
        file = paste0(
            "../../05_Cosmic/data/BaseTable/",
            "kmer_", kmer, ".csv"
        ),
        select = "ID",
        showProgress = FALSE
    )

    #' Remove the sample ID
    df <- as_tibble(df.parsed) %>%
        tidyr::separate(
            col = ID,
            into = c(
                "Chr", "Start", 
                "Tissue", "Cancer", 
                "BreakTFBS", 
                "SampleID"
            ),
            sep = "_",
            remove = TRUE
        ) %>% 
        tidyr::unite(
            TC_ID,
            c("Tissue", "Cancer"), 
            remove = FALSE
        ) %>% 
        tidyr::unite(
            Chr_Start_ID, 
            c("Chr", "Start"),
            remove = FALSE
        ) %>% 
        dplyr::select(-c(SampleID, Tissue, Cancer)) %>% 
        dplyr::mutate(Start = as.numeric(Start)) %>% 
        dplyr::select(Chr, Start) %>% 
        dplyr::rename_with(~c("seqnames", "start"))

    # liftover breakpoints to the telomere-to-telomere genome version
    chain <- import.chain("../../05_COSMIC/data/liftover/hg38-chm13v2.over.chain")
    df <- df %>% dplyr::mutate(width = 1, strand = "+")
    df <- plyranges::as_granges(df)
    df <- liftOver(df, chain)
    df <- unlist(as(df, "GRangesList"))
    df <- df %>% 
        dplyr::arrange(seqnames, start) %>% 
        as_tibble() %>% 
        dplyr::select(seqnames, start)
    
    return(df)
}

# global vars
kmer <- 8
df_true <- get_cosmic_breaks(kmer = kmer)

# get thresholds and column names
thresholds <- fread("../data/models/python/lightgbm/best_LGBM_model_thresholds.csv")
col_names <- thresholds$FPR
col_names <- gsub("\\.", "_", paste0("Pred_", col_names))

# import predictions
df_pred_all <- lapply(col_names, function(x){
    return(arrow::read_parquet(paste0("../data/Whole_genome_pred/", x, ".parquet")))
})
names(df_pred_all) <- col_names

# prepare GRanges objects
true_breaks_granges <- df_true %>% 
    dplyr::mutate(width = 1) %>% 
    plyranges::as_granges() %>% 
    unique()
pred_breaks_granges <- pbapply::pblapply(df_pred_all, function(df){
    return(
        df %>% 
            dplyr::mutate(width = 1) %>% 
            plyranges::as_granges() %>% 
            unique()
    )
})

get_norm_overlap_counts <- function(breaks_table, feature_table, overlap_vector){
    # Overlaps per genomic feature
    overlap_counts <- table(queryHits(overlap_vector))
    overlap_counts <- as.data.table(overlap_counts)
    setnames(overlap_counts, c("ID", "Count"))
    overlap_counts[, ID := as.integer(ID)]

    # get labels and total width of each TFBS
    TFBS_names <- TFBS_table$TF[match(feature_table[overlap_counts$ID]$label, TFBS_table$TF)]
    overlap_counts[, TFBS := TFBS_names]

    # total overlap grouped by TFBS
    overlap_counts <- overlap_counts[, .(Count = sum(Count)), by = TFBS]

    # normalise by length of TFBS
    TFBS_len <- TFBS_table$len[match(overlap_counts$TFBS, TFBS_table$TF)]
    overlap_counts[, TFBS_len := TFBS_len]
    overlap_counts[, norm_counts := Count / TFBS_len]

    # TFBS not present will be zero
    TFBS_names <- TFBS_table$TF
    missing_TFBS <- which(is.na(match(TFBS_names, overlap_counts$TFBS)))
    missing_TFBS <- data.table(
        TFBS_len = width(feature_table[missing_TFBS]),
        TFBS = TFBS_names[missing_TFBS]
    )
    missing_TFBS[, `:=`(Count = 0, norm_counts = 0)]
    setcolorder(missing_TFBS, c("TFBS", "Count", "TFBS_len", "norm_counts"))

    # combine table
    res <- rbind(overlap_counts, missing_TFBS)
    setorder(res, TFBS)

    return(res)
}

# count cosmic breaks in non-overlapping windows and return granges 
true_overlaps <- GenomicRanges::findOverlaps(
    liftover_granges_all, true_breaks_granges
)
res_true <- get_norm_overlap_counts(
    breaks_table = true_breaks_granges, 
    feature_table = liftover_granges_all, 
    overlap_vector = true_overlaps
)

# predicted breaks
res_pred <- pbapply::pblapply(1:length(col_names), function(x){
    pred_overlaps <- GenomicRanges::findOverlaps(
        liftover_granges_all, pred_breaks_granges[[x]]
    )
    res_pred <- get_norm_overlap_counts(
        breaks_table = pred_breaks_granges[[x]],
        feature_table = liftover_granges_all, 
        overlap_vector = pred_overlaps
    )
    return(res_pred)
})
names(res_pred) <- col_names
ppb <- c(list("true" = res_true), res_pred)

# only keep norm_counts
ppb <- lapply(1:length(ppb), function(x){
    return(ppb[[x]]$norm_counts)
})
ppb <- do.call(cbind, ppb)
colnames(ppb) <- c("True", col_names)
ppb_df <- as.data.frame(ppb)
rownames(ppb_df) <- TFBS_table$TF

# Normalise between 0 and 1 to generate a relative importance plot
ppb_df <- apply(ppb_df, 2, function(x) x / max(x))

hex <- c("#990000", "#fd952c", "#669e85","#394d7e", "#8A2BE2")
hex <- setNames(hex, c("True", col_names))

df_ppb <- as_tibble(ppb_df) %>% 
    dplyr::mutate(TFBS = TFBS_table$TF) %>% 
    dplyr::arrange(desc(ppb_df))

dir.create(path = "../figures/TFBS/", showWarnings = FALSE)
plots_ppb <- lapply(c("True", col_names), function(column_name){
    temp <- df_ppb[c(column_name, "TFBS")] %>% 
        dplyr::rename_with(~c("ppb", "TFBS")) %>% 
        dplyr::slice_max(prop = 0.05, order_by = ppb) %>%
        dplyr::arrange(ppb) %>% 
        dplyr::mutate(
            TFBS = forcats::fct_inorder(TFBS),
            hex = hex[[column_name]]
        )

    plot_title <- ifelse(
        column_name == "True", "True",
        stringr::str_replace(
            string = column_name, 
            pattern = "^(Pred)_0_(.*)$", 
            replacement = "\\1: 0.\\2 (FPR)"
        )
    )

    p1 <- temp %>% 
        ggplot2::ggplot(aes(x = ppb, y = TFBS)) +
        ggplot2::geom_segment(aes(
            x = 0, xend = ppb, 
            y = TFBS, yend = TFBS),
            linewidth = 1.1,
            col = 'grey80'
        ) +
        ggplot2::geom_point(size = 2, col = 'black') +
        ggplot2::scale_color_identity() + 
        ggplot2::theme_bw() + 
        ggplot2::theme_classic() +
        suppressWarnings(ggplot2::theme(
            text = element_text(size = 15),
            axis.text.y = element_text(colour = temp$hex)
        )) +
        ggplot2::coord_cartesian(xlim = c(0, NA)) + 
        ggplot2::labs(
            subtitle = plot_title,
            x = "", y = ""
        )

    # ggsave(
    #     filename = paste0(
    #         "../figures/TFBS/",
    #         "PPB_", column_name, ".pdf"
    #     ),
    #     plot = p1,
    #     height = 11, width = 8
    # )

    return(p1)
})

pdf(
    file = "../figures/TFBS/PPB_All.pdf",
    height = 6, width = 18
)
do.call(gridExtra::grid.arrange, c(plots_ppb, nrow = 1))
plots.saved <- dev.off()