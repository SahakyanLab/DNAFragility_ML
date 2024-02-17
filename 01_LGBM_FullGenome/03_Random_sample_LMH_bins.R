# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))
pbapply::pboptions(char = "=", type = "txt")

# source functions
setwd("/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/01_LGBM_FullGenome")

cat(paste0("Processing random breakpoints for controls...\n"))

# liftover to T2T genome version.
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
), ]
refseq.table <- data.table(
    seqnames = rownames(refseq.table),
    end = refseq.table$seqlengths
)
human_genome_len <- sum(refseq.table$end)

# get thresholds and column names
thresholds <- fread("../data/models/python/lightgbm/best_LGBM_model_thresholds.csv")
col_names <- thresholds$FPR
col_names <- gsub("\\.", "_", paste0("Pred_", col_names))

# import predicted breakpoints of the full human genome
path_to_pred <- "../data/Whole_genome_pred"
df_pred_all <- pbapply::pblapply(col_names, function(x){
    return(arrow::read_parquet(paste0(path_to_pred, "/", x, ".parquet")))
})
names(df_pred_all) <- col_names

#' 1. Bin the genome and extract the low/medium/high fragile regions
get_breaks_in_bins <- function(break_table, bin_width, chr_len_table, filter_for_chr){
    # create bins
    temp_chr <- chr_len_table[chr == filter_for_chr]
    bins <- seq(0, temp_chr$end, by = bin_width)
    # final bin to include last bases
    if(tail(bins, n = 1) < temp_chr$end) bins <- c(bins, temp_chr$end)

    # Bin genome by non-overlapping bins
    df_bp <- break_table %>% 
        dplyr::filter(seqnames == filter_for_chr) %>% 
        dplyr::distinct()

    # count breaks per bin
    df_bp_granges <- df_bp %>% 
        dplyr::mutate(bins = cut(start, breaks = bins)) %>% 
        dplyr::group_by(bins) %>% 
        dplyr::summarise(Breaks = dplyr::n()) %>% 
        dplyr::arrange(bins) %>% 
        dplyr::mutate(seqnames = filter_for_chr) %>% 
        dplyr::filter(!is.na(bins)) %>% 
        dplyr::select(bins, Breaks) %>% 
        dplyr::mutate(
            seqnames = filter_for_chr,
            start = stringr::str_extract(
                string = bins, pattern = "(?<=\\()([^,]+)"
            ),
            start = as.numeric(start)+1,
            end = pmin(start+bin_width-1, temp_chr$end)
        ) %>% 
        dplyr::select(-bins) %>% 
        plyranges::as_granges()

    return(df_bp_granges)
}

calculate_overlaps <- function(df_breaks, chr, refseq_table, bw){
    # create bins
    temp_chr <- refseq_table[chr, ] 
    bins <- seq(0, temp_chr$end, by = bw)
    # final bin to include last bases
    if(tail(bins, n = 1) < temp_chr$end) bins <- c(bins, temp_chr$end)

    # Bin genome by non-overlapping bins
    df_bp_keep_all <- df_breaks %>% 
        dplyr::filter(seqnames == paste0("chr", chr)) %>% 
        dplyr::distinct() %>% 
        dplyr::mutate(bins = cut(start, bins)) 

    df_bp <- df_breaks %>% 
        dplyr::filter(seqnames == paste0("chr", chr)) %>% 
        dplyr::distinct() %>% 
        dplyr::mutate(width = 1) %>% 
        plyranges::as_granges()

    # bin genome into non-overlapping bins
    num_bins <- ceiling(temp_chr$end/bw)
    df_bp_granges <- GRanges(
        seqnames = paste0("chr", chr),
        ranges = IRanges(
            start = seq(1, by = bw, length.out = num_bins),
            end = c(seq(bw, by = bw, length.out = num_bins-1), temp_chr$end)
        )
    )
    df_bp_granges <- df_bp_granges %>% 
        dplyr::mutate(
            Breaks = countOverlaps(df_bp_granges, df_bp),
            bins = cut(start, bins)
        )
    
    return(list(df_bp_granges, df_bp_keep_all))
}

get_lmh_groups <- function(x, bw){
    df_bp_granges <- pbapply::pblapply(1:22, function(chr){
        pred_bins <- calculate_overlaps(
            df_breaks = df_pred_all[[x]],
            chr = chr, 
            refseq_table = refseq.table,
            bw = bw
        )
        return(pred_bins)
    })
    df_bp_granges_all <- suppressWarnings(unlist(as(sapply(df_bp_granges, `[`, 1), "GRangesList")))
    df_bp_granges_keep_all <- do.call(dplyr::bind_rows, sapply(df_bp_granges, `[`, 2))

    df_breaks <- tibble(
            seqnames = as.character(seqnames(df_bp_granges_all)),
            Breaks = mcols(df_bp_granges_all)$Breaks
        ) %>% 
        dplyr::mutate(ID = 1:nrow(.))

    # top 5%
    top_5 <- df_breaks %>% 
        dplyr::slice_max(prop = 0.05, order_by = Breaks) %>% 
        dplyr::pull(ID)
    # bottom 5%
    bottom_5 <- df_breaks %>% 
        dplyr::slice_min(prop = 0.05, order_by = Breaks) %>% 
        dplyr::pull(ID)
    # rest
    df_breaks <- df_breaks %>% 
        dplyr::mutate(
            group = dplyr::case_when(
                ID %in% top_5 ~ "High",
                ID %in% bottom_5 ~ "Low",
                .default = "Mid"
            ),
            group = as.factor(group)
        ) %>% 
        dplyr::select(-ID) 

    df_bp_granges_all <- df_bp_granges_all %>% 
        dplyr::mutate(group = df_breaks$group)

    #################################################################################
    #' What percent of the genome is a low, mid, and high fragile zone?
    # df_breaks %>% 
    #     dplyr::group_by(group) %>% 
    #     dplyr::summarise(count = dplyr::n()) %>% 
    #     dplyr::ungroup() %>% 
    #     dplyr::mutate(frac = (count / human_genome_len) * bw * 100)
    #################################################################################

    #' What percent of random sampled breakpoints are found in low, mid, and high fragile zones?
    overlaps <- GenomicRanges::findOverlaps(df_random_granges, df_bp_granges_all)
    overlaps <- df_random_granges[queryHits(overlaps)] %>%
        dplyr::mutate(group = mcols(df_bp_granges_all)$group[subjectHits(overlaps)])
        
    # group_counts <- group_fracs %>% 
    #     dplyr::select(group, count) %>% 
    #     tidyr::pivot_wider(
    #         names_from = group, 
    #         values_from = count
    #     )

    # group_frac <- group_fracs %>% 
    #     dplyr::select(group, frac) %>% 
    #     tidyr::pivot_wider(
    #         names_from = group, 
    #         values_from = frac
    #     )

    # round to 2 significant figures and keep sum to 100%. 
    round_percent <- function(values){ 
        scale_factor <- 100 / sum(values)
        res <- round(values * scale_factor, digits = 2)
        diff <- 100 - sum(res)
        res[which.max(res)] <- res[which.max(res)] + diff
        return(res)
    }

    group_fracs <- as_tibble(overlaps) %>% 
        dplyr::group_by(group) %>% 
        dplyr::summarise(count = dplyr::n()) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(
            frac = count / sum(count) * 100,
            frac = round_percent(frac),
            bin_width = bw,
            threshold = x
        )
    
    return(group_fracs)
}

#' control: randomly sample positions of the genome and calculate the overlap in each group
to_sample_max <- 10000000
set.seed(1234)

# sample sequences
df_random <- refseq.table[sample(
    x = .N, size = to_sample_max, replace = TRUE
)]
df_random[, `:=`(
    start = pbapply::pbsapply(end, function(l){
        sample(x = 1:l, size = 1)
    }),
    end = NULL
)]
setorder(df_random, seqnames, start)
df_random <- distinct(df_random)

dir.create(
    paste0(path_to_pred, "/random_sampling"),
    showWarnings = FALSE
)

fwrite(
    df_random,
    paste0(path_to_pred, "/random_sampling/random_pos.csv"),
    showProgress = FALSE
)

# import range of sequence influence to determine the max bin width
range_effets <- fread("../data/range_effects/MaxValuesFromClustersByType.csv")
bws <- c(ceiling(max(range_effets$long)), 10000, 20000)
bws <- as.integer(bws)

df_random_granges <- as_tibble(df_random) %>%
    dplyr::mutate(width = 1) %>% 
    plyranges::as_granges()

df_all_categories <- lapply(bws, function(bw){
    pred_categories <- lapply(col_names, function(x){
        return(get_lmh_groups(x = x, bw = bw))
    })
    df_all_cat <- do.call(dplyr::bind_rows, pred_categories)
    df_all_cat <- as.data.table(df_all_cat)
    return(df_all_cat)
})
df_all_cat <- rbindlist(df_all_categories)

fwrite(
    df_all_cat,
    paste0(path_to_pred, "/random_sampling/random_frac.csv"),
    showProgress = FALSE
)