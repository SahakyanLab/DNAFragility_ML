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

cat(paste0("Processing overlaps with low/mid/high fragile regions...\n"))

args <- commandArgs(trailingOnly = TRUE)
bw <- as.integer(args[1])

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

# get thresholds and column names
thresholds <- fread("../data/models/python/lightgbm/best_LGBM_model_thresholds.csv")
col_names <- thresholds$FPR
col_names <- gsub("\\.", "_", paste0("Pred_", col_names))

# import predictions
path_to_pred <- "../data/Whole_genome_pred"

get_genic_feat <- function(){
    # import all genomic features for analysis
    files_to_load <- list.files(
        path = "../../05_Cosmic/data/annotations",
        pattern = "group_*",
        full.names = TRUE
    )
    group_names <- stringr::str_extract(
        string = files_to_load,
        pattern = "(?<=group_)[^.]+(?=\\.csv)"
    )

    all_groups <- lapply(files_to_load, fread, showProgress = FALSE)
    names(all_groups) <- group_names

    return(list(all_groups, group_names))
}

# global vars
df_genic <- get_genic_feat()
group_names <- df_genic[[2]]
df_genic <- df_genic[[1]]
groups_combined <- rbindlist(df_genic, fill = TRUE)
groups_combined <- groups_combined[, gene_id := NULL]
groups_split <- split(groups_combined, by = "type")
group_map <- lapply(1:length(df_genic), function(ind){
    return(data.table(
        type = unique(df_genic[[ind]]$type),
        group = names(df_genic)[[ind]]
    ))
})
group_map <- rbindlist(group_map)

hex <- c("#990000", "#fd952c", "#669e85", "#394d7e")
hex <- setNames(hex, group_names)
group_map[, hex := hex[group]]

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

# count cosmic breaks in non-overlapping windows and return granges 
df_ppb_all <- pbapply::pblapply(1:length(col_names), function(x){
    df_pred <- arrow::read_parquet(paste0(path_to_pred, "/", col_names[x], ".parquet"))

    df_bp_granges <- lapply(1:22, function(chr){
        pred_bins <- calculate_overlaps(
            df_breaks = df_pred, 
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

    #' If we bin the genome, and get the low/mid/high fragile zones, 
    #' then what is the fraction of the whole genome with these regions?
    # as_tibble(df_bp_granges_all) %>% 
    #     dplyr::group_by(group) %>%
    #     dplyr::summarise(count = dplyr::n() * bw) %>%
    #     dplyr::ungroup() %>% 
    #     dplyr::mutate(count = count / sum(refseq.table$end) * 100)

    # map classification back onto table of individual breakpoints
    mcols_all <- as.data.table(mcols(df_bp_granges_all))
    map_bins <- match(df_bp_granges_keep_all$bins, mcols_all$bins)
    df_bp_granges_keep_all[, group := mcols_all$group[map_bins]]
    df_bp_granges_keep_all[, bins := NULL]

    # get the overlap of breaks, per genic element, per L/M/H group.
    df_bp_granges_all_split <- split(df_bp_granges_keep_all, df_bp_granges_keep_all$group)
    ppb_by_group <- lapply(1:length(df_bp_granges_all_split), function(group){
        df_bp_granges_group <- df_bp_granges_all_split[[group]]
        df_bp_granges_group[, width := 1]
        df_bp_granges_group <- plyranges::as_granges(df_bp_granges_group)

        # perform same calculations for all genomic features
        ppb <- lapply(1:length(df_genic), function(g){
            # split group into chunks
            df_annot_split <- split(x = df_genic[[g]], by = "type")

            ppb_for_feat <- lapply(1:length(df_annot_split), function(ind){
                this_type <- as.character(df_annot_split[[ind]]$type[1])

                temp_annot <- as_tibble(df_annot_split[[ind]]) %>%
                    dplyr::arrange(seqnames, start) %>% 
                    plyranges::as_granges() %>%
                    reduce(.) 
                    
                temp_annot <- temp_annot %>% 
                    dplyr::mutate(type = this_type)

                # calculate overlaps count
                feature_widths <- sum(width(temp_annot))
                # calculate probability per base
                overlaps <- GenomicRanges::findOverlaps(temp_annot, df_bp_granges_group)

                if(length(overlaps) == 0){
                    overlap_counts <- 0 
                } else {
                    overlap_counts <- sum(table(queryHits(overlaps)))
                }
                ppb <- overlap_counts / feature_widths
                
                return(ppb)
            })
            names(ppb_for_feat) <- names(df_annot_split)
            return(ppb_for_feat)
        })
        ppb <- unlist(ppb)
        ppb <- ppb[order(names(ppb))]
        return(ppb)
    })
    ppb_by_group_all <- do.call(rbind, ppb_by_group)
    ppb_by_group_all <- apply(ppb_by_group_all, 2, function(x) x / sum(x))
    rownames(ppb_by_group_all) <- names(df_bp_granges_all_split)

    df_ppb <- as.data.frame(t(ppb_by_group_all)) %>% 
        tibble::rownames_to_column() %>%
        as_tibble() %>%
        dplyr::rename(feature = rowname) %>%
        dplyr::arrange(High)

    df_ppb <- df_ppb %>%
        tidyr::pivot_longer(
            cols = -feature,
            names_to = "Key",
            values_to = "Value"
        ) %>% 
        dplyr::mutate(
            feature = factor(feature, levels = df_ppb$feature),
            type = group_map$group[match(feature, group_map$type)],
            feat_hex = group_map$hex[match(feature, group_map$type)],
            Key_hex = dplyr::case_when(
                Key == "High" ~ "#b2182b",
                Key == "Mid" ~ "#006400",
                Key == "Low" ~ "#ff9e27"
            ),
            pred = col_names[x]
        )

    return(df_ppb)
})
df_ppb_rbind <- do.call(dplyr::bind_rows, df_ppb_all)

# plot stacked barplot, ranked by decreasing order of H group.
ppb_plots <- lapply(1:length(df_ppb_all), function(x){
    df_ppb_hex <- df_ppb_all[[x]] %>% 
        dplyr::group_by(feature) %>%
        dplyr::slice(1)

    plot_title <- stringr::str_replace(
        string = col_names[[x]], 
        pattern = "^(Pred)_0_(.*)$", 
        replacement = "\\1: 0.\\2 (FPR)"
    )

    p1 <- df_ppb_all[[x]] %>%
        dplyr::mutate(Key = factor(Key, levels = c("High", "Mid", "Low"))) %>%
        ggplot(aes(x = Value, y = feature, fill = Key)) + 
        geom_bar(position = "stack", stat = "identity") +
        theme_bw() + 
        theme_classic() + 
        scale_fill_manual(values = df_ppb_all[[x]]$Key_hex) + 
        suppressWarnings(theme(
            text = element_text(size = 15),
            legend.position = "none",
            axis.text.y = element_text(colour = df_ppb_hex$feat_hex)
        )) +
        labs(subtitle = plot_title, x = "Fraction", y = "")

    return(p1)
})

pdf(
    file = paste0(
        "../figures/Whole_genome_pred/",
        "feature_breakability_by_category_",
        "bw_", bw, ".pdf"
    ),
    height = 10, width = 24
)
do.call(gridExtra::grid.arrange, c(ppb_plots, nrow = 1))
plot.saved <- dev.off()