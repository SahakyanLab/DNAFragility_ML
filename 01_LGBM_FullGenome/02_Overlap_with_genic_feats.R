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
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

cat(paste0("Processing overlaps with genomic features...\n"))

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

#' 1. Get all predicted breaks files
pred_files <- list.files(
    path = "../data/models/python/lightgbm/final_predictions",
    full.names = TRUE,
    recursive = TRUE
)

# move files from chromosome folders to final_predictions folder
if(length(pred_files) < 22){
    path_to_chrs <- paste0(
        "/media/hert6114/Paddy_6TB/",
        "ProjectBoard_Patrick/04_DNAFragility/",
        "data/models/python/lightgbm"
    )

    pred_files <- list.files(
        path = path_to_chrs,
        pattern = "final_pred_chr([0-9]|1[0-9]|2[0-2])",
        full.names = TRUE,
        recursive = TRUE
    )
    pred_files <- stringr::str_sort(pred_files, numeric = TRUE)

    pred_files <<- pred_files[!grepl(paste0(
        path_to_chrs, "/final_predictions/"
    ), pred_files)]

    file.mv <- pbapply::pblapply(1:length(pred_files), function(f){
        chr_num <- stringr::str_extract(
            string = basename(pred_files[f]),
            pattern = "(\\d)+"
        )

        pos_files <- list.files(
            path = paste0("../data/models/python/lightgbm/chr", chr_num),
            pattern = "FullGenome_Ranges_1",
            full.names = TRUE,
            recursive = TRUE
        )
        pos_files <- stringr::str_sort(pos_files, numeric = TRUE)
        
        pos_files <- sapply(pos_files, function(p){
            out <- arrow::read_parquet(p)
            return(out$start)
        })
        start_pos <- unlist(pos_files, use.names = FALSE)

        preds <- arrow::read_parquet(pred_files[f])
        preds <- preds %>% dplyr::mutate(start = start_pos, .before = 1)

        arrow::write_parquet(
            preds,
            paste0(
                "../data/models/python/lightgbm/final_predictions/", 
                "final_pred_chr_", chr_num, ".parquet"
            )
        )
    })

    pred_files <- list.files(
        path = "../data/models/python/lightgbm/final_predictions",
        full.names = TRUE,
        recursive = TRUE
    )
}
pred_files <- stringr::str_sort(pred_files, numeric = TRUE)

# get thresholds and column names
thresholds <- fread("../data/models/python/lightgbm/best_LGBM_model_thresholds.csv")
col_names <- thresholds$FPR
col_names <- gsub("\\.", "_", paste0("Pred_", col_names))

#' 2. Import true cosmic breakpoint data
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
                "BreakType", 
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
kmer <- 8
df_true <- get_cosmic_breaks(kmer = kmer)

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

group_map[, ontop := ifelse(type == "Genes", 1, 0)]
setorder(group_map, -ontop, group, type)
group_map[, ontop := NULL]
groups_split <- groups_split[match(group_map$type, names(groups_split))]

# bin each chromosome based on length of each bin size
df_bw <- tibble(
    bw = as.integer(c(500, 5000, 10000, 20000, 30000, 50000, 90000)),
    x_range = c(15, 15, 30, 40, 40, 60, 80)
)

cur_path <- "../figures/Whole_genome_pred"
dir.create(path = cur_path, showWarnings = FALSE, recursive = TRUE)

#' 3. Plot distribution of break bins for true and predicted ones 
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
        dplyr::mutate(Breaks = countOverlaps(df_bp_granges, df_bp))

    return(df_bp_granges)
}

# count cosmic breaks in non-overlapping windows and return granges 
bw <- 20000
df_bp_granges <- pbapply::pblapply(1:length(pred_files), function(chr){
    chr_name <- stringr::str_extract(
        string = basename(pred_files[[chr]]),
        pattern = "(\\d)+"
    ) %>% as.integer()

    true_bins <- calculate_overlaps(
        df_breaks = df_true, 
        chr = chr_name, 
        refseq_table = refseq.table,
        bw = bw
    )

    df_pred <- arrow::read_parquet(pred_files[[chr]])
    colnames(df_pred) <- c("start", "True", col_names)
    df_pred$seqnames <- as.factor(paste0("chr", chr))
    cols_of_int <- col_names

    pred_bins <- lapply(1:length(cols_of_int), function(column){
        temp <- df_pred[, c("seqnames", "start", cols_of_int[column])]
        colnames(temp) <- c("seqnames", "start", "pred")
        temp <- temp %>% 
            dplyr::filter(pred != 0) %>% 
            dplyr::select(seqnames, start)

        pred_bins <- calculate_overlaps(
            df_breaks = temp,
            chr = chr_name, 
            refseq_table = refseq.table,
            bw = bw
        )
        pred_bins <- as.data.table(mcols(pred_bins))
        setnames(pred_bins, cols_of_int[column])
        return(pred_bins)
    })
    pred_bins <- do.call(cbind, pred_bins)
    
    all_bins <- as_tibble(true_bins) %>% 
        dplyr::rename(True = Breaks) %>% 
        dplyr::bind_cols(pred_bins) %>% 
        plyranges::as_granges()

    return(all_bins)
})
df_bp_granges_all <- suppressWarnings(unlist(as(df_bp_granges, "GRangesList")))

plot_overlaps <- as_tibble(mcols(df_bp_granges_all)) %>% 
    tidyr::pivot_longer(
        cols = everything(),
        names_to = "Key",
        values_to = "Value"
    ) %>% 
    dplyr::mutate(
        hex = dplyr::case_when(
            Key == 'True' ~ 'black',
            Key == col_names[1] ~ '#990000',
            Key == col_names[2] ~ '#fd952c',
            Key == col_names[3] ~ '#669e85',
            Key == col_names[4] ~ '#394d7e'
        ),
        Key = factor(Key, levels = c('True', col_names))
    )
plot_overlaps <- dplyr::group_split(plot_overlaps, Key)

plots <- lapply(1:length(plot_overlaps), function(x){
    cutoff <- as_tibble(plot_overlaps[[x]]) %>% 
        dplyr::group_by(Value) %>% 
        dplyr::summarise(count = dplyr::n()) %>% 
        dplyr::mutate(
            frac = count/sum(count)*100,
            cum_frac = cumsum(frac)
        ) %>% 
        dplyr::filter(cum_frac <= 99) %>%
        tail(n = 1) %>% 
        dplyr::pull(Value)

    p1 <- plot_overlaps[[x]] %>% 
        ggplot(aes(x = Value, y = after_stat(density), col = hex)) + 
        geom_line(stat = "density", linewidth = 1.5) + 
        facet_wrap(vars(Key), nrow = 1, scales = "free") +  
        scale_color_identity() + 
        theme_bw() + 
        theme_classic() + 
        theme(text = element_text(size = 15)) + 
        coord_cartesian(xlim = c(0, cutoff)) + 
        labs(x = "", y = "")

    return(p1)
})

# plotting the 99% percentile
pdf(
    file = paste0(cur_path, "/Breaks_in_bins_density.pdf"),
    height = 5, width = 16
)
do.call(gridExtra::grid.arrange, c(plots, nrow = 1))
plot.saved <- dev.off()

#' 4. Calculate the probability per base for true/pred breaks
get_pred_data <- function(chr){
    cols_to_keep <- c("seqnames", "start", col_names)
    chr_name <- basename(pred_files[chr])
    chr_name <- as.integer(gsub(
        "final_pred_chr_|\\.parquet", "", chr_name
    ))

    df <- arrow::read_parquet(pred_files[chr])
    df <- as.data.table(df)
    setnames(df, c("start", "True", col_names))
    df[, seqnames := as.factor(paste0("chr", chr_name))]
    df <- df[, ..cols_to_keep]

    df_pos <- lapply(1:length(col_names), function(x){
        filter_cols <- c("seqnames", "start", col_names[x])
        temp <- df[, ..filter_cols]
        setnames(temp, c("seqnames", "start", "pred"))
        temp <- temp[pred != 0]
        temp[, pred := NULL]
        return(temp)
    })
    names(df_pos) <- col_names
    return(df_pos)
}

df_pred_all <- list.files(
    path = "../data/Whole_genome_pred",
    pattern = "Pred_|\\.parquet"
)
if(length(df_pred_all) != length(col_names)){
    df_predictions <- lapply(1:length(pred_files), function(chr){
        return(get_pred_data(chr = chr))
    })
    df_pred_one <- rbindlist(sapply(df_predictions, `[`, 1))
    df_pred_two <- rbindlist(sapply(df_predictions, `[`, 2))
    df_pred_three <- rbindlist(sapply(df_predictions, `[`, 3))
    df_pred_four <- rbindlist(sapply(df_predictions, `[`, 4))
    df_pred_all <- list(df_pred_one, df_pred_two, df_pred_three, df_pred_four)
    names(df_pred_all) <- col_names

    dir.create("../data/Whole_genome_pred/", showWarnings = FALSE)
    arrow::write_parquet(df_pred_one, paste0("../data/Whole_genome_pred/", col_names[1], ".parquet"))
    arrow::write_parquet(df_pred_two, paste0("../data/Whole_genome_pred/", col_names[2], ".parquet"))
    arrow::write_parquet(df_pred_three, paste0("../data/Whole_genome_pred/", col_names[3], ".parquet"))
    arrow::write_parquet(df_pred_four, paste0("../data/Whole_genome_pred/", col_names[4], ".parquet"))
} else {
    df_pred_all <- lapply(col_names, function(x){
        return(arrow::read_parquet(paste0("../data/Whole_genome_pred/", x, ".parquet")))
    })
    names(df_pred_all) <- col_names
}

# prepare GRanges objects
granges_list <- lapply(groups_split, plyranges::as_granges)
true_breaks_granges <- df_true %>% 
    dplyr::mutate(width = 1) %>% 
    plyranges::as_granges() %>% 
    unique()
pred_breaks_granges <- lapply(df_pred_all, function(df){
    return(
        df %>% 
            dplyr::mutate(width = 1) %>% 
            plyranges::as_granges() %>% 
            unique()
    )
})

# #' Thought:
# #' 1. Sum Count column of Genes = Sum_Count_Genes
# #' 2. Sum width column of Genes = Sum_width_Genes
# #' 3. Gens PPB is Sum_Count_Genes / Sum_width_Genes
# #' 4. For each type, get subset of gene_id, sum Count, sum width, divide Count by width.
# genes_ppb_filtered_split <- genes_ppb %>% 
#     dplyr::select(-norm_counts, -hex, -type) %>% 
#     tidyr::pivot_longer(
#         cols = -c(gene_id, Count, width, Pred_0_001, Pred_0_005, Pred_0_01, Pred_0_05),
#         names_to = "Key",
#         values_to = "Value"
#     ) %>%
#     dplyr::group_split(., Key)

# group_map_subset <- group_map[group == "genes"]
# group_map_subset[, other_label := c(
#     "Genes", "CFS", "CDG", "Fusion", 
#     "Housekeeping", "Oncogene", "TSG"
# )]

# genes_ppb_filtered <- lapply(1:length(genes_ppb_filtered_split), function(x){
#     feat_name <- genes_ppb_filtered_split[[x]]$Key[1]
#     temp_genes_subset <- genes_ppb_filtered_split[[x]] %>% 
#         dplyr::filter(Value) %>%
#         dplyr::select(-Value, -Key, -gene_id) %>% 
#         dplyr::rename(True = Count)

#     gene_widths <- sum(temp_genes_subset$width, na.rm = TRUE)

#     temp_genes_subset <- temp_genes_subset %>% 
#         dplyr::select(-width)
#     count_all_breaks <- apply(temp_genes_subset, 2, sum, na.rm = TRUE)
#     count_all_breaks <- count_all_breaks / gene_widths

#     count_all_breaks <- as.data.frame(t(count_all_breaks))
#     count_all_breaks$type <- feat_name
#     return(count_all_breaks)
# })
# genes_ppb_filtered <- do.call(rbind, genes_ppb_filtered)
# rownames(genes_ppb_filtered) <- group_map_subset$type[match(
#     genes_ppb_filtered$type, group_map_subset$other_label
# )]
# genes_ppb_filtered$type <- NULL

calc_prob_per_base <- function(df_breaks, df_genic_feat){
    overlaps <- GenomicRanges::findOverlaps(df_genic_feat, df_breaks)

    if(length(overlaps) == 0){
        overlap_counts <- 0 
    } else {
        overlap_counts <- sum(table(queryHits(overlaps)))
    }
    return(overlap_counts)
}

ppb <- pbapply::pbsapply(1:length(granges_list), function(ind){
    df_genic_feat <- granges_list[[ind]]
    genic_feat_len <- sum(width(reduce(df_genic_feat)))

    # Calculate for true breaks
    res_true <- calc_prob_per_base(
        df_breaks = true_breaks_granges, 
        df_genic_feat = df_genic_feat
    )
    res_true <- res_true / genic_feat_len

    # predicted breaks
    res_pred <- lapply(1:length(col_names), function(x){
        res_pred <- calc_prob_per_base(
            df_breaks = pred_breaks_granges[[x]], 
            df_genic_feat = df_genic_feat
        )
        res_pred <- res_pred / genic_feat_len
        return(res_pred)
    })
    names(res_pred) <- col_names

    return(c("True" = res_true, res_pred))
})

ppb_t <- t(ppb)
ppb_df <- as_tibble(ppb_t) %>% 
    dplyr::mutate(across(where(is.list), as.numeric)) %>% 
    as.data.frame()
rownames(ppb_df) <- names(groups_split)

# rownames(ppb_df) <- names(groups_split)[-(2:7)]
# temp = data.frame(
#     type = group_map[group == "genes"]$type,
#     True = 0,
#     Pred_0_001 = 0,
#     Pred_0_005 = 0,
#     Pred_0_01 = 0,
#     Pred_0_05 = 0
#     )
# rownames(temp) <- group_map[group == "genes"]$type
# temp <- temp[,-1]
# ppb_df <- rbind(temp, ppb_df[-1,])

# # replace the gene values
# match_ind <- match(rownames(genes_ppb_filtered), rownames(ppb_df))
# ppb_df[match_ind,] <- genes_ppb_filtered

# Normalise between 0 and 1 to generate a relative importance plot
ppb_df <- apply(ppb_df, 2, function(x) x / max(x))

hex <- c("#990000", "#fd952c", "#669e85","#394d7e", "#8A2BE2")
hex <- setNames(hex, group_names)

df_ppb <- as_tibble(ppb_df) %>% 
    dplyr::mutate(type = names(groups_split)) %>% 
    dplyr::arrange(desc(ppb_df)) %>% 
    dplyr::left_join(group_map, by = "type") %>% 
    dplyr::mutate(
        type = forcats::fct_inorder(type),
        group = as.factor(group),
        hex = hex[group]
    )

fwrite(
    as.data.table(df_ppb),
    "../data/Whole_genome_pred/Probability_per_base.csv"
)

plots_ppb <- lapply(c("True", col_names), function(column_name){
    temp <- df_ppb[c(column_name, "type", "group", "hex")] %>% 
        dplyr::rename_with(~c("ppb", "type", "group", "hex")) %>% 
        dplyr::arrange(ppb) %>% 
        dplyr::mutate(type = forcats::fct_inorder(type))


    plot_title <- ifelse(
        column_name == "True", "True",
        stringr::str_replace(
            string = column_name, 
            pattern = "^(Pred)_0_(.*)$", 
            replacement = "\\1: 0.\\2 (FPR)"
        )
    )

    p1 <- temp %>% 
        ggplot2::ggplot(aes(x = ppb, y = type)) +
        ggplot2::geom_segment(aes(
            x = 0, xend = ppb, 
            y = type, yend = type),
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
    #         "../figures/Whole_genome_pred/",
    #         "All_Groups_PPB_", column_name, ".pdf"
    #     ),
    #     plot = p1,
    #     height = 11, width = 8
    # )

    return(p1)
})

pdf(
    file = paste0(
        "../figures/Whole_genome_pred/",
        "All_Groups_PPB_All.pdf"
    ),
    height = 9, width = 25
)
do.call(gridExtra::grid.arrange, c(plots_ppb, nrow = 1))
plots.saved <- dev.off()

#' 5. Get overlaps of genic features with true and predicted breaks
max_unique_hex <- as_tibble(group_map) %>% 
    dplyr::group_by(group) %>% 
    dplyr::summarise(count = dplyr::n()) %>% 
    dplyr::pull(count) %>% 
    max()

# Generate two sets of contrasting colors
feat_hex <- c(
    "#720000",
    "#b81132",
    "#A52A2A",
    "#9d5002",
    "#e93c3e",
    "#d86b00",
    "#B8860B",
    "#fd9f3f",
    "#1a6b1a",
    "#228B22",
    "#4DAF4A",
    "#669e85",
    "#7edf7e",
    "#2d3c63",
    "#455e99",
    "#546fb2",
    "#377EB8",
    "#1E90FF",
    "#984EA3",
    "#761cca",
    "#8A2BE2",
    "#a75fe9"
)
pred_hex <- c("#b2b2b2", "#7f7f7f", "#4c4c4c", "#000000")

#' Top 5% most fragile:
#'      1) Cancer census genes (cancer driver genes)
#'      2) Oncogenes
#'      3) Tumour suppressor genes

genes_subset_list <- c(
    "Genes",
    "Cancer_driver_genes",
    "Oncogenes",
    "Tumour_Suppressor_Genes"
)
genesgenes_subset_list_hex = c("#8B8000", "#FFA435", "#BFABCB", "#BB357E")
genes_subset_list_hex <- setNames(genesgenes_subset_list_hex, genes_subset_list)

genes_subset <- df_genic[["genes"]] %>% 
    dplyr::filter(type %in% genes_subset_list)

ppb_fragile_genes <- pbapply::pblapply(genes_subset_list, function(g){
    df_genic_feat <- dplyr::filter(genes_subset, type == g)
    df_genic_feat <- plyranges::as_granges(df_genic_feat)

    get_norm_overlap_counts <- function(breaks_table, feature_table, overlap_vector){
        # Overlaps per genomic feature
        overlap_counts <- table(queryHits(overlap_vector))
        overlap_counts <- as.data.table(overlap_counts)
        setnames(overlap_counts, c("ID", "Count"))
        overlap_counts[, ID := as.integer(ID)]

        # For each gene, get the total length (width)
        gene_lengths <- width(feature_table[overlap_counts$ID])
        overlap_counts[, width := gene_lengths]

        # get gene names
        overlap_counts[, gene_id := mcols(feature_table[ID,])$gene_id]

        # get total length of the gene
        overlap_counts <- overlap_counts[, .(Count = sum(Count), width = sum(width)), by = gene_id]

        # get normalised count for each gene
        overlap_counts[, norm_counts := Count / width]

        # genes not present will be zero
        gene_names <- mcols(feature_table)$gene_id
        missing_genes <- which(is.na(match(gene_names, overlap_counts$gene_id)))
        missing_genes <- data.table(
            width = width(feature_table[missing_genes]),
            gene_id = gene_names[missing_genes]
        )
        missing_genes <- missing_genes[, .(width = sum(width)), by = gene_id]
        missing_genes[, `:=`(Count = 0, norm_counts = 0)]
        setcolorder(missing_genes, c("gene_id", "Count", "width", "norm_counts"))

        # combine table
        res <- rbind(overlap_counts, missing_genes)
        setorder(res, gene_id)
        
        return(res)
    }

    true_overlaps <- GenomicRanges::findOverlaps(
        df_genic_feat, true_breaks_granges
    )
    res_true <- get_norm_overlap_counts(
        breaks_table = true_breaks_granges, 
        feature_table = df_genic_feat, 
        overlap_vector = true_overlaps
    )

    # predicted breaks
    res_pred <- sapply(col_names, function(x){
        pred_overlaps <- GenomicRanges::findOverlaps(
            df_genic_feat, pred_breaks_granges[[x]]
        )
        res_pred <- get_norm_overlap_counts(
            breaks_table = pred_breaks_granges[[x]], 
            feature_table = df_genic_feat, 
            overlap_vector = pred_overlaps
        )
        return(res_pred$Count)
    })
    res_pred <- as.data.table(res_pred)

    res_true <- cbind(res_true, res_pred)
    res_true[, type := g]

    # pred_overlaps <- GenomicRanges::findOverlaps(
    #     df_genic_feat, pred_breaks_granges[[1]]
    # )
    # res_pred <- get_norm_overlap_counts(
    #     breaks_table = pred_breaks_granges[[1]], 
    #     feature_table = df_genic_feat, 
    #     overlap_vector = pred_overlaps
    # )
    # res_true[, `:=`(
    #     pred = res_pred$norm_counts[match(res_pred$gene_id, res_true$gene_id)],
    #     type = g
    # )]
    return(res_true)
})
ppb_fragile_genes_dt <- rbindlist(ppb_fragile_genes)
ppb_fragile_genes_dt <- as_tibble(ppb_fragile_genes_dt) %>% 
    dplyr::mutate(hex = genes_subset_list_hex[type])

##########################################
#' Note:
#'  norm_counts: ppb from true dataset.
#'  pred: ppb from pred dataset.
##########################################

# Plot distribution of all genes vs. distribution of above
p1 <- ppb_fragile_genes_dt %>%
    dplyr::select(
        gene_id, Count, width, norm_counts, 
        pred = Pred_0_001, type, hex
    ) %>% 
    dplyr::mutate(pred = pred / width) %>% 
    ggplot(aes(x = pred, y = after_stat(density), group = type, col = type)) + 
    geom_line(stat = "density", linewidth = 1.2, adjust = 3) + 
    theme_bw() + 
    theme_classic() + 
    scale_color_manual(values = ppb_fragile_genes_dt$hex) + 
    theme(text = element_text(size = 15)) + 
    coord_cartesian(xlim = c(0, 0.04)) + 
    labs(
        x = "Probability per base",
        y = "Density"
    )

ggsave(
    filename = paste0(
        "../figures/Whole_genome_pred/",
        "relative_fragility_density.pdf"
    ),
    plot = p1,
    height = 8, width = 9
)

# Get list of top 5% most fragile cancer driver genes, oncogenes and tumour suppressor genes
plots_ppb <- lapply(setdiff(genes_subset_list, "Genes"), function(g){
    temp <- ppb_fragile_genes_dt %>% 
        dplyr::filter(type == g) %>% 
        dplyr::select(
            gene_id, Count, width, norm_counts, 
            pred = Pred_0_001, type, hex
        ) %>% 
        dplyr::mutate(pred = pred / width) %>% 
        dplyr::slice_max(prop = 0.05, order_by = pred) %>%
        dplyr::mutate(pred = pred / max(pred)) %>% 
        dplyr::arrange(pred) %>% 
        dplyr::mutate(gene_id = forcats::fct_inorder(gene_id))

    p1 <- temp %>% 
        ggplot2::ggplot(aes(x = pred, y = gene_id)) +
        ggplot2::geom_segment(aes(
            x = 0, xend = pred, 
            y = gene_id, yend = gene_id),
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
            subtitle = gsub("_", " ", unique(temp$type)),
            x = "", y = ""
        )

    return(p1)
})

pdf(
    file = paste0(
        "../figures/Whole_genome_pred/",
        "relative_fragility_density_top5.pdf"
    ),
    height = 7, width = 14
)
do.call(gridExtra::grid.arrange, c(plots_ppb, nrow = 1))
plots.saved <- dev.off()

# Get list of all cancer driver genes, oncogenes and tumour suppressor genes
genes_ppb <- ppb_fragile_genes_dt %>% 
    dplyr::filter(type == "Genes") %>%
    dplyr::mutate(
        Genes = TRUE,
        Housekeeping = FALSE,
        Oncogene = FALSE,
        TSG = FALSE,
        Fusion = FALSE,
        CDG = FALSE,
        CFS = FALSE
    )

genes_ppb$Housekeeping[match(
    df_genic$genes[type == "Housekeeping_genes", gene_id], 
    genes_ppb$gene_id
)] <- TRUE
genes_ppb$Oncogene[match(
    df_genic$genes[type == "Oncogenes",gene_id], 
    genes_ppb$gene_id
)] <- TRUE
genes_ppb$TSG[match(
    df_genic$genes[type == "Tumour_Suppressor_Genes",gene_id], 
    genes_ppb$gene_id
)] <- TRUE
genes_ppb$Fusion[match(
    df_genic$genes[type == "Fusion_Genes",gene_id], 
    genes_ppb$gene_id
)] <- TRUE
genes_ppb$CDG[match(
    df_genic$genes[type == "Cancer_driver_genes",gene_id], 
    genes_ppb$gene_id
)] <- TRUE
genes_ppb$CFS[match(
    df_genic$genes[type == "CFS_genes",gene_id], 
    genes_ppb$gene_id
)] <- TRUE

genes_split <- genes_ppb %>% 
    dplyr::select(
        gene_id, Count, width, Pred_0_001, 
        Genes, Housekeeping, Oncogene, TSG,
        Fusion, CDG, CFS
    ) %>%
    tidyr::pivot_longer(
        cols = -c(gene_id, Count, width, Pred_0_001),
        names_to = "Key",
        values_to = "Value"
    ) %>%
    dplyr::group_split(., Key)

genes_split <- lapply(1:length(genes_split), function(x){
    return(
        genes_split[[x]] %>% 
            dplyr::filter(Value) %>%
            dplyr::select(-Value)
    )
})
genes_to_keep <- do.call(rbind, genes_split)

genes_split <- genes_to_keep %>%
    dplyr::filter(
        Key != "Genes" & 
        Key != "Fusion" & 
        Key != "CDG" & 
        Key != "CFS"
    ) %>% 
    dplyr::mutate(
        Key = forcats::fct_inorder(Key),
        hex = dplyr::case_when(
            Key == "Housekeeping" ~ "#e7e7e7",
            Key == "Oncogene" ~ "#BFABCB",
            Key == "TSG" ~ "#BB357E"
        )
    ) %>% 
    dplyr::rename(pred = Pred_0_001) 

labels <- as.character(unique(genes_split$Key))

genes_split %>% 
    dplyr::group_by(Key) %>% 
    dplyr::summarise(count = dplyr::n()) %>% 
    dplyr::ungroup()
 
p1 <- genes_split %>%
    ggplot(aes(x = Key, y = pred, fill = hex)) + 
    geom_boxplot(alpha = 0.75, outlier.shape = NA) + 
    ggsignif::geom_signif(
        comparisons = list(
            c(labels[1], labels[2]),
            c(labels[1], labels[3])
        ),
        test = "wilcox.test",
        map_signif_level = TRUE,
        y_position = c(1, 500),
        textsize = 5,
        margin_top = 0.0025,
        tip_length = 0
    ) +
    scale_fill_identity() + 
    scale_y_continuous(labels = scales::label_number(scale = 1e-3)) +
    coord_cartesian(ylim = c(0, 7000)) + 
    theme_bw() + 
    theme_classic() +
    theme(text = element_text(size = 15)) + 
    labs(
        x = "", 
        y = expression("Absolute number of breaks, x10"^3*"")
    )

ggsave(
    filename = paste0(
        "../figures/Whole_genome_pred/",
        "relative_fragility_absolute.pdf"
    ),
    plot = p1,
    height = 8, width = 8
)

p2 <- genes_split %>%
    dplyr::mutate(pred = pred / width) %>% 
    ggplot(aes(x = Key, y = pred, fill = hex)) + 
    geom_boxplot(alpha = 0.75, outlier.shape = NA) + 
    ggsignif::geom_signif(
        comparisons = list(
            c(labels[1], labels[2]),
            c(labels[1], labels[3])
        ),
        map_signif_level = TRUE,
        y_position = c(0.011, 0.0143),
        test = "wilcox.test",
        textsize = 5,
        vjust = -0.1,
        margin_top = 0.03,
        tip_length = 0
    ) +
    scale_fill_identity() + 
    coord_cartesian(ylim = c(0, 0.05)) + 
    theme_bw() + 
    theme_classic() +
    theme(text = element_text(size = 15)) + 
    labs(x = "", y = "Probability per base")

ggsave(
    filename = paste0(
        "../figures/Whole_genome_pred/",
        "relative_fragility_normalised.pdf"
    ),
    plot = p2,
    height = 8, width = 8
)

pdf(
    file = paste0(
        "../figures/Whole_genome_pred/",
        "relative_fragility_genes.pdf"
    ),
    height = 8, width = 16
)
gridExtra::grid.arrange(p1, p2, nrow = 1)
plot.saved <- dev.off()