# Order of installing packages for Ubuntu to get ggpattern
# sudo apt-get install libmagick++-dev
# sudo apt install libgdal-dev
# sudo apt-get install -y libudunits2-dev
# install.packages("units")
# install.packages("sf")
# install.packages("gridpattern")
# install.packages("ggpattern")

# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(ggpattern)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))
suppressPackageStartupMessages(suppressWarnings(library(arrow)))

setwd("/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/05_DeltaFragility")

args <- commandArgs(trailingOnly = TRUE)
bw <- as.numeric(args[1])

t1 <- Sys.time()
cur.msg <- "Importing datasets"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

# get thresholds and column names
thresholds <- fread("../data/models/python/lightgbm/best_LGBM_model_thresholds.csv")
col_names <- thresholds$FPR
col_names <- gsub("\\.", "_", paste0("Pred_", col_names))

plot_titles <- stringr::str_replace(
    string = col_names, 
    pattern = "^(Pred)_0_(.*)$", 
    replacement = "\\1: 0.\\2 (FPR)"
)
plot_titles <- setNames(plot_titles, col_names)

path_to_pred <- "../data/Whole_genome_pred"
df_random_frac <- fread(paste0(path_to_pred, "/random_sampling/random_frac.csv"))

df_pred_all <- lapply(col_names, function(x){
    return(arrow::read_parquet(paste0(path_to_pred, "/", x, ".parquet")))
})
names(df_pred_all) <- col_names

# import labels and predictions
df_original <- fread("./data/all_SV_predictions.csv")
labels <- arrow::read_parquet("./data/all_SV_labels.parquet")
Assembly <- labels$Assembly[match(df_original$labels, labels$label)]
df_original[, Assembly := Assembly]

# get short range of influence
range_cutoff <- fread("../data/Ranges_cutoffs.csv", select = c("Cluster", "kmer_8"))
LR <- range_cutoff[Cluster == "long.range", "kmer_8"][[1]]

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")

#' 1. Do pathogenic / benign mutations increase or decrease fragility?
#' Answer: pathogenic SNVs have a higher probability per base than benign.
t1 <- Sys.time()
cur.msg <- "Calculating change in SV fragility"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

df_all_split <- as_tibble(df_original) %>% 
    dplyr::group_split(before_after) %>%
    purrr::set_names(purrr::map(., ~ unique(.x$before_after)[1]))

# calculate delta fragilities
df_all <- df_all_split[['after']] %>%
    dplyr::select(dplyr::starts_with('Pred')) %>%
    dplyr::mutate(dplyr::across(
        dplyr::starts_with('Pred'), ~ . - df_all_split[['before']][[dplyr::cur_column()]]
    )) %>% 
    dplyr::bind_cols(df_all_split[['after']] %>% 
    dplyr::select(-starts_with('Pred')))

res_summary <- df_all %>% 
    dplyr::select(dplyr::starts_with("Pred_"), ClinicalSignificance) %>%
    tidyr::pivot_longer(
        cols = -c("ClinicalSignificance"),
        names_to = "Key",
        values_to = "Value"
    ) %>%
    dplyr::rename(`Clinical Sig.` = ClinicalSignificance) %>%
    dplyr::mutate(
        Key = plot_titles[Key],
        Key = forcats::fct_inorder(Key)
    )

res_signif <- res_summary %>% 
    split(.$Key) %>%
        purrr::map_dfr(~{
        d_split <- split(.x, .x$`Clinical Sig.`)
        t_test <- t.test(d_split$Pathogenic$Value, d_split$Benign$Value)
        
        tibble(
            Key = unique(.x$Key),
            p_value = t_test$p.value,
            Benign_Value = max(d_split$Benign$Value),
            Pathogenic_Value = max(d_split$Pathogenic$Value),
            Value = mean(Benign_Value, Pathogenic_Value) / 25
        )
    }) %>%
    dplyr::mutate(
        Key = plot_titles[Key],
        Key = forcats::fct_inorder(Key),
        Signif_level = dplyr::case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            TRUE ~ "N.S."
        ),
        # y_position = c(0.006, 0.013, 0.016, 0.032),
        y_position = c(10,20,30,60),
        xmin = (1:length(col_names)) - 0.188,
        xmax = (1:length(col_names)) + 0.188,
        annotation = Signif_level
    )

totals <- res_summary %>% 
    dplyr::group_by(`Clinical Sig.`) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::pull(Count, `Clinical Sig.`)

ylim_range <- mean(res_summary$Value) + sd(res_summary$Value) * 1.4

p1 <- res_summary %>%
    ggplot(aes(
        x = Key, 
        y = Value, 
        fill = `Clinical Sig.`,
        pattern = `Clinical Sig.`
    )) + 
    geom_boxplot_pattern(
        outlier.shape = NA,
        col = "black",
        alpha = 0.75,
        pattern_colour = "white",
        pattern_spacing = 0.03,
        pattern_angle = 45
    ) +
    ggsignif::geom_signif(
        y_position = res_signif$y_position,
        xmin = res_signif$xmin,
        xmax = res_signif$xmax,
        annotation = res_signif$annotation,
        textsize = 5,
        vjust = -0.1,
        margin_top = 0.1,
        tip_length = 0.001
    ) +
    coord_cartesian(ylim = c(-ylim_range, ylim_range)) + 
    scale_fill_manual(values = rep("black", 2)) + 
    scale_pattern_manual(values = c(
        'crosshatch', 'none'
    )) +
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) +
    labs(
        title = paste0(
            "Delta fragility prediction of ClinVar's ",
            nrow(res_summary), " SVs"
        ),
        subtitle = paste0(
            "Total Benign entries: ", totals[["Benign"]], ". ", 
            "Total Pathogenic entries: ", totals[["Pathogenic"]], "."
        ),
        x = "", y = "Change in fragility upon SV"
    )

dir.create("../figures/deltafragility", showWarnings = FALSE)
ggsave(
    filename = "../figures/deltafragility/benign_vs_pathogenic_SV_boxplot.pdf",
    # filename = "demo.pdf",
    plot = p1,
    height = 6, width = 9
)

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")

#' Focusing on the biggest delta, do they happen in highly fragile regions?
t1 <- Sys.time()
cur.msg <- "Lifting over breaks to T2T version"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

hg38_chain <- import.chain("../../05_Cosmic/data/liftover/hg38-chm13v2.over.chain")
hg19_chain <- import.chain("../../05_Cosmic/data/liftover/hg19ToHg38.over.chain")

all_df_delta_liftover <- lapply(1:length(col_names), function(x){
    df_delta <- as_tibble(df_all) %>% 
        dplyr::select(-setdiff(col_names, col_names[x]), -labels) %>%
        dplyr::rename(Value = col_names[x]) %>% 
        dplyr::arrange(before_after, GeneSymbol, ClinicalSignificance)

    # liftover positions to t2t genome version.
    df_delta_split <- df_delta %>% 
        dplyr::group_split(Assembly) %>% 
        purrr::set_names(purrr::map(., ~ unique(.x$Assembly)[1]))

    df_delta_liftover <- lapply(df_delta_split, function(x){
        x <- x %>% 
            dplyr::mutate(width = 1) %>%
            plyranges::as_granges()

        genome_assembly <- mcols(x)$Assembly[1]

        # liftover to t2t genome
        if(grepl("GRCh37", genome_assembly)){
            liftover_hg19 <- suppressMessages(liftOver(x, hg19_chain))
            x <- unlist(as(liftover_hg19, "GRangesList"))
        } 
        liftover_hg38 <- suppressMessages(liftOver(x, hg38_chain))
        liftover_hg38 <- unlist(as(liftover_hg38, "GRangesList"))
        liftover_hg38 <- liftover_hg38 %>% 
            dplyr::select(-Assembly, -before_after)

        return(liftover_hg38)
    })
    df_delta_liftover <- do.call(plyranges::bind_ranges, df_delta_liftover)  

    df_delta_liftover <- df_delta_liftover %>% 
        dplyr::mutate(key = col_names[x])

    return(df_delta_liftover)
})
names(all_df_delta_liftover) <- col_names

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")

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

# import randomly sampled breakpoints coinciding with low/mid/high fragile zones
t1 <- Sys.time()
cur.msg <- "Overlap of breaks with hot/cold break zones"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

pred_categories <- lapply(1:length(col_names), function(x){
    df_bp_granges <- lapply(1:22, function(chr){
        pred_bins <- calculate_overlaps(
            df_breaks = df_pred_all[[col_names[x]]], 
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
                TRUE ~ "Mid"
            ),
            group = as.factor(group)
        ) %>% 
        dplyr::select(-ID) 

    df_bp_granges_all <- df_bp_granges_all %>% 
        dplyr::mutate(group = df_breaks$group)

    # Overlap of SV prediction with low/medium/high fragile regions of the human genome
    overlaps <- GenomicRanges::findOverlaps(all_df_delta_liftover[[col_names[x]]], df_bp_granges_all)
    overlaps <- all_df_delta_liftover[[col_names[x]]][queryHits(overlaps)] %>%
        dplyr::mutate(group = mcols(df_bp_granges_all)$group[subjectHits(overlaps)])

    dt_overlaps <- as.data.table(mcols(overlaps))
    dt_overlaps[, key := col_names[x]]
    return(dt_overlaps)
})
df_all_cat <- rbindlist(pred_categories)

fwrite(
    df_all_cat,
    paste0("./data/df_all_cat_bw_", bw, ".csv")
)

#' perform t-test between:
#'      A) mean change in fragility per clinical label
#'      B) mean change in fragility per fragile zone per clinical label
clin_p <- lapply(c("Benign", "Pathogenic"), function(clin_label){
    groups <- unique(df_all_cat$group)
    threshold_p <- sapply(col_names, function(model_threshold){
        population_df <- as_tibble(df_all_cat) %>% 
            dplyr::filter(
                ClinicalSignificance == clin_label & 
                key == model_threshold
            ) %>% 
            dplyr::pull(Value)

        group_p <- sapply(groups, function(g){
                sample_df <- as_tibble(df_all_cat) %>% 
                    dplyr::filter(
                        ClinicalSignificance == clin_label & 
                        key == model_threshold & 
                        group == g
                    ) %>% 
                    dplyr::pull(Value)

                t_test <- t.test(sample_df, population_df)
                return(t_test$p.value)
        })
        return(group_p)
    })
    threshold_p <- as.data.frame(threshold_p)
    threshold_p$ClinicalSignificance <- clin_label
    threshold_p$group <- groups
    return(threshold_p)
})
p_values <- do.call(rbind, clin_p)

fwrite(
    as.data.table(p_values),
    "./data/Difference_in_means_pvalues.csv"
)

only_sig_pvalues <- as_tibble(p_values) %>% 
    tidyr::pivot_longer(
        cols = -c(ClinicalSignificance, group),
        names_to = "key",
        values_to = "p_value"
    ) %>% 
    dplyr::mutate(
        Signif_level = dplyr::case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            TRUE ~ "N.S."
        )
    ) %>% 
    dplyr::filter(Signif_level != "N.S.")

fwrite(
    as.data.table(only_sig_pvalues),
    "./data/Difference_in_means_pvalues.csv"
)

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")

#' Focusing on the biggest delta, do they happen in highly fragile regions?
temp_all_sv <- df_all_cat %>%
    dplyr::group_by(key, ClinicalSignificance, group) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::group_by(key, ClinicalSignificance) %>%
    dplyr::mutate(frac = count / sum(count) * 100) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        ClinicalSignificance = as.factor(ClinicalSignificance),
        group = factor(group, levels = c("High", "Mid", "Low")),
        hex = dplyr::case_when(
            group == "High" ~ "#b2182b",
            group == "Mid" ~ "#006400",
            group == "Low" ~ "#ff9e27"
        )
    ) %>%
    dplyr::rename(`Fragile zone` = group) %>%
    dplyr::mutate(
        key = plot_titles[key],
        key = forcats::fct_inorder(key)
    ) %>%
    suppressMessages()

# get the random sampled breakpoint positions and merge with predicted ones
df_random_frac_subset <- as_tibble(df_random_frac) %>% 
    dplyr::filter(bin_width == bw) %>% 
    dplyr::mutate(key = plot_titles[threshold]) %>% 
    dplyr::select(-bin_width, -threshold) %>% 
    dplyr::mutate(ClinicalSignificance = c("Random")) %>% 
    dplyr::arrange(key, ClinicalSignificance) %>% 
    dplyr::mutate(
        group = factor(group, levels = c("High", "Mid", "Low")),
        hex = dplyr::case_when(
            group == "High" ~ "#b2182b",
            group == "Mid" ~ "#006400",
            group == "Low" ~ "#ff9e27"
        )
    ) %>% 
    dplyr::rename(`Fragile zone` = group)

temp_all_sv_rand <- dplyr::bind_rows(
    x = temp_all_sv, 
    y = df_random_frac_subset
    ) %>% 
    dplyr::arrange(key, ClinicalSignificance) %>% 
    dplyr::mutate(ClinicalSignificance = as.factor(ClinicalSignificance))

p1 <- temp_all_sv_rand %>%
    ggplot(aes(x = ClinicalSignificance, y = frac, fill = `Fragile zone`)) + 
    facet_wrap(vars(key), nrow = 1) + 
    geom_bar(
        position = "stack", 
        stat = "identity",
        col = "black",
        alpha = 0.75
    ) +
    geom_text(
        # data = temp_all_sv %>%
        #     dplyr::filter(`Fragile zone` %in% c("High", "Mid")), 
        aes(
            label = sprintf("%.1f%%", frac), 
            y = frac
        ), 
        position = position_stack(vjust = 0.5),
        col = "black",
        size = 5
    ) +
    theme_bw() + 
    theme_classic() + 
    scale_fill_manual(values = temp_all_sv$hex) + 
    theme(text = element_text(size = 15)) +
    labs(
        title = "ClinVar SV fragility prediction",
        subtitle = "Fraction of all SVs in low/mid/high fragile zones",
        x = "", y = "Fraction, %"
    )

ggsave(
    filename = paste0(
        "../figures/deltafragility/All_Delta_",
        "bw_", bw, ".pdf"
    ),
    plot = p1,
    height = 5, width = 12
)

# calculate p-values for all proportions 
# bw=20000
# df_all_cat <- fread(paste0("./data/df_all_cat_bw_", bw, ".csv"))

true_frac_clean <- df_all_cat %>%
    dplyr::group_by(key, ClinicalSignificance, group) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::group_by(key, ClinicalSignificance) %>%
    dplyr::mutate(
        total = sum(count),
        frac = count / total * 100
    ) %>%
    dplyr::ungroup() %>% 
    dplyr::select(key, ClinicalSignificance, group, count, total, frac) %>% 
    suppressMessages()

random_frac_clean <- as_tibble(df_random_frac) %>% 
    dplyr::filter(bin_width == bw) %>% 
    dplyr::group_by(bin_width, threshold) %>% 
    dplyr::mutate(
        total = sum(count),
        frac_check = count / total * 100,
        ClinicalSignificance = "Random"
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(key = threshold) %>% 
    dplyr::select(key, ClinicalSignificance, group, count, total, frac, -bin_width)

frac_clean <- dplyr::bind_rows(
        true_frac_clean, 
        random_frac_clean
    ) %>% 
    dplyr::arrange(key, ClinicalSignificance) %>% 
    dplyr::mutate(ClinicalSignificance = as.factor(ClinicalSignificance))

# calculate prop.tests
frac_clean_pvalues <- frac_clean %>% 
    dplyr::group_split(key, group) %>% 
    purrr::map_dfr(~{
        # Extract 'rand' and 'pathogenic' data
        random_data <- .x %>% filter(ClinicalSignificance == "Random")
        pathogenic_data <- .x %>% filter(ClinicalSignificance == "Pathogenic")
        benign_data <- .x %>% filter(ClinicalSignificance == "Benign")

        # Perform z-test of difference in proportions
        prop_test_pathogenic <- prop.test(
            x = c(random_data$count, pathogenic_data$count), 
            n = c(random_data$total, pathogenic_data$total)
        )

        prop_test_benign <- prop.test(
            x = c(random_data$count, benign_data$count), 
            n = c(random_data$total, benign_data$total)
        )

        # Create a summary tibble
        tibble(
            key = .x$key[1],
            group = unique(.x$group),
            p_value = c(prop_test_pathogenic$p.value, prop_test_benign$p.value),
            ClinicalSignificance = c("Pathogenic", "Benign")
        )
    }) %>% 
    dplyr::mutate(
        key = plot_titles[key],
        key = forcats::fct_inorder(key),
        signif_level = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            TRUE ~ "N.S."
        )
    )

fwrite(
    as.data.table(frac_clean_pvalues),
    paste0("./data/Difference_in_means_pvalues_bw_", bw, ".csv")
)

#################################################################################
#' SVs in highly fragile regions: do they become more or less fragile after SV?
sv_by_clin <- as_tibble(df_all_cat) %>% 
    dplyr::rename(type = ClinicalSignificance) %>% 
    dplyr::mutate(group = factor(group, levels = c("Low", "Mid", "High"))) %>% 
    dplyr::select(-Origin)

SV_high <- as_tibble(sv_by_clin) %>% dplyr::filter(group == "High")
SV_mid <- as_tibble(sv_by_clin) %>% dplyr::filter(group == "Mid")
SV_low <- as_tibble(sv_by_clin) %>% dplyr::filter(group == "Low")

group_cols <- c(
    # green - yellow - red
    "#006400", "#ff9e27", "#b2182b",

    # green - yellow - red
    "#006400", "#ff9e27", "#b2182b"
)

p_sv_by_clin <- dplyr::bind_rows(SV_low, SV_mid, SV_high) %>% 
    dplyr::group_by(type, group, key) %>% 
    dplyr::summarise(
        mean = mean(Value),
        sd = sd(Value),
        n = dplyr::n(),  
        # se = sd / sqrt(n),
        se = qnorm(0.975) * sd / sqrt(n),
        lower_se = mean - se,
        upper_se = mean + se
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
        combined = interaction(group, type),
        key = plot_titles[key],
        key = forcats::fct_inorder(key)
    ) %>% 
    suppressMessages() %>%
    ggplot(aes(
        x = group, y = mean, 
        fill = combined,
        pattern = combined
    )) + 
    geom_bar_pattern(
        stat = "identity", 
        position = position_dodge(), 
        color = "black",
        alpha = 0.75,
        pattern_colour = "white",
        pattern_spacing = 0.06,
        pattern_frequency = 10, 
        pattern_angle = 45
    ) +
    geom_errorbar(
        aes(ymin = lower_se, ymax = upper_se),
        position = position_dodge(.9), 
        width = 0.25
    ) + 
    facet_wrap(vars(key), nrow = 1) + 
    theme_bw() + 
    theme_classic() +
    theme(
        text = element_text(size = 15),
        legend.position = "none"
    ) + 
    scale_fill_manual(values = group_cols) +
    scale_pattern_manual(values = c(
        rep('crosshatch', 3),
        rep('none', 3)
    )) +
    labs(x = "Fragile zones", y = "Change in fragility upon SV")
    
ggsave(
    filename = paste0(
        "../figures/deltafragility/Avg_delta_fragility_by_clin_",
        "bw_", bw, ".pdf"
    ),
    plot = p_sv_by_clin,
    height = 4, width = 8
)

p_sv_all <- dplyr::bind_rows(SV_low, SV_mid, SV_high) %>% 
    dplyr::group_by(group, key) %>% 
    dplyr::summarise(
        mean = mean(Value),
        sd = sd(Value),
        n = dplyr::n(),  
        # se = sd / sqrt(n),
        se = qnorm(0.975) * sd / sqrt(n),
        lower_se = mean - se,
        upper_se = mean + se
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
        key = plot_titles[key],
        key = forcats::fct_inorder(key),
        hex = dplyr::case_when(
            group == "High" ~ "#b2182b",
            group == "Mid" ~ "#ff9e27",
            group == "Low" ~ "#006400"
        )
    ) %>% 
    suppressMessages() %>% 
    ggplot(aes(x = group, y = mean, fill = hex)) + 
    geom_bar(
        position = position_dodge(),
        stat = "identity",
        col = "black",
        alpha = 0.75
    ) +
    geom_errorbar(
        aes(ymin = lower_se, ymax = upper_se),
        position = position_dodge(.9), 
        width = 0.25
    ) + 
    facet_wrap(vars(key), nrow = 1) + 
    theme_bw() + 
    theme_classic() +
    theme(text = element_text(size = 15)) + 
    scale_fill_identity() + 
    labs(x = "Fragile zones", y = "Change in fragility upon SV")
    
ggsave(
    filename = paste0(
        "../figures/deltafragility/Avg_delta_fragility_comb_",
        "bw_", bw, ".pdf"
    ),
    plot = p_sv_all,
    height = 4, width = 8
)

sv_by_clin <- dplyr::bind_rows(SV_low, SV_mid, SV_high) %>% 
    dplyr::mutate(
        group = factor(group, levels = c("Low", "Mid", "High", "All"))
    ) %>% 
    dplyr::group_by(type, group, key) %>% 
    dplyr::summarise(
        mean = mean(Value),
        sd = sd(Value),
        n = dplyr::n(),  
        # se = sd / sqrt(n),
        se = qnorm(0.975) * sd / sqrt(n),
        lower_se = mean - se,
        upper_se = mean + se
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
        key = plot_titles[key],
        key = forcats::fct_inorder(key),
        hex = dplyr::case_when(
            group == "High" ~ "#b2182b",
            group == "Mid" ~ "#ff9e27",
            group == "Low" ~ "#006400"
        )
    ) %>% 
    dplyr::mutate(group = factor(group, levels = c("Low", "Mid", "High"))) %>%
    dplyr::filter(key == plot_titles[3]) %>% 
    dplyr::select(-key) %>%
    suppressMessages()

combines_labels <- interaction(sv_by_clin$group, sv_by_clin$type)
levels(combines_labels) <- as.character(combines_labels)
sv_by_clin$combined <- combines_labels

p1_sv_by_clin <- sv_by_clin %>%
    ggplot(aes(
        x = group, y = mean, 
        fill = combined,
        pattern = combined
        # pattern_type = combined
    )) + 
    geom_bar_pattern(
        stat = "identity", 
        position = position_dodge(), 
        color = "black",
        alpha = 0.75,
        pattern_colour = "white",
        pattern_spacing = 0.02,
        pattern_angle = 45
    ) +
    # geom_errorbar(
    #     aes(ymin = lower_se, ymax = upper_se),
    #     position = position_dodge(.9), 
    #     stat = "identity", 
    #     color = "black",
    #     alpha = 0.75
    # ) +
    scale_fill_manual(values = group_cols) +
    scale_pattern_manual(values = c(
        rep('crosshatch', 3),
        rep('none', 3)
    )) +
    coord_cartesian(ylim = c(NA, 3)) + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        plot.subtitle = element_text(size = 12),
        legend.position = "none"
    ) +
    labs(x = "", y = "Change in regional fragility upon SV")
    
ggsave(
    filename = paste0(
        "../figures/deltafragility/For_Paper_Avg_delta_fragility_by_clin_",
        "bw_", bw, ".pdf"
    ),
    plot = p1_sv_by_clin,
    height = 6, width = 6
)

# #################################################################################
# df_top_5_perc <- as_tibble(df_all_cat) %>% 
#     dplyr::group_by(key, ClinicalSignificance) %>%
#     dplyr::slice_max(order_by = delta, prop = 0.05) %>%
#     dplyr::ungroup()

# df_bot_5_perc <- as_tibble(df_all_cat) %>% 
#     dplyr::group_by(key, ClinicalSignificance) %>%
#     dplyr::slice_min(order_by = delta, prop = 0.05) %>%
#     dplyr::ungroup()

# # plot stacked barplot, ranked by decreasing order of H group.
# df_all_cat <- dplyr::bind_rows(df_top_5_perc, df_bot_5_perc)
# temp_top_sv <- df_all_cat %>%
#     dplyr::group_by(key, ClinicalSignificance, group) %>%
#     dplyr::summarise(count = dplyr::n()) %>%
#     dplyr::group_by(key, ClinicalSignificance) %>%
#     dplyr::mutate(frac = count / sum(count) * 100) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(
#         ClinicalSignificance = as.factor(ClinicalSignificance),
#         group = factor(group, levels = c("High", "Mid", "Low")),
#         hex = dplyr::case_when(
#             group == "High" ~ "#b2182b",
#             group == "Mid" ~ "#006400",
#             group == "Low" ~ "#ff9e27"
#         )
#     ) %>%
#     dplyr::rename(`Fragile zone` = group) %>%
#     dplyr::mutate(
#         key = plot_titles[key],
#         key = forcats::fct_inorder(key)
#     ) %>%
#     suppressMessages() 

# p1 <- temp_top_sv %>%
#     ggplot(aes(x = ClinicalSignificance, y = frac, fill = `Fragile zone`)) + 
#     facet_wrap(vars(key), nrow = 1) + 
#     geom_bar(
#         position = "stack", 
#         stat = "identity",
#         col = "black",
#         alpha = 0.75
#     ) +
#     geom_text(
#         # data = temp_top_sv %>%
#         #     dplyr::filter(`Fragile zone` %in% c("High", "Mid")), 
#         aes(
#             label = sprintf("%.1f%%", frac), 
#             y = frac
#         ), 
#         position = position_stack(vjust = 0.5),
#         col = "black",
#         size = 5
#     ) +
#     theme_bw() + 
#     theme_classic() + 
#     scale_fill_manual(values = temp_top_sv$hex) + 
#     theme(text = element_text(size = 15)) +
#     labs(
#         title = "ClinVar SV fragility prediction",
#         subtitle = paste0(
#             "Fraction of the top 5% highest delta ",
#             "fragility in high/mid/low fragile zones"
#         ),
#         x = "", y = "Fraction, %"
#     )

# ggsave(
#     filename = paste0(
#         "../figures/deltafragility/Top5Per_Delta_",
#         "bw_", bw, ".pdf"
#     ),
#     plot = p1,
#     height = 5, width = 11
# )

# # temp_all_and_top <- dplyr::left_join(
# #         x = temp_all_sv,
# #         y = temp_top_sv,
# #         by = join_by(
# #             key, ClinicalSignificance, `Fragile zone`, hex
# #         ),
# #         suffix = c("_all", "_top_5")
# #     ) %>%
# #     dplyr::filter(`Fragile zone` == "High") %>%
# #     tidyr::pivot_longer(
# #         cols = -c(key, ClinicalSignificance, `Fragile zone`, hex),
# #         names_to = "Key",
# #         values_to = "Value" 
# #     )

# # round_to_nearest_even <- function(x) round(x/2)*2

# # temp_all_and_top_nocounts <- temp_all_and_top %>%
# #     dplyr::filter(!grepl("count", Key)) %>%
# #     dplyr::mutate(
# #         key = plot_titles[key],
# #         key = forcats::fct_inorder(key)
# #     )

# # # calculate proportion test
# # res_signif <- temp_top_sv %>% 
# #     dplyr::group_by(key, ClinicalSignificance) %>%
# #     dplyr::mutate(total = sum(count)) %>%
# #     dplyr::ungroup() %>%
# #     dplyr::filter(`Fragile zone` == "High") %>%
# #     split(.$key) %>%
# #     purrr::map_dfr(~{
# #         prop_test <- prop.test(.$count, .$total)
# #         tibble(
# #             ClinicalSignificance = .x$ClinicalSignificance,
# #             key = unique(.x$key),
# #             p_value = prop_test$p.value,
# #             Value = prop_test$estimate
# #         )
# #     }) %>%
# #     dplyr::filter(ClinicalSignificance == "Pathogenic") %>%
# #     dplyr::mutate(
# #         key = plot_titles[key],
# #         key = forcats::fct_inorder(key),
# #         Signif_level = dplyr::case_when(
# #             p_value < 0.001 ~ "***",
# #             p_value < 0.01 ~ "**",
# #             p_value < 0.05 ~ "*",
# #             .default = "N.S."
# #         ),
# #         y_position = c(61, 49, 43, 22),
# #         xmin = 1.188,
# #         xmax = 2.21,
# #         annotation = Signif_level
# #     )

# # unique_keys <- as.character(unique(temp_all_and_top_nocounts$key))

# # plots <- lapply(unique_keys, function(unique_id){
# #     filtered_df <- temp_all_and_top_nocounts %>% 
# #         dplyr::filter(key == unique_id)
# #     filtered_sig <- res_signif %>% 
# #         dplyr::filter(key == unique_id)
    
# #     p1 <- filtered_df %>%
# #         ggplot(aes(x = ClinicalSignificance, y = Value, fill = Key)) + 
# #         geom_bar(
# #             position = position_dodge(), 
# #             stat = "identity",
# #             col = "black",
# #             alpha = 0.75
# #         ) +
# #         geom_text(
# #             data = filtered_df,
# #             aes(label = round_to_nearest_even(Value)), 
# #             position = position_dodge(width = 0.9), 
# #             vjust = -0.25,
# #             size = 5
# #         ) +
# #         ggsignif::geom_signif(
# #             y_position = filtered_sig$y_position,
# #             xmin = res_signif$xmin,
# #             xmax = res_signif$xmax,
# #             annotation = filtered_sig$annotation,
# #             textsize = 5,
# #             vjust = -0.1,
# #             margin_top = 0.1,
# #             tip_length = 0.01
# #         ) +
# #         coord_cartesian(ylim = c(0, max(temp_all_and_top_nocounts$Value) * 1.15)) + 
# #         theme_bw() + 
# #         theme_classic() + 
# #         scale_fill_manual(
# #             values = c("#2166ac", "#b2182b"),
# #             labels = c("All", "Top 5%")
# #         ) + 
# #         theme(
# #             text = element_text(size = 15),
# #             legend.position = "none"
# #         ) + 
# #         labs(x = "", y = "")

# #     if(unique_id != unique_keys[1]){
# #         p1 <- p1 +
# #             theme(
# #                 axis.title.y = element_blank(),
# #                 axis.text.y = element_blank(),
# #                 axis.ticks.y = element_blank(),
# #                 axis.line.y = element_blank()
# #             )
# #     }

# #     if(unique_id == unique_keys[1]){
# #         p1 <- p1 + 
# #             labs(x = "", y = "Fraction, %")
# #     }

# #     return(p1)
# # })

# # main_title <- grid::textGrob(
# #     "ClinVar SV fragility prediction", 
# #     gp = grid::gpar(fontsize = 14)
# # )
# # sub_title <- grid::textGrob(
# #     paste0(
# #         "Fraction of all (left) and top 5% highest (right) delta ",
# #         "fragility within the high fragile zones"
# #     ),
# #     gp = grid::gpar(fontsize = 11)
# # )

# # combined_plots <- do.call(gridExtra::grid.arrange, c(plots, nrow = 1))
# # pdf(
# #     file = paste0(
# #         "../figures/deltafragility/All_and_Top5Per_Delta_",
# #         "bw_", bw, ".pdf"
# #     ),
# #     height = 5, width = 12
# # )
# # gridExtra::grid.arrange(
# #     main_title, sub_title, 
# #     combined_plots,
# #     heights = c(0.9, 0.5, 12)
# # )
# # plot.saved <- dev.off()

# ########################################################################################
# #' Focusing on the biggest delta, is there a pattern with genomic features?
# get_genic_feat <- function(){
#     # import all genomic features for analysis
#     files_to_load <- list.files(
#         path = "../../05_Cosmic/data/annotations",
#         pattern = "group_*",
#         full.names = TRUE
#     )
#     group_names <- stringr::str_extract(
#         string = files_to_load,
#         pattern = "(?<=group_)[^.]+(?=\\.csv)"
#     )

#     all_groups <- lapply(files_to_load, fread, showProgress = FALSE)
#     names(all_groups) <- group_names

#     return(list(all_groups, group_names))
# }

# df_genic <- get_genic_feat()
# group_names <- df_genic[[2]]
# df_genic <- df_genic[[1]]
# groups_combined <- rbindlist(df_genic, fill = TRUE)
# groups_combined <- groups_combined[, gene_id := NULL]
# groups_split <- split(groups_combined, by = "type")
# group_map <- lapply(1:length(df_genic), function(ind){
#     return(data.table(
#         type = unique(df_genic[[ind]]$type),
#         group = names(df_genic)[[ind]]
#     ))
# })
# group_map <- rbindlist(group_map)
# granges_list <- lapply(groups_split, plyranges::as_granges)

# calc_prob_per_base <- function(df_breaks, df_genic_feat){
#     overlaps <- GenomicRanges::findOverlaps(df_genic_feat, df_breaks)

#     if(length(overlaps) == 0){
#         overlap_counts <- 0 
#     } else {
#         overlap_counts <- sum(table(queryHits(overlaps)))
#     }
#     return(overlap_counts)
# }

# ppb <- pbapply::pbsapply(1:length(granges_list), function(ind){
#     df_genic_feat <- granges_list[[ind]]
#     genic_feat_len <- sum(width(reduce(df_genic_feat)))

#     # predicted breaks
#     res_pred <- lapply(1:length(col_names), function(x){
#         df_top_5_perc <- as_tibble(all_df_delta_liftover[[col_names[x]]]) %>%
#             dplyr::group_by(ClinicalSignificance) %>%
#             dplyr::slice_max(order_by = delta, prop = 0.05) %>%
#             dplyr::ungroup()

#         df_bot_5_perc <- as_tibble(all_df_delta_liftover[[col_names[x]]]) %>%
#             dplyr::group_by(ClinicalSignificance) %>%
#             dplyr::slice_min(order_by = delta, prop = 0.05) %>%
#             dplyr::ungroup()

#         temp_top_sv_granges <- dplyr::bind_rows(df_top_5_perc, df_bot_5_perc) %>% 
#             plyranges::as_granges()

#         res_pred <- calc_prob_per_base(
#             df_breaks = temp_top_sv_granges, 
#             df_genic_feat = df_genic_feat
#         )
#         res_pred <- res_pred / genic_feat_len
#         return(res_pred)
#     })
#     names(res_pred) <- col_names
#     return(res_pred)
# })

# ppb_t <- t(ppb)
# ppb_df <- as_tibble(ppb_t) %>% 
#     dplyr::mutate(across(where(is.list), as.numeric)) %>% 
#     as.data.frame()
# rownames(ppb_df) <- names(groups_split)

# # Normalise between 0 and 1 to generate a relative importance plot
# ppb_df <- apply(ppb_df, 2, function(x) x / max(x))

# hex <- c("#990000", "#fd952c", "#669e85","#394d7e", "#8A2BE2")
# hex <- setNames(hex, group_names)

# df_ppb <- as_tibble(ppb_df) %>% 
#     dplyr::mutate(type = names(groups_split)) %>% 
#     dplyr::arrange(desc(ppb_df)) %>% 
#     dplyr::left_join(group_map, by = "type") %>% 
#     dplyr::mutate(
#         type = forcats::fct_inorder(type),
#         group = as.factor(group),
#         hex = hex[group]
#     )

# plots_ppb <- lapply(col_names, function(column_name){
#     temp <- df_ppb[c(column_name, "type", "group", "hex")] %>% 
#         dplyr::rename_with(~c("ppb", "type", "group", "hex")) %>% 
#         dplyr::arrange(ppb) %>% 
#         dplyr::mutate(type = forcats::fct_inorder(type))

#     p1 <- temp %>% 
#         ggplot2::ggplot(aes(x = ppb, y = type)) +
#         ggplot2::geom_segment(aes(
#             x = 0, xend = ppb, 
#             y = type, yend = type),
#             linewidth = 1.1,
#             col = 'grey80'
#         ) +
#         ggplot2::geom_point(size = 2, col = 'black') +
#         ggplot2::scale_color_identity() + 
#         ggplot2::theme_bw() + 
#         ggplot2::theme_classic() +
#         suppressWarnings(ggplot2::theme(
#             text = element_text(size = 15),
#             axis.text.y = element_text(colour = temp$hex)
#         )) +
#         ggplot2::coord_cartesian(xlim = c(0, NA)) + 
#         ggplot2::labs(
#             subtitle = plot_titles[[column_name]],
#             x = "", y = ""
#         )

#     return(p1)
# })

# pdf(
#     file = paste0(
#         "../figures/deltafragility/",
#         "Top5Per_Delta_Overlap_with_Genic_feats.pdf"
#     ),
#     height = 9, width = 20
# )
# do.call(gridExtra::grid.arrange, c(plots_ppb, nrow = 1))
# plots.saved <- dev.off()