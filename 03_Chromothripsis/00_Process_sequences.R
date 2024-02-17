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
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))
pbapply::pboptions(char = "=", type = "txt")

# source functions
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)
only_breaks <- TRUE

cat(paste0("Processing chromothripsis breakpoints...\n"))

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
    chr = rownames(refseq.table),
    end = refseq.table$seqlengths
)
human_genome_len <- sum(refseq.table$end)

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

# import predicted breakpoints of the full human genome
path_to_pred <- "../data/Whole_genome_pred"

# import chromothripsis breakpoint junctions from ChromothripsisDB
# paper: https://doi.org/10.1093/bioinformatics/btv757
# download: http://cailab.labshare.cn/ChromothripsisDB/download.html
casedata <- fread('./data/CaseData.txt')

if(only_breaks){
    casedata <- casedata[Connection != ""]
} else {
    casedata <- casedata[CNA != "" & Connection != ""]
}

all_breaks <- pbapply::pblapply(1:nrow(casedata), function(x){    
    # extract the copy number alteration coordinates
    assembly <- casedata$GenomeAssembly[x]
    label <- gsub(" ", "_", casedata$Disease[x])
    cna <- casedata$CNA[x]
    cna <- stringr::str_split(string = cna, pattern = ";")[[1]]
    bp_connections <- casedata$Connection[x]
    bp_connections <- stringr::str_split(string = bp_connections, pattern = ";")[[1]]

    cna <- stringr::str_split(string = cna, pattern = ":")
    cna <- cna[which(grepl("^chr([1-9]|1[0-9]|2[0-2])$", sapply(cna, `[[`, 1)))]
    extracted_cna <- stringr::str_split(string = sapply(cna, `[[`, 2), pattern = "-")
    start_bp <- sapply(extracted_cna, `[`, 1) %>% as.integer()
    end_bp <- sapply(extracted_cna, `[`, 2) %>% as.integer()
    chrs <- as.integer(gsub("chr", "", sapply(cna, `[`, 1)))
    copy <- sapply(cna, `[`, 3) %>% as.integer()
    
    if(!only_breaks){
        dt_all <- dt_cna <- tibble(
                seqnames = chrs, 
                start = start_bp,
                end = end_bp,
                copy = copy
            ) %>% 
            tidyr::pivot_longer(
                cols = -c(seqnames, copy),
                names_to = "Key", 
                values_to = "start"
            ) %>% 
            dplyr::select(seqnames, start, copy) %>%
            as.data.table()
    }

    # extract breakpoint junctions
    if(any(bp_connections != "")){
        bps <- stringr::str_split(string = bp_connections, pattern = ",")
        extract_unique_numbers <- function(x){
            chromosome_match <- regexpr("hs(\\d+):", x)
            chromosome_number <- ifelse(chromosome_match > -1, 
                as.numeric(substring(x, chromosome_match + 2, regexpr(":", x) - 1)), 
                NA
            )

            number_matches <- regmatches(x, gregexpr("(?<=:)[0-9]+|(?<=-)[0-9]+", x, perl=TRUE))
            numbers <- unique(as.numeric(unlist(number_matches)))

            list(chromosome_number, numbers)
        }
        extracted_bps <- lapply(bps, function(x){
            c(extract_unique_numbers(x[1]), extract_unique_numbers(x[2]))
        })
        extracted_bps <- extracted_bps[which(!is.na(sapply(extracted_bps, `[[`, 1)))]

        expand_entry <- function(entry){
            entry <- entry[lengths(entry) > 0]

            # seqnames and start
            starts <- entry[[2]]
            seqnames_start <- rep(entry[[1]], length(start))
            dt_all <- dt_start <- data.table(
                seqnames = seqnames_start,
                start = starts,
                copy = NA
            )

            # seqnames and end
            if(length(entry) == 4){
                ends <- entry[[4]]
                seqnames_end <- rep(entry[[1]], length(ends))

                dt_end <- data.table(
                    seqnames = seqnames_end,
                    start = ends,
                    copy = NA
                )
                dt_all <- rbind(dt_start, dt_end)
            }

            return(dt_all)
        }

        dt_bps <- rbindlist(lapply(extracted_bps, expand_entry))

        if(only_breaks){
            dt_all <- dt_bps
        } else {
            dt_all <- rbind(dt_bps, dt_cna)
        }
        dt_all[, copy := NULL]
    }

    # file format to run analysis on.
    df_bps <- as_tibble(dt_all) %>% 
        dplyr::arrange(seqnames, start) %>%
        dplyr::mutate(seqnames = paste0("chr", seqnames)) %>%
        dplyr::mutate(width = 1) %>%
        plyranges::as_granges()

    base_path <- "../../05_Cosmic/data/liftover"
    if(grepl("hg17", assembly)){
        # hg17 -> hg38
        liftover_chain <- import.chain(paste0(base_path, "/hg17ToHg38.over.chain"))
        df_bps <- suppressMessages(liftOver(df_bps, liftover_chain))
        df_bps <- unlist(as(df_bps, "GRangesList"))
    } else if(grepl("hg18", assembly)){
        # hg18 -> hg38
        liftover_chain <- import.chain(paste0(base_path, "/hg18ToHg38.over.chain"))
        df_bps <- suppressMessages(liftOver(df_bps, liftover_chain))
        df_bps <- unlist(as(df_bps, "GRangesList"))
    } else if(grepl("hg19", assembly)){
        # hg19 -> hg38
        liftover_chain <- import.chain(paste0(base_path, "/hg19ToHg38.over.chain"))
        df_bps <- suppressMessages(liftOver(df_bps, liftover_chain))
        df_bps <- unlist(as(df_bps, "GRangesList"))
    } 
    # hg38 -> t2t version
    liftover_chain <- import.chain(paste0(base_path, "/hg38-chm13v2.over.chain"))
    df_bps <- suppressMessages(liftOver(df_bps, liftover_chain))
    df_bps <- unlist(as(df_bps, "GRangesList"))
    df_bps <- unique(df_bps)
    seqlevels(df_bps, pruning.mode = "coarse") <- paste0("chr", 1:22)
    df_bps <- df_bps %>% dplyr::mutate(label = label)
    
    return(df_bps)
})
dt_all_breaks <- do.call(plyranges::bind_ranges, all_breaks)

hex_codes <- c("#000000", "#990000", "#fd952c", "#669e85", "#394d7e")
names(hex_codes) <- c("True", col_names)

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

map_group_to_predictions <- function(data, x, bw){
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
                .default = "Mid"
            ),
            group = as.factor(group)
        ) %>% 
        dplyr::select(-ID) 

    df_bp_granges_all <- df_bp_granges_all %>% 
        dplyr::mutate(group = df_breaks$group)

    # Where are the chromothripsis events taking place?
    overlaps <- GenomicRanges::findOverlaps(data, df_bp_granges_all)
    overlaps <- data[queryHits(overlaps)] %>%
        dplyr::mutate(group = mcols(df_bp_granges_all)$group[subjectHits(overlaps)])

    dt_overlaps <- data.table(
        group = mcols(overlaps)$group,
        label = mcols(overlaps)$label
    )
    dt_overlaps[, `:=`(key = col_names[x], bw = bw)]
    return(dt_overlaps)
}

df_pred_all <- pbapply::pblapply(col_names, function(x){
    return(arrow::read_parquet(paste0("../data/Whole_genome_pred/", x, ".parquet")))
})
names(df_pred_all) <- col_names

#' true: predicted overlaps
bws <- c(1960, 10000, 20000)
df_all_categories <- lapply(bws, function(bw){
    pred_categories <- pbapply::pblapply(1:length(col_names), function(x){
        return(map_group_to_predictions(
            data = dt_all_breaks,
            x = x, bw = bw
        ))
    })
    return(rbindlist(pred_categories))
})
df_all_cat_all <- rbindlist(df_all_categories)

#' 2. Proportion of low/medium/high fragile regions within chromothripsis sites
dir.create(path = "../figures/Chromothripsis/", showWarnings = FALSE)

group_cols <- c(
    # green - yellow - red
    "#006400", "#ff9e27", "#b2182b",

    # green - yellow - red
    "#006400", "#ff9e27", "#b2182b"
)

fwrite(
    df_all_cat_all,
    "./data/overlap_by_groups.csv"
)

for(bw in bws){
    df_all_cat <- as_tibble(df_all_cat_all) %>% 
        dplyr::filter(bw == bw) %>% 
        dplyr::group_by(key, label, group) %>%
        dplyr::summarise(count = dplyr::n()) %>%
        dplyr::group_by(key, label) %>%
        dplyr::mutate(frac = count / sum(count) * 100) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(key, label, group) %>%
        suppressMessages()

    df_all_cat_split <- dplyr::group_split(df_all_cat, key)
    names(df_all_cat_split) <- col_names

    # plot stacked barplot, ranked by decreasing order of H group.
    stacked_category_plots <- lapply(1:length(df_all_cat_split), function(x){
        temp_df <- df_all_cat_split[[x]] %>% 
            dplyr::group_by(label) %>%
            dplyr::mutate(frac = frac / sum(frac)) %>%
            dplyr::ungroup()

        # add high group if any missing. 
        labels_with_high <- temp_df %>%
            dplyr::filter(group == "High") %>%
            dplyr::distinct(label)
        
        missing_high_labels <- setdiff(unique(temp_df$label), labels_with_high$label)

        if(length(missing_high_labels) > 0){
            missing_high_rows <- temp_df[match(missing_high_labels, temp_df$label),]
            missing_high_rows <- missing_high_rows %>% 
                dplyr::mutate(
                    group = "High",
                    count = 0,
                    frac = 0,
                    hex = "#006400"
                )

            temp_df <- dplyr::bind_rows(temp_df, missing_high_rows)
        }

        ranked_unique_labels <- temp_df %>% 
            dplyr::filter(group == "High") %>%
            dplyr::arrange(frac) %>%
            dplyr::pull(label) %>%
            unique()

        plot_title <- stringr::str_replace(
            string = col_names[[x]], 
            pattern = "^(Pred)_0_(.*)$", 
            replacement = "\\1: 0.\\2 (FPR)"
        )

        temp_df <- temp_df %>%
            dplyr::mutate(
                hex = dplyr::case_when(
                    group == "High" ~ "#b2182b",
                    group == "Mid" ~ "#ff9e27",
                    group == "Low" ~ "#006400"
                ),
                label = factor(label, levels = ranked_unique_labels),
                group = factor(group, levels = c("High", "Mid", "Low"))
            ) %>% 
            dplyr::arrange(label, group) 

        p1 <- temp_df %>%
            ggplot(aes(x = frac, y = label, fill = group)) + 
            geom_bar(position = "stack", stat = "identity") +
            theme_bw() + 
            theme_classic() + 
            scale_fill_manual(values = unique(temp_df$hex)) + 
            suppressWarnings(theme(
                text = element_text(size = 15),
                legend.position = "none"
            )) +
            labs(subtitle = plot_title, x = "", y = "")

        return(p1)
    })

    pdf(
        file = paste0(
            "../figures/Chromothripsis/",
            ifelse(only_breaks,
                "Chromothripsis_by_category_onlybreaks",
                "Chromothripsis_by_category"
            ),
            "_bw_", bw, ".pdf"
        ),
        height = 10, width = 24
    )
    do.call(gridExtra::grid.arrange, c(stacked_category_plots, nrow = 1))
    plot.saved <- dev.off()

    #' control: randomly sample positions of the genome and calculate the overlap in each group
    df_random_frac <- fread(paste0(path_to_pred, "/random_sampling/random_frac.csv"))
    group_frac_random <- as_tibble(df_random_frac) %>% 
        dplyr::filter(bin_width == bw) %>% 
        dplyr::select(key = threshold, group, count, frac) %>% 
        dplyr::group_by(key) %>% 
        dplyr::mutate(total = sum(count), .before = frac) %>% 
        dplyr::ungroup()

    #' 3. Distribution of low/mid/high fragile region proportions for true and control samples
    get_group_frac <- function(data){
        return(
            as_tibble(df_all_cat) %>% 
                dplyr::select(group, key) %>%
                dplyr::group_by(key, group) %>%
                dplyr::summarise(count = dplyr::n()) %>%
                dplyr::group_by(key) %>%
                dplyr::mutate(
                    total = sum(count),
                    frac = count / total * 100
                ) %>% 
                dplyr::ungroup() %>%
                dplyr::mutate(group = factor(group, levels = c("Low", "Mid", "High"))) %>%
                dplyr::arrange(key, group) %>%
                suppressMessages()
        )
    }

    group_frac_pred <- as_tibble(df_all_cat) %>% 
        dplyr::select(-label) %>%
        dplyr::group_by(key, group) %>%
        dplyr::summarise(count = sum(count)) %>%
        dplyr::group_by(key) %>%
        dplyr::mutate(
            total = sum(count),
            frac = count / total * 100
        ) %>% 
        dplyr::ungroup() %>%
        dplyr::mutate(group = factor(group, levels = c("Low", "Mid", "High"))) %>%
        dplyr::arrange(key, group) %>%
        suppressMessages()

    totals <- paste0(
        "Total number of entries. ",
        "Model: ", group_frac_pred$total[1], ". ",
        "Rand. sample: ", group_frac_random$total[1], ". "
    )

    group_frac_all <- dplyr::bind_rows(
            group_frac_random %>% 
                dplyr::mutate(type = "rand"),
            group_frac_pred %>% 
                dplyr::mutate(type = "pred")
        ) %>% 
        dplyr::mutate(
            group = factor(group, levels = c("Low", "Mid", "High")),
            hex = dplyr::case_when(
                group == "High" ~ "#b2182b",
                group == "Mid" ~ "#006400",
                group == "Low" ~ "#ff9e27"
            ),
            combined = interaction(group, type)
        ) %>% 
        dplyr::arrange(key, group)

    fwrite(
        group_frac_all,
        paste0(
            "./data/fractional_values_by_group_",
            "bw_", bw, ".csv"
        )
    )

    group_frac_all_split <- dplyr::group_split(group_frac_all, key)
    names(group_frac_all_split) <- col_names

    # plot stacked barplot, ranked by decreasing order of H group.
    plots <- lapply(1:length(group_frac_all_split), function(x){
        # z-test of difference in proportions
        res_signif <- group_frac_all_split[[x]] %>% 
            split(.$group) %>%
            purrr::map_dfr(~{
                # Extract 'rand' and 'pred' data
                rand_data <- .x %>% filter(type == "rand")
                pred_data <- .x %>% filter(type == "pred")

                # Perform z-test of difference in proportions
                prop_test <- prop.test(
                    c(rand_data$count, pred_data$count), 
                    c(rand_data$total, pred_data$total)
                )

                # Create a summary tibble
                tibble(
                    key = .x$key,
                    group = unique(.x$group),
                    hex = unique(.x$hex),
                    p_value = prop_test$p.value,
                    estimate_rand = prop_test$estimate[1],
                    estimate_pred = prop_test$estimate[2]
                )
            }) %>% 
            dplyr::mutate(
                key = plot_titles[key],
                key = forcats::fct_inorder(key),
                Signif_level = dplyr::case_when(
                    p_value < 0.001 ~ "***",
                    p_value < 0.01 ~ "**",
                    p_value < 0.05 ~ "*",
                    TRUE ~ "N.S."
                ),
                y_position = ifelse(
                    estimate_rand < estimate_pred,
                    estimate_pred,
                    estimate_rand
                ),
                y_position = y_position + 0.05,
                annotation = Signif_level
            )

        res_signif <- res_signif %>% 
            dplyr::distinct() %>% 
            dplyr::select(-key) %>% 
            tidyr::pivot_longer(
                cols = -c(
                    group, p_value, 
                    Signif_level, y_position, 
                    annotation, hex
                ),
                names_to = "Key",
                values_to = "Value"
            ) %>% 
            dplyr::mutate(
                type = ifelse(grepl("rand", Key), "Rand. Sample", "Model"),
                type = forcats::fct_inorder(type),
                Value = Value * 100,
                combined = interaction(group, type),
                subtitle = plot_titles[[x]]
            )

        fwrite(
            res_signif,
            paste0(
                "./data/prop_test_values_",
                "bw_", bw, ".csv"
            )
        )

        res_annot <- res_signif %>% 
            dplyr::select(-type, -combined, -Key, -Value) %>% 
            dplyr::distinct() %>% 
            dplyr::mutate(
                y_position = y_position * 100,
                xmin = (1:3) - 0.225,
                xmax = (1:3) + 0.225,
            )

        p1 <- res_signif %>%
            ggplot(aes(
                x = group, y = Value, 
                fill = combined,
                pattern = combined
            )) + 
            geom_bar_pattern(
                stat = "identity", 
                position = position_dodge(), 
                color = "black",
                alpha = 0.75,
                pattern_colour = "black",
                pattern_spacing = 0.07,
                pattern_angle = 45
            ) +
            ggsignif::geom_signif(
                y_position = res_annot$y_position,
                xmin = res_annot$xmin,
                xmax = res_annot$xmax,
                annotation = res_annot$annotation,
                textsize = 5,
                vjust = -0.1,
                margin_top = 0.1,
                tip_length = 0.01
            ) +
            coord_cartesian(ylim = c(0, 100)) + 
            facet_wrap(vars(subtitle), nrow = 1) + 
            scale_fill_manual(values = group_cols) +
            scale_pattern_manual(values = c(
                rep('stripe', 3),
                rep('none', 3)
            )) +
            theme_bw() + 
            theme_classic() + 
            theme(
                text = element_text(size = 15),
                legend.position = "none"
            ) +
            labs(x = "", y = "")

        if(x == 1 & bw < 10000){
            ggsave(
                filename = paste0(
                    "../figures/Chromothripsis/",
                    "Chromothripsis_random_vs_pred_bar_sample.pdf"
                ),
                plot = p1,
                height = 5, width = 6
            )
        }

        # if(x > 1){
            p1 <- p1 +
                theme(
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.line.y = element_blank()
                )
        # }

        if(x == 3 & bw == 20000){
            combines_labels <- interaction(res_signif$group, res_signif$type)
            levels(combines_labels) <- as.character(combines_labels)
            res_signif$combined <- combines_labels

            p2 <- res_signif %>%
                # dplyr::mutate(
                #     Value = c(5,5,45,45,50,50)
                # ) %>% 
                ggplot(aes(
                    x = group, y = Value, 
                    fill = combined,
                    pattern = combined
                    # pattern_type = combined
                )) + 
                geom_bar(
                    stat = "identity", 
                    position = position_dodge(), 
                    color = "black",
                    alpha = 0.75
                ) +
                geom_bar_pattern(
                    stat = "identity", 
                    position = position_dodge(), 
                    color = "black",
                    alpha = 0.75,
                    pattern_colour = "white",
                    pattern_spacing = 0.04,
                    pattern_angle = 45
                ) +
                scale_fill_manual(values = group_cols) +
                scale_pattern_manual(values = c(
                    rep('stripe', 3),
                    rep('none', 3)
                )) +
                ggsignif::geom_signif(
                    y_position = c(8,94,10),
                    xmin = res_annot$xmin,
                    xmax = res_annot$xmax,
                    annotation = res_annot$annotation,
                    textsize = 5,
                    vjust = -0.1,
                    tip_length = 0
                ) +
                coord_cartesian(ylim = c(0, 120)) + 
                theme_bw() + 
                theme_classic() + 
                theme(
                    text = element_text(size = 15),
                    plot.subtitle = element_text(size = 12),
                    legend.position = "none"
                ) +
                labs(x = "", y = "Fraction, %")

            ggsave(
                filename = paste0(
                    "../figures/Chromothripsis/",
                    "For_Paper_Chromothripsis_random_vs_pred_bar_",
                    "bw_", bw, ".pdf"
                ),
                plot = p2,
                height = 6, width = 6
            )
        }
        
        return(p1)
    })

    pdf(
        file = paste0(
            "../figures/Chromothripsis/",
            "Chromothripsis_random_vs_pred_bar_",
            "bw_", bw, ".pdf"
        ),
        height = 5, width = 12
    )
    gridExtra::grid.arrange(
        plots[[1]], plots[[2]],
        plots[[3]], plots[[4]],
        nrow = 1
    )
    plot.saved <- dev.off()
}

# ####################################################
# # where the thyroid breaks located? All chromosome 9
# dt_all_breaks %>% 
#     dplyr::filter(grepl("Thyroid", label)) %>% 
#     dplyr::as_tibble() %>% 
#     dplyr::group_by(seqnames) %>% 
#     dplyr::summarise(count = dplyr::n()) %>% 
#     dplyr::ungroup() %>% 
#     dplyr::mutate(frac = count / sum(count) * 100)

# # where do the rest of the chromothripsis cases locate? Varies.
# dt_all_breaks %>% 
#     dplyr::as_tibble() %>% 
#     dplyr::group_by(seqnames) %>% 
#     dplyr::summarise(count = dplyr::n()) %>% 
#     dplyr::ungroup() %>% 
#     dplyr::mutate(frac = count / sum(count) * 100) %>% 
#     dplyr::arrange(desc(frac))

# # Proportion of cancers on chr 9?
# dt_all_breaks %>% 
#     dplyr::as_tibble() %>% 
#     dplyr::filter(seqnames == "chr9") %>% 
#     dplyr::group_by(label) %>% 
#     dplyr::summarise(count = dplyr::n()) %>% 
#     dplyr::ungroup() %>% 
#     dplyr::mutate(frac = count / sum(count) * 100) %>% 
#     dplyr::arrange(desc(frac))

# # Are the breakpoints clustered? NOPE.
# dt_all_breaks %>% 
#     dplyr::filter(grepl("Thyroid", label)) %>% 
#     dplyr::as_tibble() %>% 
#     dplyr::arrange(seqnames, start) %>% 
#     dplyr::group_by(seqnames) %>% 
#     dplyr::mutate(gap = lead(start) - start) %>% 
#     dplyr::ungroup() %>% 
#     dplyr::select(gap) %>% 
#     dplyr::arrange(desc(gap)) %>% 
#     tidyr::drop_na() %>% 
#     dplyr::summarise(
#         mean_gap = mean(gap),
#         sd_gap = sd(gap)
#     ) %>% 
#     dplyr::ungroup()

# # What is the general sequence composition of chromothripsis break points?
# temp_seq <- getSeq(
#     refseq,
#     dt_all_breaks %>% 
#         dplyr::filter(grepl("Thyroid", label)) %>% 
#         dplyr::mutate(start = start - 200, end = end + 200)
# )

# # count base frequencies
# letter.counts <- Biostrings::letterFrequency(
#     temp_seq, letters = "ACGT", OR = 0
# )
# letter.counts.norm <- letter.counts / width(temp_seq)
# gc.content <- letter.counts.norm[, "G"]+letter.counts.norm[, "C"]

# pdf("demo.pdf")
# hist(gc.content, breaks = 50)
# dev.off()

# # count triplet kmers
# kmer.counts <- Biostrings::oligonucleotideFrequency(
#     x = temp_seq,
#     width = 4,
#     step = 1
# )
# kmer.counts <- as_tibble(as.data.frame(kmer.counts))

# kmer.counts %>% 
#     tidyr::pivot_longer(cols = everything()) %>% 
#     dplyr::group_by(name) %>% 
#     dplyr::summarise(count = sum(value)) %>% 
#     dplyr::ungroup() %>% 
#     dplyr::mutate(frac = count / sum(count) * 100) %>% 
#     dplyr::arrange(desc(frac)) %>% 
#     print(n = 30)

# ####################################################