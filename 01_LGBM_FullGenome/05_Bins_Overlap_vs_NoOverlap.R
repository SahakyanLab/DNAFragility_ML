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

cat(paste0("Processing different overlap methods...\n"))

Rcpp::sourceCpp('./lib/countIntervals.cpp')

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

# import one of the predicted breaks as reference
df_pred <- arrow::read_parquet(paste0(
    path_to_pred, "/", col_names[1], ".parquet"
))

# bin each chromosome based on length of each bin size
df_bw <- tibble(
    bw = as.integer(c(
        1960, 2500, 5000, 7500, 
        10000, 15000, 20000, 25000, 
        30000, 40000, 45000, 50000
    ))
)

path_to_save <- "Overlapping_NonOverlappingBreaks"
dir.create(path = paste0("../figures/", path_to_save), showWarnings = FALSE)
dir.create(path = paste0("../data/", path_to_save), showWarnings = FALSE)

#####################################################################################
#' Compare binning the genome into non-overlapping and overlapping regions. 
#####################################################################################
df_break_bins <- pbapply::pblapply(1:nrow(df_bw), function(row){
    # comparing non-overlapping vs. overlapping process
    df_break_bins <- lapply(1:22, function(x){
        temp <- refseq.table[chr == paste0("chr", x)]

        # create bins
        bins <- seq(1, temp$end, by = df_bw$bw[row])
        
        # final bin to include last bases
        if(tail(bins, n = 1) < temp$end) bins <- c(bins, temp$end)

        # Bin data by non-overlapping bins
        df_bps <- df_pred %>% 
            dplyr::filter(seqnames == paste0("chr", x)) %>% 
            dplyr::distinct()

        # summary stats
        t1 <- Sys.time()
        NonOverlapping_bins <- df_bps %>%
            dplyr::mutate(bins = cut(start, breaks = bins)) %>% 
            dplyr::group_by(bins) %>% 
            dplyr::summarise(count = dplyr::n()) %>%
            dplyr::arrange(bins) %>% 
            dplyr::mutate(seqnames = paste0("chr", x)) %>% 
            dplyr::filter(!is.na(bins)) %>% 
            dplyr::select(count) %>% 
            dplyr::rename_with(~c("ID")) %>% 
            dplyr::group_by(ID) %>%
            dplyr::summarise(Total = dplyr::n()) %>%
            dplyr::ungroup() %>% 
            dplyr::mutate(type = "NonOverlapping_bins")
        
        total.time <- Sys.time() - t1
        time_diff_nooverlaps <- signif(total.time[[1]])
        unit_nooverlaps <- attr(total.time, "units")

        NonOverlapping_bins <- NonOverlapping_bins %>% 
            dplyr::mutate(
                time = time_diff_nooverlaps, 
                units = unit_nooverlaps
            )

        # Bin data by rolling window bins
        t1 <- Sys.time()
        bins_start <- 1:(temp$end-df_bw$bw[row])
        bins_end <- bins_start+df_bw$bw[row]
        df_bins <- data.table(
            start = bins_start, 
            end = bins_end
        )

        df_break_counts <- countIntervalForBreaks(
            bins_start = df_bins$start,
            bins_end = df_bins$end,
            start_bp = df_bps$start
        )
        df_bins[, Count := df_break_counts]
        df_bins <- df_bins[Count != 0]

        # summary of breaks
        Overlapping_bins <- df_bins$Count
        Overlapping_bins <- Overlapping_bins[Overlapping_bins != 0]
        Overlapping_bins <- as.data.table(table(Overlapping_bins))
        setnames(Overlapping_bins, c("ID", "Total"))
        Overlapping_bins[, `:=`(ID = as.integer(ID), type = "Overlapping_bins")]
        
        total.time <- Sys.time() - t1
        time_diff_overlaps <- signif(total.time[[1]])
        unit_overlaps <- attr(total.time, "units")

        Overlapping_bins[, `:=`(
            time = time_diff_overlaps, 
            units = unit_overlaps
        )]

        # combine stats
        df_break_bins <- rbind(Overlapping_bins, NonOverlapping_bins)
        df_break_bins[, seqnames := paste0("chr", x)]

        # non_overlapping_breaks
        gc(verbose = FALSE)
        
        return(df_break_bins)
    })
    df_break_bins <- rbindlist(df_break_bins)

    time_diff <- as_tibble(df_break_bins) %>% 
        dplyr::select(type, time, units) %>% 
        dplyr::distinct() %>% 
        dplyr::group_by(type) %>% 
        dplyr::summarise(time = sum(time)) %>% 
        dplyr::ungroup() %>% 
        dplyr::arrange(desc(type)) %>% 
        dplyr::mutate(bw = df_bw$bw[row])

    fwrite(
        as.data.table(time_diff),
        paste0(
            "../data/", path_to_save, 
            "/TimeDiff_bw_", df_bw$bw[row], ".csv"
        )
    )

    df_break_bins <- as_tibble(df_break_bins) %>% 
        dplyr::group_by(type, ID) %>%  
        dplyr::summarise(Count = sum(Total)) %>% 
        dplyr::mutate(Total_norm = Count/sum(Count)) %>% 
        dplyr::ungroup() %>% 
        dplyr::filter(ID > 0) %>% 
        dplyr::mutate(
            method = as.factor(type),
            hex = dplyr::case_when(
                method == "NonOverlapping_bins" ~ "#BB357E",
                TRUE ~ "#4592C8"
            ),
            bw = df_bw$bw[row]
        ) %>% 
        dplyr::ungroup() %>% 
        suppressMessages()

    fwrite(
        df_break_bins,
        paste0(
            "../data/",  path_to_save, 
            "/BreakBins_bw_", df_bw$bw[row], ".csv"
        )
    )

    # # plot results
    # p1 <- df_break_bins %>% 
    #     ggplot2::ggplot(aes(x = ID, y = Total_norm, col = method)) + 
    #     ggplot2::geom_point(alpha = 1) +
    #     ggplot2::geom_line(linewidth = 1.1, alpha = 0.8) + 
    #     ggplot2::scale_x_continuous(minor_breaks = seq(1, df_bw$x_range[row], 1)) + 
    #     ggplot2::coord_cartesian(
    #         xlim = c(1, df_bw$x_range[row]),
    #         ylim = c(0, NA)
    #     ) + 
    #     # ggplot2::scale_color_identity() + 
    #     ggplot2::scale_color_manual(values = c(
    #         "NonOverlapping_bins" = "#BB357E", 
    #         "Overlapping_bins" = "#4592C8"
    #     )) +
    #     ggplot2::theme_bw() + 
    #     ggplot2::theme_classic() +
    #     ggplot2::theme(text = element_text(size = 15)) +
    #     ggplot2::labs(
    #         subtitle = paste0(
    #             "Count of COSMIC breaks per bin of size: ", df_bw$bw[row]
    #         ),
    #         x = "Number of breaks in window",
    #         y = "Normalised count"
    #     )

    # ggsave(
    #     filename = paste0(
    #         "../figures/",  path_to_save, "/", 
    #         "bw_", df_bw$bw[row], ".pdf"
    #     ),
    #     plot = p1,
    #     height = 7, width = 10
    # )

    gc(verbose = FALSE)

    return(df_break_bins)
})
df_break_bins <- do.call(rbind, df_break_bins)

# smoothen adjacent two bins to reduce noise in plot
df_break_bins_smooth <- as_tibble(df_break_bins) %>% 
    dplyr::group_by(method, bw) %>% 
    dplyr::mutate(group = (dplyr::row_number()+1) %/% 2) %>% 
    dplyr::group_by(method, group, bw) %>% 
    dplyr::mutate(Mean_count = mean(Count)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(-group) %>% 
    dplyr::group_by(method, bw) %>% 
    dplyr::mutate(Total_norm = Mean_count / sum(Mean_count)) %>% 
    dplyr::arrange(method, bw) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        hex = ifelse(
            method == "NonOverlapping_bins",
            "#b2182b", "#2166ac"
        )
    )

fwrite(
    as.data.table(df_break_bins_smooth),
    paste0(
        "../data/",  path_to_save, "/All.csv"
    )
)

p1 <- df_break_bins_smooth %>% 
    ggplot2::ggplot(aes(x = ID, y = Total_norm, col = method)) + 
    ggplot2::geom_line(linewidth = 1, alpha = 0.5) + 
    ggplot2::facet_wrap(vars(bw), ncol = 4, scales = "free_y") + 
    ggplot2::coord_cartesian(
        xlim = c(1, 1500),
        ylim = c(0, NA)
    ) + 
    ggplot2::scale_color_manual(values = c(
        "NonOverlapping_bins" = "#BB357E", 
        "Overlapping_bins" = "#4592C8"
    )) +
    ggplot2::theme_bw() + 
    ggplot2::theme_classic() +
    ggplot2::theme(text = element_text(size = 15)) +
    ggplot2::labs(
        subtitle = "Most conservative fragility predictions",
        x = "Number of breaks in window",
        y = "Normalised count"
    )

ggsave(
    filename = paste0(
        "../figures/",  path_to_save, 
        "/Overlapping_and_NonOverlappingBins.pdf"
    ),
    plot = p1,
    height = 10, width = 15
)

# time difference between processing non-overlapping and overlapping bins.
time_diff <- sapply(1:nrow(df_bw), function(row){
    time_diff <- fread(
        paste0(
            "../data/", path_to_save, 
            "/TimeDiff_bw_", df_bw$bw[row], ".csv"
        )
    )
    return(time_diff$time[1] / time_diff$time[2])
})
mean(time_diff, na.rm = TRUE)