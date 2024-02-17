# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))
suppressPackageStartupMessages(suppressWarnings(library(arrow)))
suppressPackageStartupMessages(suppressWarnings(library(rtracklayer)))
suppressPackageStartupMessages(suppressWarnings(library(moments)))
suppressPackageStartupMessages(suppressWarnings(library(minpack.lm)))
pbapply::pboptions(char = "=", type = "txt")

# source functions
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

cat(paste0("Processing the optimal bin size...\n"))

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

#' for biggest bin-width calculated (90,000), the majority of the 
#' distribution seems to be within 200 breaks. Hence, that will be 
#' the cut-off used for the skewness calculation.
#' 
#' 1. cut into bins
#' 2. filter for number of breaks between 0-20
#' 3. perform normality stat test
# count cosmic breaks in non-overlapping windows and return granges 

all_bws <- seq(500, 100000, by = 1000)
df_skewness <- pbapply::pblapply(all_bws, function(bw){
    df_bp_granges <- lapply(1:22, function(chr){
        return(get_breaks_in_bins(
            break_table = df_pred, 
            bin_width = bw, 
            chr_len_table = refseq.table, 
            filter_for_chr = paste0("chr", chr)
        ))
    })
    df_bp_granges_all <- suppressWarnings(unlist(as(df_bp_granges, "GRangesList")))
    
    max_bw_val <- as_tibble(df_bp_granges_all$Breaks) %>% 
        dplyr::group_by(value) %>% 
        dplyr::summarise(count = dplyr::n()) %>% 
        dplyr::mutate(
            frac = count/sum(count)*100,
            cum_frac = cumsum(frac)
        ) %>% 
        dplyr::filter(cum_frac <= 95) %>%
        tail(n = 1) %>% 
        dplyr::pull(value)
    
    df_bp_bw <- df_bp_granges_all$Breaks
    df_bp_bw <- df_bp_bw[df_bp_bw <= max_bw_val]

    return(list(
        skew = moments::skewness(df_bp_bw),
        kurt = moments::kurtosis(df_bp_bw),
        avg = mean(df_bp_bw),
        stdev = sd(df_bp_bw)
    ))
})

df_plateau <- tibble(
        skew  = sapply(df_skewness, `[[`, 1),
        kurt  = sapply(df_skewness, `[[`, 2),
        avg   = sapply(df_skewness, `[[`, 3),
        stdev = sapply(df_skewness, `[[`, 4), 
        bw = all_bws
    ) %>% 
    dplyr::filter(bw <= 100000) %>%
    dplyr::mutate(delta_param = 1 / (lead(avg) - avg)) %>% 
    tidyr::drop_na()

# plot results
# Define parameter grid
param_grid <- expand.grid(
    a = seq(0.001, 1, by = 0.05), 
    b = seq(0.001, 1, by = 0.05),
    c = seq(0.001, 1, by = 0.05)
)

df_plateau <- df_plateau %>% 
    dplyr::do({            
        # iterate through the parameter grid
        best_rss <- Inf
        results <- list()

        for(i in 1:nrow(param_grid)){
            params <- param_grid[i, ]

            # exponential decay
            fit_exp <- tryCatch({
                minpack.lm::nlsLM(
                    delta_param ~ a * exp(-b * bw) + c,
                    data = .,
                    start = params,
                    control = nls.lm.control(maxiter = 100),
                    trace = FALSE
                )
            }, error = function(e) NULL) %>% 
            suppressWarnings()
            
            if(!is.null(fit_exp)){
                rss_exp <- sqrt(sum(residuals(fit_exp)^2)/length(residuals(fit_exp)))
                rss_exp <- ifelse(is.na(rss_exp), Inf, rss_exp)

                # only accept improved model
                if(rss_exp < best_rss){
                    best_rss <- rss_exp
                    results[[i]] <- list(
                        params = params, 
                        rss = rss_exp, 
                        model = fit_exp
                    )
                }
            }
        }

        # extract RSS values and find best model
        results[sapply(results, is.null)] <- NULL
        rss_values <- sapply(results, function(r) r$rss)
        best_result <- results[[which.min(rss_values)]]

        # use best model for predictions
        predicted_values <- predict(
            best_result$model, 
            newdata = data.frame(bw = .$bw)
        )
        tibble(., model_fit = predicted_values)
    }) %>% 
    dplyr::mutate(
        model_fit_lower = signif(model_fit, 4), 
        model_fit_upper = signif(model_fit, 8), 
        lower_bound = lead(model_fit_lower) - model_fit_lower,
        upper_bound = lead(model_fit_upper) - model_fit_upper
    ) %>% 
    tidyr::drop_na()

min_cutoff <- df_plateau %>% 
    dplyr::filter(lower_bound == 0) %>% 
    head(1) %>% 
    dplyr::pull(bw)

max_cutoff <- df_plateau %>% 
    dplyr::filter(upper_bound == 0) %>% 
    head(1) %>% 
    dplyr::pull(bw)

plot_plateau <- df_plateau %>% 
    ggplot2::ggplot() +
    ggplot2::geom_rect(aes(
        xmin = min_cutoff, xmax = max_cutoff, 
        ymin = -Inf, ymax = Inf
        ), 
        alpha = 0.1, 
        fill = "#ffc4c4"
    ) + 
    ggplot2::geom_line(
        aes(x = bw, y = delta_param), 
        alpha = 0.5,
        col = "black"
    ) + 
    ggplot2::geom_line(
        aes(x = bw, y = model_fit), 
        linewidth = 1.5, 
        col = "darkred"
    ) + 
    ggplot2::theme_bw() + 
    ggplot2::theme_classic() + 
    ggplot2::theme(text = element_text(size = 15)) + 
    ggplot2::coord_cartesian(
        xlim = c(0, 100000)
    ) + 
    ggplot2::labs(
        title = "Rate of change of average number of breaks.",
        subtitle = paste0(
            "Convergence between ", min_cutoff, " ",
            "and ", max_cutoff, "."
        ),
        x = "Bin width",
        y = "Rate of change"
    )

path_to_save <- "Overlapping_NonOverlappingBreaks"
ggsave(
    filename = paste0(
        "../figures/",
        path_to_save,
        "/AvgBreaks_RateOfChange_Cutoff.pdf"
    ),
    plot = plot_plateau,
    height = 6, width = 8
)

fwrite(
    as.data.table(df_plateau),
    paste0(
        "../data/",
        path_to_save,
        "/AvgBreaks_RateOfChange_Cutoff.csv"
    )
)