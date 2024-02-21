# import libraries
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
pbapply::pboptions(char = "=")

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

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

# import predictions
path_to_pred <- "../data/Whole_genome_pred"

sequence <- Biostrings::readDNAStringSet(
    filepath = "./data/SNP_nullomer_sequence.fasta.gz",
    format = "fasta"
)
seq_lens <- width(sequence)

# import predictions
pred_files <- list.files(
    path = "../data/SNP_nullomer_sequence_SNPs/SNPs",
    pattern = "final_pred",
    full.names = TRUE,
    recursive = TRUE
)
pred_names <- stringr::str_extract(
    string = pred_files,
    pattern = "(Nullomer|Present)"
)

res_list <- pbapply::pblapply(1:length(pred_names), function(x){
    prediction <- arrow::read_parquet(pred_files[x])
    prediction <- dplyr::select(prediction, -True)
    colnames(prediction) <- gsub("\\.", "_", colnames(prediction))
    
    preds <- prediction %>% 
        dplyr::mutate(Sequence = pred_names[x])

    return(preds)
})
names(res_list) <- pred_names
res <- do.call(dplyr::bind_rows, res_list)

res_summary <- res %>% 
    tidyr::pivot_longer(
        cols = -Sequence,
        names_to = "Key",
        values_to = "Value"
    ) %>%
    dplyr::mutate(
        Key = plot_titles[Key],
        Key = forcats::fct_inorder(Key)
    ) %>% 
    dplyr::group_by(Key, Sequence, Value) %>% 
    dplyr::summarise(count = dplyr::n()) %>% 
    dplyr::mutate(frac = count / sum(count) * 100) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
        Condition = ifelse(Value == 0, "No break", "Break"),
        Condition = factor(Condition, levels = c(
            "Break", "No break"
        )),
        Type = "Change in frag."
    ) %>%
    suppressMessages()

p1 <- res_summary %>% 
    ggplot(aes(x = Sequence, y = frac, fill = Condition)) +
    geom_bar(
        position = "stack",
        stat = "identity",
        col = "black",
        alpha = 0.75
    ) +
    facet_wrap(vars(Key), nrow = 1) +
    geom_text(
        aes(label = sprintf("%.1f%%", frac), y = frac), 
        position = position_stack(vjust = 0.5),
        col = "black",
        size = 5
    ) +
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15)
    ) + 
    scale_fill_manual(values = c("skyblue", "grey80")) + 
    labs(x = "", y = "Fraction, %")

dir.create("../figures/Nullomers/", showWarnings = FALSE)
ggsave(
    filename = "../figures/Nullomers/SNP_True_vs_nullomers.pdf",
    plot = p1,
    height = 5, width = 11
)

# calculate statistical significance of proportions
p_values <- res_summary %>%
    split(.$Key) %>%
    purrr::map_dfr(~{
        d_split <- split(.x, .x$Sequence)
        
        # Create a contingency table
        contingency_table <- matrix(c(
            d_split$Nullomer$count[d_split$Nullomer$Value == 1],
            d_split$Present$count[d_split$Present$Value == 1],
            d_split$Nullomer$count[d_split$Nullomer$Value == 0],
            d_split$Present$count[d_split$Present$Value == 0]),
            nrow = 2
        )

        # Perform chi-squared test
        chi_test <- chisq.test(contingency_table)
        
        # Create a summary table
        tibble(
            Key = unique(.x$Key),
            p_value = chi_test$p.value,
            Nullomer_Break = d_split$Nullomer$count[d_split$Nullomer$Value == 1],
            Present_Break = d_split$Present$count[d_split$Present$Value == 1],
            Nullomer_NoBreak = d_split$Nullomer$count[d_split$Nullomer$Value == 0],
            Present_NoBreak = d_split$Present$count[d_split$Present$Value == 0]
        )
    }) %>% 
    dplyr::mutate(
        Signif_level = dplyr::case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            TRUE ~ "N.S."
        )
    ) %>% 
    dplyr::select(Key, p_value, Signif_level)

fwrite(
    as.data.table(p_values), 
    "./SNP_pvalues.csv"
)