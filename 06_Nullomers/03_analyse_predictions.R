# import libraries
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
pbapply::pboptions(char = "=")

# source functions
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
# human_pred <- arrow::read_parquet(paste0(path_to_pred, "/", col_names[1], ".parquet"))

sequence <- Biostrings::readDNAStringSet(
    filepath = "./data/nullomer_and_shuffle.fasta.gz",
    format = "fasta"
)
seq_lens <- width(sequence)

# import predictions
pred_files <- list.files(
    path = "../data/nullomer_and_shuffle",
    pattern = "final_pred",
    full.names = TRUE,
    recursive = TRUE
)
pred_names <- stringr::str_extract(
    string = pred_files,
    pattern = "(Nullomer|Shuffle)"
)

res <- pbapply::pblapply(1:length(pred_names), function(x){
    prediction <- arrow::read_parquet(pred_files[x])
    prediction <- dplyr::select(prediction, -True)
    colnames(prediction) <- gsub("\\.", "_", colnames(prediction))
    
    preds <- colSums(prediction, na.rm = TRUE)
    ppb <- preds / seq_lens[x]
    return(ppb)
})
nullomer_ind <- which("Nullomer" == pred_names)
shuffle_ind <- which("Shuffle" == pred_names)

res_nullomer <- do.call(rbind, res[nullomer_ind])
res_nullomer <- as.data.table(res_nullomer)
res_nullomer[, Sequence := "Nullomer"]

res_shuffle <- do.call(rbind, res[shuffle_ind])
res_shuffle <- as.data.table(res_shuffle)
res_shuffle[, Sequence := "Shuffled"]

res_all <- rbind(res_nullomer, res_shuffle)
res_all <- as_tibble(res_all) %>% 
    tidyr::pivot_longer(
        cols = -Sequence,
        names_to = "Key",
        values_to = "Value"
    ) %>% 
    dplyr::mutate(
        Key = plot_titles[Key],
        Key = forcats::fct_inorder(Key)
    )

p1 <- res_all %>% 
    ggplot(aes(x = Sequence, y = Value, fill = Sequence)) + 
    geom_boxplot(
        outlier.shape = NA,
        col = "black",
        alpha = 0.75
    ) + 
    geom_jitter(
        width = 0.1, 
        alpha = 0.3
    ) + 
    ggsignif::geom_signif(
        comparisons = list(c("Nullomer", "Shuffled")),
        map_signif_level = TRUE,
        textsize = 5,
        vjust = -0.1,
        margin_top = -0.05,
        tip_length = 0
    ) +
    facet_wrap(vars(Key), nrow = 1) + 
    coord_cartesian(ylim = c(0,1)) + 
    scale_fill_manual(values = c("#e07b00", "#004ead")) + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    labs(
        subtitle = paste0(
            "Fragility of nullomer-generated sequences vs. ",
            "sampled from human genome"
        ),
        x = "", y = "Probability per base"
    )

dir.create("../figures/Nullomers", showWarnings = FALSE)
ggsave(
    filename = "../figures/Nullomers/True_vs_nullomers.pdf",
    plot = p1,
    height = 5, width = 11
)