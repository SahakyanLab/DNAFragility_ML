# import libraries
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
pbapply::pboptions(char = "=", type = "txt")

# source functions
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
shuffle <- as.logical(args[2])
setwd(my.path)

name <- ifelse(shuffle, "shuffle", "genome")
upper_name <- stringr::str_to_title(name)

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
sequence <- Biostrings::readDNAStringSet(
    filepath = paste0("./data/nullomer_and_", name, ".fasta.gz"),
    format = "fasta"
)
seq_lens <- width(sequence)

# import predictions
# pred_files <- list.files(
#     path = paste0("../data/nullomer_and_", name),
#     pattern = "final_pred",
#     full.names = TRUE,
#     recursive = TRUE
# )
pred_files <- ifelse(
    shuffle, 
    "./data/nullomer_and_shuffle.csv",
    "./data/nullomer_and_genome.csv"
)
pred_files <- fread(pred_files, header = FALSE)
pred_files <- pred_files$V1
pred_names <- stringr::str_extract(
    string = pred_files,
    pattern = paste0("(Nullomer|", upper_name, ")")
)

res <- pbapply::pblapply(1:length(pred_names), function(x){
    prediction <- arrow::read_parquet(paste0(
        paste0(base_folder, pred_files[x])
    ))
    prediction <- dplyr::select(prediction, -True)
    colnames(prediction) <- gsub("\\.", "_", colnames(prediction))
    
    preds <- colSums(prediction, na.rm = TRUE)
    ppb <- preds / seq_lens[x]
    return(ppb)
})
nullomer_ind <- which("Nullomer" == pred_names)
genome_ind <- setdiff(1:length(pred_names), nullomer_ind)

res_nullomer <- do.call(rbind, res[nullomer_ind])
res_nullomer <- as.data.table(res_nullomer)
res_nullomer[, Sequence := "Nullomer"]

res_genome <- do.call(rbind, res[genome_ind])
res_genome <- as.data.table(res_genome)
res_genome[, Sequence := upper_name]

res_all <- rbind(res_nullomer, res_genome)
res_all <- as_tibble(res_all) %>% 
    tidyr::pivot_longer(
        cols = -Sequence,
        names_to = "Key",
        values_to = "Value"
    ) %>% 
    dplyr::mutate(
        Key = plot_titles[Key],
        Key = forcats::fct_inorder(Key)
    ) %>% 
    tidyr::drop_na()

p1 <- res_all %>% 
    dplyr::mutate(Sequence = factor(Sequence, 
        levels = c("Nullomer", upper_name)
    )) %>%
    ggplot(aes(x = Sequence, y = Value, fill = Sequence)) + 
    geom_boxplot(
        # outlier.shape = NA,
        col = "black",
        alpha = 0.75,
        outlier.color = "#414a4c",
        outlier.alpha = 0.1
    ) + 
    # geom_jitter(
    #     width = 0.1, 
    #     alpha = 0.3
    # ) + 
    ggsignif::geom_signif(
        comparisons = list(c("Nullomer", upper_name)),
        map_signif_level = TRUE,
        textsize = 5,
        vjust = -0.1,
        y_position = 1,
        tip_length = 0
    ) +
    facet_wrap(vars(Key), nrow = 1) + 
    coord_cartesian(ylim = c(0,1.1)) +
    scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
    scale_fill_manual(values = c("#e07b00", "#004ead")) + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    labs(
        subtitle = paste0(
            "Fragility of nullomer-generated sequences vs. ",
            ifelse(shuffle, 
                "against its own shuffled version",
                "sampled from human genome"
            )
        ),
        x = "", y = "Probability per base"
    )

dir.create("../figures/Nullomers", showWarnings = FALSE)
ggsave(
    filename = paste0(
        "../figures/Nullomers/",
        upper_name, "_vs_nullomers.pdf"
    ),
    plot = p1,
    height = 5, width = 11
)