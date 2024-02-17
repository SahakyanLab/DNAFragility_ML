suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(arrow)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))

setwd("/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/01_LGBM_FullGenome")

#' 1. Import AUROC and feature importance data
df_auroc <- fread("../data/models/python/lightgbm/best_LGBM_model_ROC_curve.csv")
df_fi <- fread("../data/models/python/lightgbm/best_LGBM_model_feature_importance.csv")

#' 2. Plot both as vertical plots
p1 <- as_tibble(df_auroc) %>% 
    ggplot(aes(x = FPR, y = TPR)) + 
    geom_line(col = 'orange', linewidth = 1.3) + 
    geom_abline(
        slope = 1, intercept = 0, 
        linewidth = 1,
        linetype = "dashed", col = "grey"
    ) + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
    theme_bw() + 
    theme_classic() +
    theme(text = element_text(size = 25)) + 
    labs(x = "False Positive Rate", y = "True Positive Rate") 

ggsave(
    filename = "auroc.pdf",
    plot = p1,
    height = 6, width = 6
)

p2 <- as_tibble(df_fi) %>% 
    dplyr::arrange(importance) %>% 
    dplyr::mutate(feature = forcats::fct_inorder(feature)) %>% 
    ggplot(aes(x = importance, y = feature)) + 
    geom_segment(aes(x=0, xend=importance, y=feature, yend=feature), color="grey") +
    geom_point(color="orange", size = 5) +
    theme_bw() + 
    theme_classic() + 
    labs(x = "Importance", y = "") +
    theme(text = element_text(size = 10))

ggsave(
    filename = "fi.pdf",
    plot = p2,
    height = 5, width = 5
)
