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
suppressPackageStartupMessages(suppressWarnings(library(arrow)))
pbapply::pboptions(char = "=", type = "txt")

setwd("/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/05_DeltaFragility")

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

# import labels and predictions
path_to_delta_pred <- paste0(
    "/media/hert6114/Paddy_6TB/ProjectBoard_Patrick/",
    "04_DNAFragility/data/deltafragility_SNPs"
)
labels <- arrow::read_parquet("./data/all_SNP.parquet")

pred_files <- list.files(
    path = paste0(path_to_delta_pred, "/SNPs"),
    pattern = "final_pred",
    full.names = TRUE,
    recursive = TRUE
)
pred_names <- stringr::str_extract(
    string = pred_files, 
    pattern = paste0("before|after")
)

res_list <- pbapply::pblapply(1:length(pred_names), function(x){
    prediction <- arrow::read_parquet(pred_files[x])
    prediction <- dplyr::select(prediction, -True)
    colnames(prediction) <- gsub("\\.", "_", colnames(prediction))

    dt <- prediction %>% 
        dplyr::mutate(
            before_after = pred_names[x],
            labels
        )
    return(dt)
})
names(res_list) <- pred_names
res <- do.call(dplyr::bind_rows, res_list)

#' 1. Do pathogenic / benign mutations increase or decrease fragility?
#' Answer: pathogenic SNVs have a higher probability per base than benign.
res_split <- res %>% 
    dplyr::group_split(before_after) %>% 
    purrr::set_names(purrr::map(., ~ unique(.x$before_after)[1]))

res_all <- res_split[["after"]] %>% 
    dplyr::select(dplyr::starts_with('Pred')) %>% 
    dplyr::mutate(dplyr::across(
        dplyr::starts_with('Pred'),
        ~ . -res_split[["before"]][[dplyr::cur_column()]]
    )) %>% 
    dplyr::bind_cols(
        res_split[['after']] %>% 
            dplyr::select(-dplyr::starts_with('Pred'))
    )

res_summary <- res_all %>% 
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
    ) %>% 
    dplyr::group_by(Key, `Clinical Sig.`, Value) %>% 
    dplyr::summarise(count = dplyr::n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(Key, Value) %>% 
    dplyr::mutate(
        total = sum(count),
        frac = count / total * 100
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
        `Clinical Sig.` = as.factor(`Clinical Sig.`),
        Value = dplyr::case_when(
            Value == -1 ~ "Lower",
            Value == 0 ~ "Unch",
            Value == 1 ~ "Higher"
        ),
        Value = factor(Value, levels = c(
            "Lower", "Unch", "Higher"
        )),
        label_col = ifelse(
            `Clinical Sig.` == "Benign", 
            "black", "white"
        )
    ) %>% 
    suppressMessages()

combines_labels <- interaction(res_summary$Value, res_summary$`Clinical Sig.`)
levels(combines_labels) <- as.character(combines_labels)
res_summary$combined <- combines_labels

p1 <- res_summary %>%
    ggplot(aes(
        x = Value, y = frac, 
        fill = combined,
        pattern = combined
    )) + 
    # ggplot(aes(x = Value, y = frac, fill = `Clinical Sig.`)) +
    geom_bar_pattern(
        position = "stack",
        stat = "identity",
        col = "black",
        alpha = 0.75,
        pattern_colour = "white",
        pattern_spacing = 0.13,
        pattern_angle = 45
    ) +
    facet_wrap(vars(Key), nrow = 1) +
    geom_text(
        aes(
            label = sprintf("%.1f%%", frac),
            y = frac
        ), 
        position = position_stack(vjust = 0.5),
        col = res_summary$label_col,
        size = 5
    ) +
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        legend.position = "none"
    ) + 
    scale_fill_manual(values = c(
        rep("grey80", 3), 
        rep("#141414", 3)
    )) + 
    scale_pattern_manual(values = c(
        rep('crosshatch', 3), 
        rep('none', 3)
    )) +
    labs(x = "", y = "Fraction, %")

dir.create("../figures/deltafragility", showWarnings = FALSE)
ggsave(
    filename = "../figures/deltafragility/benign_vs_pathogenic_SNPs_bar.pdf",
    plot = p1,
    height = 5, width = 11
)

# calculate statistical significance
stat_test <- res_summary %>% 
    dplyr::select(-c(label_col, combined)) %>% 
    dplyr::group_split(Key, Value) %>% 
    purrr::map_dfr(~{
        # Extract 'rand' and 'pred' data
        pathogenic_data <- .x %>% filter(`Clinical Sig.` == "Pathogenic")
        benign_data <- .x %>% filter(`Clinical Sig.` == "Benign")

        # Perform z-test of difference in proportions
        prop_test <- prop.test(
            x = c(pathogenic_data$count, benign_data$count), 
            n = c(pathogenic_data$total, benign_data$total)
        )

        tibble(
            Key = unique(.x$Key),
            p_value = prop_test$p.value,
            delta = unique(.x$Value)
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
        )
    )

fwrite(
    as.data.table(stat_test),
    "./data/SNP_ClinSig_byGroup_pvalues.csv"
)

frac_by_group <- res_summary %>% 
    dplyr::select(-c(label_col, combined, frac, `Clinical Sig.`)) %>%
    dplyr::group_by(Key, Value) %>%
    dplyr::summarise(total = sum(count)) %>% 
    dplyr::ungroup() %>%
    dplyr::group_by(Key) %>% 
    dplyr::mutate(frac = total / sum(total) * 100) %>% 
    dplyr::arrange(Key, desc(frac)) %>% 
    dplyr::ungroup() %>% 
    suppressMessages()

frac_by_clin <- res_summary %>% 
    dplyr::select(-c(label_col, combined, frac)) %>%
    dplyr::group_by(Key, `Clinical Sig.`, Value) %>%
    dplyr::summarise(total = sum(count)) %>% 
    dplyr::ungroup() %>%
    dplyr::group_by(Key, `Clinical Sig.`) %>% 
    dplyr::mutate(frac = total / sum(total) * 100) %>% 
    dplyr::arrange(Key, `Clinical Sig.`, desc(frac)) %>% 
    dplyr::ungroup() %>% 
    suppressMessages()

frac_by_clin_pathogenic <- frac_by_clin %>% 
    dplyr::filter(`Clinical Sig.` == "Pathogenic") %>% 
    dplyr::filter(Value != "Unch") %>% 
    dplyr::mutate(diff = total / lag(total) * 100) %>% 
    dplyr::filter(Value == "Higher")

frac_by_clin_benign <- frac_by_clin %>% 
    dplyr::filter(`Clinical Sig.` == "Pathogenic") %>% 
    dplyr::filter(Value != "Unch") %>% 
    dplyr::mutate(diff = total / lag(total) * 100) %>% 
    dplyr::filter(Value == "Higher")