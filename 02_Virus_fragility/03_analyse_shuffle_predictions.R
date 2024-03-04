# import libraries
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))
pbapply::pboptions(char = "=", type = "txt")

# source functions
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

# get thresholds and column names
thresholds <- fread("../data/models/python/lightgbm/best_LGBM_model_thresholds.csv")
col_names <- thresholds$FPR
col_names <- gsub("\\.", "_", paste0("Pred_", col_names))
target_fpr_ind <- 3 # false positive rate of 0.05

# import predictions
path_to_pred <- "../data/Whole_genome_pred"
human_pred <- arrow::read_parquet(paste0(path_to_pred, "/", col_names[1], ".parquet"))

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
human_genome_len <- sum(refseq.table$end)

# import labels and predictions
labels <- fread("./data/virus_labels_shuffle.csv")
pred_files <- list.files(
    path = "../data/human_viruses_shuffle",
    pattern = "final_pred",
    full.names = TRUE,
    recursive = TRUE
)
pred_names <- stringr::str_remove_all(
    string = pred_files, 
    pattern = "../data/human_viruses_shuffle/|/final_pred.parquet"
)

res_list <- pbapply::pblapply(1:length(pred_names), function(x){
    prediction <- arrow::read_parquet(pred_files[x])
    prediction <- dplyr::select(prediction, -True)
    colnames(prediction) <- gsub("\\.", "_", colnames(prediction))

    ind_of_int <- match(pred_names[x], labels$label)
    dt <- data.table(
        pred = unname(colSums(prediction)[target_fpr_ind]),
        label = labels$label[ind_of_int],
        len = labels$length[ind_of_int]
    )

    return(dt)
})

res <- as_tibble(rbindlist(res_list)) %>% 
    dplyr::mutate(
        species = dplyr::case_when(
            grepl("Hepatitis_B_virus", label, ignore.case = TRUE) ~ "Hepatitis_B_virus",
            grepl("Hepatitis_c_virus", label, ignore.case = TRUE) ~ "Hepatitis_C_virus",
            grepl("Hepacivirus", label, ignore.case = TRUE) ~ "Hepatitis_C_virus",
            grepl("Epstein_Barr_virus", label, ignore.case = TRUE) ~ "Epstein_Barr_virus",
            grepl("Gammaherpesvirus", label, ignore.case = TRUE) ~ "Epstein_Barr_virus",
            grepl("Epstein", label, ignore.case = TRUE) ~ "Epstein_Barr_virus",
            grepl("Human_papillomavirus", label, ignore.case = TRUE) ~ "Human_papillomavirus"
        )
        # species = ifelse(
        #     grepl("_0", label, ignore.case = TRUE),
        #     stringr::str_replace(string = label, pattern = "_0", replacement = "_True"),
        #     species
        # )
    ) %>% 
    dplyr::bind_rows(
        tibble(
            pred = nrow(human_pred),
            label = "Homo_sapiens",
            species = "Homo_sapiens",
            len = human_genome_len
        ),
    ) %>% 
    dplyr::mutate(
        ppb = pred / len
    ) %>% 
    dplyr::arrange(desc(ppb)) %>% 
    dplyr::group_split(species) %>% 
    purrr::set_names(purrr::map(., ~ unique(.x$species)[1]))

res_human <- res$Homo_sapiens
res$Homo_sapiens <- NULL

# plot prob. per base distribution per species
plots <- lapply(names(res), function(label){
    xmax <- 1.0001
    temp <- res[[label]]
    
    # get density values
    get_density_ppb <- density(temp$ppb)
    df_density_ppb <- tibble(
        x = get_density_ppb$x,
        y = get_density_ppb$y
    )

    if(min(df_density_ppb$x, na.rm = TRUE) > 0){
        df_density_ppb <- df_density_ppb %>% 
            dplyr::add_row(x = 0, y = 0)   
    }

    get_nearest_val <- function(val, data){
        return(data[["y"]][which.min(abs(data[["x"]] - val))])
    }

    # map density values to ppb
    temp <- temp %>% 
        dplyr::mutate(
            nearest_y_ppb = sapply(.$ppb, get_nearest_val, df_density_ppb)
        )
    
    # actual species prob. per base
    true_df <- temp %>% 
        dplyr::filter(grepl("_0", label, ignore.case = TRUE)) %>% 
        dplyr::rename(x = ppb) %>% 
        dplyr::mutate(
            y = nearest_y_ppb,
            label = stringr::str_replace_all(
                string = species, 
                pattern = "_",
                replacement = " "
            ), 
            hex = "#b2182b"
        ) %>% 
        dplyr::select(-c(pred, len, species), x, y, label, hex)

    true_human <- res_human %>% 
        dplyr::rename(x = ppb) %>% 
        dplyr::mutate(
            nearest_y_ppb = true_df$nearest_y_ppb,
            y = nearest_y_ppb, 
            label = stringr::str_replace_all(
                string = species, 
                pattern = "_",
                replacement = " "
            ), 
            hex = "#2166ac"
        )  %>% 
        dplyr::select(-c(pred, len, species), x, y, label, hex)
    
    p1 <- df_density_ppb %>%
        ggplot(aes(x = x, y = y)) + 
        geom_line(linewidth = 1, col = "#636363") + 
        geom_segment(
            data = true_human,
            aes(
                x = x, y = 0, 
                xend = x, yend = max(df_density_ppb$y, na.rm = TRUE),
                col = hex
            ), 
            linetype = "dashed",
            linewidth = 0.7
        ) +
        ggrepel::geom_text_repel(
            data = true_human, 
            aes(x = x, y = nearest_y_ppb, label = label, col = hex),
            box.padding = 4, 
            max.overlaps = Inf,
            min.segment.length = unit(0.01, "lines"),
            segment.color = true_human$hex,
            segment.size = 0.3,
            size = 3
        ) + 
        geom_segment(
            data = true_df,
            aes(
                x = x, y = 0, 
                xend = x, yend = nearest_y_ppb,
                col = hex
            ), 
            linewidth = 0.7
        ) +
        ggrepel::geom_text_repel(
            data = true_df, 
            aes(x = x, y = nearest_y_ppb, label = label, col = hex),
            box.padding = 4, 
            max.overlaps = Inf,
            min.segment.length = unit(0.01, "lines"),
            segment.color = true_df$hex,
            segment.size = 0.3,
            size = 3
        ) + 
        scale_color_identity() + 
        coord_cartesian(xlim = c(0, xmax)) + 
        theme_bw() + 
        theme_classic() + 
        theme(text = element_text(size = 15)) + 
        labs(
            subtitle = true_df$label,
            x = "", y = ""
        )

    return(p1)    
})

pdf(
    file = paste0(
        "../figures/Virus_fragility/",
        "VirusShuffled.pdf"
    ),
    height = 10, width = 10
)
do.call(gridExtra::grid.arrange, c(plots, ncol = 2))
plot.saved <- dev.off()