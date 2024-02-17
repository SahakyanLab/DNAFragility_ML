# import libraries
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))
suppressPackageStartupMessages(suppressWarnings(library(taxonomizr)))
pbapply::pboptions(char = "=")

# source functions
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

# get thresholds and column names
thresholds <- fread("../data/models/python/lightgbm/best_LGBM_model_thresholds.csv")
col_names <- thresholds$FPR
col_names <- gsub("\\.", "_", paste0("Pred_", col_names))

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
labels <- fread("./virus_labels_shuffle.csv")

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
    prediction <- dplyr::select(prediction, -start)
    colnames(prediction) <- gsub("\\.", "_", colnames(prediction))

    ind_of_int <- match(pred_names[x], labels$title)
    dt <- data.table(
        pred = unname(colSums(prediction)[2]),
        label = labels$title[ind_of_int],
        organism = labels$organism[ind_of_int],
        len = labels$length[ind_of_int],
        accession = labels$accession[ind_of_int]
    )

    return(dt)
})
res <- as_tibble(rbindlist(res_list)) %>% 
    dplyr::filter(!grepl(paste0(
        "ORF|tissue|Ori|origin|hunan|protein|^_|esterase|domain|",
        "terminus|FMDV|junction|primer|primary|neuraminidase"
    ), label), ignore.case = TRUE) %>% 
    dplyr::mutate(
        species = gsub(" ", "_", organism),
        species = gsub("^(.*virus).*", "\\1", species),
        species = stringr::str_to_title(species)
    ) %>%  
    dplyr::select(-organism) %>%
    dplyr::mutate(
        species = dplyr::case_when(
            species == "Hepatitis_B_virus" ~ "Hepatitis_B_virus",
            species == "Hepatitis_c_virus" ~ "Hepatitis_C_virus",
            species == "Hepacivirus" ~ "Hepatitis_C_virus",
            species == "Epstein_Barr_virus" ~ "Epstein_Barr_virus",
            grepl("Epstein", label, ignore.case = TRUE) ~ "Epstein_Barr_virus",
            TRUE ~ species
        )
    ) %>% 
    dplyr::bind_rows(
        tibble(
            pred = nrow(human_pred),
            label = "Homo_sapiens",
            species = "Homo_sapiens",
            len = human_genome_len
        ),
        .
    ) %>% 
    dplyr::mutate(
        ppb = pred / len
    ) %>% 
    dplyr::arrange(desc(ppb)) %>% 
    dplyr::mutate(
        perc_rank = 1:nrow(.),
        perc_rank = (1 - (perc_rank / nrow(.))) * 100
    )

#' Before:
# get_density_ppb <- density(res$ppb)
# df_density_ppb <- tibble(
#     x = get_density_ppb$x,
#     y = get_density_ppb$y
# )
#' After:
res_by_species <- res %>% 
    dplyr::group_by(species) %>% 
    dplyr::summarise(ppb = mean(ppb, na.rm = TRUE)) %>% 
    dplyr::ungroup()
get_density_ppb <- density(res_by_species$ppb)
df_density_ppb <- tibble(
    x = get_density_ppb$x,
    y = get_density_ppb$y
)

get_density_raw <- res %>% 
    dplyr::filter(!grepl("Homo", species)) %>% 
    dplyr::pull(pred) %>%
    density()
df_density_raw <- tibble(
    x = get_density_raw$x,
    y = get_density_raw$y
)

get_nearest_val <- function(val, data){
    return(data[["y"]][which.min(abs(data[["x"]] - val))])
}

to_label <- res %>% 
    dplyr::filter(grepl(
        "^Hepatitis_B|^Hepatitis_C|Human_papillomavirus|Epstein_Barr_virus",
        species, ignore.case = TRUE
    )) %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
        mean_ppb = mean(ppb),
        mean_raw = mean(pred),
        sd_ppb = sd(ppb),
        sd_raw = sd(pred),      
        n = dplyr::n(),  
        se_ppb = sd_ppb / sqrt(n),
        se_raw = sd_raw / sqrt(n)
    ) %>%
    dplyr::mutate(
        nearest_y_ppb = sapply(.$mean_ppb, get_nearest_val, df_density_ppb),
        nearest_y_raw = sapply(.$mean_raw, get_nearest_val, get_density_raw),
        hex = "#990000",
        species = stringr::str_replace_all(
            string = species, pattern = "_", replacement = " "
        ),
        species = stringr::str_to_title(species),
        species = stringr::str_replace(
            string = species, pattern = "Virus", replacement = "virus"
        )
    )

# plot prediction distributions
vline_val <- res %>% 
    dplyr::filter(grepl("^Homo_sapiens$", label)) %>%
    dplyr::pull(ppb)
left_part <- dplyr::filter(df_density_ppb, x < vline_val)
right_part <- dplyr::filter(df_density_ppb, x >= vline_val)

p1 <- df_density_ppb %>% 
    ggplot(aes(x = x, y = y)) + 
    geom_ribbon(
        data = left_part, 
        aes(x = x, ymin = 0, ymax = y), 
        fill = "#2166ac", alpha = 0.35
    ) +
    geom_ribbon(
        data = right_part, 
        aes(x = x, ymin = 0, ymax = y), 
        fill = "#b2182b", alpha = 0.2
    ) +
    geom_line(linewidth = 1) + 
    geom_segment(
        data = data.frame(
            x = vline_val, 
            y = 0, 
            xend = vline_val, 
            yend = get_nearest_val(val = vline_val, data = df_density_ppb)
        ), 
        aes(x = x, y = y, xend = xend, yend = yend),
        col = "darkblue",
        linetype = "dashed",
        linewidth = 0.7
    )+
    geom_errorbar(
        data = to_label, aes(
            x = mean_ppb, y = nearest_y_ppb,
            xmin = mean_ppb - se_ppb,
            xmax = mean_ppb + se_ppb,
            col = hex
        ),
        width = 0.1,
        linewidth = 1.3
    ) + 
    ggrepel::geom_text_repel(
        data = to_label, aes(
            x = mean_ppb, y = nearest_y_ppb, 
            label = species, col = hex
        ),
        box.padding = 6, 
        max.overlaps = Inf,
        min.segment.length = unit(0.01, "lines"),
        segment.color = to_label$hex,
        segment.size = 0.3,
        size = 6
    ) + 
    scale_color_identity() + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    labs(
        title = "DNA virus sequence fragility averaged by species",
        subtitle = paste0(
            nrow(res), " viruses and ",
            length(unique(res$species)), " unique species."
        ),
        x = "Probability per base",
        y = "Density"
    )

path_to_save <- "../figures/Virus_fragility_shuffle/"
dir.create(path =path_to_save, showWarnings = FALSE)
ggsave(
    filename = paste0(
       path_to_save,
        "Norm_sequence_fragility.pdf"
    ),
    plot = p1,
    height = 7, width = 8
)

p1 <- df_density_ppb %>% 
    ggplot(aes(x = x, y = y)) + 
    geom_ribbon(
        data = left_part, 
        aes(x = x, ymin = 0, ymax = y), 
        fill = "#2166ac", alpha = 0.35
    ) +
    geom_ribbon(
        data = right_part, 
        aes(x = x, ymin = 0, ymax = y), 
        fill = "#b2182b", alpha = 0.25
    ) +
    geom_line(linewidth = 1) + 
    geom_segment(
        data = data.frame(
            x = vline_val, 
            y = 0, 
            xend = vline_val, 
            yend = get_nearest_val(val = vline_val, data = df_density_ppb)
        ), 
        aes(x = x, y = y, xend = xend, yend = yend),
        col = "darkblue",
        linetype = "dashed",
        linewidth = 1
    )+
    geom_errorbar(
        data = to_label, aes(
            x = mean_ppb, y = nearest_y_ppb,
            xmin = mean_ppb - se_ppb,
            xmax = mean_ppb + se_ppb,
            col = hex
        ),
        width = 0.1,
        linewidth = 1.1
    ) + 
    scale_color_identity() + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    coord_cartesian(ylim = c(0, 1.7)) + 
    labs(
        x = "Probability per base",
        y = "Density"
    )

ggsave(
    filename = paste0(
       path_to_save,
        "For_Paper_Norm_sequence_fragility.pdf"
    ),
    plot = p1,
    height = 6, width = 7
)

# Density plot of virus fragility, one per species  
res_by_organism <- res %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
        sum_pred = sum(pred),
        sum_len = sum(len),
        ppb = sum_pred / sum_len
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(ppb))

get_density_ppb_by_oragnism <- density(res_by_organism$ppb)
df_density_ppb_by_organism <- tibble(
    x = get_density_ppb_by_oragnism$x,
    y = get_density_ppb_by_oragnism$y
)

to_label_by_organism <- res_by_organism %>% 
    dplyr::filter(grepl(
        "^Hepatitis_B|^Hepatitis_C|Human_papillomavirus|Epstein_Barr_virus",
        species, ignore.case = TRUE
    )) %>%
    dplyr::mutate(
        nearest_y_ppb = sapply(
            .$ppb, 
            get_nearest_val, 
            df_density_ppb_by_organism
        ),
        hex = "#990000",
        species = stringr::str_replace_all(
            string = species, pattern = "_", replacement = " "
        ),
        species = stringr::str_to_title(species),
        species = stringr::str_replace(
            string = species, pattern = "Virus", replacement = "virus"
        )
    )

p1 <- df_density_ppb_by_organism %>% 
    ggplot(aes(x = x, y = y)) + 
    geom_line(linewidth = 1) + 
    geom_vline(
        xintercept = res %>% 
            dplyr::filter(grepl("^Homo_sapiens$", label)) %>%
            dplyr::pull(ppb),
        col = "darkblue",
        linetype = "dashed"
    ) + 
    ggrepel::geom_text_repel(
        data = to_label_by_organism, aes(
            x = ppb, y = nearest_y_ppb, 
            label = species, col = hex
        ),
        box.padding = 5, 
        max.overlaps = Inf,
        min.segment.length = unit(0.01, "lines"),
        segment.color = to_label_by_organism$hex,
        segment.size = 0.3,
        size = 5
    ) + 
    scale_color_identity() + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    labs(
        title = "DNA virus sequence fragility averaged by species",
        subtitle = paste0(
            "Total number of sequences: ",
            nrow(df_density_ppb_by_organism)
        ),
        x = "Probability per base",
        y = "Density"
    )

# ggsave(
#     filename = paste0(
#        path_to_save,
#         "Norm_sequence_fragility_by_organism.pdf"
#     ),
#     plot = p1,
#     height = 7, width = 8
# )