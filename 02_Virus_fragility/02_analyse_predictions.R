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

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

# get thresholds and column names
thresholds <- fread("../data/models/python/lightgbm/best_LGBM_model_thresholds.csv")
col_names <- thresholds$FPR
col_names <- gsub("\\.", "_", paste0("Pred_", col_names))
target_fpr_ind <- 3 # false positive rate of 0.05

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
labels <- fread("./data/virus_labels.csv")
pred_files <- fread("./data/filenames.csv", header = FALSE)
pred_files <- pred_files$V1
pred_files <- paste0("../data/", pred_files)
pred_names <- stringr::str_remove_all(
    string = pred_files, 
    pattern = "../data/human_viruses/|/final_pred.parquet"
)

# import predictions
path_to_pred <- "../data/Whole_genome_pred"
human_pred <- arrow::read_parquet(paste0(
    path_to_pred, "/", col_names[target_fpr_ind], ".parquet"
))

res_list <- pbapply::pblapply(1:length(pred_names), function(x){
    prediction <- arrow::read_parquet(pred_files[x])
    prediction <- dplyr::select(prediction, -True)
    colnames(prediction) <- gsub("\\.", "_", colnames(prediction))

    ind_of_int <- match(pred_names[x], labels$title)
    dt <- data.table(
        pred = unname(colSums(prediction)[target_fpr_ind]),
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
        "terminus|FMDV|junction|primer|primary|neuraminidase|intron|",
        "iteration|carcinoma|for|delet"
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
            grepl("Gammaherpesvirus", label, ignore.case = TRUE) ~ "Epstein_Barr_virus",
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
        x = "Breakage prob. per base",
        y = "Density"
    )

dir.create(path = "../figures/Virus_fragility/", showWarnings = FALSE)
ggsave(
    filename = paste0(
        "../figures/Virus_fragility/",
        "Norm_sequence_fragility.pdf"
    ),
    plot = p1,
    height = 7, width = 8
)

res_by_humans <- res %>% 
    dplyr::filter(
        grepl("^Human", species, ignore.case = TRUE) | 
        grepl(
            "^Hepatitis_B|^Hepatitis_C|Human_papillomavirus|Epstein_Barr_virus",
            species, ignore.case = TRUE
        )
    ) %>% 
    dplyr::group_by(species) %>% 
    dplyr::summarise(ppb = median(ppb, na.rm = TRUE)) %>% 
    dplyr::ungroup() 

get_density_ppb <- density(res_by_humans$ppb)
df_density_human_ppb <- tibble(
    x = get_density_ppb$x,
    y = get_density_ppb$y
    ) %>% 
    dplyr::filter(x >= 0)
    # dplyr::filter(x >= 0 & x <= 1)

to_label <- res %>% 
    dplyr::filter(grepl(
        "^Hepatitis_B|^Hepatitis_C|Human_papillomavirus|Epstein_Barr_virus",
        species, ignore.case = TRUE
    )) %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
        mean_ppb = mean(ppb),
        sd_ppb = sd(ppb),
        n = dplyr::n(),  
        se_ppb = sd_ppb / sqrt(n)
    ) %>%
    dplyr::mutate(
        nearest_y_ppb = sapply(.$mean_ppb, get_nearest_val, df_density_human_ppb),
        hex = "#990000",
        species = stringr::str_replace_all(
            string = species, pattern = "_", replacement = " "
        ),
        species = stringr::str_to_title(species),
        species = stringr::str_replace(
            string = species, pattern = "Virus", replacement = "virus"
        )
    )

left_part_human <- dplyr::filter(
    df_density_human_ppb, x < vline_val
)
right_part_human <- dplyr::filter(
    df_density_human_ppb, x >= vline_val
)

p1_humans <- df_density_human_ppb %>% 
    ggplot(aes(x = x, y = y)) + 
    geom_ribbon(
        data = left_part_human, 
        aes(x = x, ymin = 0, ymax = y), 
        fill = "#2166ac", alpha = 0.35
    ) +
    geom_ribbon(
        data = right_part_human, 
        aes(x = x, ymin = 0, ymax = y), 
        fill = "#b2182b", alpha = 0.25
    ) +
    geom_line(linewidth = 1) + 
    geom_segment(
        data = data.frame(
            x = vline_val, 
            y = 0, 
            xend = vline_val, 
            yend = get_nearest_val(val = vline_val, data = df_density_human_ppb)
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
    scale_color_identity() + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    coord_cartesian(xlim = c(0, 1.3), ylim = c(0, 1.7)) +  
    labs(
        x = "Breakage prob. per base",
        y = "Density"
    )

ggsave(
    filename = paste0(
        "../figures/Virus_fragility/",
        "For_Paper_Norm_sequence_human_fragility.pdf"
    ),
    plot = p1_humans,
    height = 6, width = 7
)

p1_humans_labels <- p1_humans + 
    ggrepel::geom_text_repel(
        data = to_label, aes(
            x = mean_ppb, y = nearest_y_ppb, 
            label = species, col = hex
        ),
        box.padding = 4, 
        max.overlaps = Inf,
        min.segment.length = unit(0.01, "lines"),
        segment.color = to_label$hex,
        segment.size = 0.3,
        size = 6
    ) + 
    coord_cartesian(xlim = c(0, 1))

ggsave(
    filename = paste0(
        "../figures/Virus_fragility/",
        "Norm_sequence_human_fragility_labels.pdf"
    ),
    plot = p1_humans_labels,
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
        x = "Breakage prob. per base",
        y = "Density"
    )

ggsave(
    filename = paste0(
        "../figures/Virus_fragility/",
        "Norm_sequence_fragility_by_organism.pdf"
    ),
    plot = p1,
    height = 7, width = 8
)

#' Thought: fragility scores should be normalised by the length of the sequence 
#' to cross-compare, but we should also look at the absolute number of breaks. 
#' If a host genome gets infected by an enormously large virus and it's predicted
#' to be more fragile, then it will bring in more breakage potential than a tiny virus.
res_raw_by_organism <- res %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
        sum_pred = sum(pred),
        sum_len = sum(len),
        ppb = sum_pred / sum_len,
        mean_pred = mean(pred)
    ) %>%
    dplyr::ungroup() %>%
    # dplyr::filter(!grepl("Homo", species)) %>% 
    dplyr::arrange(desc(mean_pred))

to_label_raw_by_organism_cancer <- res_raw_by_organism %>% 
    dplyr::filter(grepl(
        "^Hepatitis_B|^Hepatitis_C|Human_papillomavirus|Epstein_Barr_virus",
        species, ignore.case = TRUE
    )) %>%
    dplyr::mutate(
        hex = "#990000",
        species = stringr::str_replace_all(
            string = species, pattern = "_", replacement = " "
        ),
        species = stringr::str_to_title(species),
        species = stringr::str_replace(
            string = species, pattern = "Virus", replacement = "virus"
        )
    )

to_label_raw_by_organism_top <- res_raw_by_organism %>% 
    dplyr::arrange(desc(mean_pred), desc(ppb)) %>%
    head(10) %>%
    dplyr::mutate(
        hex = dplyr::case_when(
            grepl("^Homo_sapiens", species) ~ "#000000",
            .default = "#006400"
        ),
        species = stringr::str_replace_all(
            string = species, pattern = "_", replacement = " "
        ),
        species = stringr::str_to_title(species)
    )

to_label_raw_by_organism_sars <- res_raw_by_organism %>% 
    dplyr::filter(grepl(
        "SARS|^Bat_sars|^Bat_corona|^Human_coronavirus",
        species, ignore.case = TRUE
    )) %>%
    dplyr::mutate(
        hex = "#FFA435",
        species = stringr::str_replace_all(
            string = species, pattern = "_", replacement = " "
        ),
        species = stringr::str_to_title(species),
        species = stringr::str_replace(
            string = species, pattern = "Virus", replacement = "virus"
        )
    )

to_label_raw_by_organism <- dplyr::bind_rows(
    to_label_raw_by_organism_top, 
    to_label_raw_by_organism_cancer,
    to_label_raw_by_organism_sars
)

p1 <- res_raw_by_organism %>% 
    ggplot(aes(x = ppb, y = mean_pred)) + 
    geom_point(
        data = res_raw_by_organism %>% 
            dplyr::filter(!species %in% to_label_raw_by_organism$species),
        size = 2, alpha = 0.7, col = "#e7e7e7"
    ) + 
    geom_point(
        data = to_label_raw_by_organism,
        aes(col = hex),
        size = 2, alpha = 1
    ) + 
    ggrepel::geom_text_repel(
        data = to_label_raw_by_organism, 
        aes(label = species, col = hex),
        box.padding = 2.5, 
        max.overlaps = Inf,
        min.segment.length = unit(0.01, "lines"),
        segment.color = to_label_raw_by_organism$hex,
        segment.size = 0.3,
        size = 3
    ) + 
    coord_cartesian(xlim = c(0, 1)) + 
    scale_color_identity() + 
    scale_y_continuous(labels = scales::label_number(scale = 1e-6)) +
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    labs(
        title = "DNA virus sequence fragility averaged by organisms",
        subtitle = "Cancer-associated (red), most fragile (green), and SARS (orange) viruses",
        x = "Breakage prob. per base",
        y = expression("Absolute number of breaks, x10"^6*"")
    )

ggsave(
    filename = paste0(
        "../figures/Virus_fragility/",
        "Raw_sequence_fragility_by_organism.pdf"
    ),
    plot = p1,
    height = 7, width = 9
)

p2 <- res_raw_by_organism %>% 
    ggplot(aes(x = ppb, y = log2(mean_pred))) + 
    geom_point(
        data = res_raw_by_organism %>% 
            dplyr::filter(!species %in% to_label_raw_by_organism$species),
        size = 2, alpha = 0.7, col = "#e7e7e7"
    ) + 
    geom_point(
        data = to_label_raw_by_organism,
        aes(col = hex),
        size = 2, alpha = 1
    ) + 
    ggrepel::geom_text_repel(
        data = to_label_raw_by_organism, 
        aes(label = species, col = hex),
        box.padding = 2, 
        max.overlaps = Inf,
        min.segment.length = unit(0.01, "lines"),
        segment.color = to_label_raw_by_organism$hex,
        segment.size = 0.3,
        size = 3
    ) + 
    coord_cartesian(xlim = c(0, 1)) + 
    scale_color_identity() + 
    # scale_y_continuous(labels = scales::label_number(scale = 1e-6)) +
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    labs(
        title = "DNA virus sequence fragility averaged by organisms",
        subtitle = "Cancer-associated (red), most fragile (green), and SARS (orange) viruses",
        x = "Breakage prob. per base",
        y = "Log2-scaled number of breaks"
    )

ggsave(
    filename = paste0(
        "../figures/Virus_fragility/",
        "Log2_scaled_sequence_fragility_by_organism.pdf"
    ),
    plot = p2,
    height = 7, width = 9
)

# combine the plots
p1 <- res_raw_by_organism %>% 
    dplyr::filter(!grepl("Homo", species)) %>% 
    ggplot(aes(x = ppb, y = mean_pred)) + 
    geom_point(
        data = res_raw_by_organism %>% 
            dplyr::filter(!grepl("Homo", species)) %>% 
            dplyr::filter(!species %in% to_label_raw_by_organism$species),
        size = 2, alpha = 0.7, col = "#e7e7e7"
    ) + 
    geom_point(
        data = to_label_raw_by_organism %>%
            dplyr::filter(!grepl("Homo", species)),
        aes(col = hex),
        size = 2, alpha = 1
    ) + 
    ggrepel::geom_text_repel(
        data = to_label_raw_by_organism %>% 
            dplyr::filter(!grepl("Homo", species)), 
        aes(label = species, col = hex),
        box.padding = 3, 
        max.overlaps = Inf,
        min.segment.length = unit(0.01, "lines"),
        segment.color = to_label_raw_by_organism %>% 
            dplyr::filter(!grepl("Homo", species)) %>% 
            dplyr::pull(hex),
        segment.size = 0.3,
        size = 3
    ) + 
    coord_cartesian(xlim = c(0, 1)) + 
    scale_color_identity() + 
    scale_y_continuous(labels = scales::label_number(scale = 1e-6)) +
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    labs(x = "", y = expression("Absolute number of breaks, x10"^6*""))
    
p2 <- p2 + labs(
    title = "", subtitle = "", 
    x = "Breakage prob. per base"
)
main_title <- grid::textGrob(
    "DNA virus sequence fragility averaged by organisms", 
    gp = grid::gpar(fontsize = 16)
)
sub_title <- grid::textGrob(
    "Cancer-associated (red), most fragile (green), and SARS (orange) viruses",
    gp = grid::gpar(fontsize = 13)
)

combined_plots <- gridExtra::grid.arrange(p1, p2, ncol = 1)
pdf(
    file = paste0(
        "../figures/Virus_fragility/",
        "Raw_and_Log2_scaled_sequence_fragility_by_organism.pdf"
    ),
    height = 12, width = 9
)
gridExtra::grid.arrange(
    main_title, sub_title, 
    combined_plots,
    heights = c(0.9, 0.4, 12)
)
plot.saved <- dev.off()

########################################################################
#' Human viruses vs. hot blooded viruses vs. cold blooded viruses.
#' 
#' 1. Download Supplement Table S1 Dataset from paper: 10.1186/s12862-014-0178-z.
#' 2. Download taxonomizr 
taxonomizr::prepareDatabase('accessionTaxa.sql')

warm_cold_labels <- as_tibble(fread("./data/Endotherm_Ectotherm_animals.csv"))
warm_cold_labels <- warm_cold_labels %>% 
    dplyr::mutate(
        taxaId = getId(Species, "accessionTaxa.sql"),
        class = as.matrix(getTaxonomy(
                taxaId, 'accessionTaxa.sql'
            ))[, "class"]
    ) %>% 
    suppressWarnings()

get_genbank_common_name <- function(df){
    genbank_row <- which(df$type == "genbank common name")[1]
    common_row <- which(df$type == "common name")[1]
    if(length(genbank_row) > 0){
        return(df$name[genbank_row])
    } else if(length(common_row) > 0){
        return(df$name[common_row])
    } 
    return("")
}
df_ID <- getCommon(warm_cold_labels$taxaId, 'accessionTaxa.sql')
df_ID <- lapply(1:length(df_ID), function(x){
    if(!is.null(df_ID[[x]])){
        output <- get_genbank_common_name(df_ID[[x]])
        if(is.na(output)) output <- ""
        return(output)
    } 
    return("")
})
warm_cold_labels <- warm_cold_labels %>% 
    dplyr::mutate(common_name = unlist(df_ID)) %>% 
    dplyr::filter(common_name != "")

this_common_label <- stringr::str_split(
    string = warm_cold_labels$common_name,
    pattern = " "
)

warm_cold_labels <- warm_cold_labels %>% 
    dplyr::mutate(
        common_name = sapply(this_common_label, tail, n = 1),
        common_name = stringr::str_to_title(common_name)
    )

# match names
warm_cold_labels <- warm_cold_labels %>% 
    dplyr::select(Class, Blood, common_name) %>% 
    dplyr::distinct()

# Manual cleaning
warm_cold_labels_modified <- warm_cold_labels %>% 
    dplyr::filter(!grepl(paste0(
            "Tit"
        ), common_name, ignore.case = TRUE
    ))

warm_cold_labels_modified <- warm_cold_labels_modified %>% 
    dplyr::select(Class, Blood, common_name) %>% 
    dplyr::distinct()

res_modified <- res %>% 
    dplyr::select(-label, -perc_rank) %>% 
    dplyr::mutate(Class = "", Blood = "", common_name = "")

for(x in 1:nrow(warm_cold_labels_modified)){
    this_common_name <- warm_cold_labels_modified$common_name[x]
    this_blood <- warm_cold_labels_modified$Blood[x]
    this_class <- warm_cold_labels_modified$Class[x]
    
    # find match if any
    ind_match <- which(grepl(
        this_common_name,
        res_modified$species, 
        ignore.case = TRUE
    ))

    # only replace if empty in table
    if(length(ind_match) > 0){
        ind_match <- ind_match[res_modified$common_name[ind_match] == ""]
        if(length(ind_match) > 0){
            res_modified$common_name[ind_match] <- this_common_name
            res_modified$Blood[ind_match] <- this_blood
            res_modified$Class[ind_match] <- this_class
        }
    }
}

# Manual cleaning and editing
res_modified <- res_modified %>% 
    dplyr::filter(Class != "")

man_modify_as_needed <- function(ind, Class, Blood, common_name){
    res_modified$Class[ind] <<- Class
    res_modified$Blood[ind] <<- Blood
    res_modified$common_name[ind] <<- common_name
}

# Giant pandas are warm blooded mammals
to_change <- which(grepl("Ailuropoda_melanoleuca", res_modified$species))
man_modify_as_needed(ind = to_change, Class = "Mammals", Blood = "Warm", common_name = "Giant Panda")

# Some snakes were misclassified
to_change <- which(grepl("snake", res_modified$species))
man_modify_as_needed(ind = to_change, Class = "Reptilian", Blood = "Cold", common_name = "Snake")

# Some snakes were misclassified
to_change <- which(grepl("Porcine", res_modified$species))
man_modify_as_needed(ind = to_change, Class = "Mammals", Blood = "Warm", common_name = "Pigs")

# Only keep rats if it actually is a rat
to_change <- which(grepl("^Rat|rattus", res_modified$species))
man_modify_as_needed(ind = to_change, Class = "Mammals", Blood = "Warm", common_name = "Rat")

# Sheeps
to_change <- which(grepl("Ovine", res_modified$species))
man_modify_as_needed(ind = to_change, Class = "Mammals", Blood = "Warm", common_name = "Sheep")

# Fish
to_change <- which(grepl("Sparus_aurata", res_modified$species))
man_modify_as_needed(ind = to_change, Class = "Fish", Blood = "Cold", common_name = "Fish")

# Butterfly
to_change <- which(grepl("Heliconius_erato", res_modified$species))
man_modify_as_needed(ind = to_change, Class = "Insecta", Blood = "Cold", common_name = "Butterfly")

# Canine
to_change <- which(grepl("Canine", res_modified$species))
man_modify_as_needed(ind = to_change, Class = "Mammals", Blood = "Warm", common_name = "Dog")

# Metapenaeopsis
to_change <- which(grepl("Metapenaeopsis|Metapenaeus", res_modified$species))
man_modify_as_needed(ind = to_change, Class = "Malacostraca", Blood = "Cold", common_name = "Shrimps")

# remove the following
to_discard <- which(
    grepl("Rat", res_modified$common_name) & 
    !grepl("Mammals", res_modified$Class)
)
man_modify_as_needed(ind = to_discard, Class = "", Blood = "", common_name = "")

to_discard <- which(
    grepl("Ape", res_modified$common_name) & 
    grepl("Mammals", res_modified$Class)
)
man_modify_as_needed(ind = to_discard, Class = "", Blood = "", common_name = "")

# Grapevine, sapelo
to_discard <- which(grepl(paste0(
    "Grapevine|Sapelo|Chloris|Cucurbita"
), res_modified$species))
man_modify_as_needed(ind = to_discard, Class = "", Blood = "", common_name = "")

res_modified <- res_modified %>% 
    dplyr::filter(Class != "")

# human viruses
blood_human <- res_raw_by_organism %>% 
    dplyr::filter(
        grepl("^Human", species, ignore.case = TRUE) | 
        grepl(
            "^Hepatitis_B|^Hepatitis_C|Human_papillomavirus|Epstein_Barr_virus",
            species, ignore.case = TRUE
        )
    ) %>%
    dplyr::mutate(
        Blood = "Human viruses",
        Class = "Mammals",
        hex = "#e7e7e7"
    ) %>% 
    dplyr::select(mean_pred, species, ppb, Class, Blood) %>% 
    dplyr::rename(pred = mean_pred)

#' warm blooded (endotherm): Mammals, primates, humans, dolphin, birds, sinsectivores
#' cold blooded (ectotherm): Fishes, reptiles, amphibians
blood_all <- res_modified %>% 
    dplyr::select(pred, species, ppb, Class, Blood) %>% 
    dplyr::bind_rows(., blood_human) %>% 
    dplyr::mutate(hex = "#e7e7e7") %>% 
    dplyr::group_by(Blood) %>%
    dplyr::mutate(median_ppb = median(ppb)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(median_ppb)) %>%
    dplyr::mutate(
        hex = dplyr::case_when(
            Blood == "Warm" ~ "#b2182b",
            Blood == "Cold" ~ "#2166ac",
            Blood == "Human viruses" ~ "#b2182b",
            .default = hex
        ),
        Blood = dplyr::case_when(
            Blood == "Warm" ~ "Warm blooded viruses",
            Blood == "Cold" ~ "Cold blooded viruses",
            .default = Blood
        ),
        Blood = forcats::fct_inorder(Blood)
    )

blood_shades <- c(
    rgb(178,24,43, alpha = 0.1, maxColorValue = 255),
    rgb(33,102,172, alpha = 0.1, maxColorValue = 255),
    rgb(33,102,172, alpha = 0.1, maxColorValue = 255)
)

p3 <- blood_all %>% 
    ggplot(aes(x = Blood, y = ppb, fill = hex)) +
    geom_boxplot(alpha = 0.6, col = "#000000") + 
    ggsignif::geom_signif(
        comparisons = list(
            c("Warm blooded viruses", "Human viruses"),
            c("Cold blooded viruses", "Human viruses")
        ),
        map_signif_level = TRUE,
        y_position = c(1, 1.1),
        test = "t.test",
        textsize = 5,
        vjust = -0.1,
        margin_top = 0.1,
        tip_length = c(0.05, 0.05, 0.15, 0.05)
    ) +
    geom_jitter(width = 0.1, alpha = 0.75) + 
    coord_cartesian(ylim = c(0, 1.3)) + 
    scale_fill_identity() + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        # axis.text.y = element_blank(),
        legend.position = "none"
    ) + 
    labs(
        title = "DNA virus sequence fragility by blood class",
        x = "Host species",
        y = "Breakage prob. per base"
    )

ggsave(
    filename = paste0(
        "../figures/Virus_fragility/",
        "PPB_sequence_fragility_by_blood.pdf"
    ),
    plot = p3,
    height = 7, width = 8
)

# plot by species type
p4 <- blood_all %>% 
    dplyr::group_by(Class) %>%
    dplyr::mutate(median_ppb = median(ppb)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(median_ppb)) %>%
    dplyr::mutate(Class = forcats::fct_inorder(Class)) %>% 
    ggplot(aes(x = Class, y = ppb, fill = Class)) +
    geom_boxplot(alpha = 0.75, col = "#000000") + 
    geom_jitter(width = 0.1, alpha = 0.75) + 
    coord_cartesian(ylim = c(0, 1)) + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        legend.position = "none"
    ) + 
    labs(
        title = "DNA virus sequence fragility by species class",
        x = "", y = "Breakage prob. per base"
    )

ggsave(
    filename = paste0(
        "../figures/Virus_fragility/",
        "PPB_sequence_fragility_by_species-class.pdf"
    ),
    plot = p4,
    height = 7, width = 9
)

blood_shades <- c(
    rgb(33,102,172, alpha = 0.1, maxColorValue = 255),
    rgb(178,24,43, alpha = 0.1, maxColorValue = 255),
    rgb(178,24,43, alpha = 0.1, maxColorValue = 255)
)

p3 <- blood_all %>% 
    ggplot(aes(x = ppb, y = Blood, fill = hex)) +
    ggsignif::geom_signif(
        comparisons = list(
            c("Warm blooded viruses", "Human viruses"),
            c("Cold blooded viruses", "Human viruses")
        ),
        map_signif_level = TRUE,
        y_position = c(1.4, 1.5),
        test = "t.test",
        textsize = 5,
        vjust = -0.1,
        margin_top = -0.2,
        tip_length = 0,
        size = 0.8
    ) +
    geom_boxplot(
        aes(x = ppb, y = Blood, fill = hex), 
        alpha = 0.6, col = "#000000"
    ) + 
    geom_jitter(height = 0.2, alpha = 0.75) + 
    coord_cartesian(xlim = c(0, 1.3)) + 
    scale_fill_identity() + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        axis.text.y = element_blank(),
        legend.position = "none"
    ) + 
    # scale_fill_manual(values = blood_shades) + 
    labs(
        x = "Breakage prob. per base",
        y = ""
    )

ggsave(
    filename = paste0(
        "../figures/Virus_fragility/",
        "For_Paper_PPB_sequence_fragility_by_blood.pdf"
    ),
    plot = p3,
    height = 6, width = 7
)