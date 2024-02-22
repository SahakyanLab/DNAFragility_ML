# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
data.table::setDTthreads(threads = 4)

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
exp <- as.character(args[2])
homer_path <- as.character(args[3])
setwd(my.path)

# other global vars
assembly <- "hg19"
file_name <- paste0(exp, "_kmer-8_seed-41964_zscore.csv")
true_breaks <- paste0("../data/feature_matrix/true_break_locations/", file_name)
control_breaks <- paste0("../data/feature_matrix/control_break_locations/", file_name)

# get true and negative control DNA strand breaks
true_breaks <- fread(true_breaks) %>% plyranges::as_granges()
control_breaks <- fread(control_breaks) %>% plyranges::as_granges()

# get medium range sequence effect
df_ranges <- fread(paste0(
    "../../03_Breakpoints_v2/03_FitCurves/data/ranges/",
    "kmer_8_Ranges_cutoffs_from_clustering_all-exp.csv"
))

mid_range <- as_tibble(df_ranges) %>% 
    dplyr::rename(full_exp_name = exp) %>%  
    dplyr::filter(stringr::str_detect(
        string = full_exp_name,
        pattern = paste0(
            "K562_Top2_mediated_DSBs/",
            ifelse(grepl("DMSO", exp), "DMSO", "ETO")
        )
    )) %>% 
    dplyr::filter(rowid == "ranges") %>%
    dplyr::select(dplyr::contains("range")) %>% 
    tidyr::pivot_longer(cols = everything()) %>% 
    dplyr::filter(!is.na(value)) %>% 
    dplyr::arrange(value) %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarise(value = mean(value)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(value) %>% 
    dplyr::pull(value) %>% 
    median(.)

round_to_nearest_even <- function(x) round(x/2)*2
mid_range <- round_to_nearest_even(mid_range)

# extract medium range sequence effect
true_MR <- plyranges::stretch(true_breaks, mid_range)
control_MR <- plyranges::stretch(control_breaks, mid_range)

# extract sequences for STREAM
dir.create("./data", showWarnings = FALSE)
# hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
# max_nr <- 1000000

# true_MR_seq <- getSeq(hg19, true_MR)
# true_MR_seq <- true_MR_seq[1:max_nr]
# names(true_MR_seq) <- paste0("seq_", 1:length(true_MR_seq))
# Biostrings::writeXStringSet(
#     true_MR_seq,
#     filepath = "./data/STREAM_TRUE.fasta",
#     format = "fasta"
# )

# control_MR_seq <- getSeq(hg19, control_MR)
# control_MR_seq <- control_MR_seq[1:max_nr]
# names(control_MR_seq) <- paste0("seq_", 1:length(control_MR_seq))
# Biostrings::writeXStringSet(
#     control_MR_seq,
#     filepath = "./data/STREAM_CONTROL.fasta",
#     format = "fasta"
# )

# # BiocManager::install("memes")
# library(memes)
# # options(meme_bin = "/home/imm/hert6114/Downloads/meme-5.5.5/bin/")
# check_meme_install()

# runStreme(
#   input = "./data/STREAM_TRUE.fasta",
#   control = "./data/STREAM_CONTROL.fasta",
#   outdir = "auto",
#   objfun = "de",
#   alph = "dna"
# )
# importStremeXML("./data/STREAM_TRUE_vs_STREAM_CONTROL/streme.xml")

#######################################################################################

#' download marge - the R API to Homer
#' mkdir lib
#' wget -P lib/ http://homer.ucsd.edu/homer/configureHomer.pl
#' vi ~/.bashrc
#' PATH=$PATH:/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/09_HOMER/lib/.//bin/
#' source ~/.bashrc
#' path_to_file = "/home/imm/hert6114/R/x86_64-pc-linux-gnu-library/4.3/marge-master"
#' devtools::install(
#'      pkg = path_to_file,
#'      quiet = FALSE,
#'      force = TRUE
#' )
options('homer_path' = homer_path)
suppressPackageStartupMessages(library(marge))
check_homer()
list_homer_packages()

# save files for HOMER
plyranges::write_bed(
    true_MR,
    "./data/HOMER_TRUE.bed"
)
plyranges::write_bed(
    control_MR,
    "./data/HOMER_CONTROL.bed"
)

marge::find_motifs_genome(
    x = "./data/HOMER_TRUE.bed",
    path = "./data",
    genome = "hg19",
    motif_length = c(8,10,12,14,16,18),
    scan_size = "given",
    background = "./data/HOMER_CONTROL.bed",
    only_known = FALSE,
    only_denovo = FALSE,
    cores = 16,
    overwrite = TRUE,
    keep_minimal = TRUE,
    scale_logos = TRUE
)

dir.create(
    path = paste0("./data/HOMER/", exp, "/"), 
    showWarnings = FALSE,
    recursive = TRUE
)
files_to_mv <- c(
    "HOMER_CONTROL.bed", "HOMER_TRUE.bed", 
    "homerMotifs.all.motifs", "knownResults.txt"
)
files_moved <- lapply(files_to_mv, function(f){
    system(paste0("mv ./data/", f, " ./data/HOMER/", exp))
})

# read the motifs
known <- marge::read_known_results(
        path = paste0("./data/HOMER/", exp), 
        homer_dir = TRUE
    ) %>% 
    suppressWarnings()
denovo <- marge::read_denovo_results(
        path = paste0("./data/HOMER/", exp), 
        homer_dir = TRUE
    ) %>% 
    suppressWarnings()

dir.create(
    path = "../figures/HOMER",
    showWarnings = FALSE
)
plot_denovo <- denovo %>% 
    dplyr::select(consensus, log_p_value, tgt_pct, bgd_pct) %>% 
    dplyr::filter(!grepl(pattern = "N", x = consensus)) %>%     
    dplyr::filter(is.finite(log_p_value)) %>% 
    dplyr::arrange(log_p_value) %>% 
    dplyr::mutate(ID = paste0(consensus, "_", tgt_pct*100, "%")) %>%  
    dplyr::filter(log_p_value >= 50) %>% 
    ggplot(aes(x = log_p_value, y = forcats::fct_inorder(ID))) +
    geom_point(size = 2.3) +
    geom_vline(xintercept = -log10(0.05), col = "red") + 
    theme_bw() + 
    theme(
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank()
    ) + 
    labs(
        title = "De novo motifs",
        x = "-log10 p-values",
        y = "Consensus motif w/target pct"
    )
ggsave(
    filename = paste0(
        "../figures/HOMER/", exp,
        "_denovo_motifs.pdf"),
    plot = plot_denovo,
    height = 10, width = 7
)

plot_known_seq <- known %>% 
    dplyr::select(motif_name, consensus, log_p_value, tgt_pct, bgd_pct) %>% 
    dplyr::filter(!grepl(pattern = "N", x = consensus)) %>%     
    dplyr::filter(is.finite(log_p_value)) %>% 
    dplyr::mutate(
        tgt_pct = tgt_pct*100, 
        bgd_pct = bgd_pct*100
    ) %>% 
    dplyr::filter(log_p_value >= 50) %>% 
    dplyr::arrange(log_p_value) %>% 
    dplyr::mutate(ID = paste0(motif_name, "_", tgt_pct, "%")) %>%  
    ggplot(aes(x = log_p_value, y = forcats::fct_inorder(ID))) +
    geom_point(size = 2.3) +
    geom_vline(xintercept = -log10(0.05), col = "red") + 
    theme_bw() + 
    theme_classic() + 
    labs(
        title = "Known motifs",
        x = "-log10 p-values",
        y = "Consensus motif w/target pct"
    )

ggsave(
    filename = paste0(
        "../figures/HOMER/", exp,
        "_known_motifs.pdf"),    
    plot = plot_known_seq,
    height = 10, width = 7
)