# Load in all the cosmic mutations, work out a grouping, then do mutations and
# plot deltaRTs with those groupings
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.UCSC.hg38)))
pbapply::pboptions(char = '=')

my.path = as.character(args[1])
setwd(my.path)

# global vars
kmer <- 8
TopTC_ID <- readRDS(paste0(
        "../data/break_persistence/",
        "Top_Tissue_Cancer_ID.RData"
))

#' Remove the sample ID
df.parsed <- fread(
    file = paste0(
        "../data/BaseTable/",
        "kmer_", kmer, ".csv"
    ),
    select = "ID",
    showProgress = FALSE
)

#' Remove the sample ID
df <- as_tibble(df.parsed) %>%
    tidyr::separate(
        col = ID,
        into = c(
            "Chr", "Start", 
            "Tissue", "Cancer", 
            "BreakType", 
            "SampleID"
        ),
        sep = "_",
        remove = TRUE
    ) %>% 
    tidyr::unite(
        TC_ID,
        c("Tissue", "Cancer"), 
        remove = FALSE
    ) %>% 
    tidyr::unite(
        Chr_Start_ID, 
        c("Chr", "Start"),
        remove = FALSE
    ) %>% 
    dplyr::select(-c(SampleID, Tissue, Cancer)) %>% 
    dplyr::mutate(Start = as.numeric(Start))

#' Look at only the top occurrences of tisse_cancer IDs
df <- df[which(!is.na(match(df$TC_ID, TopTC_ID))),]

# for mapping breakage type with the colour
map.mut.cols <- c(
    "del" = "#6D7997", 
    "delins" = "#FBA32D", 
    "dup" = "#619B61", 
    "ins" = "#612940", 
    "inv" = "#4592C8",
    "intrachr" = "#1d4a7d",
    "interchr" = "#606a9f"
)

#' For each tisse_cancer IDs,
#'  1) get count of cases
plot.gaps <- pbapply::pblapply(1:length(TopTC_ID), function(x){
    total.case <- df %>% 
        dplyr::filter(TC_ID == TopTC_ID[x]) %>% 
        nrow(.)

    TopTC_ID.df <- df %>% 
        dplyr::filter(TC_ID == TopTC_ID[x]) %>% 
        dplyr::mutate(
            total = nrow(.),
            Start = as.numeric(Start)
        ) %>% 
        dplyr::arrange(Chr, Start) %>% 
        dplyr::group_by(Chr, BreakType) %>% 
        dplyr::mutate(Gap = Start-lag(Start)) %>% 
        dplyr::group_by(Gap, BreakType) %>% 
        dplyr::summarise(count = dplyr::n()) %>% 
        dplyr::arrange(Gap, desc(count), .by_group = TRUE) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(frac = count/total.case*100) %>% 
        dplyr::arrange(Gap, desc(frac), .by_group = TRUE) %>% 
        dplyr::mutate(hex.col = map.mut.cols[BreakType]) %>%
        suppressMessages()

    plot.title <- stringr::str_split(
        string = TopTC_ID[x],
        pattern = "_"
    )[[1]]

    plot.title <- paste0(
        "Tissue: ", plot.title[1], "\n",
        "Cancer: ", plot.title[2]
    )

    p1 <- TopTC_ID.df %>%
        dplyr::filter(Gap <= 10) %>% 
        dplyr::group_by(Gap) %>% 
        ggplot(aes(x = Gap, y = frac, fill = hex.col)) + 
        geom_bar(stat = "identity", position = "stack") + 
        scale_fill_identity() + 
        scale_x_continuous(
            breaks = 0:10
        ) + 
        theme_bw() + 
        theme_classic() + 
        theme(
            text = element_text(size = 15),
            plot.title = element_text(size = 11),
            plot.subtitle = element_text(size = 11)
        ) + 
        labs(
            title = plot.title,
            subtitle = paste0("Total: ", total.case),
            x = "", y = ""
        )

    return(p1)
})

pdf(
    file = paste0(
        "../figures/EDA/", 
        "GapBetweenTissueCancerID.pdf"
    ),
    height = 19, width = 13
    # units = "in", res = 600
)
do.call(gridExtra::grid.arrange, c(plot.gaps, ncol = 3))
plot.saved <- dev.off()

# dfs <- pbapply::pbsapply(1:length(TopTC_ID), function(x){
#     total.case <- df %>% 
#         dplyr::filter(TC_ID == TopTC_ID[x]) %>% 
#         nrow(.)

#     TopTC_ID.df <- df %>% 
#         dplyr::filter(TC_ID == TopTC_ID[x]) %>% 
#         dplyr::mutate(
#             total = nrow(.),
#             Start = as.numeric(Start)
#         ) %>% 
#         dplyr::arrange(Chr, Start) %>% 
#         dplyr::group_by(Chr, BreakType) %>% 
#         dplyr::mutate(Gap = Start-lag(Start)) %>% 
#         dplyr::group_by(Gap, BreakType) %>% 
#         dplyr::summarise(count = dplyr::n()) %>% 
#         dplyr::arrange(Gap, desc(count), .by_group = TRUE) %>% 
#         dplyr::ungroup() %>% 
#         dplyr::mutate(frac = count/total.case*100) %>% 
#         dplyr::arrange(Gap, desc(frac), .by_group = TRUE) %>% 
#         dplyr::mutate(hex.col = map.mut.cols[BreakType]) %>%
#         suppressMessages()

#     return(
#         TopTC_ID.df %>% 
#             dplyr::filter(frac <= 10) %>% 
#             dplyr::pull(frac) %>% 
#             sum(.)
#     )
# })

# mean(dfs)
# tibble(
#     TopTC_ID = TopTC_ID,
#     frac = dfs
#     ) %>%
#     dplyr::arrange(desc(frac))