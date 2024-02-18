# Load in all the cosmic mutations, work out a grouping, then do mutations and
# plot deltaRTs with those groupings
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.UCSC.hg38)))

args = commandArgs(trailingOnly = TRUE)
my.path = as.character(args[1])
setwd(my.path)

# global vars
kmer <- 8

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

#' Explore the possibility of a persistent breakage. That is, 
#' same chromosome and start position across multiple types of  
#' tissue and cancer combinations.
UniqueTC_per_ChrStart <- df %>% 
    dplyr::group_by(Chr_Start_ID) %>% 
    dplyr::summarise(Bp = dplyr::n_distinct(TC_ID)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(desc(Bp))

BreakPersistenceTable <- UniqueTC_per_ChrStart %>%
    dplyr::group_by(Bp) %>% 
    dplyr::summarise(count = dplyr::n()) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(frac = count/sum(count)*100) %>% 
    dplyr::arrange(Bp) %>%
    dplyr::select(-count) %>%  
    dplyr::rename_with(~c("Bp", "%"))

# BreakPersistenceTable %>% 
#     dplyr::rename_with(~c("Bp", "y")) %>% 
#     dplyr::mutate(cumsum = cumsum(y))

dir.create(path = "../data/break_persistence", showWarnings = FALSE)
fwrite(
    BreakPersistenceTable, 
    file = paste0(
        "../data/break_persistence/", 
        "Fractional_table.csv"
    )
)

p1 <- BreakPersistenceTable %>% 
    dplyr::rename_with(~c("Bp", "perc")) %>% 
    ggplot(aes(x = Bp, y = perc)) + 
    geom_line() + 
    geom_point() + 
    coord_cartesian(ylim = c(0, 100)) + 
    theme_bw() + 
    theme_classic() +
    theme(text = element_text(size = 15)) + 
    labs(
        title = "Percentage of breakage persistence (Bp) in all unique IDs",
        x = "Breakage persistence, Bp",
        y = "Contribution, %"
    )

dir.create(
    path = "../figures/break_persistence/",
    showWarnings = FALSE,
    recursive = TRUE
)
ggsave(
    filename = paste0(
        "../figures/break_persistence/", 
        "PercentBPAllUniqueIDs.pdf"
    ),
    plot = p1,
    height = 6, width = 8
)

BreakPersistence <- dplyr::left_join(
    x = df,
    y = UniqueTC_per_ChrStart,
    by = "Chr_Start_ID"
)

# BreakPersistence %>% 
#     dplyr::group_by(Bp) %>% 
#     dplyr::summarise(count = dplyr::n()) %>% 
#     dplyr::arrange(desc(Bp)) %>% 
#     head()

# Define the starting, middle, and end colors
gradient_colors <- c(
    "#443679", "#534294", "#5b48a2", "#445ca2", "#3e78b1", 
    "#538cc3", "#5298a5", "#6cb89b", "#88c89a", "#aad89b", 
    "#c4e798", "#e4f592", "#effaa5","#feffb9", "#faee9e", 
    "#f8de85", "#f4c371", "#f1a95c","#ea894e", "#e46941", 
    "#d15443", "#c14148", "#a32a40", "#891a3a", "#68142c"
)

# Q: which chromosomes have the highest proportion of high Bps? 
p1 <- BreakPersistence %>% 
    dplyr::group_by(Chr, Bp) %>% 
    dplyr::summarise(count = dplyr::n()) %>% 
    dplyr::mutate(frac = count/sum(count)*100) %>%  
    dplyr::group_by(Chr) %>% 
    dplyr::mutate(
        max_Bp = max(Bp),
        frac_of_max_Bp = frac[which.max(Bp)]
    ) %>% 
    dplyr::arrange(desc(max_Bp), desc(frac_of_max_Bp)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
        Chr = forcats::fct_inorder(Chr),
        Bp = factor(Bp)
    ) %>% 
    suppressMessages() %>%
    ggplot(aes(x = Chr, y = count, fill = Bp)) + 
    geom_bar(position = "fill", stat = "identity") + 
    scale_fill_manual(values = gradient_colors) + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    ) + 
    labs(
        x = "", 
        y = "Fraction, %",
        fill = "Break persistence"
    )

ggsave(
    filename = paste0(
        "../figures/break_persistence/", 
        "BpFractionPerChromosome.pdf"
    ),
    plot = p1,
    height = 6, width = 8
)

# Q: which breakage type have the highest proportion of high Bps?
p1 <- BreakPersistence %>% 
    dplyr::group_by(BreakType, Bp) %>% 
    dplyr::summarise(count = dplyr::n()) %>% 
    dplyr::mutate(frac = count/sum(count)*100) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
        BreakType = forcats::fct_inorder(BreakType),
        Bp = factor(Bp)
    ) %>% 
    suppressMessages() %>% 
    ggplot(aes(x = BreakType, y = count, fill = Bp)) + 
    geom_bar(position = "fill", stat = "identity") + 
    scale_fill_manual(values = gradient_colors) + 
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) + 
    labs(x = "", y = "Fraction, %", fill = "Break persistence")

ggsave(
    filename = paste0(
        "../figures/break_persistence/", 
        "BpFractionPerBreaktype.pdf"
    ),
    plot = p1,
    height = 6, width = 8
)

#' filter for the top 95% of tissue cancer IDs
TopTC <- BreakPersistence %>% 
    dplyr::group_by(TC_ID) %>% 
    dplyr::summarise(count = dplyr::n()) %>% 
    dplyr::arrange(desc(count)) %>% 
    dplyr::mutate(
        frac = count/sum(count),
        cum.frac = cumsum(frac)
    ) %>% 
    dplyr::filter(cum.frac <= 0.95) %>% 
    dplyr::arrange(desc(frac))

fwrite(
    TopTC %>% 
        dplyr::select(TC_ID, count) %>% 
        as.data.table(),
    file = paste0(
        "../data/BaseTable/",
        "TopTC_Table.csv")
)

TopTC_ID <- dplyr::pull(TopTC, TC_ID)
saveRDS(
    object = TopTC_ID, 
    file = paste0(
        "../data/break_persistence/",
        "Top_Tissue_Cancer_ID.RData"
))

to.save <- TopTC %>% 
    dplyr::select(TC_ID, frac) %>% 
    dplyr::mutate(frac = frac*100) %>% 
    as.data.table()

fwrite(
    to.save, 
    file = paste0(
        "../data/break_persistence/",
        "Top_Tissue_Cancer_ID_withFrac.csv"
    )
)

# filter by the top 95% of tissue cancer IDs
BreakPersistence_TopTC <- BreakPersistence[which(
    !is.na(match(BreakPersistence$TC_ID, TopTC_ID))
),]

df.for.p1 <- BreakPersistence %>% 
    dplyr::group_by(TC_ID, Bp) %>% 
    dplyr::summarise(count = dplyr::n()) %>% 
    dplyr::group_by(TC_ID) %>%
    dplyr::mutate(
        TC_ID.total = sum(count),
        frac = count/TC_ID.total*100
    ) %>% 
    dplyr::group_by(TC_ID) %>% 
    dplyr::mutate(
        max_Bp = max(Bp),
        frac_of_max_Bp = frac[which.max(Bp)]
    ) %>% 
    dplyr::arrange(desc(frac_of_max_Bp)) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(max_Bp >= 20) %>% 
    suppressMessages()

p1 <- df.for.p1 %>% 
    dplyr::mutate(
        TC_ID = forcats::fct_inorder(TC_ID),
        Bp = factor(Bp)
    ) %>% 
    ggplot(aes(x = TC_ID, y = count, fill = Bp)) + 
    geom_bar(position = "fill", stat = "identity") + 
    scale_fill_manual(values = gradient_colors) + 
    scale_x_discrete(labels = function(x){
        stringr::str_wrap(
            string = x, 
            width = 35, 
            whitespace_only = FALSE
        )
    }) +
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)
    ) + 
    labs(
        title = paste0(
            "Sorted by decreasing order of tissue:cancer ",
            "ID with highest fraction of max. Bp"
        ),
        subtitle = "Cutoff max Bp: 20",
        x = "", y = "Fraction, %", fill = "Break persistence"
    )

ggsave(
    filename = paste0(
        "../figures/break_persistence/", 
        "BpFractionPerTissueCancer_BpCutOff20.pdf"
    ),
    plot = p1,
    height = 9, width = 18
)

p1 <- BreakPersistence_TopTC %>% 
    dplyr::group_by(TC_ID, Bp) %>% 
    dplyr::summarise(count = dplyr::n()) %>% 
    dplyr::group_by(TC_ID) %>%
    dplyr::mutate(
        TC_ID.total = sum(count),
        frac = count/TC_ID.total*100
    ) %>% 
    dplyr::group_by(TC_ID) %>% 
    dplyr::mutate(
        max_Bp = max(Bp),
        frac_of_max_Bp = frac[which.max(Bp)]
    ) %>% 
    dplyr::arrange(desc(max_Bp), desc(frac_of_max_Bp)) %>% 
    dplyr::ungroup() %>%     
    dplyr::mutate(
        TC_ID = forcats::fct_inorder(TC_ID),
        Bp = factor(Bp)
    ) %>% 
    suppressMessages() %>% 
    ggplot(aes(x = TC_ID, y = count, fill = Bp)) + 
    geom_bar(position = "fill", stat = "identity") + 
    scale_fill_manual(values = gradient_colors) + 
    scale_x_discrete(labels = function(x){
        stringr::str_wrap(
            string = x, 
            width = 35, 
            whitespace_only = FALSE
        )
    }) +
    theme_bw() + 
    theme_classic() + 
    theme(
        text = element_text(size = 15),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5)
    ) + 
    labs(x = "", y = "Fraction, %", fill = "Break persistence")

ggsave(
    filename = paste0(
        "../figures/break_persistence/", 
        "BpFractionPerTopTissueCancer.pdf"
    ),
    plot = p1,
    height = 9, width = 12
)

# Q: what percentage do the top 95% TCs contribute to the total number of TCs of a particular Bp?
TC_BP_overlaps <- sapply(1:max(BreakPersistence$Bp), function(BP){
    Unique_TCs_Per_Bp <- BreakPersistence %>% 
        dplyr::filter(Bp == BP) %>% 
        dplyr::summarise(TC_unique = unique(TC_ID)) %>% 
        dplyr::pull(TC_unique)

    Unique_TopTCs_Per_Bp <- BreakPersistence_TopTC %>% 
        dplyr::filter(Bp == BP) %>% 
        dplyr::summarise(TC_unique = unique(TC_ID)) %>% 
        dplyr::pull(TC_unique)

    shared.TCs <- intersect(Unique_TCs_Per_Bp, Unique_TopTCs_Per_Bp)
    fraction.TCs <- length(shared.TCs)/length(Unique_TCs_Per_Bp)
    return(fraction.TCs)
})

p1 <- tibble(TC_BP_overlaps) %>% 
    dplyr::mutate(BP = 1:nrow(.), .before = 1) %>% 
    tidyr::drop_na() %>% 
    ggplot(aes(x = BP, y = TC_BP_overlaps*100)) +
    geom_line(linewidth = 1.3) + 
    coord_cartesian(ylim = c(0, 100)) + 
    theme_bw() + 
    theme_classic() + 
    theme(text = element_text(size = 15)) + 
    labs(
        title = paste0(
            "Contribution of the top 95% occurring tissue:cancers (TCs) \n",
            "towards all TCs of a given Bp"
        ),
        x = "Breakage persistence, Bp",
        y = "Contribution, %"
    )

ggsave(
    filename = paste0(
        "../figures/break_persistence/", 
        "TopTCsContributionToAllTCsAtBp.pdf"
    ),
    plot = p1,
    height = 6, width = 8
)