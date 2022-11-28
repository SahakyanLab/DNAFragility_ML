# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(gtools)))
suppressPackageStartupMessages(suppressWarnings(library(matrixStats)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(mclust)))
suppressPackageStartupMessages(suppressWarnings(library(cluster)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
options(dplyr.summarise.inform = FALSE)

# source functions
setwd("/Users/paddy/Documents/DPhil/04_DNAFragility/scripts/")
source("../lib/KmerTable.R")
source("../lib/Breakpoints.R")
source("../lib/Features.R")

# global variables
k = 4
statistic = "mean"
exp = "DSBCapture"
seed = 1234
cols = "all"
break_score = "zscore"

kmer_table <- KmerTable$new(
    k = k, 
    break_score = break_score, 
    statistic = statistic, 
    exp = exp,
    all_exp = TRUE,
    group_exp = TRUE
)
kmer_table$generate_querytable()

kmer.df <- fread(paste0("../data/kmertone/QueryTable/",
                        "QueryTable-breaks_kmer-", 
                        k, "_", break_score, ".csv"))

# correlation coefficients
categories <- kmer.df$category
kmers <- colnames(kmer.df)
kmer.df <- apply(kmer.df[, -"category"], 2, scale, 
                 center = TRUE, scale = TRUE)

# # split into target and rest df
# temp <- as.data.table(kmer.df)
# temp[, category := categories]
# group.of.interest <- "Biological_NHEK_cells"
# df2 <- temp[category != group.of.interest]
# df1 <- temp[category == group.of.interest]
# df1 <- do.call(
#     rbind, 
#     replicate(nrow(df2), df1, simplify = FALSE)
# )
# cor(t(df1[, -"category"]), t(df2[, -"category"]))

# by experiment
rownames(kmer.df) <- categories
kmer.df.by.exp <- cor(t(kmer.df), use = "complete.obs", method = "pearson")
dist.df <- dist(kmer.df.by.exp, method = "euclidean")

# function to compute coefficient
hclust.co.dist <- function(x){
  #' CC compares the actual pairwise distances of 
  #' all your samples to those implied by the 
  #' hierarchical clustering. Thus, values closer
  #' to 1 suggest clustering preserves original
  #' distances.
  h.res <- hclust(dist.df, method = x)
  c.res <- cophenetic(h.res)
  return(cor(dist.df, c.res))
}
# get agglomerative coefficient for each linkage method
hclust.m <- c(
    "average" = "average",
    "single" = "single",
    "complete" = "complete",
    "ward" = "ward.D2"
)
co.dist.results <- purrr::map_dbl(hclust.m, hclust.co.dist)
co.dist.max <- names(which.max(co.dist.results))
dataset.tree <- hclust(dist.df, method = unname(hclust.m[co.dist.max]))

suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(FactoMineR))
# optimal number of clusters based on dendograms
cuts = 8
fviz_dend(
    dataset.tree,
    k = cuts,
    horiz = TRUE,
    rect = TRUE,
    rect_fill = TRUE,
    rect_border = "jco",
    k_colors = "jco",
    cex = 0.5) %>% 
    suppressWarnings()

# by k-mer
kmer.df <- as.data.table(kmer.df)
kmer.df[, category := categories]

cor(kmer.df, use = "complete.obs", method = "pearson")

kmer.df.by.exp <- cor(t(kmer.df), use = "complete.obs", method = "pearson")
dist.df <- dist(kmer.df.by.exp, method = "euclidean")