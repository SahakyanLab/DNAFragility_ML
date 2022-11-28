# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
break_score <- as.character(args[2])
k <- as.numeric(args[3])
cols <- as.character(args[4])

# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(gtools)))
suppressPackageStartupMessages(suppressWarnings(library(matrixStats)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(ggcorrplot)))

# source functions
# setwd("/Users/paddy/Documents/DPhil/04_DNAFragility/scripts/")
setwd(my.path)
source("../lib/KmerTable.R")
source("../lib/Breakpoints.R")
source("../lib/Features.R")

# global variables
statistic = "mean"
seed = 1234
# k = 8; exp = "DSBCapture"; break_score = "all"; cols = "all"

# extract features
features <- Features$new(k = k, exp = exp, 
                         seed = seed, 
                         break_score = break_score,
                         assembly = "hg19",
                         scores_with_kmers = FALSE)
features$get_features_from_csv()
features$select_columns(cols = cols)

# standardise
features$feature_matrix <- features$feature_matrix[, -"predictor"]
dataset <- apply(features$feature_matrix, 2, scale, 
                 center = TRUE, scale = TRUE)

# correlation heatmap
dataset.cor <- cor(dataset, use = "complete.obs")

# with all labels
p1 <- ggcorrplot::ggcorrplot(dataset.cor) + 
    labs(
        title = "Correlation heatmap of all features from full data set",
        subtitle = paste0("Exp: ", exp, 
                          ". Features: ", cols,
                          ". Breakage: ", k, "-mer.")
    )

dir.create("../figures/original_dataset/correlation/", showWarnings = FALSE)
ggsave(
  filename = paste0("../figures/original_dataset/correlation/",
                    exp, "_", cols, 
                    "_kmer-", k, "_", 
                    break_score,
                    "_all-features_with-labels.pdf"),
  plot = p1,
  height = 25,
  width = 40
)

# with all labels
p1 <- ggcorrplot::ggcorrplot(dataset.cor) + 
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
    ) +
    labs(
        title = "Correlation heatmap of all features from full data set",
        subtitle = paste0("Exp: ", exp, 
                          ". Features: ", cols,
                          ". Breakage: ", k, "-mer.")
    )

ggsave(
  filename = paste0("../figures/original_dataset/correlation/",
                    exp, "_", cols, 
                    "_kmer-", k, "_", 
                    break_score,
                    "_all-features_without-labels.pdf"),
  plot = p1
)

# filter for most pos/neg correlated features
dataset.cor.filtered <- dataset.cor
dataset.cor.filtered[(dataset.cor.filtered < 0.8 & dataset.cor.filtered > -0.8) |
                     (dataset.cor.filtered <= -0.99 | dataset.cor.filtered >= 0.99)] <- 0

p1 <- ggcorrplot::ggcorrplot(dataset.cor.filtered) + 
    labs(
        title = "|Correlation| >= 0.8 heatmap of all features from full data set",
        subtitle = paste0("Exp: ", exp, 
                          ". Features: ", cols,
                          ". Breakage: ", k, "-mer.")
    )

ggsave(
  filename = paste0("../figures/original_dataset/correlation/",
                    exp, "_", cols, 
                    "_kmer-", k, "_", 
                    break_score,
                    "_all-features_highest-correlations",
                    "_with-labels.pdf"),
  plot = p1,
  height = 25,
  width = 40
)

p1 <- ggcorrplot::ggcorrplot(dataset.cor.filtered) + 
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank()
    ) +
    labs(
        title = "|Correlation| >= 0.8 heatmap of all features from full data set",
        subtitle = paste0("Exp: ", exp, 
                          ". Features: ", cols,
                          ". Breakage: ", k, "-mer.")
    )

ggsave(
  filename = paste0("../figures/original_dataset/correlation/",
                    exp, "_", cols, 
                    "_kmer-", k, "_", 
                    break_score,
                    "_all-features_highest-correlations",
                    "_without-labels.pdf"),
  plot = p1
)

# as_tibble(dataset.cor.filtered) %>% 
#     dplyr::mutate(
#         feature_x = rownames(dataset.cor.filtered),
#         .before = 1) %>% 
#     tidyr::gather(feature_y, correlation, -feature_x) %>% 
#     dplyr::filter(correlation != 0) %>% 
#     dplyr::arrange(desc(correlation))

# # Compute a matrix of correlation p-values
# p.mat <- ggcorrplot::cor_pmat(dataset)

# p1 <- ggcorrplot::ggcorrplot(
#     dataset.cor, 
#     type = "lower", 
#     p.mat = p.mat) + 
#     labs(
#         title = "Correlation heatmap w/signif level of all features from full data set",
#         subtitle = paste0("Exp: ", exp, 
#                           ". Features: ", cols,
#                           ". Breakage: ", k, "-mer.")
#     )

# ggsave(
#   filename = paste0("../figures/original_dataset/correlation/",
#                     exp, "_", cols, 
#                     "_kmer-", k, "_", 
#                     break_score,
#                     "_all-features_signif-level_with-labels.pdf"),
#   plot = p1,
#   height = 25,
#   width = 40
# )

# correlation heatmap with hierarchical clustering
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer)))
suppressPackageStartupMessages(suppressWarnings(library(gplots)))
suppressPackageStartupMessages(suppressWarnings(library(cluster)))

dataset.dist <- dist(dataset.cor)
# methods to assess
m <- c(
  "average" = "average", 
  "single" = "single", 
  "complete" = "complete", 
  "ward" = "ward"
)
# function to compute coefficient
ac <- function(x){
  agnes(dataset.dist, method = x)$ac
}
# get agglomerative coefficient for each linkage method
ac.results <- purrr::map_dbl(m, ac)
# select method with highest coefficient
ac.max <- names(which.max(ac.results))
hclust.m <- c(
    "average" = "average",
    "single" = "single",
    "complete" = "complete",
    "ward" = "ward.D2"
)

# function to compute coefficient
hclust.co.dist <- function(x){
  #' CC compares the actual pairwise distances of 
  #' all your samples to those implied by the 
  #' hierarchical clustering. Thus, values closer
  #' to 1 suggest clustering preserves original
  #' distances.
  h.res <- hclust(dataset.dist, method = x)
  c.res <- cophenetic(h.res)
  return(cor(dataset.dist, c.res))
}
# get agglomerative coefficient for each linkage method
co.dist.results <- purrr::map_dbl(hclust.m, hclust.co.dist)
co.dist.max <- names(which.max(co.dist.results))
# dataset.tree <- hclust(dataset.dist, method = unname(hclust.m[ac.max]))
dataset.tree <- hclust(dataset.dist, method = unname(hclust.m[co.dist.max]))
dataset.dend <- as.dendrogram(dataset.tree)
color.scheme <- rev(brewer.pal(10, "RdBu"))
row.clust <- as.dendrogram(hclust(
    dist(1-t(dataset.cor)), method = unname(hclust.m[co.dist.max])
))

pdf(paste0("../figures/original_dataset/correlation/",
           exp, "_", cols, 
           "_kmer-", k, "_", 
           break_score,
           "_all-features_",
           "hierarchical-clustering.pdf"),
    height = 11, width = 15)
heatmap.2(
    dataset.cor, 
    Rowv = row.clust,
    Colv = dataset.dend,
    dendrogram = "both", 
    revC = TRUE, 
    trace = "none", 
    density.info = "none",
    col = color.scheme, 
    key = TRUE, 
    key.xlab = "Correlation Coefficient",
    symkey = TRUE,
    margins = c(10,10),
    cexRow = 0.3,
    cexCol = 0.3
)
plot.save <- dev.off()

# Construct dendorgram
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(FactoMineR))
# optimal number of clusters based on dendograms
cuts = 13
p1 <- fviz_dend(
  dataset.tree,
  k = cuts,
  horiz = TRUE,
  rect = TRUE,
  rect_fill = TRUE,
  rect_border = "jco",
  k_colors = "jco",
  cex = 0.5) %>% 
  suppressWarnings()

# cut trees
sub.groups <- cutree(dataset.tree, k = cuts)
table(sub.groups)

ggsave(
  filename = paste0("../figures/original_dataset/correlation/",
                    exp, "_", cols, 
                    "_kmer-", k, "_", 
                    break_score,
                    "_all-features_",
                    "hierarchical-clustering_",
                    cuts, "-clusters.pdf"),
  plot = p1,
  height = 15,
  width = 12
)

#####
# dataset.dist.opt <- dist(
#     dataset.cor, 
#     method = "euclidean", 
#     diag = FALSE
# )

# res <- NbClust::NbClust(
#     data = as.matrix(dataset.cor),
#     diss = dataset.dist.opt,
#     # diss = as.matrix(as.dist(dataset.cor)), 
#     distance = NULL, 
#     min.nc = 2, 
#     max.nc = 10, 
#     method = unname(hclust.m[ac.max])
#     # index = "all"
# )

# fviz_nbclust(
#     dataset.cor, 
#     hcut, 
#     method = "silhouette", 
#     k.max = 100) + 
#     theme_minimal() + 
#     ggtitle("The Silhouette Plot") + 
#     theme(axis.text.x = element_text(angle = 90))

# dend_plot <- fviz_dend(dataset.tree)
# dend_data <- attr(dend_plot, "dendrogram")
# dend_cuts <- cut(dend_data, h = 8)
# fviz_dend(dend_cuts$lower[[2]])