# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
break_score <- as.character(args[2])
k <- as.numeric(args[3])
cols <- as.character(args[4])
exp <- as.character(args[5])

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

# source functions
# setwd(my.path)
setwd("/Users/paddy/Documents/DPhil/04_DNAFragility/scripts/")
source("../lib/KmerTable.R")
source("../lib/Breakpoints.R")
source("../lib/Features.R")

# global variables
statistic = "mean"
seed = 1234
k = 8
exp = "DSBCapture"
break_score = "all"
cols = "all"

# extract features
dir.create("../data/kmertone/QueryTable/", showWarnings = FALSE)
for(score in c("ratio", "zscore")){
    table <- KmerTable$new(
        k = k, 
        break_score = score, 
        statistic = statistic, 
        exp = exp,
        all_exp = FALSE,
        group_exp = TRUE
    )
    table$generate_querytable()
}

# private=self=NULL
# self$k <- k
# private$bp_dir <- paste0("../data/experiments/", exp)
# self$seed <- seed
# assembly = "hg19"
# private$break_score=break_score
# private$scores_with_kmers = FALSE

features <- Features$new(k = k, exp = exp, 
                         seed = seed, 
                         break_score = break_score,
                         assembly = "hg19",
                         scores_with_kmers = FALSE)

######################################################################################################################################################
features$get_features(
    FEAT_G4_REGEX = TRUE, g4_type = "GPQS", 
    FEAT_GC_COUNT = TRUE,
    FEAT_KMER_COUNTS = TRUE, kmer_window = 3,
    FEAT_VIENNA_RNA = TRUE, sliding_window = NULL, 
    nuc_type = "DNA", 
    RNAfold.CALL = "/home/imm/hert6114/anaconda3/bin/RNAfold",
    FEAT_DNA_SHAPE = TRUE,
    FEAT_TFBS_EUCLIDEAN_DISTANCE = FALSE,
    FEAT_OCCUPANCY_SCORES = TRUE,
    FEAT_QM_PARAMETERS = TRUE,  
    SAVE_OUTPUT = TRUE
)
features$feature_matrix
training <- apply(features$feature_matrix[, -"predictor"], 2, 
                  scale, center = TRUE, scale = TRUE)
training <- as.data.table(training)
training$predictor <- features$feature_matrix[, "predictor"]
my_pca <- prcomp(
  training[, -"predictor"], 
  retx = TRUE, 
  center = FALSE, 
  scale. = FALSE
)
pc_eigen <- my_pca$sdev^2
pc_cumvar <- tibble(
  pc = 1:length(pc_eigen),
  variance = pc_eigen) %>% 
  dplyr::mutate(
    pct = variance/sum(variance),
    pct_cum = cumsum(pct))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(FactoMineR))
p1 <- fviz_pca_ind(
  my_pca,
  alpha = 0.1,
  geom.ind = "point",
  col.ind = training$predictor, 
  palette = c("#00AFBB", "#FC4E07"),
  addEllipses = TRUE,
  legend.title = "Breakpoints?") + 
  labs(subtitle = paste0("Exp: ", exp, 
                         ". Features: ", "Only breaks",
                         ". Breakage: ", k, "-mer.")
  )

ggsave(
  filename = paste0("../figures/PCA/test_all.pdf"),
  plot = p1,
  height = 8,
  width = 10
)
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################

features$get_features_from_csv()
features$select_columns(cols = cols)

# standardise
training <- apply(features$feature_matrix[, -"predictor"], 2, 
                  scale, center = TRUE, scale = TRUE)
training <- as.data.table(training)
training$predictor <- features$feature_matrix[, "predictor"]

# features$get_features_from_csv()
# features$select_columns(cols = cols)
# features$train_test_split(train.split = 0.7)
# features$standardise(cols = cols)
# training <- features$train_matrix
# testing <- features$test_matrix
# training[, predictor := as.factor(predictor)]
# testing[, predictor := as.factor(predictor)]

# training=training[c(1:10000, (nrow(training)-10000):(nrow(training))),]

# compute the PCA of Data
my_pca <- prcomp(
  training[, -"predictor"], 
  retx = TRUE, 
  center = FALSE, 
  scale. = FALSE
)

# Cumulative and proportion of variance explained (CVE, PVE)
pc_eigen <- my_pca$sdev^2
pc_cumvar <- tibble(
  pc = 1:length(pc_eigen),
  variance = pc_eigen) %>% 
  dplyr::mutate(
    pct = variance/sum(variance),
    pct_cum = cumsum(pct))
  
p1 <- pc_cumvar %>% 
  dplyr::select(-variance) %>% 
  dplyr::rename_with(~c("PC", "PVE", "CVE")) %>% 
  tidyr::gather(key, value, -PC) %>% 
  ggplot(aes(x = PC, y = value)) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(~key, ncol = 1, scales = "free_y") + 
  labs(y = "Variance Explained")

dir.create("../figures/PCA/", showWarnings = FALSE)
ggsave(
    filename = paste0("../figures/PCA/",
                      exp, "_", cols, 
                      "_kmer-", k, "_", 
                      break_score,
                      "_CVE.pdf"), 
    plot = p1,
    width = 9,
    height = 8
)

# CVE threshold
max.pcs <- which(pc_cumvar$pct_cum <= 0.8)
my_pca.training <- as.data.table(my_pca$x)[, 1:max(max.pcs)]
# my_pca.training <- cbind(training[, "predictor"], my_pca.rotated)

# Eigenvalue threshold
# max.pcs <- which(pc_eigen >= 1)
# my_pca.rotated <- as.data.table(my_pca$x)[, 1:max(max.pcs)]
# my_pca.training <- cbind(training[, "predictor"], my_pca.rotated)

###################################################################
# Feature loadings
################################################################### 
feature.loading <- as_tibble(my_pca$rotation[, 1:max(max.pcs)]) %>%
  dplyr::mutate(
    feature = rownames(my_pca$rotation[, 1:max(max.pcs)]),
    .before = 1) %>% 
  tidyr::gather(key, value, -feature) %>% 
  dplyr::group_by(key) %>% 
  dplyr::arrange(value, .by_group = TRUE)
  
loading.key <- stringr::str_sort(unique(feature.loading$key), numeric = TRUE)
p1 <- pblapply(loading.key, function(pca){
  p1 <- feature.loading %>% 
    dplyr::filter(key == pca) %>% 
    dplyr::slice_tail(n = 60) %>% 
    dplyr::mutate(feature = forcats::fct_inorder(feature)) %>% 
    ggplot(aes(x = value, y = feature)) + 
    geom_segment(aes(xend = min(value), yend = feature), col = "darkgrey") + 
    geom_point(col = "orange") + 
    theme_bw() + 
    labs(
      x = "",
      y = "",
      title = paste0("Feature Loading ", pca)
    ) + 
    theme(text = element_text(size = 10))
  
  return(p1)
})

ggsave(
    filename = paste0("../figures/PCA/",
                      exp, "_", cols, 
                      "_kmer-", k, "_", 
                      break_score,
                      "_feature_loadings_PCs.pdf"), 
  plot = do.call(gridExtra::grid.arrange, c(p1, ncol = 5)),
  height = 47,
  width = 23
)

###################################################################
# Loading plots by PCA
###################################################################
# make sure to install without source compilation
# install.packages("carData")
# install.packages("lme4",
#                  repos=c("http://lme4.r-forge.r-project.org/repos",
#                          getOption("repos")[["CRAN"]]))
# install.packages("car")
# install.packages("factoextra")
# install.packages("FactoMineR")

suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(FactoMineR))

# group by variable category
col.groups <- readRDS("../data/feature_matrix/DSBCapture_kmer-8_all-features_COLNAMES.RData")
col.hex <- c("#e81123", "#68217a", "#00bcf2", "#009e49", "#ff8c00",
             "#ea00ff", "#8800ff")

# PC1/PC2 plot. Groupings: breaks, biophysics, protein, shape, counts
only.breaks.names <- unlist(c(col.groups$BREAKS), use.names = FALSE)
only.breaks <- rep(1, length(only.breaks.names))
names(only.breaks) <- only.breaks.names

biophysics.names <- unlist(c(col.groups$QM_PARAMETERS), use.names = FALSE)
biophysics <- rep(2, length(biophysics.names))
names(biophysics) <- biophysics.names

# protein.names <- unlist(c(col.groups$OCCUPANCY_SCORES), use.names = FALSE)
# protein <- rep(3, length(protein.names))
# names(protein) <- protein

shape.names <- unlist(c(
  col.groups$DNA_SHAPE,
  col.groups$VIENNA_RNA), 
  use.names = FALSE)
shape <- rep(4, length(shape.names))
names(shape) <- shape.names

counts.names <- unlist(c(
  col.groups$G4_REGEX,
  col.groups$SINGLETON,
  col.groups$GC_COUNT,
  col.groups$KMER_COUNTS), 
  use.names = FALSE)
counts <- rep(5, length(counts.names))
names(counts) <- counts.names

smaller.groups <- c(only.breaks, biophysics, shape, counts)
smaller.groups <- as.factor(smaller.groups)

p1 <- fviz_pca_var(
  my_pca,
  repel = TRUE,
  col.var = smaller.groups,
  alpha = 0.65,
  addEllipses = TRUE, 
  ellipse.level = 0.95,
  habillage = "none",
  geom.var = c("arrow")) + 
  scale_fill_discrete(guide = "none") + 
  scale_color_manual(
    name = "Feature Groups",
    labels = c("Only breakages", "Biophysical params",
              #  "Protein occupancy", 
               "DNA shape",
               "Counts"),
    values = col.hex) + 
  labs(subtitle = paste0("Exp: ", exp, 
                         ". Features: ", cols, 
                         ". Breakage: ", k, "-mer.")
  )

ggsave(
  filename = paste0("../figures/PCA/",
                    exp, "_", cols, 
                    "_kmer-", k, "_", 
                    break_score,
                    "_variable-graph_all",
                    "-features_PC1-PC2.pdf"),
  plot = p1,
  height = 8,
  width = 12
)

# PCA by breakage scores
my_pca_breaks <- prcomp(
  training[, ..only.breaks.names], 
  retx = TRUE, 
  center = FALSE, 
  scale. = FALSE
)

breaks.labels <- ifelse(
  grepl(pattern = "Biological", x = only.breaks.names), 1,
  ifelse(grepl(pattern = "Cell_free", x = only.breaks.names), 2,
  ifelse(grepl(pattern = "Mechanical", x = only.breaks.names), 3,
  ifelse(grepl(pattern = "Natural", x = only.breaks.names), 4,
  ifelse(grepl(pattern = "Enzymatic", x = only.breaks.names), 5, 6)))))
names(breaks.labels) <- only.breaks.names
breaks.labels <- as.factor(breaks.labels)

p1 <- fviz_pca_var(
  my_pca_breaks,
  repel = TRUE,
  col.var = breaks.labels,
  alpha = 0.65,
  habillage = "none",
  geom.var = c("arrow", "text")) + 
  scale_fill_discrete(guide = "none") + 
  scale_color_manual(
    name = "Feature Groups",
    labels = c("Biological", "Cell free DNA",
               "Mechanical", "Ancient DNA",
               "Enzymatic"),
    values = col.hex) + 
  labs(subtitle = paste0("Exp: ", exp, 
                         ". Features: ", "Only breaks", 
                         ". Breakage: ", k, "-mer.")
  )

ggsave(
  filename = paste0("../figures/PCA/",
                    exp, "_", cols, 
                    "_kmer-", k, "_", 
                    break_score,
                    "_variable-graph_only-breaks",
                    "-PC1-PC2.pdf"),
  plot = p1,
  height = 12,
  width = 17
)

# PCA by biophysical parameters
my_pca_biophysics <- prcomp(
  training[, ..biophysics.names], 
  retx = TRUE, 
  center = FALSE, 
  scale. = FALSE
)

biophysics.labels <- ifelse(
  grepl(pattern = "A_", x = biophysics.names), 1,
  ifelse(grepl(pattern = "B_", x = biophysics.names), 2,
  ifelse(grepl(pattern = "Z_", x = biophysics.names), 3, 4)))
names(biophysics.labels) <- biophysics.names
biophysics.labels <- as.factor(biophysics.labels)

p1 <- fviz_pca_var(
  my_pca_biophysics,
  repel = TRUE,
  col.var = biophysics.labels,
  alpha = 0.65,
  habillage = "none",
  geom.var = c("arrow", "text")) + 
  scale_fill_discrete(guide = "none") + 
  scale_color_manual(
    name = "Feature Groups",
    labels = c("A_DNA", "B_DNA", "Z_DNA"),
    values = col.hex) + 
  labs(subtitle = paste0("Exp: ", exp, 
                         ". Features: ", "Only biophysical params", 
                         ". Breakage: ", k, "-mer.")
  )

ggsave(
  filename = paste0("../figures/PCA/",
                    exp, "_", cols, 
                    "_kmer-", k, "_", 
                    break_score,
                    "_variable-graph_biophysical-params",
                    "-PC1-PC2.pdf"),
  plot = p1,
  height = 10,
  width = 14
)

# PCA by DNA shape
my_pca_shape <- prcomp(
  training[, ..shape.names], 
  retx = TRUE, 
  center = FALSE, 
  scale. = FALSE
)

shape.labels <- ifelse(grepl(pattern = "viennaRNA", x = shape.names), 1, 2)
names(shape.labels) <- shape.names
shape.labels <- as.factor(shape.labels)

p1 <- fviz_pca_var(
  my_pca_shape,
  repel = TRUE,
  col.var = shape.labels,
  alpha = 0.65,
  habillage = "none",
  geom.var = c("arrow", "text")) + 
  scale_fill_discrete(guide = "none") + 
  scale_color_manual(
    name = "Feature Groups",
    labels = c("viennaRNA", "DNA shape"),
    values = col.hex) + 
  labs(subtitle = paste0("Exp: ", exp, 
                         ". Features: ", "Only DNA shape", 
                         ". Breakage: ", k, "-mer.")
  )

ggsave(
  filename = paste0("../figures/PCA/",
                    exp, "_", cols, 
                    "_kmer-", k, "_", 
                    break_score,
                    "_variable-graph_dna-shape",
                    "-PC1-PC2.pdf"),
  plot = p1,
  height = 10,
  width = 14
)

# # PCA by protein occupancy
# my_pca_protein <- prcomp(
#   training[, ..protein.names], 
#   retx = TRUE, 
#   center = FALSE, 
#   scale. = FALSE
# )

# protein.labels <- ifelse(
#   grepl(pattern = "DNA_", x = protein.names), 1,
#   ifelse(grepl(pattern = "Glioblastoma_", x = protein.names), 2,
#   ifelse(grepl(pattern = "PIASY|SUMO|UBC", x = protein.names), 3,
#   ifelse(grepl(pattern = "XRCC5|RAD51", x = protein.names), 4, 
#   ifelse(grepl(pattern = "RNA_", x = protein.names), 5,
#   ifelse(grepl(pattern = "TOP2B", x = protein.names), 6, 7))))))
# names(protein.labels) <- protein.names
# protein.labels <- as.factor(protein.labels)

# p1 <- fviz_pca_var(
#   my_pca_protein,
#   repel = TRUE,
#   col.var = protein.labels,
#   alpha = 0.65,
#   habillage = "none",
#   geom.var = c("arrow", "text")) + 
#   scale_fill_discrete(guide = "none") + 
#   scale_color_manual(
#     name = "Feature Groups",
#     labels = c("DNA methylation", "GBM TIC", "SUMO machinery",
#                "DNA repair", "RNA_PolII", "Topoisomerase_II"),
#     values = col.hex) + 
#   labs(subtitle = paste0("Exp: ", exp, 
#                          ". Features: ", "Only protein occupancy", 
#                          ". Breakage: ", k, "-mer.")
#   )

# ggsave(
#   filename = paste0("../figures/PCA/",
#                     exp, "_", cols, 
#                     "_kmer-", k, "_", 
#                     break_score,
#                     "_variable-graph_protein-occupancy",
#                     "-PC1-PC2.pdf"),
#   plot = p1,
#   height = 10,
#   width = 14
# )

# by coordinate points
p1 <- fviz_pca_ind(
  my_pca,
  # select.ind = list(contrib = length(training$predictor)*0.2),
  alpha = 0.1,
  geom.ind = "point",
  col.ind = training$predictor, 
  palette = c("#00AFBB", "#FC4E07"),
  addEllipses = TRUE,
  legend.title = "Breakpoints?") + 
  labs(subtitle = paste0("Exp: ", exp, 
                         ". Features: ", cols,
                         ". Breakage: ", k, "-mer.")
  )

ggsave(
  filename = paste0("../figures/PCA/",
                    exp, "_", cols, 
                    "_kmer-", k, "_", 
                    break_score,
                    "_coordinate-graph_all",
                    "-features_PC1-PC2.pdf"),
  plot = p1,
  height = 8,
  width = 10
)
###################################################################
# Exploring KMeans clustering on principal components
###################################################################
example <- my_pca.training[, -"predictor"]
kmeans_clustering <- kmeans(
  x = example[, 1:2],
  centers = 2,
  nstart = 5
)

# Extract cluster centers
kmeans_clustering_centers <- kmeans_clustering$centers

# Extract clusters
p1 <- example %>% 
  dplyr::mutate(
    predictor = my_pca.training$predictor,
    cluster = kmeans_clustering$cluster
  ) %>% 
  ggplot(aes(x = PC1, y = PC2, col = as.factor(cluster))) + 
  geom_point(alpha = 0.3) + 
  facet_wrap(~predictor)

ggsave(
  filename = paste0("../figures/PCA/h2o_PC_kmeans_clustering_plots_", 
                    cols, "_kmer_", k, "_", score,".png"), 
  plot = p1,
  width = 9,
  height = 9
)

# # plot pca pairs
# png("test.png")
# pairs(example, col = c(1:3)[kmeans_clustering$cluster])
# # with(iris, pairs(dat, col = c(1:3)[clus$cluster])) 
# plot.off <- dev.off()

###################################################################
# Exploring distribution of training data
###################################################################
p1 <- pblapply(2:(ncol(training)), function(x){
  temp <- as_tibble(training) %>% 
    dplyr::select(predictor, x)
  
  p1 <- temp %>% 
    dplyr::rename_with(~c("predictor", "col")) %>% 
    ggplot(aes(x = col, fill = predictor)) +
    geom_density(alpha = 0.7) + 
    facet_wrap(~predictor)+ 
    labs(
      x = ifelse(x == 2, "Range", ""),
      y = ifelse(x == 2, "Density", ""),
      title = colnames(temp)[2]
    ) + 
    theme(
      text = element_text(size = 10)
    )
  
  return(p1)
})

ggsave(
    filename = paste0("../figures/Original_dataset/training_data_density_plots_", 
                      cols, "_kmer_", k, "_", score,".pdf"), 
  plot = do.call(gridExtra::grid.arrange, c(p1, ncol = 5)),
  height = 45,
  width = 22
)

###################################################################
# Exploring potential cut-off
###################################################################
pc1_2.df <- as_tibble(my_pca.training) %>% 
  dplyr::select(c("predictor", "PC1", "PC2")) %>% 
  dplyr::mutate(Label = ifelse(PC2 > -20, "ABOVE", "BELOW"))

p1 <- pc1_2.df %>% 
  ggplot(aes(x = PC1, y = PC2, col = as.factor(Label))) + 
  geom_point(alpha = 0.2) + 
  facet_wrap(~predictor)

ggsave(
  filename = "test.png",
  plot = p1,
  width = 8, 
  height = 9
)

# Subset original feature matrix from PC1/PC2 separation
pca.mat <- as.matrix(my_pca$rotation[, 1:max(max.pcs)])
PC1_2 <- as.matrix(training[, -"predictor"]) %*% pca.mat

# PC1.rows <- na.omit(match(PC1_2, filter(pc1_2.df, Label == "BELOW")$PC1))
# PC2.rows <- na.omit(match(PC1_2, filter(pc1_2.df, Label == "BELOW")$PC2))
# PC1_2.rows <- c(PC1.rows, PC2.rows)


# # Comparison using box plots
# p1 <- as_tibble(training[PC1_2.rows, ]) %>% 
#   tidyr::gather(key, value, -predictor) %>% 
#   ggplot(aes(
#     x = as.factor(predictor),
#     y = value,
#     fill = predictor)) + 
#   ggplot2::geom_boxplot(alpha = 0.5) +
#   facet_wrap(~key, scales = "free")

# ggsave(
#     filename = "test.pdf",
#     # filename = "../figures/PCA/PC1_feature_matrix_comparison_boxplot.pdf", 
#     plot = p1,
#     height = 17,
#     width = 22
# )

# training.pc1_2 <- as_tibble(training[, predictor:PC2]) %>% 
#   dplyr::mutate(Label = ifelse(PC_2 >= 0.12*PC1, "ABOVE", "BELOW"))

# p1 <- training.pc1_2 %>% 
#   ggplot(aes(x = PC_1, y = PC_2, col = Label)) + 
#     geom_point(alpha = 0.1) + 
#     facet_wrap(~predictor)

# ggsave(
#     filename = "../figures/PCA/PC1_PC2_comparison.pdf", 
#     plot = p1
# )

# subset.pc1_2 <- training.pc1_2 %>% 
#   dplyr::filter(
#     (predictor == "YES" & Label == "ABOVE") | 
#     (predictor == "NO" & Label == "BELOW"))

# p1 <- subset.pc1_2 %>% 
#   ggplot(aes(x = PC_1, y = PC_2, col = Label)) + 
#     geom_point(alpha = 0.1) + 
#     facet_wrap(~predictor)

# ggsave(
#     filename = "../figures/PCA/PC1_PC2_comparison_filtered.pdf", 
#     plot = p1
# )

# # Subset original feature matrix from PC1/PC2 separation
# PC_1 <- feat.asmatrix %*% pca.asmatrix[, "pc1"]
# PC_1.rows <- na.omit(match(PC_1, filter(subset.pc1_2, predictor == "YES")$PC_1))
# # feat.mat[PC_1.rows, ]

# PC_2 <- feat.asmatrix %*% pca.asmatrix[, "pc2"]
# PC_2.rows <- na.omit(match(PC_2, filter(subset.pc1_2, predictor == "NO")$PC_2))
# # feat.mat[PC_2.rows, ]

# PC_1.subset <- as_tibble(feat.mat[PC_1.rows, ]) %>% 
#   dplyr::select(-1) %>% 
#   tidyr::gather("feature", "value") %>% 
#   dplyr::mutate(predictor = "YES")

# PC_2.subset <- as_tibble(feat.mat[PC_2.rows, ]) %>% 
#   dplyr::select(-1) %>% 
#   tidyr::gather("feature", "value") %>% 
#   dplyr::mutate(predictor = "NO")

# PC1_2.subset <- rbind(PC_1.subset, PC_2.subset)

# # Comparison of density plots
# p1 <- PC1_2.subset %>% 
#   ggplot(aes(
#     x = value, 
#     col = predictor, 
#     fill = predictor)) + 
#   ggplot2::geom_density(alpha = 0.2) +
#   facet_wrap(~feature, scales = "free")

# ggsave(
#     filename = "../figures/PCA/PC1_feature_matrix_comparison_density.pdf", 
#     plot = p1,
#     height = 17,
#     width = 22
# )

# # Comparison using box plots
# p1 <- PC1_2.subset %>% 
#   ggplot(aes(
#     x = as.factor(predictor),
#     y = value,
#     fill = predictor)) + 
#   ggplot2::geom_boxplot(alpha = 0.5) +
#   facet_wrap(~feature, scales = "free")

# ggsave(
#     filename = "../figures/PCA/PC1_feature_matrix_comparison_boxplot.pdf", 
#     plot = p1,
#     height = 17,
#     width = 22
# )

###################################################################
# Machine Learning models based on PCAs
################################################################### 
# setup h2o cluster
h2o.init(max_mem_size = "5g")
training.h2o <- as.h2o(my_pca.training)

testing <- features$test_matrix
testing[, predictor := as.factor(predictor)]
testing.mat <- as.matrix(testing[, -"predictor"])
pca.rotation.mat <- as.matrix(my_pca$rotation[, 1:max(max.pcs)])
testing.mat <- testing.mat %*% pca.rotation.mat
testing.mat <- data.table::as.data.table(testing.mat)
colnames(testing.mat) <- paste0("PC", 1:max(max.pcs))
testing <- cbind(testing[, "predictor"], testing.mat)
testing.h2o <- as.h2o(testing)

# Linear Regression
glm.fit <- h2o.glm(
    family = "binomial",
    x = colnames(training.h2o)[2:ncol(training.h2o)],
    y = "predictor",
    training_frame = training.h2o,
    lambda_search = TRUE,
    nfolds = 5,
    seed = 1234
)

glm.perf <- h2o.performance(glm.fit, newdata = testing.h2o)
glm.auc <- h2o.auc(glm.perf); glm.auc

# Naive Bayes
nb.fit <- h2o.naiveBayes(
    x = colnames(training.h2o)[2:ncol(training.h2o)],
    y = "predictor",
    training_frame = training.h2o,
    nfolds = 5,
    seed = 1234
)

nb.perf <- h2o.performance(nb.fit, newdata = testing.h2o)
nb.auc <- h2o.auc(nb.perf); nb.auc

# Random Forest
rf.fit <- h2o.randomForest(
    x = colnames(training.h2o)[2:ncol(training.h2o)],
    y = "predictor",
    training_frame = training.h2o,
    nfolds = 5,
    seed = 1234,
    stopping_rounds = 2,
    ntrees = 200,
    score_each_iteration = TRUE
)

rf.perf <- h2o.performance(rf.fit, newdata = testing.h2o)
rf.auc <- h2o.auc(rf.perf); rf.auc

h2o.shutdown(prompt = FALSE)

h2o.model.perf <- list(glm.perf, rf.perf, nb.perf) %>% 
    purrr::map(function(x) x %>% 
        .@metrics %>% 
        .$thresholds_and_metric_scores %>% 
        .[c('fpr', 'tpr')] %>% 
        dplyr::add_row(tpr = 0, fpr = 0, .before = FALSE) %>% 
        dplyr::add_row(tpr = 0, fpr = 0, .before = FALSE)) %>% 
    purrr::map2(c(
        paste0('Linear', ' (', signif(glm.auc, 3),')'), 
        paste0('Random Forest', ' (', signif(rf.auc, 3),')'), 
        paste0('Naive Bayes', ' (', signif(nb.auc, 3),')')),
        function(x, y) x %>% dplyr::mutate(model = y)
    ) %>% 
    rbindlist()

p1 <- h2o.model.perf %>% 
    ggplot(aes(x = fpr, y = tpr, col = model)) + 
    geom_line(size = 1.3) + 
    geom_segment(aes(
        x = 0, y = 0,
        xend = 1, yend = 1),
        linetype = 2,
        col = "darkgrey"
    ) + 
    labs(
        x = "False Positive Rate",
        y = "True Positive Rate",
        title = "ROC Curves"
    )

p1

ggsave(
    filename = paste0("../figures/PCA/h2o_plots_", cols, 
                      "_kmer_", k, "_", score,".pdf"), 
    plot = p1,
    width = 8,
    height = 8
)