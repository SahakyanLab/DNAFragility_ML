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
k = 8
statistic = "mean"
exp = "DSBCapture"
seed = 1234
cols = "all"

breaks.dt <- lapply(c("ratio", "zscore"), function(break_score){
    # generate query table of k-mers
    kmer_table <- KmerTable$new(
        k = k, 
        break_score = break_score, 
        statistic = statistic, 
        exp = exp,
        all_exp = TRUE,
        group_exp = TRUE
    )
    kmer_table$generate_querytable()
    
    # extract z-scores and probability ratios of breaks
    breaks <- Breakpoints$new(
        k = k, 
        exp = exp, 
        seed = seed,
        assembly = "hg19"
    )
    breaks$get_breaks(break_score = break_score, 
                      scores_with_kmers = TRUE)
    feature.mat <- rbind(
        breaks$true_breaks[[break_score]]$scores,
        breaks$control_breaks[[break_score]]$scores
    )
    colnames(feature.mat) <- stringr::str_remove(
        string = colnames(feature.mat),
        pattern = paste0("_", break_score)
    )
    predictor <- c(
        rep("YES", dim(breaks$true_breaks[[break_score]]$scores)[1]), 
        rep("NO", dim(breaks$control_breaks[[break_score]]$scores)[1])
    )
    predictor <- as.factor(predictor)
    feature.mat <- cbind(predictor, feature.mat)
    return(feature.mat[complete.cases(feature.mat)])
})

breaks.dt.ratio <- breaks.dt[[1]]
breaks.dt.zscores <- breaks.dt[[2]]
categories <- colnames(breaks.dt.ratio)
top_percent=ifelse(k == 6, 0.025, 0.005)
break_score = "all"

dir.create("../figures/original_dataset/breakage_insights/breakage_insights/", 
           recursive = TRUE, showWarnings = FALSE)
prelim_analysis <- function(x){
    to.keep <- c("predictor", "kmer", categories[x])

    dt <- as_tibble(breaks.dt.ratio[, ..to.keep]) %>% 
        dplyr::rename_with(~c("predictor", "kmer", "ratio")) %>% 
        dplyr::mutate(
            predictor = as.character(predictor),
            ratio = log2(ratio),
            z.scores = abs(breaks.dt.zscores[[x]]))

    dt.subset <- dt %>% 
        dplyr::select(predictor, kmer) %>% 
        dplyr::group_by(predictor) %>% 
        dplyr::mutate(count = dplyr::n()) %>% 
        dplyr::group_by(predictor, kmer) %>% 
        dplyr::summarise(norm.count = dplyr::n()/count) %>% 
        dplyr::distinct(kmer, .keep_all = TRUE) %>% 
        dplyr::arrange(predictor) %>% 
        dplyr::group_by(kmer) %>% 
        dplyr::summarise(diff = base::diff(norm.count)) %>% 
        dplyr::ungroup()

    dt.subset.top <- dt.subset %>% 
        dplyr::slice_max(diff, prop = top_percent) %>% 
        dplyr::arrange(dplyr::desc(diff)) %>% 
        dplyr::mutate(label = "Top")

    dt.subset.bottom <- dt.subset %>% 
        dplyr::slice_min(diff, prop = top_percent) %>% 
        dplyr::arrange(dplyr::desc(diff)) %>% 
        dplyr::mutate(label = "Bottom")

    dt.subset.all <- rbind(
        dt.subset.top, 
        dt.subset.bottom) %>% 
        dplyr::mutate(kmer = forcats::fct_inorder(kmer))

    dt.ref <- dt %>% 
        dplyr::distinct(kmer, .keep_all = TRUE)

    dt <- dt.subset.all %>% 
        dplyr::mutate(
            kmer = as.character(kmer),
            ratio = dt.ref$ratio[match(kmer, dt.ref$kmer)],
            z.scores = dt.ref$z.scores[match(kmer, dt.ref$kmer)])

    # drop inf values
    dt <- dt[is.finite(rowSums(dt[, c("ratio", "z.scores")])),]   

    p1 <- dt %>% 
        ggplot(aes(x = ratio, y = z.scores, col = label)) + 
        geom_point(alpha = 0.5) + 
        theme_bw() + 
        theme(legend.position = "none") + 
        labs(
            x = ifelse(x == 3, "Log2 probability ratio", ""),
            y = ifelse(x == 3, "Absolute z-scores", ""),
            title = categories[x]
        )    

    # k-means clustering
    set.seed(seed)
    kmeans_clustering <- kmeans(
        x = dt[c("ratio", "z.scores")],
        centers = 2,
        nstart = 10,
    )

    kmeans_clustering <- data.frame(
        ratio = dt$ratio,
        z.scores = dt$z.scores,
        cluster = kmeans_clustering$cluster
    )

    # model-based clustering
    dt.mat <- as.matrix(dt[c("ratio", "z.scores")])

    # Apply GMM model with 3 components
    set.seed(seed)
    if(k == 6){
        model.names <- switch(categories[x],
            "Biological_recombination_rates" = "EVV",
            "Cell_free_DNA_disease" = "EVV",
            NULL
        )
    } else if(k == 8){
        model.names <- switch(categories[x],
            "Biological_HCT116_cells_rep-init" = "EEV",
            "Biological_HCT116_cells_WRN_loss" = "EEV",
            "Biological_K562_cells" = "EEV",
            "Biological_KM12_cells_WRN_loss" = "EEV",
            "Biological_neuronal_cells" = "EEV",
            "Biological_recombination_rates" = "EEE",
            "Biological_RPE1_cells_WT" = "EEV",
            "Biological_TK6_cells" = "EEV",
            "Cell_free_DNA_cancer" = "EEV",
            "Cell_free_DNA_disease" = "EEV",
            "Enzymatic_Twist_library" = "EEV",
            "Mechanical_nebulisation" = "EEV",
            "Natural_fragmentation_Neandertal" = "EEV",
            NULL
        )
    }

    # models = c(
    #     "EII", # spherical, equal volume
    #     "VII", # spherical, unequal volume
    #     "EEI", # diagonal, equal volume and shape
    #     "VEI", # diagonal, varying volume, equal shape
    #     "EVI", # diagonal, equal volume, varying shape
    #     "VVI", # diagonal, varying volume and shape
    #     "EEE", # ellipsoidal, equal volume, shape, and orientation
    #     "EVE", # ellipsoidal, equal volume and orientation
    #     "VEE", # ellipsoidal, equal shape and orientation
    #     "VVE", # ellipsoidal, equal orientation
    #     "EEV", # ellipsoidal, equal volume and equal shape
    #     "VEV", # ellipsoidal, equal shape
    #     "EVV", # ellipsoidal, equal volume
    #     "VVV"  # ellipsoidal, varying volume, shape, and orientation
    # )

    # for(model in models){
    #     dt_mc <- mclust::Mclust(
    #         dt.mat, 
    #         G = 2,
    #         verbose = FALSE,
    #         modelNames = model
    #     )
    #     plot.Mclust(dt_mc, what = "classification",
    #     xlab = model)
    # }

    dt_mc <- mclust::Mclust(
        dt.mat, 
        G = 2,
        verbose = FALSE,
        modelNames = model.names
    )

    temp <- kmeans_clustering[, c("ratio", "z.scores")]
    temp <- cbind(temp, kmeans_clustering$cluster)
    colnames(temp) <- c("ratio", "z.scores", "cluster")
    temp[, "cluster"] <- ifelse(temp[, "cluster"] == 1, 
                                "#d32321", "#1c87ee")

    # Plot results
    pdf(
        paste0("../figures/original_dataset/breakage_insights/",
               exp, "_", cols, 
               "_kmer-", k, "_",
               "scores-", break_score, "_", 
               "ratio_zscores_mclust_",
               categories[x], ".pdf"),
        width = 21, height = 7
    )
    par(mfrow = c(1,3))
    plot.Mclust(dt_mc, what = "uncertainty", main = TRUE)
    plot.Mclust(dt_mc, what = "classification", 
                main = TRUE, xlab = "", ylab = "")
    base::plot(x = temp[, "ratio"], 
        y = temp[, "z.scores"], 
        col = temp[, "cluster"],
        type = "p",
        pch = 18,
        cex = 1.5,
        main = "K-Means Clustering",
        xlab = "", ylab = "")
    dev.off()

    # euclidean distance between cluster centres
    euc.dist <- dist(dt_mc$parameters$mean, method = "euclidean")

    # euclidean distance of the silhouttes of clusters
    dt_mc.sil <- cluster::silhouette(
        x = dt_mc$classification,
        dist = dist(dt.mat, method = "euclidean")
    )
    avg.sil <- mean(dt_mc.sil[, "sil_width"])
    cluster.sum <- tibble(
        exp = categories[x],
        cluster.model = dt_mc$modelName,
        euc.dist = as.numeric(euc.dist),
        avg.sil = avg.sil)

    dt <- dt %>% 
        dplyr::mutate(
            exp = categories[x],
            cluster = dt_mc$classification
        )
    return(list(p1, dt, cluster.sum))
}

results <- pbapply::pblapply(3:length(categories), function(ind){
    prelim_analysis(x = ind)
})

# clustering results
dt.clusters <- rbindlist(sapply(results, `[`, 3))
dir.create("../data/kmertone/Clustering/", showWarnings = FALSE)
fwrite(
    dt.clusters,
    file = paste0("../data/kmertone/Clustering/",
                  exp, "_", cols, 
                  "_kmer-", k, ".csv")
)

# Out of the k-mers that are intrinsically susceptible to breakages,
# how many were actually broken as a percentage
res <- rbindlist(sapply(results, `[`, 2))

res.subset <- as_tibble(res) %>% 
    dplyr::group_by(exp) %>% 
    dplyr::mutate(
        cluster = ifelse(cluster == 1, "Top", "Bottom")) %>% 
    dplyr::select(c(kmer, ratio, z.scores, exp))

unique.exp.names <- unique(res.subset$exp[res.subset$exp != "Biological_NHEK_cells"])
res.dsbcapture <- res.subset %>% 
    dplyr::filter(exp == "Biological_NHEK_cells") %>% 
    dplyr::arrange(kmer)

cor.res <- lapply(unique.exp.names, function(x){
    res.combined <- left_join(
        res.dsbcapture,
        res %>% dplyr::filter(exp == x), by = "kmer"
    )

    #' multiple correlation coefficients by multiple linear regression.
    #' Then, calculate the coefficient of correlation between the 
    #' predicted and observed values of the dependent variable.
    
    # predicting probability ratios
    res.model.ratio <- lm(ratio.y ~ z.scores.y+ratio.x+z.scores.x, data = res.combined)
    cor.coef.ratio <- stats::cor(
        res.model.ratio$model$ratio.y, res.model.ratio$fitted.values, 
        method = "pearson", 
        use = "complete.obs"
    )

    # predicting z-scores
    res.model.zscores <- lm(z.scores.y ~ ratio.y+ratio.x+z.scores.x, data = res.combined)
    cor.coef.z.scores <- stats::cor(
        res.model.zscores$model$z.scores.y, res.model.zscores$fitted.values, 
        method = "pearson", 
        use = "complete.obs"
    )

    out <- data.table(
        exp = x, 
        cor.ratio = cor.coef.ratio,
        ratio.coef.intercept = res.model.ratio$coefficients["(Intercept)"],
        ratio.coef.z.scores.y = res.model.ratio$coefficients["z.scores.y"],
        ratio.coef.z.scores.x = res.model.ratio$coefficients["z.scores.x"],
        ratio.coef.ratio.x = res.model.ratio$coefficients["ratio.x"],        
        cor.z.scores = cor.coef.z.scores,
        z.scores.coef.intercept = res.model.zscores$coefficients["(Intercept)"],
        z.scores.coef.ratio.y = res.model.zscores$coefficients["ratio.y"],
        z.scores.coef.ratio.x = res.model.zscores$coefficients["ratio.x"],
        z.scores.coef.z.scores.x = res.model.zscores$coefficients["z.scores.x"],
        cor.avg = mean(c(cor.coef.ratio, cor.coef.z.scores))
    )
    return(out)
})

cor.res <- rbindlist(cor.res)
cor.res <- as_tibble(cor.res) %>% 
    dplyr::arrange(desc(cor.avg))

fwrite(
    cor.res,
    file = paste0("../data/kmertone/Clustering/",
                  exp, "_", cols, 
                  "_kmer-", k, "_",
                  "scores-", break_score,
                  "_multiple-correlation-coefficient.csv")
)

p1 <- cor.res %>% 
    dplyr::select(c(exp, cor.avg)) %>% 
    dplyr::arrange(cor.avg) %>% 
    dplyr::mutate(exp = forcats::fct_inorder(exp)) %>% 
    tidyr::gather(key, value, -exp) %>% 
    ggplot(aes(x = value, y = exp)) + 
    geom_segment(aes(xend = 0, yend = exp), col = "darkgrey") +
    geom_point(col = "orange") + 
    coord_cartesian(xlim = c(0, 1)) + 
    theme_bw() +
    labs(
        x = "Correlation Coefficients",
        y = "",
        title = paste0("Multiple correlation coefficient of z-scores ", 
                       "and probability ratios vs. DSBCapture"),
        subtitle = paste0("Method: multiple linear regression based on ", 
                          k, "-mers")
    )

ggsave(
    filename = paste0("../figures/original_dataset/breakage_insights/",
                      exp, "_", cols, 
                      "_kmer-", k, "_",
                      "scores-", break_score,
                      "_multiple_correlation.pdf"),
    plot = p1,
    height = 8,
    width = 10
)

res.filtered <- as_tibble(res) %>% 
    dplyr::filter(exp != "Biological_NHEK_cells") %>% 
    dplyr::group_by(exp) %>% 
    dplyr::mutate(
        cluster = ifelse(cluster == 1, "Top", "Bottom")) %>% 
    dplyr::summarise(
        broken.vs.enriched = sum(label == "Top" & cluster == "Top")/
                            (sum(label == "Top" & cluster == "Top")+
                             sum(label == "Bottom" & cluster == "Top")),
        not.broken.vs.depleted = sum(label == "Bottom" & cluster == "Bottom")/
                                (sum(label == "Bottom" & cluster == "Bottom")+
                                 sum(label == "Top" & cluster == "Bottom")),
        accuracy = mean(label == cluster)) %>% 
    dplyr::arrange(desc(accuracy))

p1 <- res.filtered %>% 
    dplyr::rename_with(~c(
        "exp", "Broken vs enriched", "Not broken vs depleted", "Accuracy"
    )) %>%
    dplyr::arrange(desc(Accuracy)) %>% 
    dplyr::mutate(exp = forcats::fct_inorder(exp)) %>% 
    tidyr::gather(key, value, -exp) %>% 
    ggplot(aes(x = exp, y = value, col = key)) + 
    geom_point(aes(shape = key), size = 3) + 
    coord_cartesian(ylim = c(0, 1)) + 
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1)
    ) +
    labs(
        x = "",
        y = "Percentage",
        title = paste0("Intrinsic susceptibility of ", k, "-mers ",
                       "versus their true fragility")
    )

ggsave(
    filename = paste0("../figures/original_dataset/breakage_insights/",
                      exp, "_", cols, 
                      "_kmer-", k, "_",
                      "scores-", break_score,
                      "_fragility-vs-susceptibility.pdf"),
    plot = p1,
    height = 8,
    width = 10
)

# plot probability ratio vs. z-scores
ggsave(
    filename = paste0("../figures/original_dataset/breakage_insights/",
                      exp, "_", cols, 
                      "_kmer-", k, "_",
                      "scores-", break_score,
                      "_true_control_breaks.pdf"), 
    plot = do.call(gridExtra::grid.arrange, c(
        sapply(results, `[`, 1), 
        ncol = 3, 
        top = paste0("\nBottom ", top_percent*100, "% (red) vs. top ", 
                                top_percent*100, "% (blue) of the ", 
                               "diff=true-control breakpoint ", k, "-kmer counts\n")
    )),
    height = 25,
    width = 13
)