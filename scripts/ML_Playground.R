setwd("/Users/paddy/Documents/DPhil/04_DNAFragility/scripts/")
source("./Main.R")

# score = "ratio"
score = "z-score"
k = 8

cols = "all"
# cols = "complex"
# cols = "only_breaks"
# cols = "only_triplets"

# features$get_features(
#     G4_REGEX = FALSE, g4_type = g4_type, 
#     GC_COUNT = FALSE,
#     KMER_COUNTS = FALSE, kmer_window = kmer_window,
#     VIENNA_RNA = FALSE, sliding_window = sliding_window, 
#     nuc_type = FALSE, RNAfold.CALL = RNAfold.CALL,
#     DNA_SHAPE = FALSE,
#     TFBS_EUCLIDEAN_DISTANCE = FALSE,
#     OCCUPANCY_SCORES = TRUE
#     SAVE_OUTPUT = FALSE
# )

# col.names <- colnames(features$feature_matrix)[grepl(
#     pattern = "TFBS_EucDist_PC*",
#     x = colnames(features$feature_matrix)
# )]
# temp = features$feature_matrix[, ..col.names]; temp
# features$get_features_from_csv()

# colnames(features$feature_matrix)
# features$feature_matrix = features$feature_matrix[, 1:74]
# features$feature_matrix <- cbind(features$feature_matrix, temp)

# fwrite(
#     features$feature_matrix,
# )

# col.names <- colnames(features$feature_matrix)[which(!is.na(
#     str_extract(string = colnames(features$feature_matrix), 
#                 pattern = "HelT|MGW|ProT|Roll")
# ))]

# temp = features$feature_matrix[, ..col.names]; temp
# features$get_features_from_csv()
# features$feature_matrix[, 67:ncol(temporary)]
# features$feature_matrix[, 67:ncol(temporary)] = temp
###########################################################################
# extract features
features$get_features_from_csv()
features$feature_matrix

features$select_columns(cols = cols)
# features$select_columns(cols = "only_singleton")
features$feature_matrix

features$train_test_split(train.split = 0.7)
features$standardise(cols = cols)
training <- features$train_matrix
testing <- features$test_matrix

################################################################################
# ML models using h2o framework
################################################################################
# convert data to h2o object
training[, predictor := as.factor(predictor)]
testing[, predictor := as.factor(predictor)]

suppressPackageStartupMessages(library(h2o))

############################################
# Random Forest
############################################
h2o.init(min_mem_size = "8g")
training.h2o <- as.h2o(training)
RANDOM_SEARCH_MINUTES=60

tuneGrid <- list(
    # Maximum depth of a tree, try within [3,10]
    max_depth = seq(6, 14, 1),
    # Denotes the fraction of observations
    # to be randomly sampled for each tree, try within [0.5, 1].
    # Avoids trees becoming highly correlated
    sample_rate = seq(0.4, 1, 0.1),
    # Number of decision trees in the random forest model
    ntrees = 5000
)

search_criteria <- list(
    strategy = "RandomDiscrete",
    max_runtime_secs = 60*RANDOM_SEARCH_MINUTES,
    seed = 1234
)

h2o.rf.grid.fit <- h2o.grid(
    grid_id = "model_grid",
    algorithm = "randomForest",
    x = colnames(training.h2o)[2:ncol(training.h2o)],
    y = "predictor",
    training_frame = training.h2o,
    hyper_params = tuneGrid,
    search_criteria = search_criteria,
    fold_assignment = "Stratified",
    nfolds = 5,
    seed = 1234
)

h2o.rf.gridperf <- h2o.getGrid(
    grid_id = "model_grid",
    sort_by = "auc",
    decreasing = TRUE
)
h2o.rf.gridperf

# Grab the top random forest model, chosen by validation AUC
best.h2o.rf <- h2o.getModel(h2o.rf.gridperf@model_ids[[1]])

h2o::h2o.saveModel(
    best.h2o.rf,
    path = "../data/models/",
    force = TRUE,
    filename = paste0("random_forest", cols, "_kmer_", 
                      k, "_", score,".RData")
)

# Evaluate model performance on test set
h2o.rf.perf <- h2o.performance(best.h2o.rf, newdata = as.h2o(testing))
h2o.rf.auc <- h2o.auc(h2o.rf.perf)
h2o.rf.auc

h2o.shutdown(prompt = FALSE)
###########################################################################

# test <- as_tibble(training) %>% 
#     dplyr::filter(predictor == "YES") %>% 
#     dplyr::select(2:24)
#     # dplyr::select(top.corr$experiment)
# top.corr
# test
# fwrite(test, "./highest_corr_exp.csv")

# example = tibble(
#     true = test$Natural_fragmentation_Neandertal,
#     predicted = 0.71*test$Natural_fragmentation_Denisovan+0.54*test$Cell_free_DNA_disease-14.92)
    
# correlation = cor(example$true, example$predicted)

# p1 <- example %>% 
#     ggplot(aes(x = true, y = predicted)) + 
#     geom_point(alpha = 0.2) +
#     ggplot2::geom_abline(
#         slope = 1,
#         linetype = "dashed"
#     ) +
#     ggplot2::geom_smooth(
#         method = lm, 
#         formula = y ~ x
#     ) + 
#     labs(
#         title = paste0("R^2 = ", signif(correlation, 3), "\n"),
#     )

# ggsave(
#     p1, 
#     filename = "eureqa_true_breaks-Natural_fragmentation_Neandertal.png", 
#     width = 8, 
#     height = 8
# )

################################################################################
# ML models using LightGBM frameworks
suppressPackageStartupMessages(library(lightgbm))

# create LGB dataset
train.label <- ifelse(training$predictor == "YES", 1, 0)
train.data <- as.matrix(training[, -"predictor"])
test.label <- ifelse(testing$predictor == "YES", 1, 0)
test.data <- as.matrix(testing[, -"predictor"])

dtrain <- lgb.Dataset(data = train.data, label = train.label)
dtest <- lgb.Dataset.create.valid(dataset = dtrain, data = test.data, label = test.label)

tuneGrid <- expand.grid(
    # Size of each step in the gradient descent method
    learning_rate = c(0.01, 0.1), 
    # Threshold for classifcation
    cutoff = c(0.3, 0.4, 0.5, 0.6), 
    # Maximum number of leafs in a tree
    num_leaves = c(5, 8, 12, 15),
    # Maximum depth of a tree, try within [3,10]
    max_depth = c(6, 8, 12, 14),
    # Same as subsample of GBM. Denotes the fraction of observations
    # to be randomly sampled for each tree, try within [0.5, 1].
    # Avoids trees becoming highly correlated
    bagging_fraction = c(0.4, 0.6, 1),
    acc = 0
)
 
# validataion data
valids <- list(test = dtest)
set.seed(1234)

for (i in 1:nrow(tuneGrid)){
    print(paste0(i, "/", nrow(tuneGrid)))
    params <- list(
        objective = 'binary',
        metric = 'auc',
        learning_rate = tuneGrid$learning_rate[i],
        cutoff = tuneGrid$cutoff[i],
        num_leaves = tuneGrid$num_leaves[i],
        max_depth = tuneGrid$max_depth[i],
        bagging_fraction = tuneGrid$bagging_fraction[i]
    ) 

    fit.lgb <- lgb.train(
        params = params,
        data = dtrain,
        nrounds = 20000L,
        valids = valids,
        verbose = 0
    )

    prob.lgb <- predict(fit.lgb, test.data)
    pred.lgb <- ifelse(prob.lgb > tuneGrid$cutoff[i], 1, 0)
    confusion.mat <- table(pred.lgb, test.label)
    acc <- sum(diag(confusion.mat))/sum(confusion.mat)
    tuneGrid[i, "acc"] <- acc
}

# find best model
tuneGrid[which.max(tuneGrid$acc), ]
bestParam <- tuneGrid[which.max(tuneGrid$acc), ]

bestParam <- list(
    learning_rate = fit.lgb$params$learning_rate,
    cutoff = fit.lgb$params$cutoff,
    num_leaves = fit.lgb$params$num_leaves,
    max_depth = fit.lgb$params$max_depth,
    bagging_fraction = fit.lgb$params$bagging_fraction
)

# final model
params <- list(
    objective = 'binary',
    metric = 'auc',
    learning_rate = bestParam$learning_rate,
    cutoff = bestParam$cutoff,
    num_leaves = bestParam$num_leaves,
    max_depth = bestParam$max_depth,
    bagging_fraction = bestParam$bagging_fraction
) 

fit.lgb <- lgb.train(
    params = params,
    data = dtrain,
    nrounds = 20000L,
    valids = valids
)

saveRDS.lgb.Booster(
    fit.lgb,
    file = paste0("../data/models/lightgbm_", cols, 
                  "_kmer_", k, "_", score, 
                  "_cutoff-", 
                  gsub(
                    pattern = "\\.", 
                    replacement = "_", 
                    x = as.character(bestParam$cutoff)),
                  ".RData")
)

# prob.lgb <- predict(fit.lgb, test.data)
# pred.lgb <- ifelse(prob.lgb > bestParam$cutoff, 1, 0)
# confusion.mat <- table(pred.lgb, test.label)
# lgb.auc <- sum(diag(confusion.mat))/sum(confusion.mat)

# Evaluation Curve
pred <- ROCR::prediction(prob.lgb, test.label)
roc <- ROCR::performance(pred, "tpr", "fpr")
lgb.auc <- ROCR::performance(pred, "tpr", "fpr", measure = "auc")
lgb.auc <- lgb.auc@y.values[[1]]

# feature importance
tree.imp <- lgb.importance(fit.lgb, percentage = TRUE)
tree.imp[, feature_importance := (Gain-min(Gain))/(max(Gain)-min(Gain))*100]

p1 <- as_tibble(tree.imp) %>% 
    dplyr::arrange(feature_importance) %>% 
    dplyr::mutate(Feature = forcats::fct_inorder(Feature)) %>% 
    ggplot(aes(x = feature_importance, y = Feature)) + 
    geom_segment(aes(xend = 0, yend = Feature)) + 
    geom_point(col = "orange") + 
    theme_bw() + 
    labs(
        x = "Weights",
        y = "Features",
        title = "Variable Importance: LightGBM Model"
    )

ggsave(
    filename = paste0("../figures/lightgbm_varimp_plot_", 
                      cols, "_kmer_", k, "_", score,".pdf"),
    plot = p1,
    width = 10, 
    height = 10
)

################################################################################
# ML models using h2o framework
################################################################################
# convert data to h2o object
training[, predictor := as.factor(predictor)]
testing[, predictor := as.factor(predictor)]

suppressPackageStartupMessages(library(h2o))

############################################
# Random Forest
############################################
h2o.init(min_mem_size = "8g")
training.h2o <- as.h2o(training)
RANDOM_SEARCH_MINUTES=60

# construct hyper-parameter space
# tuneGrid <- list(
#     max_depth = c(6, 8, 14),
#     sample_rate = c(0.4, 0.8, 1),
#     ntrees = 5000
# )

tuneGrid <- list(
    # Maximum depth of a tree, try within [3,10]
    max_depth = seq(6, 14, 1),
    # Denotes the fraction of observations
    # to be randomly sampled for each tree, try within [0.5, 1].
    # Avoids trees becoming highly correlated
    sample_rate = seq(0.4, 1, 0.1),
    # Number of decision trees in the random forest model
    ntrees = 5000
)

search_criteria <- list(
    strategy = "RandomDiscrete",
    max_runtime_secs = 60*RANDOM_SEARCH_MINUTES,
    seed = 1234
)

h2o.rf.grid.fit <- h2o.grid(
    grid_id = "model_grid",
    algorithm = "randomForest",
    x = colnames(training.h2o)[2:ncol(training.h2o)],
    y = "predictor",
    training_frame = training.h2o,
    hyper_params = tuneGrid,
    search_criteria = search_criteria,
    fold_assignment = "Stratified",
    nfolds = 5,
    seed = 1234
)

h2o.rf.gridperf <- h2o.getGrid(
    grid_id = "model_grid",
    sort_by = "auc",
    decreasing = TRUE
)
h2o.rf.gridperf

# Grab the top random forest model, chosen by validation AUC
best.h2o.rf <- h2o.getModel(h2o.rf.gridperf@model_ids[[1]])

h2o::h2o.saveModel(
    best.h2o.rf,
    path = "../data/models/",
    force = TRUE,
    filename = paste0("random_forest", cols, "_kmer_", 
                      k, "_", score,".RData")
)

# Evaluate model performance on test set
h2o.rf.perf <- h2o.performance(best.h2o.rf, newdata = as.h2o(testing))
h2o.rf.auc <- h2o.auc(h2o.rf.perf)
h2o.rf.auc

h2o.shutdown(prompt = FALSE)

###################################
h2o.init(min_mem_size = "8g")
training.h2o <- as.h2o(training)

# Linear Regression
glm.fit <- h2o.glm(
    family = "binomial",
    x = colnames(training)[2:ncol(training)],
    y = "predictor",
    training_frame = training.h2o,
    fold_assignment = "Stratified",
    standardize = FALSE,
    lambda_search = TRUE,
    nfolds = 5,
    seed = 1234
)

h2o::h2o.saveModel(
    glm.fit,
    path = "../data/models/",
    force = TRUE,
    filename = paste0("linear_", cols, "_kmer_", 
                      k, "_", score,".RData")
)

glm.perf <- h2o.performance(glm.fit, newdata = as.h2o(testing))
glm.auc <- h2o.auc(glm.perf)
glm.auc
h2o.shutdown(prompt = FALSE)

###################################
h2o.init(min_mem_size = "8g")
training.h2o <- as.h2o(training)
RANDOM_SEARCH_MINUTES = 10

# Naive Bayes
tuneGrid <- list(
    laplace = seq(0, 10, 0.1)
)

search_criteria <- list(
    strategy = "RandomDiscrete",
    max_runtime_secs = 60*RANDOM_SEARCH_MINUTES,
    seed = 1234
)

h2o.nb.grid.fit <- h2o.grid(
    grid_id = "model_grid",
    algorithm = "naivebayes",
    x = colnames(training.h2o)[2:ncol(training.h2o)],
    y = "predictor",
    training_frame = training.h2o,
    hyper_params = tuneGrid,
    search_criteria = search_criteria,
    fold_assignment = "Stratified",
    nfolds = 5,
    seed = 1234
)

h2o.nb.gridperf <- h2o.getGrid(
    grid_id = "model_grid",
    sort_by = "auc",
    decreasing = TRUE
)
h2o.nb.gridperf

# Grab the top naive bayes model, chosen by validation AUC
best.h2o.nb <- h2o.getModel(h2o.nb.gridperf@model_ids[[1]])

# Evaluate model performance on test set
h2o.nb.perf <- h2o.performance(best.h2o.nb, newdata = as.h2o(testing))
h2o.nb.auc <- h2o.auc(h2o.nb.perf)
h2o.nb.auc

h2o::h2o.saveModel(
    best.h2o.nb,
    path = "../data/models/",
    force = TRUE,
    filename = paste0("naive_bayes_", cols, "_kmer_", 
                      k, "_", score,".RData")
)

#########################################################
testing.h2o <- as.h2o(testing)

glm.fit <- h2o::h2o.loadModel("../data/models/linear_all_kmer_8_ratio.RData")
glm.perf <- h2o.performance(glm.fit, newdata = testing.h2o)
glm.auc <- h2o.auc(glm.perf)

nb.fit <- h2o::h2o.loadModel("../data/models/naive_bayes_all_kmer_8_ratio.RData")
nb.perf <- h2o.performance(nb.fit, newdata = testing.h2o)
nb.auc <- h2o.auc(nb.perf)

rf.fit <- h2o::h2o.loadModel("../data/models/random_forestall_kmer_8_ratio.RData")
rf.perf <- h2o.performance(rf.fit, newdata = testing.h2o)
rf.auc <- h2o.auc(rf.perf)

# plottings
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
    purrr::map2(c(glm.auc, rf.auc, nb.auc),
        function(x, y) x %>% dplyr::mutate(auc = y)
    ) %>% 
    rbindlist()

# combine results with light gbm model
rows.per.group <- h2o.model.perf[, .(count = .N), by = model]$count[1]

lgb.dt <- data.table(
    fpr = attr(roc, "x.values")[[1]],
    tpr = attr(roc, "y.values")[[1]],
    model = paste0('Light GBM', ' (', signif(lgb.auc, 3),')'),
    auc = lgb.auc
)

h2o.model.perf <- rbind(
    h2o.model.perf, 
    lgb.dt[seq(1, nrow(lgb.dt), length.out = rows.per.group), ]
)

p1 <- as_tibble(h2o.model.perf) %>% 
    dplyr::arrange(desc(auc)) %>% 
    dplyr::mutate(model = forcats::fct_inorder(model)) %>% 
    ggplot(aes(x = fpr, y = tpr, col = model)) + 
    geom_line(size = 1.3) + 
    geom_segment(aes(
        x = 0, y = 0,
        xend = 1, yend = 1),
        linetype = 2,
        col = "darkgrey"
    ) + 
    theme_bw() + 
    labs(
        x = "False Positive Rate",
        y = "True Positive Rate",
        title = "ROC Curves"
    )

ggsave(
    filename = paste0("../figures/h2o_plots_", cols, 
                      "_kmer_", k, "_", score,".pdf"), 
    plot = p1,
    width = 10,
    height = 8
)

# Standard coef. magnitude of linear model
p1 <- tibble(
    features = names(glm.fit@model$coefficients),
    weights = unname(glm.fit@model$coefficients)) %>% 
    dplyr::filter(features != "Intercept") %>% 
    dplyr::mutate(
        coef = ifelse(weights >= 0, "Positive", "Negative"),
        weights = abs(weights)) %>% 
    dplyr::slice(1:30) %>% 
    dplyr::arrange(weights) %>% 
    dplyr::mutate(features = forcats::fct_inorder(features)) %>% 
    ggplot(aes(x = weights, y = features, fill = coef)) + 
    geom_bar(stat = "identity") + 
    theme_bw() + 
    labs(
        x = "Weights",
        y = "Features",
        title = "Standardised Coef. Magnitudes of Linear Model"
    )

ggsave(
    filename = paste0("../figures/h2o_glm_std_coef_plot_", 
                      cols, "_kmer_", k, "_", score,".pdf"),
    plot = p1,
    width = 10, 
    height = 10
)

# Variable Importance of Linear Model
p1 <- tibble(
    features = glm.fit@model$variable_importances$variable,
    weights = glm.fit@model$variable_importances$scaled_importance) %>%  
    dplyr::arrange(weights) %>% 
    dplyr::mutate(features = forcats::fct_inorder(features)) %>% 
    ggplot(aes(x = weights, y = features)) + 
    geom_segment(aes(xend = 0, yend = features)) + 
    geom_point(col = "orange") + 
    theme_bw() + 
    labs(
        x = "Weights",
        y = "Features",
        title = "Variable Importance: Linear Model"
    )

ggsave(
    filename = paste0("../figures/h2o_glm_varimp_plot_linear_", 
                      cols, "_kmer_", k, "_", score,".pdf"),
    plot = p1,
    width = 10, 
    height = 10
)

# Variable Importance of Random Forest
h2o.init(min_mem_size = "8g")
rf.fit <- h2o::h2o.loadModel("../data/models/random_forest/random_forestall_kmer_8_z-score.RData")

p1 <- tibble(
    features = rf.fit@model$variable_importances$variable,
    weights = rf.fit@model$variable_importances$scaled_importance) %>%  
    dplyr::arrange(weights) %>% 
    dplyr::mutate(features = forcats::fct_inorder(features)) %>% 
    ggplot(aes(x = weights, y = features)) + 
    geom_segment(aes(xend = 0, yend = features)) + 
    geom_point(col = "orange") + 
    theme_bw() + 
    labs(
        x = "Weights",
        y = "Features",
        title = "Variable Importance: Random Forest Model"
    )

ggsave(
    filename = paste0("../figures/h2o_glm_varimp_plot_random_forest_", 
                      cols, "_kmer_", k, "_", score,".pdf"),
    plot = p1,
    width = 10, 
    height = 20
)
################################################################################
################################################################################
################################################################################
# JASPAR Playground
suppressPackageStartupMessages(library(JASPAR2020))
suppressPackageStartupMessages(library(TFBSTools))

# extract motifs corresponding to human transcription factors
pwm_library <- getMatrixSet(
  JASPAR2020,
  opts = list(
    collection = 'CORE',
    species    = 'Homo sapiens',
    matrixtype = 'PWM'
))

# match motifs
pwm_matrix <- sapply(pwm_library, as.matrix)
pwm_ids <- unname(ID(pwm_library))
pwm_names <- unname(name(pwm_library))
pwm_nameid <- paste0(pwm_ids, "_", pwm_names)

# Extract 8-mers for PWM matching
table <- KmerTable$new(k = k, action = score, statistic = statistic, exp = exp)
kmer.table <- DNAStringSet(table$kmer_ref$kmer)

# one-hot encoding matrix
atgc.order <- rownames(as.matrix(pwm_library[[1]]))
atgc.order <- paste(atgc.order, collapse = "")

kmer.refs <- pbapply::pblapply(1:length(kmer.table), function(x){
    t(letterFrequencyInSlidingView(
        x = kmer.table[[x]],
        view.width = 1,
        letters = atgc.order,
        OR = 0,
        as.prob = TRUE
    ))
})

breaks = Breakpoints$new(k=k, exp=exp, seed=1234)
breaks

# find the most similar motif to our motif
pwm.matches <- pbapply::pblapply(1:length(kmer.refs), function(x){
    pwm_sim <- PWMSimilarity(
        pwmSubject = pwm_library,
        pwmQuery = kmer.refs[[x]],
        method = 'Euclidean'
    )
    return(unname(pwm_sim))
})

pwm.matches <- do.call(cbind, pwm.matches)
pwm.matches <- t(pwm.matches)
colnames(pwm.matches) <- pwm_nameid
pwm.matches <- as.data.table(pwm.matches)
pwm.matches <- cbind(table$kmer_ref, pwm.matches)
fwrite(
    pwm.matches,
    file = "../data/tfbs/QueryTable-kmer_8.csv"
)

# pwm.matches <- pbapply::pblapply(1:100, function(break_seq){
#     hits <- lapply(pwm_matrix, function(pwm_seq){
#         matchPWM(
#             pwm_seq, 
#             long.ranges.seq[[break_seq]], 
#             min.score = "80%"
#         )
#     })
#     hit.counts <- sapply(hits, length)
#     return(hit.counts)
# })

################################################################################
################################################################################
################################################################################
# histone markers

################################################################################
################################################################################
################################################################################