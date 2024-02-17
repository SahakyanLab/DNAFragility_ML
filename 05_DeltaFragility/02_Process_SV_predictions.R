# import libraries
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))
suppressPackageStartupMessages(suppressWarnings(library(arrow)))
pbapply::pboptions(char = "=", type = "txt")

setwd("/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/05_DeltaFragility")

# get thresholds and column names
thresholds <- fread("../data/models/python/lightgbm/best_LGBM_model_thresholds.csv")
col_names <- thresholds$FPR
col_names <- gsub("\\.", "_", paste0("Pred_", col_names))

plot_titles <- stringr::str_replace(
    string = col_names, 
    pattern = "^(Pred)_0_(.*)$", 
    replacement = "\\1: 0.\\2 (FPR)"
)
plot_titles <- setNames(plot_titles, col_names)

# import random breakpoints in low/mid/high fragile regions
df_random_frac <- fread("../data/Whole_genome_pred/random_sampling/random_frac.csv")

# import labels and predictions
path_to_delta_pred <- paste0(
    "/media/hert6114/Paddy_6TB/ProjectBoard_Patrick/",
    "04_DNAFragility/data"
)
labels <- arrow::read_parquet("./data/all_SV_labels.parquet")

############################################################################
t1 <- Sys.time()
cur.msg <- "Listing all files and preparing for processing"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

#' Run the below command in the directory outside deltafragility
pred_files <- fread("./data/filenames.csv", header = FALSE)
pred_files <- pred_files$V1

pred_names <- stringr::str_remove_all(
    string = pred_files, 
    pattern = "deltafragility/|_final_pred.parquet"
)
before_after <- stringr::str_extract(
    string = pred_names, pattern = "before|after"
)
pred_names <- stringr::str_remove(
    string = pred_names, pattern = "/before|/after"
)

# match labels
ind <- match(pred_names, labels$labels)
labels <- labels[ind,]
labels[, before_after := before_after]

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")
############################################################################

t1 <- Sys.time()
cur.msg <- "Listing all files and preparing for processing"

res_list <- pbapply::pbsapply(1:nrow(labels), function(x){
    prediction <- arrow::read_parquet(paste0(
        path_to_delta_pred, "/", pred_files[x]
    ))
    prediction <- dplyr::select(prediction, -True)
    colnames(prediction) <- gsub("\\.", "_", colnames(prediction))
    prediction <- colSums(prediction)
    return(prediction)
})

cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))
total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")

############################################################################
t1 <- Sys.time()
cur.msg <- "Saving all prediction outcomes"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

res_all <- t(res_list)

res_all <- as.data.frame(res_all)
res_all <- as.data.table(res_all)
df_all <- cbind(res_all, labels)
fwrite(df_all, "./data/all_SV_predictions.csv")

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")