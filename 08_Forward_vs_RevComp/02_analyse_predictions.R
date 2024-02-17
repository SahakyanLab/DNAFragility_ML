# import libraries
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(pbapply)))
pbapply::pboptions(char = "=")
pbapply::pboptions(char = "=")

# source functions
setwd("/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/08_Forward_vs_RevComp")

# get thresholds and column names
thresholds <- fread("../data/models/python/lightgbm/best_LGBM_model_thresholds.csv")
col_names <- thresholds$FPR
col_names <- gsub("\\.", "_", paste0("Pred_", col_names))

# import predictions
pred_files <- list.files(
    path = "../data/Forward_vs_RevComp",
    pattern = "final_pred",
    full.names = TRUE,
    recursive = TRUE
)
pred_names <- stringr::str_extract(
    string = pred_files,
    pattern = "(Forward|RevComp)(?=_chr)"
)

res <- pbapply::pblapply(1:length(pred_names), function(x){
    prediction <- arrow::read_parquet(pred_files[x])
    prediction <- dplyr::select(prediction, -True)
    colnames(prediction) <- gsub("\\.", "_", colnames(prediction))
    
    return(colSums(prediction, na.rm = TRUE))
})
fwd_ind <- which("Forward" == pred_names)
rc_ind <- which("RevComp" == pred_names)

res_fwd <- do.call(rbind, res[fwd_ind])
res_rc <- do.call(rbind, res[rc_ind])

print(identical(res_fwd, res_rc))