suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))

# load file
df <- fread("./QueryTable_kmer-6_zscore.csv") %>% as_tibble()
column_names <- list()

# TFBS
tfbs_ind <- which(grepl(pattern = "^TFBS_", df$category))
column_names$TFBS <- df$category[tfbs_ind]

# G4seq map
g4_ind <- which(grepl(pattern = "^G4seq", df$category))
column_names$G4MAP <- df$category[g4_ind]

# epigenome marks
epigenome_ind <- which(grepl(
    pattern = "Epigenome|ATACseq|Chipseq|Dnaseseq|FAIREseq", 
    df$category
))
column_names$EPIGENOME_MARKS <- df$category[epigenome_ind]            

# remaining ones are the breakage scores
breaks_ind <- which(!grepl(
    pattern = paste0(
        "^TFBS_|^[A-Z]_ds.dEhof_ds_ss|^G4seq|",
        "Epigenome|ATACseq|Chipseq|Dnaseseq|FAIREseq"
    ), 
    df$category
))
column_names$BREAKS <- df$category[breaks_ind]

# split file by category
for(name in names(column_names)){
    ind <- match(column_names[[name]], df$category)
    df_subset <- as.data.table(df[ind,])

    fwrite(
        df_subset,
        paste0(name, ".csv")
    )
}