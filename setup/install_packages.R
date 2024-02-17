reg_libs <- c(
    # data wrangling
    "data.table", "dplyr", "stringr",
    "purrr",

    # plotting
    "ggplot2", "ggsignif", "ggpattern", "ggrepel",
    "gridExtra",

    # others
    "arrow", "progressr", "furrr", "future",
    "pbapply", "Rcpp", "caret", "DNAshapeR", "RCy3"
)
to_install <- reg_libs[!reg_libs %in% installed.packages()]

# And finally we install the missing packages, including their dependency.
for(lib in to_install){
    install.packages(
        lib, 
        dependencies = TRUE,
        repos = 'http://cran.us.r-project.org'
    )
}

bioscience_libs <- c(
    # data wrangling and more
    "plyranges", "Biostrings", "rtracklayer",
    "taxonomizr",

    # reference genomes
    "BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0",
    "BSgenome.Hsapiens.1000genomes.hs37d5",
    "BSgenome.Hsapiens.UCSC.hg19"   
)
to_install <- bioscience_libs[!bioscience_libs %in% installed.packages()]
for(lib in to_install){
    BiocManager::install(lib, update = TRUE, ask = FALSE)
}