# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))

args = commandArgs(trailingOnly = TRUE)
my.path = as.character(args[1])

my.path="/Users/paddy/Documents/DPhil/DNAFragility/setup"
setwd(my.path)

# create chromosome-separated fasta files
ref_path <- "../../data/ref/hg19"
dir.create(path = ref_path, showWarnings = FALSE, recursive = TRUE)
refseq <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
ref_files <- list.files(path = ref_path, pattern = ".fasta.gz")

if(length(ref_files) < 22){
    for(chr in 1:22){
        Biostrings::writeXStringSet(
            Biostrings::DNAStringSet(refseq[[paste0("chr", chr)]]),
            filepath = paste0(ref_path, "/chr", chr, ".fasta.gz"),
            format = "fasta",
            compress = TRUE
        )
    }
}

# load GEO datasets
files <- list.files(
    path = "../data/setup",
    pattern = "bed",
    recursive = TRUE,
    full.names = TRUE
)
file_names <- stringr::str_extract(
    string = files,
    pattern = "K562_DMSO_DSBs|K562_Top2_mediated_DSBs"
)

# process start positions
for(x in 1:length(files)){
    df <- fread(files[x], header = FALSE, showProgress = FALSE)
    df <- df[df$V1 %in% paste0("chr", 1:22)]
    df[, `:=`(V3 = NULL, V4 = NULL)]
    setnames(df, c("chromosome", "start.pos"))

    for(chr in 1:22){
        # for kmertone
        temp <- df[chromosome == paste0("chr", chr)]  
        fwrite(
            temp, 
            paste0(
                "../data/setup/", file_names[x], 
                "/kmertone/chr", chr, ".csv"
            )
        )

        # for strand breaks
        temp[, chromosome := NULL]
        fwrite(
            temp, 
            paste0(
                "../data/experiments/", file_names[x], 
                "/breakpoint_positions/chr", chr, ".csv"
            )
        )
    }
}

# run kmertone
source("../lib/Kmertone/kmertone.R")
for(x in 1:length(files)){
    dir.create(
        path = paste0("../data/setup/", file_names[x], "/kmertone_scores/"),
        showWarnings = FALSE
    )

    # obtain long-range sequence context for control region cut-off
    longest_range <- fread(paste0(
            "../data/ranges/",
            "kmer_8_Ranges_cutoffs_from_clustering_all-exp.csv"
        )) %>% 
        as_tibble() %>% 
        dplyr::filter(stringr::str_detect(
            string = exp,
            pattern = paste0(
                "K562_Top2_mediated_DSBs/",
                ifelse(grepl("DMSO", file_names[x]), "DMSO", "ETO")
            )
        )) %>% 
        dplyr::filter(rowid == "ranges") %>% 
        dplyr::select(dplyr::contains(".range")) %>% 
        tidyr::pivot_longer(cols = everything()) %>% 
        dplyr::filter(!is.na(value)) %>% 
        dplyr::mutate(value = ceiling(value / 2)) %>% 
        dplyr::pull(value) %>% 
        max(.)
        
    control_regions <- c(longest_range, longest_range+1000)

    kmertone(
        pwd="../lib/Kmertone",
        case.coor.path=paste0("../data/setup/", file_names[x], "/kmertone/"),
        genome.name="unknown", 
        strand.sensitive=FALSE,
        k=8,
        ctrl.rel.pos=control_regions, 
        case.pattern=NULL,
        output.path=paste0("../data/setup/", file_names[x], "/kmertone_scores/"),
        case.coor=NULL, 
        genome=NULL,
        genome.path=paste0(ref_path, "/"),
        rm.case.kmer.overlaps=FALSE,
        case.length=2,
        merge.replicate=TRUE, 
        kmer.table=NULL,
        module="score", 
        ncpu=1
    )
}