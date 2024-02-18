# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(Rcpp)))

args = commandArgs(trailingOnly = TRUE)
my.path = as.character(args[1])
setwd(my.path)

# global vars
kmer <- 8

df.parsed <- fread(
    file = paste0(
        "../data/BaseTable/",
        "kmer_", kmer, ".csv"
    ),
    select = "ID",
    showProgress = FALSE
)

#' Remove the sample ID
df <- as_tibble(df.parsed) %>%
    tidyr::separate(
        col = ID,
        into = c(
            "Chr", "Start", 
            "Tissue", "Cancer", 
            "BreakType", 
            "SampleID"
        ),
        sep = "_",
        remove = TRUE
    ) %>% 
    tidyr::unite(
        TC_ID,
        c("Tissue", "Cancer"), 
        remove = FALSE
    ) %>% 
    tidyr::unite(
        Chr_Start_ID, 
        c("Chr", "Start"),
        remove = FALSE
    ) %>% 
    dplyr::select(-c(SampleID, Tissue, Cancer)) %>% 
    dplyr::mutate(Start = as.numeric(Start)) %>% 
    dplyr::select(Chr, Start) %>% 
    dplyr::rename_with(~c("seqnames", "start")) %>% 
    dplyr::distinct() %>% 
    as.data.table()

# liftover to telomere-to-telomere genome version
refseq <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0::BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
seqnames(refseq) <- paste0("chr", seqnames(refseq))

# liftover breakpoints to the telomere-to-telomere genome version
suppressPackageStartupMessages(suppressWarnings(library(rtracklayer)))
chain <- import.chain("../../data/liftover/hg38-chm13v2.over.chain")
df <- df %>% dplyr::mutate(width = 1, strand = "+")
df <- plyranges::as_granges(df)
df <- liftOver(df, chain)
df <- unlist(as(df, "GRangesList"))
df <- as_tibble(df) %>% 
    dplyr::select(seqnames, start) %>% 
    as.data.table()

# split by chromosome, and run kmertone
dir.create(
    "../data/experiments/COSMIC", 
    showWarnings = FALSE,
    recursive = TRUE
)
dir.create(
    "../../data/experiments/COSMIC/breakpoint_positions/", 
    showWarnings = FALSE,
    recursive = TRUE
)

df_bp <- lapply(1:22, function(chr){
    df_chr <- df[seqnames == paste0("chr", chr)]
    setorder(df_chr, start)
    setnames(df_chr, c("chromosome", "start.pos"))

    fwrite(
        df_chr, 
        paste0("../data/experiments/COSMIC/chr", chr, ".csv")
    )

    fwrite(
        df_chr, 
        paste0("../../data/experiments/COSMIC/breakpoint_positions/chr", chr, ".csv")
    )

    return(df_chr)
})
# names(df_bp) <- paste0("chr", 1:length(df_bp))
# df_bp <- lapply(df_bp, function(df) df[, chromosome := NULL])

# get max ranges from each category of breakage
max_range <- fread(paste0(
    "../../data/range_effects/",
    "MaxValuesFromClustersByType.csv"
))
max_range <- apply(max_range[, -"break_type"], 2, max)
round_to_nearest_even <- function(x) round(x/2)*2
ranges <- round_to_nearest_even(max_range)
names(ranges) <- stringr::str_to_title(names(ranges))

# run kmertone calculations
dir.create(
    path = "../data/kmertone/COSMIC/kmertone_scores/",
    showWarnings = FALSE
)
control.regions <- c(ranges[["Long"]], ranges[["Long"]]+1000)

# create chromosome-separated fasta files
ref_path <- "../../data/ref/telo_to_telo"
dir.create(path = ref_path, showWarnings = FALSE, recursive = TRUE)
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

source("../../lib/Kmertone/kmertone.R")
kmertone(
    pwd="../../lib/Kmertone",
    case.coor.path="../../data/experiments/COSMIC/",
    genome.name="unknown", 
    strand.sensitive=FALSE,
    k=8,
    ctrl.rel.pos=control.regions, 
    case.pattern=NULL,
    output.path="../data/kmertone/COSMIC/kmertone_scores/",
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