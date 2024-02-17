# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
breakpoint.path <- as.character(args[2])
kmer <- as.integer(args[3])
ref.path <- as.character(args[4])
output.path <- as.character(args[5])
cores <- as.integer(args[6])
first.idx <- as.integer(args[7])

# my.path="/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/03_Breakpoints/02_Correlations/scripts/"
# breakpoint.path="../Raw_data/21-BLISS/K562_Top2_mediated_DSBs/ETO/kmertone"
# kmer=6
# ref.path="hg19"
# output.path="../02_Correlations/data/kmertone/K562_Top2_mediated_DSBs/ETO"
# cores=1
# first.idx=0

# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
setwd(paste0(my.path, "../../Kmertone"))
source("./kmertone.R")

# threshold <- 100

# check for strand sensitivity
sample.files <- list.files(path = breakpoint.path, pattern = "csv|txt")
if(length(sample.files) > 0){
    dt <- fread(paste0(breakpoint.path, "/", sample.files[1]), nrows=1)
}
strand.sensitive <- ifelse("strand" %in% colnames(dt), TRUE, FALSE)

# obtain long-range sequence context for control region cut-off
rmsd.path <- stringr::str_replace(
    string = breakpoint.path,
    pattern = "../Raw_data/",
    replacement = "../01_RMSD/data/"
)
rmsd.path <- stringr::str_replace(
    string = rmsd.path,
    pattern = "kmertone",
    replacement = paste0("new_stats_kmer_", kmer,".csv")
)
dt.rmsd <- fread(rmsd.path)
dt.rmsd <- dt.rmsd["ranges", on = "rowid"]
to.keep <- colnames(dt.rmsd)[grepl(pattern = "range", x = colnames(dt.rmsd))]
longest.range <- ceiling(max(dt.rmsd[, ..to.keep])/2)
print(paste0("Control region range = ", c(longest.range, longest.range+400)), 
             quote = FALSE)

# run kmertone
kmertone(
    case.coor.path=breakpoint.path, 
    genome.name="unknown", 
    strand.sensitive=strand.sensitive, 
    k=kmer,
    ctrl.rel.pos=c(longest.range, longest.range+400), 
    case.pattern=NULL,
    output.path=output.path, 
    case.coor=NULL, 
    genome=NULL,
    genome.path=paste0("../Raw_data/ref/", ref.path), 
    rm.case.kmer.overlaps=FALSE,
    case.length=2, 
    merge.replicate=TRUE, 
    kmer.table=NULL,
    module="score", 
    ncpu=cores
)