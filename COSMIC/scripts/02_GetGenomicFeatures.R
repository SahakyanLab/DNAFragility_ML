# Load in all the cosmic mutations, work out a grouping, then do mutations and
# plot deltaRTs with those groupings
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)))
suppressPackageStartupMessages(suppressWarnings(library(org.Hs.eg.db)))
suppressPackageStartupMessages(suppressWarnings(library(GenomicFeatures)))
suppressPackageStartupMessages(suppressWarnings(library(rtracklayer)))

args = commandArgs(trailingOnly = TRUE)
my.path = as.character(args[1])
setwd(my.path)

# Download from https://github.com/marbl/CHM13/tree/master#gene-annotation
# Link: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz
txdb <- GenomicFeatures::makeTxDbFromGFF(
    file = '../data/annotations/chm13v2.0_RefSeq_Liftoff_v5.1.gff3',
    format = 'gff3',
    dataSource = 'UCSC',
    organism = 'Homo sapiens',
    taxonomyId = 9606
    ) %>% 
    suppressWarnings()

# import all genomic features for analysis
# Converting IDs for visuailation (mapper, id_list, TO, FROM)
possible_id_types <- columns(org.Hs.eg.db)

# liftover from hg38 to T2T genome version
refseq <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0
if(!any(grepl(pattern = "^chr", x = seqnames(refseq)))){
    chr_names <- paste0("chr", seqnames(refseq))
    seqnames(refseq@seqinfo) <- seqnames(refseq) <- chr_names
}
# get the start and end positions of each chromosome
refseq.table <- as.data.frame(refseq@seqinfo)
refseq.table <- refseq.table[grepl(
    pattern = "^chr([1-9]|1[0-9]|2[0-2])$", 
    x = rownames(refseq.table)
), ]
refseq.table <- data.table(
    chromName = rownames(refseq.table),
    chromSize = refseq.table$seqlengths
)
t2t_chain <- import.chain("../../data/liftover/hg38-chm13v2.over.chain")

# transcripts
all_transcripts <- GenomicFeatures::transcripts(txdb, columns = "gene_id") %>% 
    dplyr::mutate(
        gene_id = as.character(gene_id),
        type = "Transcript"
    ) %>% 
    as_tibble() %>% 
    tidyr::drop_na() %>% 
    plyranges::as_granges()

# all genes
all_genes <- GenomicFeatures::genes(txdb) %>%
    dplyr::mutate(type = "Genes")

# save_file <- all_genes %>% dplyr::arrange(seqnames, start)
# seqlevels(save_file, pruning.mode = "coarse") <- paste0("chr", 1:22)
# plyranges::write_bed(save_file, "../data/annotations/genes.bed")

# write.table(
#     x = refseq.table, 
#     file = "../data/annotations/genome.txt",
#     quote = FALSE,
#     sep = "\t",
#     row.names = FALSE
# )
#' gff2bed < yourfile_sorted.gff > yourfile.bed
#' sort -k1,1 -k2,2n yourfile.bed > yourfile_sorted.bed
#' cut -f1 genome.txt | tail -n +2 > valid_chromosomes.txt
#' awk 'NR==FNR{c[$1]++;next};c[$1] > 0' valid_chromosomes.txt yourfile_sorted.bed > filtered.bed
#' bedtools complement -i filtered.bed -g genome.txt > intergenic.bed

##########################################################################################
#' group 1

t1 <- Sys.time()
cur.msg <- "Processing annotations for genic elements"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

# promoter sequences
all_promoters_upstream <- GenomicFeatures::promoters(
        txdb, upstream = 2500, downstream = 0
    ) %>% suppressWarnings()
all_promoters_upstream <- unlist(as(all_promoters_upstream, "GRangesList"))
all_promoters_upstream <- all_promoters_upstream %>% 
    reduce(.) %>% 
    dplyr::mutate(type = "Promoters_Upstream")

all_promoters_downstream <- GenomicFeatures::promoters(
        txdb, upstream = 0, downstream = 1500
    ) %>% suppressWarnings()
all_promoters_downstream <- unlist(as(all_promoters_downstream, "GRangesList"))
all_promoters_downstream <- all_promoters_downstream %>% 
    reduce(.) %>% 
    dplyr::mutate(type = "Promoters_Downstream") 

all_promoters_all <- GenomicFeatures::promoters(
        txdb, upstream = 2500, downstream = 1500
    ) %>% suppressWarnings()
all_promoters_all <- unlist(as(all_promoters_all, "GRangesList"))
all_promoters_all <- all_promoters_all %>% 
    reduce(.) %>% 
    dplyr::mutate(type = "Promoters_All") 

# transcription start sites
all_tss <- resize(all_transcripts, width = 1, fix = 'start')
all_tss <- all_tss %>% 
    reduce(.) %>% 
    dplyr::mutate(type = "TSS")

# Downloaded from https://reftss.riken.jp/datafiles/4.1/human/
# all_tss <- readRDS("../data/annotations/refTSS_v4.1_human_coordinate.hg38.rds")
# all_tss <- as_tibble(all_tss) %>% 
#     dplyr::select(seqnames, start, end, width, strand) %>% 
#     dplyr::mutate(type = "TSS") %>% 
#     plyranges::as_granges()

# introns
all_introns <- GenomicFeatures::intronsByTranscript(txdb)
all_introns <- unlist(as(all_introns, "GRangesList"))
all_introns <- all_introns %>% 
    reduce(.) %>% 
    dplyr::mutate(type = "Introns")

# exons
all_exons <- GenomicFeatures::exonsBy(txdb)
all_exons <- unlist(as(all_exons, "GRangesList"))
all_exons <- all_exons %>% 
    reduce(.) %>% 
    dplyr::mutate(type = "Exons")

# intergenic region
genic <- reduce(c(all_genes, all_promoters_all, all_introns), ignore.strand = TRUE)
# genic <- reduce(all_genes, ignore.strand = TRUE)
intergenic <- gaps(genic)

# The last step is really important, as otherwise you 
# will get an additional 2 entries per chromosome (one for each of + and -)
all_intergenic <- intergenic[strand(intergenic) == "*"]

# # https://support.bioconductor.org/p/66003/
# exbygene <- GenomicFeatures::exonsBy(txdb, "gene")
# all_intergenic <- gaps(unlist(range(exbygene)))
all_intergenic <- all_intergenic %>% 
    reduce(.) %>% 
    dplyr::mutate(type = "Intergenic")

# coding sequence
all_cds <- GenomicFeatures::cdsBy(txdb)
all_cds <- unlist(as(all_cds, "GRangesList"))
all_cds <- all_cds %>% 
    reduce(.) %>% 
    dplyr::mutate(type = "CDS")

# 5' UTR
all_fiveUTR <- GenomicFeatures::fiveUTRsByTranscript(txdb)
all_fiveUTR <- unlist(as(all_fiveUTR, "GRangesList"))
all_fiveUTR <- all_fiveUTR %>% 
    reduce(.) %>% 
    dplyr::mutate(type = "Five_prime_UTR")

# 3' UTR
all_three_UTR <- GenomicFeatures::threeUTRsByTranscript(txdb)
all_three_UTR <- unlist(as(all_three_UTR, "GRangesList"))
all_three_UTR <- all_three_UTR %>% 
    reduce(.) %>% 
    dplyr::mutate(type = "Three_prime_UTR")

# combine all genomic annotations
group_one <- plyranges::bind_ranges(
        all_introns, all_exons, all_intergenic,
        all_cds, all_promoters_all, 
        all_promoters_upstream, all_promoters_downstream, 
        all_transcripts, all_fiveUTR, all_three_UTR, all_tss
    ) %>% 
    dplyr::select(type) 

seqlevels(group_one, pruning.mode = "coarse") <- paste0("chr", 1:22)

#' intergenic region
#' Transcripts minus up/downstream region of promoters. 
#' Expect close to introns fragility. Our internal control.
# transcript_temp <- dplyr::filter(group_one, type == "Transcript")
# prom_temp <- dplyr::filter(group_one, type == "Promoters_All")
# all_intergenic <- setdiff(transcript_temp, prom_temp) %>% 
#     dplyr::mutate(type = "Intergenic")

# group_one <- plyranges::bind_ranges(group_one, all_intergenic) %>% 
#     as_tibble() %>% 

group_one <- as_tibble(group_one) %>% 
    as_tibble() %>% 
    dplyr::arrange(type, seqnames, start) %>% 
    dplyr::distinct() %>% 
    as.data.table()
    
fwrite(
    group_one,
    "../data/annotations/group_genic_features.csv"
)

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")

##########################################################################################
#' group 2

t1 <- Sys.time()
cur.msg <- "Processing annotations for housekeeping and cancer driver elements"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

# extract housekeeping genes
# Paper: https://doi.org/10.1093/nar/gkaa609
# Downloaded from: https://housekeeping.unicamp.br/?download 
housekeeping_genes <- fread(
    "../data/annotations/Human_Housekeeping_Genes.csv",
    select = "Gene.name"
)
housekeeping_genes <- as_tibble(housekeeping_genes) %>% 
    dplyr::rename_with(~c("gene_id")) %>% 
    dplyr::arrange(gene_id)

housekeeping_genes <- as_tibble(all_genes) %>% 
    dplyr::arrange(gene_id) %>% 
    dplyr::filter(gene_id %in% housekeeping_genes$gene_id) %>% 
    dplyr::mutate(type = "Housekeeping_genes") %>% 
    plyranges::as_granges()

# gene centric fragile sites
# Paper: https://doi.org/10.1186/s12864-018-5330-5
# Downloaded from: https://webs.iiitd.edu.in/raghava/humcfs/download.html
# Reference genome: human genome Ensembl (GRCh38/hg38).
CFS_genes_files <- list.files(
    path = "../data/annotations/chromosome_bed",
    # path = "../data/annotations/fragile_site_bed",
    full.names = TRUE
)
CFS_genes_files <- stringr::str_sort(CFS_genes_files, numeric = TRUE)
CFS_genes <- lapply(CFS_genes_files, function(f){
    # return(fread(f))
    return(tail(fread(f), n = -3))
})
CFS_genes <- rbindlist(CFS_genes)
setnames(CFS_genes, c(
    "seqnames", "start", "end", 
    "gene_id", "score", "strand")
) 
CFS_genes[, score := NULL]

CFS_genes <- as_tibble(all_genes) %>% 
    dplyr::arrange(gene_id) %>% 
    dplyr::filter(gene_id %in% CFS_genes$gene_id) %>% 
    dplyr::mutate(type = "CFS_genes") %>% 
    plyranges::as_granges()

# cancer driver genes
all_cdg <- fread(
    "../data/COSMIC/Cosmic_CancerGeneCensus_v98_GRCh38.tsv",
    select = c(
        "CHROMOSOME", "GENOME_START", "GENOME_STOP", 
        "GENE_SYMBOL", "ROLE_IN_CANCER"
    )
)
all_cdg[, CHROMOSOME := paste0("chr", CHROMOSOME)]
setnames(all_cdg, c(
    "seqnames", "start", "end",
    "gene_id", "role_in_cancer"
))
all_cdg <- as_tibble(all_cdg) %>% 
    tidyr::drop_na() %>% 
    dplyr::mutate(type = "Cancer_driver_genes") %>% 
    dplyr::select(-role_in_cancer)

all_cdg <- as_tibble(all_genes) %>% 
    dplyr::arrange(gene_id) %>% 
    dplyr::filter(gene_id %in% all_cdg$gene_id) %>% 
    dplyr::mutate(type = "Cancer_driver_genes") %>% 
    plyranges::as_granges()

# oncogenes
census_labels <- fread(
    "../data/COSMIC/Cosmic_CancerGeneCensus_v98_GRCh38.tsv",
    select = c(
        "GENE_SYMBOL", "NAME", "ROLE_IN_CANCER"
    )
)
setnames(census_labels, c("gene_id", "name", "role"))

census_labels <- as_tibble(all_genes) %>% 
    dplyr::left_join(.,
        as_tibble(census_labels), 
        by = "gene_id") %>%
    tidyr::drop_na()

all_oncogene <- census_labels %>% 
    dplyr::filter(grepl("oncogene", role)) %>% 
    dplyr::select(-name, -role) %>% 
    dplyr::mutate(type = "Oncogenes")

all_oncogene <- as_tibble(all_genes) %>% 
    dplyr::arrange(gene_id) %>% 
    dplyr::filter(gene_id %in% all_oncogene$gene_id) %>% 
    dplyr::mutate(type = "Oncogenes") %>% 
    plyranges::as_granges()

# tumour suppressor genes
all_TSG <- census_labels %>% 
    dplyr::filter(grepl("TSG", role)) %>% 
    dplyr::select(-name, -role) %>% 
    dplyr::mutate(type = "Tumour_Suppressor_Genes")

all_TSG <- as_tibble(all_genes) %>% 
    dplyr::arrange(gene_id) %>% 
    dplyr::filter(gene_id %in% all_TSG$gene_id) %>% 
    dplyr::mutate(type = "Tumour_Suppressor_Genes") %>% 
    plyranges::as_granges()

# fusion genes
all_fusion <- census_labels %>% 
    dplyr::filter(grepl("fusion", role)) %>% 
    dplyr::select(-name, -role) %>% 
    dplyr::mutate(type = "Fusion_Genes")

all_fusion <- as_tibble(all_genes) %>% 
    dplyr::arrange(gene_id) %>% 
    dplyr::filter(gene_id %in% all_fusion$gene_id) %>% 
    dplyr::mutate(type = "Fusion_Genes") %>% 
    plyranges::as_granges()

# combine all genomic annotations
group_two <- plyranges::bind_ranges(
        all_genes, housekeeping_genes, CFS_genes, 
        all_cdg, all_oncogene, all_TSG, all_fusion
    ) %>% 
    dplyr::select(type, gene_id)

seqlevels(group_two, pruning.mode = "coarse") <- paste0("chr", 1:22)

group_two <- as_tibble(group_two) %>% 
    dplyr::arrange(type, seqnames, start) %>% 
    dplyr::distinct() %>% 
    as.data.table()

fwrite(
    group_two,
    "../data/annotations/group_genes.csv"
)

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")

##########################################################################################
#' group 3

t1 <- Sys.time()
cur.msg <- "Processing annotations for regulatory elements"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

# UCSC Genome Browser Table. 
all_cpg <- fread("../data/annotations/output_CpG_Islands.csv")
setnames(all_cpg, c(
    "seqnames", "start", "end", 
    "name", "length", "cpgNum", "gcNum", 
    "perCpg", "perGc", "obsExp"
))
all_cpg <- plyranges::as_granges(all_cpg) %>% 
    dplyr::mutate(type = "CpG") %>% 
    dplyr::select(type)

# centromere and satellites.
# From https://github.com/marbl/CHM13/tree/master#repeat-annotation.
# Downloaded from https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.0.bed
# Downloaded from UCSC Table Browser for T2T genome version.
all_centromere <- fread(
    "../data/annotations/chm13v2.0_censat_v2.0.bed", 
    header = FALSE, select = paste0("V", 1:4)
)
setnames(all_centromere, c("seqnames", "start", "end", "name"))
all_centromere[, name := stringr::str_extract(name, "^[^_]*")]
all_centromere <- plyranges::as_granges(all_centromere) %>% 
    reduce(.) %>% 
    dplyr::mutate(type = "Centromere_and_Pericentromere")

# telomeres
# From https://github.com/marbl/CHM13/tree/master#repeat-annotation.
# Download from: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_telomere.bed
all_telomere <- fread("../data/annotations/chm13v2.0_telomere.bed")
setnames(all_telomere, c("seqnames", "start", "end"))
all_telomere <- plyranges::as_granges(all_telomere) %>% 
    dplyr::mutate(type = "Telomere")

# G4 sites from PQS data - which genome version???
# Original file is based on hg19. Liftover from hg19 to T2T genome version.
all_g4 <- fread("../data/annotations/G4_PQS.csv")
all_g4 <- as_tibble(all_g4) %>% 
    dplyr::rename_with(~c("seqnames", "start", "end", "strand")) %>% 
    dplyr::mutate(type = "G4_sites") %>% 
    plyranges::as_granges()

hg38_chain <- import.chain("../../data/liftover/hg19ToHg38.over.chain")
all_g4 <- liftOver(all_g4, hg38_chain)
all_g4 <- unlist(as(all_g4, "GRangesList"))
all_g4 <- liftOver(all_g4, t2t_chain)
all_g4 <- unlist(as(all_g4, "GRangesList"))

# hg17 isochores from UCSC, then liftover to hg38, then to t2t
# isochore family definitions: https://doi.org/10.1016/S0378-1119(99)00485-0
# also referenced here: https://doi.org/10.1093/molbev/msi115
# Download from UCSC Browser Table.
all_iso <- rtracklayer::import.bb("../data/annotations/iso_hg17.bb")
all_iso <- all_iso %>% 
    dplyr::select(name) %>% 
    as_tibble() %>% 
    dplyr::mutate(
        gc = stringr::str_extract(
            string = name, 
            pattern = "(?<=\\=)[^=]+$"
        ),
        gc = as.numeric(gc)
    ) %>%
    dplyr::select(-name) %>% 
    dplyr::mutate(
        type = dplyr::case_when(
            gc < 37 ~ "Isochore_L1",
            gc >= 37 & gc <= 42 ~ "Isochore_L2",
            gc > 42 & gc <= 47 ~ "Isochore_H1",
            gc > 47 & gc <= 52 ~ "Isochore_H2",
            gc > 52 ~ "Isochore_H3"
        )
    ) %>% 
    plyranges::as_granges()

hg17_chain <- import.chain("../../data/liftover/hg17ToHg38.over.chain")
all_iso <- liftOver(all_iso, hg17_chain)
all_iso <- unlist(as(all_iso, "GRangesList"))
all_iso <- liftOver(all_iso, t2t_chain)
all_iso <- unlist(as(all_iso, "GRangesList"))

# combine all data sets
group_three <- plyranges::bind_ranges(
        all_cpg, all_centromere, all_telomere, 
        all_g4, all_iso
    ) %>% 
    dplyr::select(type)
seqlevels(group_three, pruning.mode = "coarse") <- paste0("chr", 1:22)

group_three <- as_tibble(group_three) %>% 
    dplyr::arrange(type, seqnames, start) %>% 
    dplyr::distinct() %>% 
    as.data.table()

fwrite(
    group_three,
    "../data/annotations/group_regulatory.csv"
)

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")

##########################################################################################
#' group 4
# From https://github.com/marbl/CHM13/tree/master#repeat-annotation
# Download from: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed

# otherwise, do the below.
# repeat elements: help from biostars https://www.biostars.org/p/312097/
#' wget http://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz
#' gunzip hg38.fa.out.gz
#' convert2bed --input=rmsk < hg38.fa.out | cut -f1-3,11 > hg38.fa.out.bed

t1 <- Sys.time()
cur.msg <- "Processing annotations for repeat elements"
cat(paste0(cur.msg, paste0(rep(".", 70-nchar(cur.msg)), collapse = "")))

all_repeats <- fread(
    "../data/annotations/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed", 
    select = c(paste0("V", 1:3), paste0("V", 6:8)),
    showProgress = FALSE
)
setnames(all_repeats, c(
    "seqnames", "start", "end", 
    "strand", "repeat_class", "subclass"
))
all_repeats <- all_repeats[repeat_class != "Unknown"]
all_repeats <- all_repeats[repeat_class != "Beta"]

group_four <- as_tibble(all_repeats) %>% 
    dplyr::mutate(
        subclass = dplyr::case_when(
            grepl("LINE", repeat_class) ~ stringr::str_replace_all(subclass, "[^0-9]", ""),
            TRUE ~ subclass
        ),
        repeat_class = dplyr::case_when(
            grepl("^LINE", repeat_class) & subclass == 1 ~ paste0(repeat_class, "_", subclass),
            grepl("^LINE", repeat_class) & subclass == 2 ~ paste0(repeat_class, "_", subclass),
            grepl("^LINE", repeat_class) ~ paste0(repeat_class, "_Others"),
            TRUE ~ repeat_class
        ),
        repeat_class = dplyr::case_when(
            grepl("^SINE", repeat_class) & subclass == "Alu" ~ paste0(repeat_class, "_", subclass),
            grepl("^SINE", repeat_class) & subclass == "_tRNA" ~ paste0(repeat_class, "_tRNA"),
            TRUE ~ repeat_class
        ),
        repeat_class = dplyr::case_when(
            # grepl("^DNA", repeat_class) & subclass == "Crypton" ~ paste0(repeat_class, "_Crypton"),
            # grepl("^DNA", repeat_class) & subclass == "Merlin" ~ paste0(repeat_class, "_Merlin"),
            # grepl("^DNA", repeat_class) & subclass == "Kolobok" ~ paste0(repeat_class, "_Kolobok"),
            # grepl("^DNA", repeat_class) & grepl("^hAT", subclass) ~ "DNA_hATs",
            # grepl("^DNA", repeat_class) & grepl("^TcMar", subclass) ~ "DNA_TcMars",
            grepl("^DNA", repeat_class) ~ "DNA_Others",
            grepl("srpRNA|snRNA|scRNA", repeat_class) ~ "RNA_Others",
            grepl("^RC", repeat_class) ~ "DNA_relaxed_circular",
            TRUE ~ repeat_class
        )
    ) %>% 
    dplyr::select(-subclass) %>% 
    dplyr::rename(type = repeat_class) %>% 
    plyranges::as_granges()  

# as_tibble(all_repeats) %>%
#     dplyr::group_by(repeat_class) %>%
#     dplyr::summarise(count = dplyr::n()) %>%
#     dplyr::arrange(desc(count))

# as_tibble(group_four) %>%
#     dplyr::group_by(type) %>%
#     dplyr::summarise(count = dplyr::n()) %>%
#     dplyr::arrange(desc(count))

# temp <- as_tibble(group_four) %>% 
#     dplyr::filter(type == "tRNA") %>% 
#     dplyr::filter(seqnames %in% paste0("chr", 1:22)) %>% 
#     dplyr::distinct() %>% 
#     plyranges::as_granges()     
# temp_seq <- getSeq(refseq, temp)
# temp_content <- Biostrings::letterFrequency(temp_seq, letters="ATGC", OR=0)
# temp_content <- temp_content / width(temp_seq)
# gc <- (temp_content[,"G"] + temp_content[,"C"])
# quantile(gc)

seqlevels(group_four, pruning.mode = "coarse") <- paste0("chr", 1:22)
group_four <- as_tibble(group_four) %>% 
    dplyr::arrange(type, seqnames, start) %>% 
    dplyr::distinct() %>% 
    as.data.table()

fwrite(
    group_four,
    "../data/annotations/group_repeat_elements.csv"
)

total.time <- Sys.time() - t1
cat("DONE! --", signif(total.time[[1]], 2), 
    attr(total.time, "units"), "\n")