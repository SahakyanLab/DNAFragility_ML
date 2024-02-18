# Load in all the cosmic mutations, work out a grouping, then do mutations and
# plot deltaRTs with those groupings
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))
suppressPackageStartupMessages(suppressWarnings(library(BSgenome.Hsapiens.UCSC.hg38)))

get_SVs <- function(){
    df <- fread("../data/COSMIC/Cosmic_Breakpoints_v98_GRCh38.tsv", showProgress = FALSE)
    df <- as_tibble(df)
    
    class.table <- fread("../data/COSMIC/Cosmic_Classification_v98_GRCh38.tsv", showProgress = FALSE)
    class.table <- as_tibble(class.table)

    # join and table
    df <- dplyr::left_join(
        x = class.table,
        y = df,
        by = "COSMIC_PHENOTYPE_ID"
    )
    df <- df %>% 
        dplyr::mutate(
            CHROM_FROM = as.numeric(CHROM_FROM),
            LOCATION_FROM_MIN = as.numeric(LOCATION_FROM_MIN),
            LOCATION_FROM_MAX = as.numeric(LOCATION_FROM_MAX),
            CHROM_TO = as.numeric(CHROM_TO),
            LOCATION_TO_MIN = as.numeric(LOCATION_TO_MIN),
            LOCATION_TO_MAX = as.numeric(LOCATION_TO_MAX)
        ) %>% 
        suppressWarnings()
    df <- as.data.table(df)
    df <- df[!grepl(pattern = "MT|^[X|Y]", x = CHROM_FROM)]
    df <- df[!grepl(pattern = "MT|^[X|Y]", x = CHROM_TO)]
    df <- df[!is.na(CHROM_FROM)]
    df <- df[!is.na(CHROM_TO)]

    # change name
    df[, MUTATION_TYPE := ifelse(
        grepl(pattern = "intra", x = MUTATION_TYPE), 
        "intrachr", "interchr"
    )]

    # every location is a unique breakpoint
    df.locationfrom <- df[, .(
        PRIMARY_SITE, PRIMARY_HISTOLOGY,
        MUTATION_TYPE, CHROM_FROM,
        LOCATION_FROM_MIN, LOCATION_FROM_MAX
    )]
    setnames(df.locationfrom, c(
        "PRIMARY_SITE", "PRIMARY_HISTOLOGY", "MUTATION_TYPE",
        "CHROMOSOME", "GENOME_START", "GENOME_STOP"
    ))
    df.locationfrom[, GENOME_STOP := ifelse(GENOME_START == GENOME_STOP, 0, GENOME_STOP+1)]
    df.locationfrom <- df.locationfrom[GENOME_START <= GENOME_STOP]
    df.locationfrom <- df.locationfrom[GENOME_STOP != 0]

    # if MutationEndPos is not zero, append that data point as an additional row
    to.append.rows <- copy(df.locationfrom[, -"GENOME_START"])
    setnames(to.append.rows, gsub("GENOME_STOP", "GENOME_START", colnames(to.append.rows)))
    df.locationfrom[, GENOME_STOP := NULL]
    df.locationfrom <- rbind(df.locationfrom, to.append.rows)

    # repeat with destination breakpoints
    df.locationto <- df[, .(
        PRIMARY_SITE, PRIMARY_HISTOLOGY,
        MUTATION_TYPE, CHROM_TO,
        LOCATION_TO_MIN, LOCATION_TO_MAX
    )]
    setnames(df.locationto, c(
        "PRIMARY_SITE", "PRIMARY_HISTOLOGY", "MUTATION_TYPE",
        "CHROMOSOME", "GENOME_START", "GENOME_STOP"
    ))
    df.locationto[, GENOME_STOP := ifelse(GENOME_START == GENOME_STOP, 0, GENOME_STOP+1)]
    df.locationto <- df.locationto[GENOME_START <= GENOME_STOP]
    df.locationto <- df.locationto[GENOME_STOP != 0]

    # if MutationEndPos is not zero, append that data point as an additional row
    to.append.rows <- copy(df.locationto[, -"GENOME_START"])
    setnames(to.append.rows, gsub("GENOME_STOP", "GENOME_START", colnames(to.append.rows)))
    df.locationto[, GENOME_STOP := NULL]
    df.locationto <- rbind(df.locationto, to.append.rows)

    # combine both tables
    df <- rbind(df.locationfrom, df.locationto)
    df[, CHROMOSOME := paste0("chr", CHROMOSOME)]

     # replace any dash with underscores
    df[, `:=`(
        PRIMARY_SITE = gsub(
            pattern = "_", 
            replacement = "-", 
            x = PRIMARY_SITE
        ),
        PRIMARY_HISTOLOGY = gsub(
            pattern = "_", 
            replacement = "-", 
            x = PRIMARY_HISTOLOGY
        )
    )]

    #' get a table of unique IDs.
    #'      1. Chromosome
    #'      2. Position
    #'      3. Tissue
    #'      4. Cancer
    #'      5. Breakage type
    #'      6. Sample ID
    df[, ID := paste0(
        CHROMOSOME, "_", GENOME_START, "_",
        PRIMARY_SITE, "_", PRIMARY_HISTOLOGY, "_",
        MUTATION_TYPE
    )]
    df <- distinct(df, ID, .keep_all = TRUE)

    # # extract the short medium and long ranges before mutation
    # df[, `:=`(
    #     before_short = paste0(
    #         CHROMOSOME, ":", GENOME_START-kmer/2, "-", 
    #         GENOME_START-kmer/2+kmer-1
    #     ),
    #     before_medium = paste0(
    #         CHROMOSOME, ":", GENOME_START-ceiling(ranges["Medium"]/2), "-", 
    #         GENOME_START+ceiling(ranges["Medium"]/2)
    #     ),
    #     before_long = paste0(
    #         CHROMOSOME, ":", GENOME_START-ceiling(ranges["Long"]/2), "-", 
    #         GENOME_START+ceiling(ranges["Long"]/2)
    #     )
    # )] %>% suppressWarnings()

    # # get the start and end positions of each chromosome
    # refseq.table <- as.data.frame(refseq@seqinfo)
    # refseq.table <- refseq.table[grepl(
    #     pattern = "^chr([1-9]|1[0-9]|2[0-2])$", 
    #     x = rownames(refseq.table)
    # ),]
    # refseq.table <- data.table(
    #     seqnames = rownames(refseq.table),
    #     start = refseq.table$seqlengths
    # )
    
    # df[, CHECK_BOUNDARY := GENOME_START+ceiling(ranges["Long"]/2)]
    # df.split <- split(df, by = "CHROMOSOME")
    # df.split.filtered <- lapply(1:22, function(x){
    #     OOB.lim <- refseq.table[seqnames == paste0("chr", x)]$start
    #     df.split[[paste0("chr", x)]][, CHECK_BOUNDARY := CHECK_BOUNDARY <= OOB.lim]
    #     temp <- df.split[[paste0("chr", x)]][CHECK_BOUNDARY == TRUE]
    #     temp[, CHECK_BOUNDARY := NULL]
    #     return(temp)
    # })
    # df <- rbindlist(df.split.filtered) 

    # # extract sequences
    # df[, `:=`(
    #     before_short = paste0(getSeq(refseq, GRanges(before_short))),
    #     before_medium = paste0(getSeq(refseq, GRanges(before_medium))),
    #     before_long = paste0(getSeq(refseq, GRanges(before_long)))
    # )]

    #' The above sequence attraction processes are all correct, however, 
    #' for works requiring a single nucleotide breakpoint position, we need to 
    #' shift by one position upstream.
    #' DEL: GENOME_START-1
    #' INS: GENOME_START-1
    #' DUP: GENOME_START-1
    #' INV: GENOME_START-1
    # #' DELINS: GENOME_START-1
    df[, GENOME_START := GENOME_START-1]

    df <- df[, .(
        ID, GENOME_START
    )]

    # # extract the short medium and long ranges after mutation
    # # focusing only on deletions
    # del.delins.ind <- which(grepl(pattern = "del$|delins", x = df$MUTATION_TYPE))
    # if(length(del.delins.ind) > 0){
    #     df[del.delins.ind, MUTATION_LENGTH := nchar(GENOMIC_WT_SEQ)]
    #     df[del.delins.ind, after_short := paste0(
    #         substring(before_short, 1, kmer/2),
    #         paste0(getSeq(refseq, GRanges(paste0(
    #             CHROMOSOME, ":", GENOME_START+MUTATION_LENGTH, "-",
    #             (GENOME_START+MUTATION_LENGTH)+(kmer/2-1)
    #         ))))
    #     )]
    #     df[del.delins.ind, `:=`(
    #         after_medium = paste0(
    #             substring(before_medium, 1, ceiling(ranges["Medium"]/2)-(kmer/2)),
    #             after_short,
    #             paste0(getSeq(refseq, GRanges(paste0(
    #                 CHROMOSOME, ":", GENOME_START+MUTATION_LENGTH, "-",
    #                 (GENOME_START+MUTATION_LENGTH)+(nchar(before_medium)-ceiling(ranges["Medium"]/2)-(kmer/2)-1)
    #             ))))
    #         ),            
    #         after_long = paste0(
    #                 substring(before_long, 1, ceiling(ranges["Long"]/2)-(kmer/2)),
    #             after_short,
    #             paste0(getSeq(refseq, GRanges(paste0(
    #                 CHROMOSOME, ":", GENOME_START+MUTATION_LENGTH, "-",
    #                 (GENOME_START+MUTATION_LENGTH)+(nchar(before_long)-ceiling(ranges["Long"]/2)-(kmer/2)-1)
    #             ))))
    #         )  
    #     )]
    # }

    return(df)
}

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

# global vars
kmer <- 8
refseq <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

max_range <- fread(paste0(
    "../../data/",
    "ranges/MaxValuesFromClustersByType.csv"
))
max_range <- apply(max_range[, -"break_type"], 2, max)
round_to_nearest_even <- function(x) round(x/2)*2
ranges <- round_to_nearest_even(max_range)
names(ranges) <- stringr::str_to_title(names(ranges))

df.sv <- get_SVs()

# temporarily add a unique SAMPLE_ID
df.sv[, ID := paste0(ID, "_COSS", 1:nrow(df.sv))]

# save table
fwrite(
    df.sv,
    file = paste0(
        "../data/BaseTable/SV_kmer_", kmer, ".csv"
    )
)