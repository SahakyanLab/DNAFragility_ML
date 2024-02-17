# Load in all the cosmic mutations, work out a grouping, then do mutations and
# plot deltaRTs with those groupings
suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(suppressWarnings(library(plyranges)))

get_ranges_from_breaks <- function(coding){
    # extract the short, medium and long range sequences before a mutation
    file.name <- ifelse(coding, "Cosmic_MutantCensus_PROCESSED", "Cosmic_NonCodingVariants_PROCESSED")

    df <- fread(
        file = paste0("../data/COSMIC/", file.name, ".tsv"),
        select = {
            if(coding){
                c(
                    "PRIMARY_SITE",
                    "PRIMARY_HISTOLOGY",
                    "CHROMOSOME",
                    "GENOME_START",
                    "GENOME_STOP",
                    "GENOMIC_WT_ALLELE",
                    "GENOMIC_MUT_ALLELE",
                    "HGVSG",
                    "COSMIC_SAMPLE_ID",
                    "GENE_SYMBOL",
                    "MUTATION_DESCRIPTION"
                )
            } else {
                c(
                    "PRIMARY_SITE",
                    "PRIMARY_HISTOLOGY",
                    "CHROMOSOME",
                    "GENOME_START",
                    "GENOME_STOP",
                    "GENOMIC_WT_SEQ",
                    "GENOMIC_MUT_SEQ",
                    "HGVSG",
                    "COSMIC_SAMPLE_ID"
                )
            }
        },
        showProgress = FALSE
    )

    if(coding){
        setnames(df, c(
            c(
                "PRIMARY_SITE",
                "PRIMARY_HISTOLOGY",
                "CHROMOSOME",
                "GENOME_START",
                "GENOME_STOP",
                "GENOMIC_WT_SEQ",
                "GENOMIC_MUT_SEQ",
                "HGVSG",
                "COSMIC_SAMPLE_ID",
                "GENE_SYMBOL",
                "MUTATION_DESCRIPTION"
            )
        ))
    }

    # clean HGVSG column
    df[, HGVSG := gsub(pattern = "[0-9:]+g\\.", replacement = "", x = HGVSG)]
    # extract the type of mutation and mutation length
    df[, `:=`(
        MUTATION_TYPE = sapply(strsplit(x = HGVSG, "[^a-z]+"), `[[`, 2),
        MUTATION_LENGTH = GENOME_STOP-GENOME_START+1,
        MUTATION_START = GENOME_START,
        CHROMOSOME = paste0("chr", CHROMOSOME)
    )]

    ###########################################################################
    #' genome_position: always plus strand. Thus, if working on 
    #' the minus strand, specify "-" in the getSeq function, it gets
    #' the reverse complemented range of genome_position and outputs
    #' the results in the 5' -> 3' direction of the minus strand.
    ###########################################################################
    # focusing only on deletions and partially inversions and delins
    inv.del.dup.ind <- which(grepl(pattern = "inv|delins|del$", x = df$MUTATION_TYPE))
    if(length(inv.del.dup.ind)){
        # extract the right breakpoint start position
        df[inv.del.dup.ind, GENOME_STOP := ifelse(MUTATION_LENGTH > 1, GENOME_STOP+1, 0)]
    }

    # focusing only on insertions and duplications
    # only one bond is broken during the insertion of a base or multiple bases.
    dup.ins.ind <- which(grepl(pattern = "dup|^ins", x = df$MUTATION_TYPE))
    if(length(dup.ins.ind) > 0){
        df[dup.ins.ind, `:=`(
            GENOME_START = GENOME_START+1,
            GENOME_STOP = 0
        )]
    }

    # if MutationEndPos is not zero, append that data point as an additional row
    to.append.rows <- copy(df[GENOME_STOP != 0, -"GENOME_START"])
    setnames(to.append.rows, gsub("GENOME_STOP", "GENOME_START", colnames(to.append.rows)))
    df[, GENOME_STOP := NULL]
    df <- rbind(df, to.append.rows)

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
    if(coding){
        df[, MUTATION_DESCRIPTION := ifelse(
            MUTATION_DESCRIPTION == "", 
            "NS", MUTATION_DESCRIPTION
        )]
        df[, MUTATION_DESCRIPTION := gsub(" ", "-", MUTATION_DESCRIPTION)]
        df[, MUTATION_DESCRIPTION := gsub("_", "-", MUTATION_DESCRIPTION)]

        df[, ID := paste0(
            CHROMOSOME, "_", GENOME_START, "_",
            PRIMARY_SITE, "_", PRIMARY_HISTOLOGY, "_",
            MUTATION_TYPE, "_", 
            COSMIC_SAMPLE_ID
        )]
    } else {
        df[, ID := paste0(
            CHROMOSOME, "_", GENOME_START, "_",
            PRIMARY_SITE, "_", PRIMARY_HISTOLOGY, "_",
            MUTATION_TYPE, "_", COSMIC_SAMPLE_ID
        )]
    }
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
    #     # before_left_control = paste0(
    #     #     CHROMOSOME, ":", GENOME_START-ceiling(ranges["Long"]/2)-1000, "-", 
    #     #     GENOME_START-ceiling(ranges["Long"]/2)
    #     # ),
    #     # before_right_control = paste0(
    #     #     CHROMOSOME, ":", GENOME_START+ceiling(ranges["Long"]/2), "-", 
    #     #     GENOME_START+ceiling(ranges["Long"]/2)+1000
    #     # )
    # )] %>% suppressWarnings()

    # # extract sequences
    # df[, `:=`(
    #     before_short = paste0(getSeq(refseq, GRanges(before_short))),
    #     before_medium = paste0(getSeq(refseq, GRanges(before_medium))),
    #     before_long = paste0(getSeq(refseq, GRanges(before_long)))
    #     # before_left_control = paste0(getSeq(refseq, GRanges(before_left_control))),
    #     # before_right_control = paste0(getSeq(refseq, GRanges(before_right_control)))
    # )]

    # # extract the short medium and long ranges after mutation
    # # focusing only on deletions
    # del.delins.ind <- which(grepl(pattern = "del$|delins", x = df$MUTATION_TYPE))
    # if(length(del.delins.ind) > 0){
    #     df[del.delins.ind, MUTATION_LENGTH := nchar(GENOMIC_WT_SEQ)]
    #     df[del.delins.ind, after_short := ifelse(GENOME_START == MUTATION_START,
    #         paste0(
    #             substring(before_short, 1, kmer/2),
    #             paste0(getSeq(refseq, GRanges(paste0(
    #                 CHROMOSOME, ":", GENOME_START+MUTATION_LENGTH, "-",
    #                 (GENOME_START+MUTATION_LENGTH)+(kmer/2-1)
    #             ))))
    #         ),
    #         paste0(
    #             paste0(getSeq(refseq, GRanges(paste0(
    #                 CHROMOSOME, ":", MUTATION_START-kmer/2, "-", 
    #                 MUTATION_START-1
    #             )))),
    #             substring(before_short, kmer/2+1, kmer)
    #         )
    #     )]
    #     df[del.delins.ind, `:=`(
    #         after_medium = ifelse(GENOME_START == MUTATION_START,
    #             paste0(
    #                 substring(before_medium, 1, ceiling(ranges["Medium"]/2)-(kmer/2)),
    #                 after_short,
    #                 paste0(getSeq(refseq, GRanges(paste0(
    #                     CHROMOSOME, ":", GENOME_START+MUTATION_LENGTH, "-",
    #                     (GENOME_START+MUTATION_LENGTH)+(nchar(before_medium)-ceiling(ranges["Medium"]/2)-(kmer/2)-1)
    #                 ))))
    #             ),
    #             paste0(
    #                 paste0(getSeq(refseq, GRanges(paste0(
    #                     CHROMOSOME, ":", 
    #                     (MUTATION_START-(nchar(before_medium)-ceiling(ranges["Medium"]/2)-1)),
    #                     "-",
    #                     MUTATION_START-(kmer/2)-1
    #                 )))),
    #                 after_short,
    #                 substring(before_medium, ceiling(ranges["Medium"]/2)+(kmer/2)+1)
    #             )                
    #         ),
    #         after_long = ifelse(GENOME_START == MUTATION_START,
    #             paste0(
    #                 substring(before_long, 1, ceiling(ranges["Long"]/2)-(kmer/2)),
    #                 after_short,
    #                 paste0(getSeq(refseq, GRanges(paste0(
    #                     CHROMOSOME, ":", GENOME_START+MUTATION_LENGTH, "-",
    #                     (GENOME_START+MUTATION_LENGTH)+(nchar(before_long)-ceiling(ranges["Long"]/2)-(kmer/2)-1)
    #                 ))))
    #             ),
    #             paste0(
    #                 paste0(getSeq(refseq, GRanges(paste0(
    #                     CHROMOSOME, ":", 
    #                     (MUTATION_START-(nchar(before_long)-ceiling(ranges["Long"]/2)-1)),
    #                     "-",
    #                     MUTATION_START-(kmer/2)-1
    #                 )))),
    #                 after_short,
    #                 substring(before_long, ceiling(ranges["Long"]/2)+(kmer/2)+1)
    #             )
    #         )
    #     )]
    # }

    # # test=as.data.table(table(df$MUTATION_LENGTH))
    # # test[, `:=`(V1 = as.numeric(V1), frac = N/sum(N)*100, cum.frac = cumsum(N/sum(N)*100))]
    # # setnames(test, c("len", "count", "frac", "cum.frac"))
    # # setorder(test, -cum.frac)
    # # test[len >= 500]

    # # extract the short medium and long ranges after mutation
    # # focusing only on insertions and duplications
    # ins.dup.ind <- which(grepl(pattern = "^ins|dup", x = df$MUTATION_TYPE))
    # if(length(ins.dup.ind) > 0){
    #     df[ins.dup.ind, MUTATION_LENGTH := nchar(GENOMIC_MUT_SEQ)]
    #     df[ins.dup.ind, `:=`(
    #         after_short = paste0(
    #             substring(before_short, 1, kmer/2),
    #             ifelse(MUTATION_LENGTH >= (kmer/2), 
    #                 substring(GENOMIC_MUT_SEQ, 1, kmer/2), 
    #                 paste0(
    #                     GENOMIC_MUT_SEQ, 
    #                     substring(before_short, kmer/2+1, (kmer/2+1)+(kmer-kmer/2-MUTATION_LENGTH-1))
    #                 )
    #             )
    #         ),
    #         after_medium = paste0(
    #             substring(before_medium, 1, ceiling(ranges["Medium"]/2)),
    #             ifelse(MUTATION_LENGTH >= ceiling(ranges["Medium"]/2), 
    #                 substring(GENOMIC_MUT_SEQ, 1, ceiling(ranges["Medium"]/2)),
    #                 paste0(
    #                     GENOMIC_MUT_SEQ,
    #                     substring(
    #                         before_medium, 
    #                         ceiling(ranges["Medium"]/2)+1,
    #                         nchar(before_medium)-MUTATION_LENGTH
    #                     )
    #                 )
    #             )
    #         ),
    #         after_long = paste0(
    #             substring(before_long, 1, ceiling(ranges["Long"]/2)),
    #             ifelse(MUTATION_LENGTH >= ceiling(ranges["Long"]/2), 
    #                 substring(GENOMIC_MUT_SEQ, 1, ceiling(ranges["Long"]/2)),
    #                 paste0(
    #                     GENOMIC_MUT_SEQ,
    #                     substring(
    #                         before_long, 
    #                         ceiling(ranges["Long"]/2)+1,
    #                         nchar(before_long)-MUTATION_LENGTH
    #                     )
    #                 )
    #             )
    #         )
    #     )]
    # }

    # # focusing only on inversions
    # inv.ind <- which(grepl(pattern = "inv", x = df$MUTATION_TYPE))
    # if(length(inv.ind) > 0){
    #     df[inv.ind, MUTATION_LENGTH := nchar(GENOMIC_MUT_SEQ)]
    #     df[inv.ind, `:=`(
    #         after_short = ifelse(GENOME_START == MUTATION_START,
    #             paste0(
    #                 substring(before_short, 1, kmer/2),
    #                 ifelse(MUTATION_LENGTH >= (kmer/2), 
    #                     substring(GENOMIC_MUT_SEQ, 1, kmer/2), 
    #                     paste0(
    #                         GENOMIC_MUT_SEQ, 
    #                         substring(before_short, nchar(before_short)-MUTATION_LENGTH+1)
    #                     )
    #                 )
    #             ),
    #             paste0(
    #                 ifelse(MUTATION_LENGTH >= (kmer/2), 
    #                     substring(GENOMIC_MUT_SEQ, nchar(GENOMIC_MUT_SEQ)-(kmer/2)+1),
    #                     paste0(
    #                         substring(before_short, 1, nchar(before_short)-(kmer/2)-MUTATION_LENGTH),
    #                         GENOMIC_MUT_SEQ
    #                     )
    #                 ),
    #                 substring(before_short, kmer/2+1)
    #             )
    #         ),
    #         after_medium = ifelse(GENOME_START == MUTATION_START, 
    #             paste0(
    #                 substring(before_medium, 1, ceiling(ranges["Medium"]/2)),
    #                 ifelse(MUTATION_LENGTH >= ceiling(ranges["Medium"]/2), 
    #                     substring(GENOMIC_MUT_SEQ, 1, ceiling(ranges["Medium"]/2)),
    #                     paste0(
    #                         GENOMIC_MUT_SEQ,
    #                         substring(
    #                             before_medium, 
    #                             ceiling(ranges["Medium"]/2)+(kmer/2)-1
    #                         )
    #                     )
    #                 )
    #             ),
    #             paste0(
    #                 substring(
    #                     before_medium, 1, ceiling(ranges["Medium"]/2)-(MUTATION_LENGTH)
    #                 ),
    #                 ifelse(MUTATION_LENGTH >= ceiling(ranges["Medium"]/2), 
    #                     substring(GENOMIC_MUT_SEQ, 1, ceiling(ranges["Medium"]/2)),
    #                     paste0(
    #                         GENOMIC_MUT_SEQ,
    #                         substring(
    #                             before_medium, 
    #                             ceiling(ranges["Medium"]/2)+1
    #                         )
    #                     )
    #                 )
    #             )
    #         ),
    #         after_long = ifelse(GENOME_START == MUTATION_START, 
    #             paste0(
    #                 substring(before_long, 1, ceiling(ranges["Long"]/2)),
    #                 ifelse(MUTATION_LENGTH >= ceiling(ranges["Long"]/2), 
    #                     substring(GENOMIC_MUT_SEQ, 1, ceiling(ranges["Long"]/2)),
    #                     paste0(
    #                         GENOMIC_MUT_SEQ,
    #                         substring(
    #                             before_long, 
    #                             ceiling(ranges["Long"]/2)+(kmer/2)-1
    #                         )
    #                     )
    #                 )
    #             ),
    #             paste0(
    #                 substring(
    #                     before_long, 1, ceiling(ranges["Long"]/2)-(MUTATION_LENGTH)
    #                 ),
    #                 ifelse(MUTATION_LENGTH >= ceiling(ranges["Long"]/2), 
    #                     substring(GENOMIC_MUT_SEQ, 1, ceiling(ranges["Long"]/2)),
    #                     paste0(
    #                         GENOMIC_MUT_SEQ,
    #                         substring(
    #                             before_long, 
    #                             ceiling(ranges["Long"]/2)+1
    #                         )
    #                     )
    #                 )
    #             )
    #         )
    #     )]
    # }

    # # focusing only on the insertion process after a deletion
    # delins.ind <- which(grepl(pattern = "delins", x = df$MUTATION_TYPE))
    # if(length(delins.ind) > 0){
    #     df[delins.ind, MUTATION_LENGTH := nchar(GENOMIC_MUT_SEQ)]
    #     df[delins.ind, `:=`(
    #         after_short = ifelse(GENOME_START == MUTATION_START,
    #             paste0(
    #                 substring(after_short, 1, kmer/2),
    #                 ifelse(MUTATION_LENGTH >= (kmer/2), 
    #                     substring(GENOMIC_MUT_SEQ, 1, kmer/2), 
    #                     paste0(
    #                         GENOMIC_MUT_SEQ, 
    #                         substring(after_short, nchar(after_short)-MUTATION_LENGTH+1)
    #                     )
    #                 )
    #             ),
    #             paste0(
    #                 ifelse(MUTATION_LENGTH >= (kmer/2), 
    #                     substring(GENOMIC_MUT_SEQ, nchar(GENOMIC_MUT_SEQ)-(kmer/2)+1),
    #                     paste0(
    #                         substring(after_short, 1, nchar(after_short)-(kmer/2)-MUTATION_LENGTH),
    #                         GENOMIC_MUT_SEQ
    #                     )
    #                 ),
    #                 substring(after_short, kmer/2+1)
    #             )
    #         ),
    #         after_medium = ifelse(GENOME_START == MUTATION_START, 
    #             paste0(
    #                 substring(after_medium, 1, ceiling(ranges["Medium"]/2)),
    #                 ifelse(MUTATION_LENGTH >= ceiling(ranges["Medium"]/2), 
    #                     substring(GENOMIC_MUT_SEQ, 1, ceiling(ranges["Medium"]/2)),
    #                     paste0(
    #                         GENOMIC_MUT_SEQ,
    #                         substring(
    #                             after_medium, 
    #                             ceiling(ranges["Medium"]/2)+1,
    #                             nchar(after_medium)-MUTATION_LENGTH
    #                         )
    #                     )
    #                 )
    #             ),
    #             paste0(
    #                 substring(
    #                     after_medium, 
    #                     1+MUTATION_LENGTH,
    #                     ceiling(ranges["Medium"]/2)
    #                 ),
    #                 ifelse(MUTATION_LENGTH >= ceiling(ranges["Medium"]/2), 
    #                     substring(GENOMIC_MUT_SEQ, 1, ceiling(ranges["Medium"]/2)),
    #                     paste0(
    #                         GENOMIC_MUT_SEQ,
    #                         substring(
    #                             after_medium, 
    #                             ceiling(ranges["Medium"]/2)+1
    #                         )
    #                     )
    #                 )
    #             )
    #         ),
    #         after_long = ifelse(GENOME_START == MUTATION_START, 
    #             paste0(
    #                 substring(after_long, 1, ceiling(ranges["Long"]/2)),
    #                 ifelse(MUTATION_LENGTH >= ceiling(ranges["Long"]/2), 
    #                     substring(GENOMIC_MUT_SEQ, 1, ceiling(ranges["Long"]/2)),
    #                     paste0(
    #                         GENOMIC_MUT_SEQ,
    #                         substring(
    #                             after_long, 
    #                             ceiling(ranges["Long"]/2)+1,
    #                             nchar(after_long)-MUTATION_LENGTH
    #                         )
    #                     )
    #                 )
    #             ),
    #             paste0(
    #                 substring(
    #                     after_long, 
    #                     1+MUTATION_LENGTH,
    #                     ceiling(ranges["Long"]/2)
    #                 ),
    #                 ifelse(MUTATION_LENGTH >= ceiling(ranges["Long"]/2), 
    #                     substring(GENOMIC_MUT_SEQ, 1, ceiling(ranges["Long"]/2)),
    #                     paste0(
    #                         GENOMIC_MUT_SEQ,
    #                         substring(
    #                             after_long, 
    #                             ceiling(ranges["Long"]/2)+1
    #                         )
    #                     )
    #                 )
    #             )
    #         )
    #     )]
    # }

    #' The above sequence attraction processes are all correct, however, 
    #' for works requiring a single nucleotide breakpoint position, we need to 
    #' shift by one position upstream.
    #' DEL: GENOME_START-1
    #' INS: GENOME_START-1
    #' DUP: GENOME_START-1
    #' INV: GENOME_START-1
    #' DELINS: GENOME_START-1
    df[, GENOME_START := GENOME_START-1]

    df <- df[, .(
        ID, GENOME_START
    )]
    return(df)
}

args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
setwd(my.path)

# global vars
kmer <- 8
refseq <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

max_range <- fread(paste0(
    "../../03_Breakpoints_v2/03_FitCurves/data/",
    "ranges/MaxValuesFromClustersByType.csv"
))
max_range <- apply(max_range[, -"break_type"], 2, max)
round_to_nearest_even <- function(x) round(x/2)*2
ranges <- round_to_nearest_even(max_range)
names(ranges) <- stringr::str_to_title(names(ranges))

df.coding <- get_ranges_from_breaks(coding = TRUE)
df.noncoding <- get_ranges_from_breaks(coding = FALSE)
df.all <- rbind(df.noncoding, df.coding)
df.all <- distinct(df.all, ID, .keep_all = TRUE)

dir.create(
    path = "../data/BaseTable/",
    showWarnings = FALSE,
    recursive = TRUE
)

# save tables
fwrite(
    df.coding,
    file = paste0(
        "../data/BaseTable/CODING_kmer_", kmer, ".csv"
    )
)

fwrite(
    df.noncoding,
    file = paste0(
        "../data/BaseTable/NONCODING_kmer_", kmer, ".csv"
    )
)

fwrite(
    df.all,
    file = paste0(
        "../data/BaseTable/kmer_", kmer, ".csv"
    )
)

df.all <- df.all[, .(ID)]
source("./00_GetSVs.R")
df.sv <- fread(
    paste0(
        "../data/BaseTable/SV_kmer_", kmer, ".csv"
    )
)
df.sv <- df.sv[, .(ID)]
df.combined <- rbind(df.all, df.sv)
fwrite(
    df.combined,
    file = paste0(
        "../data/BaseTable/kmer_", kmer, ".csv"
    )
)