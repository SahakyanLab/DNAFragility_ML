# Load in all the cosmic mutations, work out a grouping, then do mutations and
# plot deltaRTs with those groupings
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(Biostrings))
pbapply::pboptions(char = "=", type = "txt")
data.table::setDTthreads(threads = 2)

# global vars
args = commandArgs(trailingOnly = TRUE)
bin_width = bw = as.integer(args[1])
overlapping_bins = as.logical(args[2])
overlap_factor = as.numeric(args[3])
remove_zero_breaks = as.logical(args[4])
my.path = as.character(args[5])
setwd(my.path)
constant = 1
kmer = 8

get_data <- function(){
    df.parsed <- fread(
        file = paste0(
            "../05_Cosmic/data/BaseTable/",
            "kmer_", kmer, ".csv"
        ),
        select = "ID",
        showProgress = TRUE
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
        dplyr::rename_with(~c("seqnames", "start"))

    df <- df %>%
        dplyr::mutate(width = 1) %>%
        dplyr::distinct() %>% 
        plyranges::as_granges() %>%
        sort()

    df <- as_tibble(df) %>% 
        dplyr::select(seqnames, start)

    return(df)
}
df <- get_data()
refseq <- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0::BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0

# liftover breakpoints to the telomere-to-telomere genome version
suppressPackageStartupMessages(suppressWarnings(library(rtracklayer)))
chain <- import.chain("../05_Cosmic/data/liftover/hg38-chm13v2.over.chain")
df <- df %>% dplyr::mutate(width = 1, strand = "+")
df <- plyranges::as_granges(df)
df <- liftOver(df, chain)
df <- unlist(as(df, "GRangesList"))
df <- as_tibble(df) %>% dplyr::select(seqnames, start)

# get the start and end positions of each chromosome
refseq.table <- as.data.frame(refseq@seqinfo)
refseq.table <- refseq.table[grepl(
    pattern = "^([1-9]|1[0-9]|2[0-2])$", 
    x = rownames(refseq.table)
),]
refseq.table <- data.table(
    chr = rownames(refseq.table),
    end = refseq.table$seqlengths
)

# find the start positions of non-'N' regions
find_non_N_positions <- function(chr){
    sequence <- getSeq(refseq, names = chr)
    ind_of_int <- match(chr, refseq.table$chr)
    chr_len <- refseq.table$end[ind_of_int]
    chunk_size <- 1000

    for(i in seq(1, chr_len, by = chunk_size)){
        chunk <- subseq(sequence, start = i, end = min(i+chunk_size-1, chr_len))
        any_non_N_positions <- letterFrequency(chunk, letters = "ACGT", OR = 0)
        if(sum(any_non_N_positions) > 0) return(i)
    }
}
refseq.table[, start := pbapply::pbsapply(chr, find_non_N_positions)]
setcolorder(refseq.table, c("chr", "start", "end"))

# count cosmic breaks in non-overlapping windows and return granges 
df_bp_granges <- pbapply::pblapply(1:22, function(chr){
    # create bins
    temp_chr <- refseq.table[chr, ] 
    bins <- seq(temp_chr$start, temp_chr$end, by = bw)
    
    # final bin to include last bases
    if(tail(bins, n = 1) < temp_chr$end) bins <- c(bins, temp_chr$end)

    # Bin genome by non-overlapping bins
    df_bp <- df %>% 
        dplyr::filter(seqnames == paste0("chr", chr)) %>% 
        dplyr::distinct() %>% 
        dplyr::mutate(width = 1) %>% 
        plyranges::as_granges()

    if(overlapping_bins){
        # bin genome into overlapping bins
        bw_half <- ceiling(bw/overlap_factor)
        bin_starts <- seq(
            temp_chr$start, 
            by = bw_half, 
            length.out = ceiling((temp_chr$end-bw+1)/bw_half)
        )
        bin_ends <- bin_starts+bw-1

        bin_ends[length(bin_ends)] <- pmax(bin_ends[length(bin_ends)], temp_chr$end)

        df_bp_granges <- GRanges(
            seqnames = paste0("chr", chr),
            ranges = IRanges(
                start = bin_starts,
                end = bin_ends
            )
        )
    } else {
        # bin genome into non-overlapping bins
        num_bins <- ceiling((temp_chr$end-temp_chr$start+1)/bw)
        df_bp_granges <- GRanges(
            seqnames = paste0("chr", chr),
            ranges = IRanges(
                start = seq(
                    temp_chr$start, 
                    by = bw, 
                    length.out = num_bins
                ),
                end = c(
                    seq(temp_chr$start+bw-1, 
                    by = bw, 
                    length.out = num_bins-1
                ), temp_chr$end)
            )
        )
    }

    df_bp_granges <- df_bp_granges %>% 
        dplyr::mutate(Breaks = countOverlaps(df_bp_granges, df_bp))

    if(remove_zero_breaks){
        df_bp_granges <- df_bp_granges %>% 
            dplyr::filter(Breaks > 0)
    }

    return(df_bp_granges)
})

if(bw == 1){
    out <- pbapply::pblapply(1:22, function(chr){
        dt <- as.data.table(df_bp_granges[[chr]])
        dt[, `:=`(seqnames = NULL, end = NULL, width = NULL, strand = NULL)]
        setnames(dt, c("start.pos", "Breaks"))
        setcolorder(dt, c("Breaks", "start.pos"))

        # split into chunks of sizes 2 million rows
        nr_rows <- nrow(dt)
        chunk_size <- 2000000
        num_chunks <- ceiling(nr_rows/chunk_size)
        dt[, chunks := rep(1:num_chunks, each = chunk_size)[1:nr_rows]]
        
        path_to_dir <- paste0(
            "../data/experiments/FullGenomeChunks/",
            "breakpoint_positions/chr", chr
        )
        dir.create(path = path_to_dir, showWarnings = FALSE, recursive = TRUE)

        for(chunk in 1:num_chunks){
            temp <- dt[chunks == chunk]
            temp[, chunks := NULL]
            
            fwrite(
                temp,
                paste0(path_to_dir, "/chunk_", chunk, ".csv"),
                showProgress = FALSE
            )
        }
    })
} else {
    df_bp_granges_all <- suppressWarnings(unlist(as(df_bp_granges, "GRangesList")))

    base_dir <- "../figures"
    df_breaks <- tibble(Breaks = log2(constant+mcols(df_bp_granges_all)$Breaks))
    p1 <- df_breaks %>% 
        ggplot(aes(x = Breaks, y = after_stat(density))) + 
        geom_line(stat = "density", linewidth = 1.5, col = "darkred", adjust = 2) + 
        theme_bw() + 
        theme_classic() + 
        labs(x = "Number of breaks (log)", y = "Density")
        
    dir.create(path = paste0(base_dir, "/DeepLearning_Trials/"), showWarnings = FALSE)
    ggsave(
        filename = paste0(
            base_dir, "/DeepLearning_Trials/", 
            "Density_bw_", bw, ".pdf"
        ),
        plot = p1
    )

    p2 <- df_breaks %>% 
        dplyr::mutate(Breaks = exp(Breaks)-constant) %>% 
        ggplot(aes(x = Breaks, y = after_stat(density))) + 
        geom_line(stat = "density", linewidth = 1.5, col = "darkred") + 
        theme_bw() + 
        theme_classic() + 
        coord_cartesian(xlim = c(0,500)) + 
        labs(x = "Number of breaks", y = "Density")

    # save for regional breakpoint prediction ML
    base_dir <- "../data/experiments/FullGenome_DeepLearning/breakpoint_positions"
    dir.create(path = base_dir, showWarnings = FALSE, recursive = TRUE)
    out <- lapply(1:22, function(chr){
        dt <- as_tibble(df_bp_granges_all) %>% 
            dplyr::filter(seqnames == paste0("chr", chr)) %>% 
            dplyr::mutate(start = start+ceiling((end-start)/2)-1) %>% 
            dplyr::select(Breaks, start) %>% 
            dplyr::rename_with(~c("Breaks", "start.pos")) %>% 
            as.data.table()

        fwrite(
            dt,
            paste0(base_dir, "/chr", chr, ".csv"),
            showProgress = FALSE
        )
    })
}