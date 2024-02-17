buildCoordinate <- function(path, case.length=NULL, k, genome.name,
                            strand.sensitive=TRUE, merge.replicate=TRUE,
                            case.pattern=NULL, chr.name=NULL, is.kmer=FALSE) {

  # Given a folder with chromosome-separated coordinate OR bedfile. This
  # function will create an S3-class object named genomic.coordinate containing
  # helpful psuedo-attributes (not in order):
  # 1. case.length         : useful for extracting k-mer in a direct or sliding
  #                          method
  # 2. k                   : size of k-mer (single digit)
  # 3. status              : A data.table of a single column: is.kmer. Reason
  #                          for data.table is to assign by reference. The
  #                          coordinate reflect the status. A function called
  #                          kmerize is written to convert inter-changebly
  #                          between is.kmer = FALSE and TRUE.A coordinate with
  #                          is.kmer=FALSE status will result in slinding k-mer
  #                          extraction.
  # 4. genome.name         : For error checking - must be the same with the
  #                          genome object name.
  # 5. rep.name            : For multiple replicate, the replicate name is
  #                          shown.
  # 6. is.strand.sensitive : A boolean whether strand polarity matters.
  #
  # path            <string>      A path to a folder containing either
  #                               (1) chromosome-separated coordinate files
  #                                   (assume replicates for subfolder) OR
  #                               (2) bedfile. (assume replicates for bedfiles)
  # case.length      <int>        Case length e.g. damage length. Default is
  #                               NULL for automatic detection or varied length.
  # k                <int>        K-mer length
  # is.kmer          <bool>       Is the coordinate refers to kmer? case or
  #                               kmer. Default is FALSE
  # chr.name        <string>      Build only a single chromosome coordinate.
  #                               Useful for looping.
  # genome.name     <string>      Genome name associated with the coordinates.
  #                               The name should be the same with the genome
  #                               object.
  # strand.sensitive <bool>       A boolean whether strand polarity matters
  # merge.replicate  <bool>       Merge replicates by removing duplicates.
  #                               If not merging, duplicates will give weight
  #                               to the kmer counting.

  suppressPackageStartupMessages( library(data.table)  )

  # Convert bedfile to chromosome-separated csv folder setting
  bedfile <- grepl("\\.bed|BED|narrowpeak|Narrowpeak|broadPeak.*$", list.files(path)[1])
  if (bedfile) {

    # Make a temporary folder containing chromosome-separated coordinates
    temp.root <- paste0("temp_", floor(runif(1, max = 100000)), "/")

    for (rep in seq_along(list.files(path))) {
      temp.path <- paste0(temp.root, "rep_", rep, "/")
      dir.create(temp.path, recursive = TRUE, showWarnings = FALSE)
      bedToCoor(path = paste0(path, "/chr", rep, ".bed"), temp.path)
    }

    # Update path to temp.path
    path <- temp.root
  }

  # List all the chromosome-separated filenames
  coor.filenames <-
    list.files(path, ifelse(!is.null(chr.name),
                            paste0("^", chr.name, "(\\.|$)"),
                            '^chr'),
               recursive = TRUE)

  # Looping by chromosome
  chr.coor <- lapply(coor.filenames, function(coor.filename) {

    # data.table::fread is blazing fast but calling it repeatedly for different
    # chromosomes causes segfault. Revert to base::read.csv for now...
    coor <- fread(paste0(path, "/", coor.filename), showProgress = FALSE)
    #coor <- read.csv(paste0(path, "/", coor.filename))
    #setDT(coor)

    # Resolve start coordinate
    if (any(idx <- grepl("start|Start|START", names(coor)))) {
      setnames(coor, names(coor)[idx], "start")
    } else {
      setnames(coor, # expect first column is start
               coor[1:10][, names(.SD)[sapply(.SD, is.integer)]][1],
               "start")
    }

    # Resolve chromosome name
    chr.name <- gsub("(\\..*$)|(^.*/)", "", coor.filename)

    # Resolve replicate name if any (sub-folder)
    rep.name <- gsub("/.*$", "", coor.filename)
    if (rep.name == coor.filename) rep.name <- NULL

    # Resolve end coordinate
    if (is.null(case.length)) { # expect column end
      if (any(idx <- grepl("end|End|END", names(coor)))) {
        setnames(coor, names(coor)[idx], "end")
      } else { # expect second column is end
        setnames(coor,
                 coor[1:10][, names(.SD)[sapply(.SD, is.integer)]][2],
                 "end")
      }

      # Resolve case.length if NULL
      #case.length <- coor[, unique(end - start + 1)]

    } else {
      if ("end" %in% names(coor)) coor[, end := NULL]
    }

    # Resolve strand polarity
    if (strand.sensitive) {
      if (any(idx <- grepl("strand|Strand|STRAND", names(coor)))) {
        setnames(coor, names(coor)[idx], "strand")
      } else {
        setnames(coor, coor[1:10][, names(.SD)[
          sapply(.SD, function(col.elm) {
            any(col.elm %in% c("+", "-", "*"))})
        ]], "strand")
      }

    }

    # Resolve necessary columns
    col.names <- c("start",
                   if (is.null(case.length)) "end",
                   if (strand.sensitive) "strand")

    for (col.num in which(!names(coor) %in% col.names)) {
      coor[, names(coor)[col.num] := NULL ]
    }
    setcolorder(coor, col.names)

    # Remove duplicates
    coor <- unique(coor)
    setkeyv(coor, col.names)

    # Initial genomic.coordinate object
    coor <- list(coor = coor, chr.name = chr.name, rep.name = rep.name)

    return(coor)
  })

  # Build genomic.coordinate object
  coor <- lapply(chr.coor, `[[`, "coor")
  names(coor) <- sapply(chr.coor, `[[`, "chr.name", USE.NAMES = FALSE)

  chr.names <- unique(names(coor))

  # Merge or combine replicates if any
  if (!is.null(chr.coor[[1]]$rep.name)) {
    coor <- lapply(chr.names, function(chr.name) {
      if (merge.replicate) {
        Reduce(function(x, y) merge.data.table(x, y, all = TRUE),
               coor[names(coor) == chr.name])
      } else {
        rbindlist(coor[names(coor) == chr.name])
      }
    })
    names(coor) <- chr.names
  }

  # Assigning useful "attributes"
  coor$rep.name <- unique(unlist(lapply(chr.coor, `[[`, "rep.name")))
  coor$chr.names <- chr.names
  coor$case.pattern <- case.pattern
  coor$status <- data.table(is.kmer = is.kmer)
  coor$k = k
  coor$genome.name <- genome.name
  coor$is.strand.sensitive <- strand.sensitive
  coor$case.length <- case.length

  # Delete the temporary folder
  if (bedfile) unlink(temp.root, recursive = TRUE)

  class(coor) <- "genomic.coordinate"

  return(coor)
}
