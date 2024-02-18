bedToCoor <- function(path, output.path="coordinate") {
  
  # Convert bedfile to normal 1-index csv file, separated by chromosome in a
  # a folder.
  
  suppressPackageStartupMessages( library(data.table)  )
  
  dir.create(output.path, showWarnings = FALSE, recursive = TRUE)
  
  # By convention, these are BedFile columns in a strict order.
  bed.cols <- c("chrom", "chromStart", "chromEnd", "name", "score",
                "strand", "thickStart", "thickEnd", "itemRgb", "blockCount",
                "blockSizes", "blockStarts")
    bed <- fread(path, showProgress = FALSE)
    setnames(bed, names(bed), bed.cols[1:ncol(bed)])
    bed[, chromStart := chromStart + 1]
    
    # Just to set it apart from bedfile
    setnames(bed, c("chromStart", "chromEnd"), c("start", "end"))
    
    # Check for strand polarity
    has.strand <- "strand" %in% names(bed)
    if (has.strand) {
      strand.type <- bed[, unique(strand)]
      strand.sensitive <- ifelse(length(strand.type) == 1, FALSE, TRUE)
    } else {
      strand.sensitive <- FALSE
    }
    
    # Calculate unique length(s)
    len <- bed[, unique(end - start + 1)]
    
    dir.create(output.path, recursive = TRUE, showWarnings = FALSE)
    
    bed[, fwrite(.SD,  paste0(output.path, "/", chrom, ".csv")),
        
        .SDcols = c("start",
                    if (length(len) > 1) "end",
                    if (strand.sensitive) "strand"),
        by = chrom]
    
    # Save metadata if any column drop
    if (has.strand & !strand.sensitive)
      writeLines(paste0("Unique strand: ", strand.type),
                 paste0(output.path, "/note.txt"))
    if (length(len) == 1)
      writeLines(paste0("Unique length: ", len),
                 paste0(output.path, "/note.txt"))
    
  
}