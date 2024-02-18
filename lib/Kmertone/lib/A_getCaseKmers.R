getCaseKmers <- function(case.coor, genome, rm.case.kmer.overlaps,
                         kmer.table=NULL, add.rev.kmers=TRUE) {
  
  t1 <- Sys.time()
  
  # Print some message
  msg <- paste0("Extracting ", case.coor$k, "-mer from the case regions",
                if(length(case.coor$chr.names) == 1)
                  paste0(" of ", case.coor$chr.names))
  l <- paste(rep(".", 60 - nchar(msg)), collapse = "")
  cat(msg, l, sep = "")
  
  # Expand case to k-mer
  if (!is.null(case.coor$case.length)) {
    if (case.coor$k > case.coor$case.length) {
      kmerize(case.coor)
    } else if (case.coor$k == case.coor$case.length) {
      case.coor$status[, is.kmer := TRUE]
    }
  }
  
  if (rm.case.kmer.overlaps) {
    cat("\n", sep = "")
  }
  
  # Extract case k-mers
  case.kmers <- extractKmers(case.coor, genome, rm.case.kmer.overlaps,
                             kmer.table, add.rev.kmers)
  
  if (rm.case.kmer.overlaps) {
    msg <- paste0("...", msg)
    l <- paste(rep(".", 60 - nchar(msg)), collapse = "")
    cat(msg, l, sep = "")
  }
  
  t <- Sys.time() - t1
  cat("DONE! --", round(t[[1]], 2), attr(t, "units"), "\n")
  
  return(case.kmers)
}