checkCasePattern <- function(case.coor, genome, rm.diff.pattern=FALSE) {
  
  # Utility function for genomic.coordinate class object
  # Check for case.pattern occurence. If different pattern found, report
  # percentage of occurrence and/or remove them from the table. By default
  # kmertone module "score" will exclude them in case k-mers but selection
  # of control regions will include them. If user want them to be removed,
  # use rm.diff.pattern = TRUE
  
  
  initial.state.is.kmer <- case.coor$status$is.kmer
  
  if (case.coor$status$is.kmer) {
    kmerize(case.coor, revert = TRUE)
  }
  
  if (! "seq" %in% names(case.coor)) {
    getCoorSeq(case.coor, genome)
  }
  
  total.pos <- rep(0, length(case.coor$chr.names))
  diff.num <- rep(0, length(case.coor$chr.names))
  for (i in seq_along(case.coor$chr.names)) {
    
    total.pos[i] <- case.coor[[i]][, .N]
    idx.diff <- !case.coor[[i]]$seq %in% case.coor$case.pattern
    diff.num[1] <- sum(idx.diff)
    
    if (rm.diff.pattern) {
      case.coor[[i]] <- case.coor[[i]][!idx.diff]
    }

    case.coor[[i]][, seq := NULL]
  }
  
  names(total.pos) <- case.coor$chr.names
  names(diff.num) <- case.coor$chr.names
  
  # ----------------------------------------------------------------------------
  # Print stats
  
  cat("Total cases", paste0(c(rep(" ", 17), ": "), collapse = ""),
      sum(total.pos),
      
      "\nDifferent DNA pattern",
      paste0(c(rep(" ", 7), ": "), collapse = ""),
      diff.num, " (",
      signif(sum(diff.num) / sum(total.pos) * 100, 2), "%)\n", sep = "")
  
  cat("The difference by chromosome:\n")
  print(signif(diff.num/total.pos * 100, 2))
  
  if(rm.diff.pattern) message("--- Different patterns are removed! ---")
  
  if (initial.state.is.kmer) {
    kmerize(case.coor)
  }
  
  if (rm.diff.pattern) {
    return(case.coor)
  }
}