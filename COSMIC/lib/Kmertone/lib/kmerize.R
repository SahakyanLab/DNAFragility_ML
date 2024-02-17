kmerize <- function(coor, revert=FALSE) {
  
  # Resolve kmer coordinates
  # - change status is.kmer from FALSE to TRUE
  #   ("revert" change the other way around)
  # - only applicable to single case.length
  #
  # coor    <genomic.coordinate>
  # revert       <bool>
  
  # Error checking
  if (class(coor) != "genomic.coordinate") {
    stop("Incompatible genomic coordinate")
  } else if (!coor$status$is.kmer & revert) {
    stop("Cannot revert. Status is.kmer is already FALSE")
  } else if (coor$status$is.kmer & !revert) {
    stop("The coordinate is already kmer coordinates")
  }
  
  if (length(coor$case.length) == 1) {
    flank <- (coor$k - coor$case.length) / 2
  }
  
  # Expand to kmer coordinate
  if (!coor$status$is.kmer & length(coor$case.length) == 1 &&
      coor$case.length < coor$k) {
    
    for (chr.name in coor$chr.names) {
      
      coor[[chr.name]][, start := start - flank]
    }
    coor$status[, is.kmer := TRUE]
    
    # Contract to case coordinate
  } else if (coor$status$is.kmer & revert) {
    for (chr.name in coor$chr.names) {
      coor[[chr.name]][, start := start + flank]
    }
    coor$status[, is.kmer := FALSE]
  } 
  # else if (coor$case.length >= coor$k) {
  #   message("The case length is either the same as k or longer.")
  # }
  
  invisible(coor)
}