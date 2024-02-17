getCoorSeq <- function(coor, genome) {
  
  # Add a column with sequence mapped to coordinate
  # Dependencies: reverseComplement.R
  
  for (chr.name in coor$chr.names) {
    
    coor[[chr.name]][, seq := {
      
      if (coor$status$is.kmer) {
        seq <- stri_sub(genome$seq[chr.name], start, length = coor$k)
      } else if (!is.null(coor$case.length)) {
        seq <- stri_sub(genome$seq[chr.name], start, length = coor$case.length)
      } else {
        seq <- stri_sub(genome$seq[chr.name], start, end)
      }
      
      if (coor$is.strand.sensitive && strand == "-") {
        seq <- reverseComplement(seq)
      }
      
      seq
    }, by = c(if(coor$is.strand.sensitive) "strand")]
  }
  
  invisible(coor)
}