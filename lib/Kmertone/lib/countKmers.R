countKmers <- function(chr.name, start, end=NULL, len=NULL, strand=NULL,
                       start2=NULL, end2=NULL, len2=NULL, genome) {
  
  # Function helper for extractKmers.R
  # Either end or len must be given
  # Output is <data.table> with two columns: kmer and N
  
  mapKmers <- function(start, end, len, strand) {
    if (is.null(end)) {
      kmer.count <- stri_sub(genome$seq[[chr.name]], start, length = len)
    } else {
      kmer.count <- stri_sub(genome$seq[[chr.name]], start, end)
    }
    if (!is.null(strand) && strand == "-") {
      kmer.count <- reverseComplement(kmer.count)
    }
    
    return(kmer.count)
  }
  
  kmer.count <- mapKmers(start, end, len, strand)
  if (!is.null(start2)) {
    kmer.count2 <- mapKmers(start2, end2, len2, strand)
    kmer.count <- data.table(kmer_part1 = kmer.count,
                             kmer_part2 = kmer.count2,
                             key = c("kmer_part1" , "kmer_part2"))
  } else {
    kmer.count <- data.table(kmer = kmer.count, key = "kmer")
  }
  
  kmer.column <- names(kmer.count)[grepl("kmer", names(kmer.count))]
  
  kmer.count <- kmer.count[, .N, by = kmer.column]
  setkeyv(kmer.count, kmer.column)
  
  return(kmer.count)
}