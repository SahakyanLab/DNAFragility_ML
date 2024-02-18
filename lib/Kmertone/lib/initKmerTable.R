initKmerTable <- function (k, case.pattern=NULL) {
  
  # Initialise k-mer table with the following columns: kmer, pos_strand, and
  # neg_strand
  # For 9- to 15-mer, the kmer sequence is separated to two columns to
  # reduce memory usage significantly.
  # k is limited to 15 because vector size is limited to .Machine$integer.max
  # R is still using 32-bit architecture. Even in 64-bit, an enormous amount
  # of memory is needed.
  #
  # Update: skewness count is considered, so count has to be separated to
  #         positive and negative strand count.
  #         column count is separated to pos_strand and neg_strand
  
  if (k > 15) {
    stop("k is limited to 15 because vector size is limited to ",
         ".Machine$integer.max")
  }
  
  if (k < 9) {
    if (is.null(case.pattern)) {
      possible.kmers <- do.call(CJ, rep(list(c("A", "C", "G", "T")), k))
    } else {
      # calculate expansion factor and pattern position
      expansion.factor <- (k - nchar(case.pattern[1])) / 2
      
      expansion.kmers <- do.call(CJ, rep(list(c("A", "C", "G", "T")),
                                         expansion.factor))
      expansion.kmers <- expansion.kmers[, do.call(paste0, .SD)]
      expansion.kmers.combn <- do.call(CJ, rep(list(expansion.kmers), 2))
      
      possible.kmers <-
        lapply(case.pattern, function(dna.pattern) {
          
          possible.kmers <- expansion.kmers.combn[, .(V1, dna.pattern, V2)]
          
          return(possible.kmers)
        })
      
      possible.kmers <- rbindlist(possible.kmers)
    }
    
    possible.kmers <- possible.kmers[, do.call(paste0,.SD)]
    kmer.table <- data.table(kmer = possible.kmers, pos_strand = 0,
                             neg_strand = 0, key = "kmer")
    
    # Separate k-mer sequence to two parts (sweet spot)
  } else {
    
    k1 <- floor(k / 2)
    k2 <- k - k1
    
    part1 <- initKmerTable(k1)$kmer
    if (k2 == k1) {
      part2 <- part1
    } else {
      part2 <- initKmerTable(k2)$kmer
    }
    
    kmer.table <- do.call(CJ, list(kmer_part1 = part1, kmer_part2 = part2))
    kmer.table[, c("pos_strand", "neg_strand") := 0]
    setkey(kmer.table, kmer_part1, kmer_part2)
  }
  
  
  return(kmer.table)
}