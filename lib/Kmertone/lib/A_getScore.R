getScore <- function(case.kmers, control.kmers, case.pattern=NULL, k,
                     output.path) {
  
  # Calculate k-mer skew
  #if (!strand.sensitive) {
    calKmerSkew(case.kmers)
    calKmerSkew(control.kmers)
  #}
  
  # Sum counts from both strands
  case.kmers[, case := pos_strand + neg_strand]
  control.kmers[, control := pos_strand + neg_strand]
  
  # Remove separate strand counts
  case.kmers[, c("pos_strand", "neg_strand") := NULL]
  control.kmers[, c("pos_strand", "neg_strand") := NULL]
  
  # Rename kmer_skew
  setnames(case.kmers, "kmer_skew", "case_skew")
  setnames(control.kmers, "kmer_skew", "control_skew")
  
  setkeyv(case.kmers, names(case.kmers)[grepl("kmer", names(case.kmers))])
  setkeyv(case.kmers, names(control.kmers)[grepl("kmer", names(case.kmers))])
  
  if (is.null(case.pattern)) {
    kmer.table <- merge(case.kmers, control.kmers, all = TRUE)
  } else {
    kmer.table <- merge(case.kmers, control.kmers, all.x = TRUE)
  }
  
  scoreKmers(kmer.table)
  
  setcolorder(kmer.table,
              c(names(kmer.table)[grepl("kmer", names(kmer.table))],
                "case", "control", "case_skew", "control_skew", "z"))
  
  fwrite(kmer.table, paste0(output.path, "/score_", k, "-mers.csv"),
         showProgress = FALSE)
  
  cat(paste(c(rep("-", 80), "\n"), collapse = ""))
  message("The ", k, "-mer scores are saved at ", output.path, "/score_",
          k,"-mer.csv")
  
  return(kmer.table)
}