prepCoordinate <- function(case.coor, case.coor.path, case.length, k,
                           genome.name, strand.sensitive, merge.replicate,
                           case.pattern, chr.name, rm.diff.pattern) {
  
  t1 <- Sys.time()
  
  # Print some message
  msg <- paste0("Building genomic coordinate",
                if(!is.null(chr.name)) paste(" of", chr.name))
  l <- paste(rep(".", 60 - nchar(msg)), collapse = "")
  cat(msg, l, sep = "")
  
  if (!is.null(case.coor)) {
    
    if (class(case.coor) != "genomic.coordinate") {
      stop("Please input genomic.coordinate class object only.")
    }
  
  # Build genomic.coordinate-class object.
  } else if (is.null(case.coor)) {
    case.coor <- buildCoordinate(path = case.coor.path, case.length, k,
                                 genome.name, strand.sensitive, merge.replicate,
                                 case.pattern, chr.name)
  }
  
  # Remove different case pattern?
  if (rm.diff.pattern) {
    cat("\n")
    case.coor <- checkCasePattern(case.coor, genome, rm.diff.pattern = TRUE)
    cat(msg, l, sep = "")
  }
  
  t <- Sys.time() - t1
  cat("DONE! --", round(t[[1]], 2), attr(t, "units"), "\n")
  
  return(case.coor)
}