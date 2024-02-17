getControlKmers <- function(case.coor, ctrl.rel.pos, genome,
                            kmer.table=NULL, add.rev.kmers=TRUE, output.path) {
  
  t1 <- Sys.time()
  
  # Print some message
  msg <- paste0("Resolving control regions",
                if(length(case.coor$chr.names) == 1)
                  paste0(" of ", case.coor$chr.names))
  l <- paste(rep(".", 60 - nchar(msg)), collapse = "")
  cat(msg, l, sep = "")

  # Resolve control regions
  control.coor <- buildControl(case.coor, ctrl.rel.pos, genome)
  
  # Save control coordinates
  saveCoor(control.coor, paste0(output.path, "/control_coordinates/"))

  t <- Sys.time() - t1
  cat("DONE! --", round(t[[1]], 2), attr(t, "units"), "\n")
  
  
  t1 <- Sys.time()
  
  # Print some message
  msg <- paste0("Extracting ", case.coor$k, "-mer from the control regions",
                if(length(case.coor$chr.names) == 1)
                  paste0(" of ", case.coor$chr.names))
  l <- paste(rep(".", 60 - nchar(msg)), collapse = "")
  cat(msg, l, sep = "")
  
  control.kmers <- extractKmers(control.coor, genome, rm.overlaps = FALSE,
                                kmer.table, add.rev.kmers)
  
  t <- Sys.time() - t1
  cat("DONE! --", round(t[[1]], 2), attr(t, "units"), "\n")
  
  return(control.kmers)
}