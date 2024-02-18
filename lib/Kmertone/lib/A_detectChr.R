detectChr <- function(case.coor, case.coor.path, genome.path) {
  
  if (!is.null(case.coor)) {
    return(case.coor$chr.names)
  }
  
  filenames <- list.files(case.coor.path, "^chr", recursive = TRUE)
  chr.names <- unique(gsub("^.*/|\\..*$", "", filenames))
  
  # If bedfile, rely on genome filename
  if(length(filenames) == 0 |
     sum(grepl("^chr", chr.names)) != length(chr.names)) {
    filenames <- list.files(genome.path, "\\.((fasta)|(FASTA)|(fa)|(FA)).*$")
    chr.names <- gsub("\\..*$", "", filenames)
  }
  
  return(chr.names)
}