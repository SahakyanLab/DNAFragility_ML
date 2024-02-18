buildGenome <- function(path, genome.name, chr.name=NULL){
  
  # Given a folder with chromosome-separated fasta sequence
  # The fasta files should be named exactly like chromosome names in
  # genomic coordinate table.
  #
  # path         <string>     A path to a folder containing chromosome-separated
  #                          fasta files. Compressed file is supported. The file
  #                          name must be similar to coordinate table with no
  #                          prefix and suffix.
  # genome.name  <string>    User-defined genome name. It is a good idea to name
  #                          according to its version e.g. hg38
  # chr.name     <string>    A single chromosome name. Useful for looping.
  #
  # Dependencies: data.table, stringi
  
  suppressPackageStartupMessages( library(data.table)  )
  suppressPackageStartupMessages( library(  stringi )  )
  
  chr.filenames <-
    list.files(path, ifelse(is.null(chr.name), ".+\\.(fa|fna|fasta)",
                            paste0(chr.name, "\\.(fa|fna|fasta)")))
  
  genome.seq <- lapply(chr.filenames, function(chr.filename) {
    
    genome.seq <- fread(paste0(path, "/", chr.filename), showProgress = FALSE)
    setnames(genome.seq, names(genome.seq), "sequence")
    genome.seq[, sequence := stri_trans_toupper(sequence)]
    genome.seq <- paste(genome.seq[[1]], collapse = "")
  })
  
  names(genome.seq) <- gsub("\\.(fa|fna|fasta).*$", "", chr.filenames)
  
  genome <- structure(list(seq = genome.seq, chr.names = names(genome.seq),
                           len = nchar(genome.seq), name = genome.name),
                      class = "genome")
  
  # For printing
  print.genome <- function(obj) {
    cat("Genome:", genome.name, "\n")
    print(obj$len)
  }
  assign("print.genome", print.genome, envir = .GlobalEnv)
  
  return(genome)
}