prepGenome <- function(genome, genome.path, genome.name, chr.name) {
  
  t1 <- Sys.time()
  
  # Print some message
  msg <- paste0("Building ", if(!is.null(chr.name)) paste(chr.name, "of "),
                "genome ", genome.name)
  l <- paste(rep(".", 60 - nchar(msg)), collapse = "")
  cat(msg, l, sep = "")
  
  if (!is.null(genome)) {
    
    if (class(genome) != "genome") {
      stop("Please input genome-class object only.")
    }
    
    # Build user's genome
  } else if (!is.null(genome.path)) {
    
    genome <- buildGenome(genome.path, genome.name, chr.name)
    
    # Build genome-class object from database
  } else if (genome.name %in% c("hg19", "hg38")) {
    
    genome.path <- paste0("data/genome/human/", genome.name, "/")
    
    # Check if genome available locally to download
    if (!dir.exists(genome.path) |
        length(list.files(genome.path, "^chr.+\\.fa")) == 0) {
      
      dir.create(genome.path, showWarnings = FALSE, recursive = TRUE)
      
      genome.url <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/",
                           genome.name, "/chromosomes/")
      chr.filenames <- unlist(strsplit(getURL(genome.url), '>|<|\n'))
      rgx <- grepl("^chr([0-9XYxy]+|[Mm].{0,3})\\.fa", chr.filenames)
      chr.filenames <- chr.filenames[rgx]
      
      dir.create(genome.path, recursive = TRUE, showWarnings = FALSE)
      for (chr.filename in chr.filenames) {
        cat("\nFetching", chr.filename, "\n")
        download.file(paste0(genome.url, chr.filename),
                      paste0(genome.path, chr.filename),
                      method = "curl")
      }
      genome <- buildGenome(genome.path, genome.name)
      saveRDS(genome, paste0(genome.path, "/genome.RDS"))
    }
    
    # Reprint message
    msg <- paste("Building genome", genome.name)
    l <- paste(rep(".", 60 - nchar(msg)), collapse = "")
    cat(msg, l, sep = "")
    
    if (!is.null(chr.name)) {
      genome <- buildGenome(genome.path, genome.name, chr.name)
    } else {
      genome <- readRDS(paste0(genome.path, "/genome.RDS"))
    }
  }
  
  t <- Sys.time() - t1
  cat("DONE! --", round(t[[1]], 2), attr(t, "units"), "\n")
  
  return(genome)
}