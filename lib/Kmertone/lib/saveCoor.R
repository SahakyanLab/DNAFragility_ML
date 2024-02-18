saveCoor <- function(coor, output.path) {
  
  dir.create(output.path, showWarnings = FALSE, recursive = TRUE)
  
  for (chr.name in coor$chr.names) {
    
    fwrite(coor[[chr.name]], paste0(output.path, "/", chr.name, ".csv"))
    
    
  }
  
  note <- paste0("Case length: ", coor$case.length, "\nGenome: ",
                 coor$genome.name)
  
  writeLines(note, paste0(output.path, "/note.txt"))
  
  invisible(NULL)
}