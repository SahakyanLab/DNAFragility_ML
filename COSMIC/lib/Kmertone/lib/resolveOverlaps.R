resolveOverlaps <- function(coor, action="mark") {
  
  # action    <string>     mark - mark the overlaps
  #                      merge - merge the overlaps

  # Assign idx.overlap to coor
  
  for (chr.name in coor$chr.names) {
    
    # Locate continuous overlapping regions, group them together ,and
    # assign group number
    setkeyv(coor[[chr.name]], c(if(coor$is.strand.sensitive) "strand", "start",
                                if(length(coor$case.length) != 1) "end"))
    
    if (coor$status$is.kmer & length(coor$case.length) == 1) {
      coor[[chr.name]][, group := cumsum(c(1, cummax(head(start + coor$k - 1,
                                                          -1)) -
                                             tail(start, -1) < 0)),
                       by = eval(if(coor$is.strand.sensitive) "strand")]
    } else {
      coor[[chr.name]][, group := cumsum(c(1, cummax(head(end, -1)) -
                                             tail(start, -1) < 0)),
                       by = eval(if(coor$is.strand.sensitive) "strand")]
    }
    
    if (action == "mark") {
      coor[[chr.name]][, is_overlap := ifelse(.N > 1, TRUE, FALSE),
                       by = c(if(coor$is.strand.sensitive) "strand", "group")]
      
      coor[[chr.name]][, group := NULL]
      
    } else if (action == "merge" & any(duplicated(coor[[chr.name]]$group))) {
      coor[[chr.name]] <-
        coor[[chr.name]][, .(start = min(start), end = max(end)),
                         by = c(if(coor$is.strand.sensitive) "strand", "group")]
      
      # remove column group
      coor[[chr.name]][, group := NULL]
      
      # rearrange the columns
      setcolorder(coor[[chr.name]], c("start", "end",
                                      if(coor$is.strand.sensitive) "strand"))
      
    } else {
      # remove column group
      coor[[chr.name]][, group := NULL]
      
      # rearrange the columns
      setcolorder(coor[[chr.name]], c("start", "end",
                                      if(coor$is.strand.sensitive) "strand"))
    }
  }
  
  return(coor)
}