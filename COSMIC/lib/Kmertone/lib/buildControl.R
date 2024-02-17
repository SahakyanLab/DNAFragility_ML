buildControl <- function(coor, ctrl.rel.pos=c(50,500), genome) {

  # control region distance is measured from case position, so...
  initial.coordinate.status.is.kmer <- coor$status$is.kmer
  if (coor$status$is.kmer) {
    kmerize(coor, revert = TRUE)
  }
  
  control.coor <- lapply(coor$chr.names, function(chr.name) {

    # Initial control region location
    control.coor <- coor[[chr.name]][, .(
      start = c(start - ctrl.rel.pos[2],
                (if(!is.null(coor$case.length)) (start + coor$k - 1) else
                  end) + ctrl.rel.pos[1]),
      end = c(start - ctrl.rel.pos[1],
              (if(!is.null(coor$case.length)) (start + coor$k - 1) else
                  end) + ctrl.rel.pos[2])
    )]

    # Build genomic.coordinate object of initial control regions
    control.coor <- list(coor = control.coor,
                         status = data.table(is.kmer = FALSE),
                         is.strand.sensitive = FALSE)
    names(control.coor)[1] <- chr.name

    # Merge overlapping regions
    control.coor <- resolveOverlaps(control.coor, action = "merge")

    # Trim out-of-range regions
    control.coor[[chr.name]][start < 1, start := 1]
    control.coor[[chr.name]][end > genome$len[[chr.name]],
                             end := genome$len[[chr.name]]]

    setkeyv(control.coor[[chr.name]], c("start", "end"))

    # Build case zone i.e. case region + buffer
    case.zone <- list(
      coor = coor[[chr.name]][, .(
        start = start - ctrl.rel.pos[1],
        end = (if(!is.null(coor$case.length)) (start + coor$k - 1) else
                  end) + ctrl.rel.pos[1])],

      chr.names = "coor", is.strand.sensitive = FALSE,
      status = data.table(is.kmer = FALSE))
    case.zone <- resolveOverlaps(case.zone, action = "merge")
    setkey(case.zone[["coor"]], start, end)

    # Remove control region portions that overlap with case region
    control.coor <- removeCaseZone(control.coor[[chr.name]], case.zone$coor)

    return(control.coor$coor)
  })

  names(control.coor) <- coor$chr.names
  control.coor$status <- data.table(is.kmer = FALSE)

  for (attrb in c("chr.names", "genome.name", "k")) {
    control.coor[[attrb]] <- coor[[attrb]]
  }

  control.coor$is.strand.sensitive <- FALSE

  class(control.coor) <- "genomic.coordinate"

  if(initial.coordinate.status.is.kmer) {
    kmerize(coor)
  }

  return(control.coor)
}
