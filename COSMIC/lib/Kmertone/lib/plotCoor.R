plotCoor <- function(case.coor, ctrl.coor, chr.name, pos.range,
                     rm.case.overlaps=FALSE) {
  
  # Input coordinates are <genomic.coordinate> object. Use buildCoordinate.R to
  # build <genomic.coordinate> object.
  # pos.range    <character>  c(start, end)
  # rm.overlaps    <bool>     Remove case k-mer overlaps only if the case.coor
  #                           the contains overlapping info,
  # this plot is not strand sensitive ,so + and - strand with the same position
  # will be shown as overlap
  
  # Error checking
  if (!chr.name %in% case.coor$chr.names) {
    stop("No selected chromosome in the coordinates.")
  }
  
  # Select only query region
  ctrl.reg <-
    ctrl.coor[[chr.name]][start >= pos.range[1] & end <= pos.range[2]]
  ctrl.reg[, tag := "control"]
  
  case.reg <-
    case.coor[[chr.name]][
      start >= pos.range[1] & 
        if(!is.null(case.coor$case.length))
          (start + case.coor$case.length - 1) <= pos.range[2] else
            end <= pos.range[2] &
        if(rm.case.overlaps) !is_overlap else TRUE]
  case.reg[, tag := "case"]
  if (!is.null(case.coor$case.length)) {
    case.reg[, end := start + case.coor$case.length - 1]
  }
  
  if (!"strand" %in% names(case.reg)) {
    case.reg[, strand := "*"]
    ctrl.reg[, strand := "*"]
  } else {
    ctrl.reg <- rbind(ctrl.reg[, strand := "+"],
                      copy(ctrl.reg)[, strand := "-"])
  }
  
  # Only take necessary columns
  unwanted.columns <- names(case.reg)[!names(case.reg) %in%
                                        c("start", "end", "strand", "tag")]
  for (column in unwanted.columns) {
    set(case.reg, j = column, value = NULL)
  }
  
  unwanted.columns <- names(ctrl.reg)[!names(ctrl.reg) %in%
                                        c("start", "end", "strand", "tag")]
  for (column in unwanted.columns) {
    set(ctrl.reg, j = column, value = NULL)
  }
  
  # combine data
  dt <- rbind(ctrl.reg, case.reg)
  
  if (nrow(dt) == 0) {
    msg <- paste0("There is no case and control regions that fit within",
                  " the query region")
    stop(msg)
  }
  
  # Locate any overlaps among case and control regions
  setkey(dt, strand, start, end)
  dt <- dt[, group := cumsum(c(1, cummax(head(end, -1)) - tail(start, -1) < 0)),
           by = strand]
  
  # Assign bin number for overlaps 
  assignBin <- function(start, end) {
    # a wrapper function for assigning bin number for merge-able regions
    
    N <- length(start)
    
    if (N == 1) {
      bin <- 1
    } else {
      
      bin <- c(1, rep(NA, N-1))
      
      for (i in 2:N) {
        
        end.before <- end[1:(i-1)]
        bin.before <- bin[1:(i-1)]
        
        if (sum(start[i] <= end.before) >= max(bin.before)) {
          
          bin[i] <- max(bin.before) + 1
          
        } else if (sum(start[i] <= end.before) == 0) {
          
          bin[i] <- 1
          
        } else {
          
          # get bin number for overlapping region
          bin.overlap <- bin.before[start[i] <= end.before]
          
          # assign bin number that does not belong to the overlaps
          bin[i] <- bin.before[!bin.before %in% bin.overlap][1]
          
        }
      }
    }
    return(bin)
  }
  
  # Assign bin number for every continuous regions i.e. non-overlapping
  dt[, bin := assignBin(start, end), by = c("strand", "group")]
  if (any(unique(dt$strand) %in% c("+", "-"))) {
    max.bin.plus <- dt[strand == "+", if(.N > 0) max(bin) else 0]
    dt[strand == "-", bin := bin + max.bin.plus]
  }

  height <- 1
  xlim <- c(min(dt$start), max(dt$end))
  bins <- dt$bin # to detect overlap and draw overlap above
  plot.new()
  plot.window(xlim = xlim, ylim =  c(0, max(bins) * (height + 0.5)))
  ybottom <- bins * (0.5 + height) - height
  ybottom1 <- ybottom[dt$tag == "case"]
  ybottom2 <- ybottom[dt$tag == "control"]
  
  bor = 0.001
  if (nrow(dt[tag == "case"]) == 0) {
    message("There is no case region in this area.")
  } else {
    rect(xleft = dt[tag == "case"]$start - bor,
         ybottom = ybottom1,
         xright = dt[tag == "case"]$end + bor,
         ytop = ybottom1 + height,
         col = "red", border=NA)
  }
  
  if (nrow(dt[tag == "control"]) == 0) {
    message("There is no control region in this area.")
  } else {
    rect(xleft = dt[tag == "control"]$start - bor,
         ybottom = ybottom2,
         xright = dt[tag == "control"]$end + bor,
         ytop = ybottom2 + height,
         col = "blue", border = NA)
  }
  
  par(xpd = TRUE)
  legend("topright", inset = c(0.05,-0.25), legend=c("Case", "Control"),
         fill=c("red", "blue"), cex=0.8, bty = "n")
  
  title(xlab="Genomic coordinate")
  
  # Draw a line to separate between plus and minus strand
  if (sum(unique(dt$strand) %in% c("+", "-")) == 2) {
    strand.line <- ( max(ybottom[dt$strand == "+"] + height) +
                       min(ybottom[dt$strand == "-"]) ) / 2
    abline(h = strand.line, lty = 2, xpd = FALSE)
    text(x = par("usr")[1], y = strand.line, pos = 1, labels = "+",
         family = "mono", font = 2)
    text(x = par("usr")[1], y = strand.line, pos = 3, labels = "-",
         family = "mono", font = 2)
  } else {
    text(x = par("usr")[1], y = par("usr")[4], pos = 3,
         labels = paste("strand", if(unique(dt$strand) == "*") "\U2731" else
           unique(dt$strand)))
  }
  
  
  axis(1)
  
  return(dt)
}

