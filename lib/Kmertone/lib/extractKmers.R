extractKmers <- function(coor, genome, rm.overlaps=TRUE, kmer.table=NULL,
                         add.rev.kmers=TRUE) {

  # Out-of-range coordinate i.e. defect kmers are removed.
  #
  # kmer.table     <data.table>    If given, update the table.
  #
  # Dependencies: countKmers.R

  if (genome$name != coor$genome.name) {
    stop("Genome name mismatch between coordinate and genome.")
  } else if (class(coor) != "genomic.coordinate") {
    stop("Incompatible genomic coordinate")
  } else if (class(genome) != "genome") {
    stop("Incompatible genome")
  }

  # Create table with all possible k-mers
  if (is.null(kmer.table) & coor$k < 16) {
    kmer.table <- initKmerTable(k = coor$k, case.pattern = coor$case.pattern)
  } else if (is.null(kmer.table) & coor$k > 15) {
    # Sample table
    kmer.table <- data.table(kmer = character(), pos_strand = numeric(),
                             neg_strand = numeric(), key = "kmer")
  }

  # environment for k > 15. This is needed because data.table is not capable of
  # updating row by reference
  if (coor$k > 15) env <- environment()

  # If remove overlaps, get idx.overlap
  if (rm.overlaps) {
    # Overlapping case k-mer regardless of strand polarity
    strand.sensitive.initial.state <- coor$is.strand.sensitive
    coor$is.strand.sensitive <- FALSE
    coor <- resolveOverlaps(coor, action = "mark")
    coor$is.strand.sensitive <- strand.sensitive.initial.state
  }

  if (coor$status$is.kmer) {
    for (chr.name in coor$chr.names) {

      coor[[chr.name]][
        # remove out-of-range coordinates and/or overlaps
        (start > 0 & start + coor$k - 1 <= genome$len[chr.name]) &
          (if(rm.overlaps) (!is_overlap) else TRUE), {

            if (coor$k < 9) {

              # count k-mer
              kmer.count <-
                countKmers(chr.name, start, len = coor$k,
                           strand = if(coor$is.strand.sensitive) strand,
                           genome = genome)

              # update kmer.table
              if (coor$is.strand.sensitive && strand == "-") {
                kmer.table[kmer.count , neg_strand := neg_strand + N]
              } else {
                kmer.table[kmer.count , pos_strand := pos_strand + N]
              }

            } else if (coor$k > 8 & coor$k < 16) {

              if (coor$is.strand.sensitive && strand == "-") {
                k2 <- floor(coor$k / 2)
                k1 <- coor$k - k2
              } else {
                k1 <- floor(coor$k / 2)
                k2 <- coor$k - k1
              }

              # count k-mers
              kmer.count <-
                countKmers(chr.name, start, len = k1,
                           strand = if(coor$is.strand.sensitive) strand,
                           start2 = start + k1, len2 = k2, genome = genome)

              # update kmer.table
              if (coor$is.strand.sensitive && strand == "-") {
                kmer.table[kmer.count , neg_strand := neg_strand + N]
              } else {
                kmer.table[kmer.count , pos_strand := pos_strand + N]
              }

            } else if (coor$k > 15) {
              # Growing kmer.table instead of updating

              # Count k-mers
              kmer.count <-
                countKmers(chr.name, start, len = coor$k,
                           strand = if(coor$is.strand.sensitive) strand,
                           genome = genome)

              # Update kmer.table for existing k-mers in the table
              if (coor$is.strand.sensitive && strand == "-") {
                kmer.table[kmer.count , neg_strand := neg_strand + N]
              } else {
                kmer.table[kmer.count , pos_strand := pos_strand + N]
              }

              # Add additional column to kmer.count
              if (coor$is.strand.sensitive && strand == "-") {
                setnames(kmer.count, "N", "neg_strand")
                kmer.count[, pos_strand := 0]
              } else {
                setnames(kmer.count, "N", "pos_strand")
                kmer.count[, neg_strand := 0]
              }

              # Bind new k-mers found
              kmer.table <- rbind(kmer.table, kmer.count[!kmer.table])
              setkey(kmer.table, kmer)

              assign("kmer.table", kmer.table, envir = env)
            }

          }, by = eval(if(coor$is.strand.sensitive) "strand")]

      # Print percent removed overlapping k-mers
      if (rm.overlaps) {
        msg <- paste0("...", signif(coor[[chr.name]][, sum(is_overlap)] /
                                      coor[[chr.name]][, .N] * 100,
                                    2),
                      "% overlapping k-mers are removed")
        l <- paste(rep(".", 60 - nchar(msg)), collapse = "")
        cat(msg, l, "\n", sep = "")
      }
    }

    if (!coor$is.strand.sensitive & add.rev.kmers) {
      kmer.table <- countRevCompKmers(kmer.table)
    }

    # Sliding kmer extraction for case or control with varied length
  } else if (!coor$status$is.kmer & is.null(coor$case.length)) {
    for (chr.name in coor$chr.names) {

      # remove region with lower than size k
      coor[[chr.name]][end - start + 1 >= coor$k, {

        # Resolve k-mer start positions
        k.start <- unlist(Map(function(start, end) {
          start:(end - coor$k + 1)
        }, start, end))

        if (coor$k < 9) {

          kmer.count <-
            countKmers(chr.name, k.start, len = coor$k,
                       strand = if(coor$is.strand.sensitive) strand,
                       genome = genome)

          # Update kmer.table
          if (coor$is.strand.sensitive && strand == "-") {
            kmer.table[kmer.count , neg_strand := neg_strand + N]
          } else {
            kmer.table[kmer.count , pos_strand := pos_strand + N]
          }

        } else if (coor$k > 8 & coor$k < 16) {

          if (coor$is.strand.sensitive && strand == "-") {
            k2 <- floor(coor$k / 2)
            k1 <- coor$k - k2
          } else {
            k1 <- floor(coor$k / 2)
            k2 <- coor$k - k1
          }

          kmer.count <-
            countKmers(chr.name, k.start, len = k1,
                       strand = if(coor$is.strand.sensitive) strand,
                       start2 = k.start + k1, len2 = k2,
                       genome = genome)

          # Update kmer.table
          if (coor$is.strand.sensitive && strand == "-") {
            kmer.table[kmer.count , neg_strand := neg_strand + N]
          } else {
            kmer.table[kmer.count , pos_strand := pos_strand + N]
          }

        } else if (coor$k > 15) {

          kmer.count <-
            countKmers(chr.name, k.start, len = coor$k,
                       strand = if(coor$is.strand.sensitive) strand,
                       genome = genome)

          # Update kmer.table for existing k-mers in the table
          if (coor$is.strand.sensitive && strand == "-") {
            kmer.table[kmer.count , neg_strand := neg_strand + N]
          } else {
            kmer.table[kmer.count , pos_strand := pos_strand + N]
          }

          # Add additional column to kmer.count
          if (coor$is.strand.sensitive && strand == "-") {
            setnames(kmer.count, "N", "neg_strand")
            kmer.count[, pos_strand := 0]
          } else {
            setnames(kmer.count, "N", "pos_strand")
            kmer.count[, neg_strand := 0]
          }

          # Bind new k-mers found
          kmer.table <- rbind(kmer.table, kmer.count[!kmer.table])
          setkey(kmer.table, kmer)

          assign("kmer.table", kmer.table, envir = env)
        }

      }, by = eval(if(coor$is.strand.sensitive) "strand")]
    }

    if (!coor$is.strand.sensitive & add.rev.kmers) {
      kmer.table <- countRevCompKmers(kmer.table)
    }

  } else if (coor$k == coor$case.length) {
    coor$status[, kmer := TRUE]
    kmer.table <- extractKmers(coor, genome, rm.overlaps, kmer.table)

    # Change to kmer coordinate first!
  } else if (!coor$status$is.kmer & length(coor$case.length) == 1) {
    stop("Please change to kmer coordinate first!")
  }

  return(kmer.table)
}
