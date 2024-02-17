kmertone <- function(case.coor.path, genome.name="unknown", strand.sensitive, k,
                     ctrl.rel.pos, case.pattern=NULL,
                     output.path="data", case.coor=NULL, genome=NULL,
                     genome.path=NULL, rm.case.kmer.overlaps=TRUE,
                     case.length=NULL, merge.replicates=TRUE, kmer.table=NULL,
                     module = "score", rm.diff.pattern=FALSE, ncpu=1) {
  
  # ----------------------------------------------------------------------------
  # Kmertone program has multiple modules:
  # A) score
  # B) explore 
  # C) 
  #
  # There are two S3-class object introduced by kmertone:
  # 1. <genome> object contains chromosome names, length and sequences.
  # 2. <genomic.coordinate> contains memory-optimised coordinate information
  #    (no redundancy), coordinate status (whether the coordinate refer to
  #    original case or kmer position), case length, k, strand sensitive status,
  #    overlapping k-mer index if applicable, replicate names (if any),
  #    chromosome name(s), and genome name.
  
  # Flag                   Format     Description
  # case.coor.path       <string>     A path to a folder containing either
  #                                   (1) chromosome-separated coordinate files
  #                                       (assume replicates for subfolder) OR
  #                                   (2) bedfile. (assume replicates for
  #                                       bedfiles)
  #
  # genome.name           <string>    Name of available genome: "hg19" or "hg38"
  #                                   Optional: user genome name.
  #
  # genome.path           <string>    A folder path to user own genome. The
  #                                   folder must contain separated chromosome
  #                                   fasta files. The name of the fasta files
  #                                   must be similar to chromosome name in
  #                                   the genomic coordinate table.
  #
  # genome                <genome>    Optional for pre-built <genome> object.
  #                                   Option genome.name and genome.path will be
  #                                   ignored.
  #
  # case.coor         <genomic.coordinate>
  #                                   Optional pre-built <genomic.coordinate>
  #                                   object. Option case.coor.path,
  #                                   case.pattern, strand.sensitive,
  #                                   case.length, k, etc. will be ignored.
  #
  # case.pattern          <string>    A single or multiple DNA pattern.
  #                                   e.g. "G", "TT", "TC", etc.
  #                                   It can be set as NULL for no pattern.
  #
  # k                     <numeric>   The size of a kmer.
  #
  # ctrl.rel.pos          <vector>    A relative position of control regions in
  #                                   a vector format i.e. c(start, end).
  #                                   For example c(80,500) means control
  #                                   regions are 80-500 bases away from the
  #                                   the case site (upstream and downstream)
  #
  # strand.sensitive       <bool>     Mode of strand. In sensitive mode, case
  #                                   kmers are extracted from the sense strand
  #                                   only. In contrast, in insensitive mode,
  #                                   the case kmers are extracted from both.
  #
  # case.length            <int>      Length of case (useful if the length is
  #                                   uniform and only start positions are
  #                                   defined in the genomic coordinates file.)
  #
  # rm.case.kmer.overlaps   <bool>     To remove overlapping case kmers or not
  #
  # merge.replicate        <bool>     If merge.replicate is TRUE, all duplicate
  #                                   positions resulting from replicates will
  #                                   be treated as one unique individual point.
  #
  # kmer.table          <data.table>  A k-mer table for other module calculation
  #
  # mode                   <string>   Mode of action - can be multiple.
  #                                   "score"  : calculate z score
  #                                   "explore": perform EDA
  #                                   "tune"   : find the best length of k-mer
  #
  # output.path            <string>   An output folder name.
  #
  # ncpu                  <numeric>   Number of cpu core to use. Only one core
  #                                   is supported at the moment.
  
  # location of the TrantoR library
  TrantoRLib = "lib/TrantoRext/"
  
  ## Dependant libraries #######################################################
  # order is important due to namespace masking effect
  suppressPackageStartupMessages( library(  RCurl   )  )
  suppressPackageStartupMessages( library(  stringi )  )
  suppressPackageStartupMessages( library(data.table)  )

  ## Dependant functions #######################################################
  source("../lib/Kmertone/lib/buildGenome.R")
  source("../lib/Kmertone/lib/reverseComplement.R")
  source("../lib/Kmertone/lib/buildCoordinate.R")
  source("../lib/Kmertone/lib/bedToCoor.R")
  source("../lib/Kmertone/lib/resolveOverlaps.R")
  source("../lib/Kmertone/lib/kmerize.R")
  source("../lib/Kmertone/lib/buildControl.R")
  source("../lib/Kmertone/lib/removeCaseZone.R")
  source("../lib/Kmertone/lib/initKmerTable.R")
  source("../lib/Kmertone/lib/countRevCompKmers.R")
  source("../lib/Kmertone/lib/extractKmers.R")
  source("../lib/Kmertone/lib/countKmers.R")
  source("../lib/Kmertone/lib/scoreKmers.R")
  source("../lib/Kmertone/lib/saveCoor.R")
  source("../lib/Kmertone/lib/calKmerSkew.R")
  source("../lib/Kmertone/lib/A_detectChr.R")
  source("../lib/Kmertone/lib/A_prepGenome.R")
  source("../lib/Kmertone/lib/A_prepCoordinate.R")
  source("../lib/Kmertone/lib/A_getCaseKmers.R")
  source("../lib/Kmertone/lib/A_getControlKmers.R")
  source("../lib/Kmertone/lib/A_getScore.R")

  dir.create(output.path, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------------------------------------------------------
  # [A] module: score
  # prepGenome + prepCoordinate --> getCaseKmers + getControlKmers --> getScore
  
  if ("score" %in% module) {
    
    # Run per chromosome when genome and/or case.coor are not given
    # Only load necessary data
    if (is.null(genome) | is.null(case.coor)) {
      
      chr.names <- detectChr(case.coor, case.coor.path, genome.path)

      # Initiate k-mer table for both case k-mers and control
      if (k < 16) {
        case.kmers <- initKmerTable(k, case.pattern)
        control.kmers <-
          initKmerTable(k, c(case.pattern,
                             if(!is.null(case.pattern))
                               reverseComplement(case.pattern)))
      } else if (k > 15) {
        # Sample table
        case.kmers <- data.table(kmer = character(), pos_strand = numeric(),
                                 neg_strand = numeric(), key = "kmer")
        control.kmers <- copy(case.kmers)
      }
      
      genome.state <- genome
      case.coor.state <- case.coor
      t1 <- Sys.time()
      
      for (chr.name in chr.names) {
        
        prepGenome      (genome = genome.state, genome.path,
                         genome.name, chr.name)                        -> genome
        prepCoordinate  (case.coor = case.coor.state, case.coor.path,
                         case.length, k, genome.name, strand.sensitive,
                         merge.replicate, case.pattern, chr.name,
                         rm.diff.pattern)                           -> case.coor

        getCaseKmers    (case.coor, genome, rm.case.kmer.overlaps,
                         case.kmers, add.rev.kmers = FALSE)        -> case.kmers

        getControlKmers (case.coor, ctrl.rel.pos, genome,
                         control.kmers, add.rev.kmers = FALSE,
                         output.path)                           -> control.kmers
        
      }
      
      if (!case.coor$is.strand.sensitive) {
        countRevCompKmers(case.kmers)
      }
      countRevCompKmers(control.kmers)
      
      getScore        (case.kmers, control.kmers, case.pattern,
                       k, output.path)                             -> kmer.table
        
      t <- Sys.time() - t1
      cat("Final time taken:", t[[1]], attr(t, "units"), "\n")
      cat(paste(c(rep("-", 80), "\n"), collapse = ""))
      
    } else {
      
      t1 <- Sys.time()
      
      prepGenome      (genome, genome.path, genome.name, chr.name)     -> genome
      prepCoordinate  (case.coor, case.coor.path, case.length, k,
                       genome.name, strand.sensitive, merge.replicate,
                       case.pattern, chr.name)                      -> case.coor
      getCaseKmers    (case.coor, genome, rm.case.kmer.overlaps)   -> case.kmers
      getControlKmers (case.coor, ctrl.rel.pos, genome,
                       output.path = output.path)               -> control.kmers
      getScore        (case.kmers, control.kmers, case.pattern,
                       k, output.path)                             -> kmer.table
      
      t <- Sys.time() - t1
      cat("Final time taken:", t[[1]], attr(t, "units"), "\n")
      cat(paste(c(rep("-", 80), "\n"), collapse = ""))
      
    }
    
    if (length(mode) == 1) {
      return(kmer.table)
    }
  }
}