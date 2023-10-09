Breakpoints <- R6::R6Class(
    classname = "Breakpoints",
    public = list(
        #' @field k Numeric vector of k-mer size.
        k = 8,

        #' @field bp_dir Numeric vector to fix the seed.
        seed = 1234,

        #' @field true_breaks GRanges of true breakpoints.
        true_breaks = NULL,

        #' @field true_breaks_expanded GRanges list of break ranges.
        true_breaks_expanded = NULL,

        #' @field control_breaks GRanges of control breakpoints.
        control_breaks = NULL,

        #' @field control_breaks_expanded GRanges list of break ranges.
        control_breaks_expanded = NULL,

        initialize = function(k, exp, seed, scores_with_kmers, true_prop,
                              assembly, break_score, break_score_all){
            if(!missing(k)) self$k <- k
            if(!missing(exp)) private$bp_dir <- paste0("../data/experiments/", exp)
            if(!missing(seed)) self$seed <- seed
            if(!missing(scores_with_kmers)) private$scores_with_kmers <- scores_with_kmers
            if(!missing(true_prop)) private$true_prop <- true_prop
            if(private$true_prop <= 0.1) stop("Too few true breakpoint samples!")
            
            # get reference sequence
            assembly <- ifelse(is.null(assembly), "hg19", assembly)
            private$get_ref_seq(assembly)

            if(!missing(break_score)) private$break_score <- break_score
            private$break_score_all <- ifelse(private$break_score == "all", TRUE, FALSE)
        },

        #' @description
        #' Function extracts true and control breakpoints.
        #' @param break_type Character vector of any of c("biological", "enzymatic", "high_frequency"). 
        #'  If user knows the type of breakage but does not have any range effect values, the maximum
        #'  range effect cross similar types of breakages will be taken. If no type is given, user 
        #'  has to input 3 range effects into the ranges parameter.
        #' @param ranges Numeric vector of length 3, where position 1 is the short range effect, 
        #'  position 2 is the medium range effect, and position 3 is the long range effect.
        #' @return None.
        get_breaks = function(break_type = "biological", ranges = NULL){
            # get range cutoffs
            private$get_ranges(break_type = break_type, ranges = ranges)

            # extend query table
            private$get_extended_querytable()

            #' @description
            #' Wrapper function that extracts the true and control breakpoint positions
            #' and extract the associated intrinsic susceptibility scores.
            #' @return None.
            run_functions = function(){
                msg <- paste0("Extracting for breakage score: ", private$break_score)
                l <- paste0(rep(".", 70-nchar(msg)))
                cat(msg, l, "\n", sep = "")

                # extract the true breakpoint locations
                private$get_true_breaks()

                # sample control breakpoints from control breakpoint regions
                if(!private$regression & private$get_controls){
                    private$get_control_breaks()
                    private$get_break_scores(type = "controls")
                }

                # as size of control breakpoints is much smaller than the true breakpoints,
                # after sub-sampling true breakpoints to reach an approximate balance bewteen true
                # and control breakpoint populations, only then we calculate the breakpoint scores
                private$get_break_scores(type = "breaks")
            }

            if(private$break_score == "all"){
                for(break_score in c("zscore", "ratio")){
                    private$break_score <- break_score
                    run_functions()
                }
            } else if(private$break_score == "zscore" || private$break_score == "ratio"){
                run_functions()
            } else{
                stop("Break_score must be all, zscore or ratio!")
            }
            
            private$expand_ranges()
        }
        
    ),
    private = list(
        #' @field bp_dir Character vector of directory containing breakpoints.
        bp_dir = NULL,

        #' @field ranges Numeric vector of upper limit of the short, mid and long ranges.
        ranges = NULL,

        #' @field break_score Character vector of "zscore" or "ratios".
        break_score = NULL,

        #' @field true_prop Numeric vector. The proportion of true breakpoint samples
        #'  compared to the control breakpoint samples. By default, will sample 
        #'  110% of the control breakpoints to have a slight shift 
        #'  towards true breakpoints. 
        true_prop = 1,       

        #' @field break_score_all Boolean. If TRUE, extract both "zscore" and "ratios".
        break_score_all = FALSE,

        #' @field scores_with_kmers Boolean. If TRUE, breakage scores return with k-mer col.
        scores_with_kmers = FALSE,

        #' @field ref BSgenome reference genome assembly.
        ref = NULL,

        #' @field chrs_len Numeric vector of the length of each chromosome of the reference sequence.
        chrs_len = NULL,

        #' @field which_chr Numeric vector. Specify which chromosome to work with, default are autosomes.
        which_chr = 1:22,

        #' @field kmer_window Numeric vector of the window size for counting k-mers.
        kmer_window = NULL,

        #' @field kmer_ref Data.Table of forward and reverse complement k-mers.
        kmer_ref = NULL,

        #' @field kmer_list Character vector of all k-mers.
        kmer_list = NULL,

        #' @field extended_preparam_table Data Table of the extended query table.
        extended_preparam_table = NULL,

        #' @description 
        #' Extract reference sequence.
        #' @param assembly Character Vector of genome assembly versions. 
        #' @return None.
        get_ref_seq = function(assembly){
            private$ref <- switch(assembly,
                "hg18" = BSgenome.Hsapiens.UCSC.hg18::BSgenome.Hsapiens.UCSC.hg18,
                "hg19" = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                "hg38" = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                "hs37d5" = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
            )

            # get the start and end positions of each chromosome
            refseq.table <- as.data.frame(private$ref@seqinfo)
            refseq.table <- refseq.table[grepl(
                pattern = "^chr([1-9]|1[0-9]|2[0-2])$", 
                x = rownames(refseq.table)
            ),]
            private$chrs_len <- refseq.table$seqlengths            
            private$chrs_len <- setNames(private$chrs_len, rownames(refseq.table))
        },

        #' @description 
        #' Extracts the short, mid and long-range upper limits around the central breakpoint.
        #' @param break_type Character vector of any of c("biological", "enzymatic", "high_frequency"). 
        #'  If user knows the type of breakage but does not have any range effect values, the maximum
        #'  range effect cross similar types of breakages will be taken. If no type is given, user 
        #'  has to input 3 range effects into the ranges parameter.
        #' @param ranges Numeric vector of length 3, where position 1 is the short range effect, 
        #'  position 2 is the medium range effect, and position 3 is the long range effect.
        #' @return None.
        get_ranges = function(break_type = "biological", ranges = NULL){
            if(!is.null(ranges)){
                all.ranges <- list(
                    short.range = ranges[1],
                    mid.range = ranges[2],
                    long.range = ranges[3]                    
                )
            } else {
                datatable <- fread(paste0(
                    "../data/range_effects/", 
                    "MaxValuesFromClustersByType.csv"
                ))

                # extract each range
                row.id <- which(datatable$type == break_type)
                all.ranges <- as.list(datatable[row.id,-"type"])
            }

            private$ranges$short_range <- all.ranges$short.range
            private$ranges$mid_range <- all.ranges$mid.range
            private$ranges$long_range <- all.ranges$long.range            
        },

        #' @description 
        #' Generate k-mer table for reference when counting k-mer frequencies.
        #' @param k Numeric vector. Size of the k-mer. 
        #' @return None.
        generate_kmer_table = function(k){
            if(missing(k)) k <- private$kmer_window            
            all.kmers <- do.call(
                data.table::CJ, 
                rep(list(c("A", "C", "G", "T")), k)
            )
            private$kmer_list <- all.kmers[, do.call(paste0, .SD)]
            kmer_ref <- data.table(
                'kmer' = private$kmer_list,
                'rev.comp' = as.character(
                    Biostrings::reverseComplement(
                        Biostrings::DNAStringSet(private$kmer_list)
                    )
                )
            )
            kmer_ref[, `:=`(cond = ifelse(seq(1:nrow(.SD)) < match(kmer, rev.comp), 
            TRUE, ifelse(kmer == rev.comp, TRUE, FALSE)))]
            private$kmer_ref <- kmer_ref[cond == TRUE, .(kmer, rev.comp)]
        },

        #' @description
        #'  Extend table to have all k-mers not just lexicographically occurring ones.
        #' @return Data table of the query table.
        get_extended_querytable = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Retrieving table of breakage parameters"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            extend_querytable <- function(dat, bp.of.interest){
                dat <- dat[category == bp.of.interest, -"category"]
                dat <- t(dat)
                df <- as.data.table(dat)
                df[, kmer := rownames(dat)]
                setnames(df, c("prob", "kmer"))
                setcolorder(df, c("kmer", "prob"))

                # extend table to have all k-mers, not just 
                # lexicographically in the correct order
                rev.comp <- as.character(
                    Biostrings::reverseComplement(
                        Biostrings::DNAStringSet(df$kmer)
                    )
                )
                df[, rev_comp := rev.comp]

                k.mers <- do.call(
                    data.table::CJ,
                    rep(list(c("A", "C", "G", "T")), self$k)
                )
                kmer_list <- k.mers[, do.call(paste0, .SD)]
                kmer_ref <- data.table(kmer=kmer_list)
                kmer_ref[, prob := 0]

                # forward k-mers
                ind_fwd <- match(df$kmer, kmer_ref$kmer)    
                kmer_ref$prob[ind_fwd] <- df$prob

                # reverse complement k-mers
                ind_rc <- match(df$rev_comp, kmer_ref$kmer)
                ind_rc[is.na(ind_rc)] <- (nrow(kmer_ref)-nrow(df)+1):(nrow(kmer_ref))
                kmer_ref$prob[ind_rc] <- df$prob

                # clean table
                kmers <- kmer_ref$kmer
                kmer_ref <- kmer_ref[, "prob"]
                setnames(kmer_ref, bp.of.interest)
                kmer_ref <- as.data.frame(kmer_ref)
                rownames(kmer_ref) <- kmers
                kmer_ref <- t(kmer_ref)
                df[, rev_comp := NULL]

                return(kmer_ref)
            }

            private$kmer_window <- 8
            private$generate_kmer_table()

            # import all values from pre-parameterised table
            preparam.table <- fread(
                paste0("../data/kmertone/QueryTable/QueryTable_",
                    "kmer-", private$kmer_window, "_",
                    private$break_score, ".csv"), 
                showProgress = FALSE
            )
            category_col <- preparam.table[, "category"]
            res <- lapply(preparam.table$category, function(b){
                extend_querytable(
                    dat = preparam.table, 
                    bp.of.interest = b
                )
            })
            res <- do.call(rbind, res)
            private$extended_preparam_table <- cbind(category_col, res)

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description 
        #' Extracts true breakpoint locations.
        #' @return None.
        get_true_breaks = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Extracting true breakpoints"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            bp.folder <- paste0(private$bp_dir, "/breakpoint_positions/chr")

            # get breakpoints for each chromosome
            datatable <- lapply(private$which_chr, function(chr){
                datatable <- fread(
                    paste0(bp.folder, chr, ".csv"), 
                    showProgress = FALSE
                )
                if("freq" %in% colnames(datatable)){
                    datatable[, `:=`(lev.dist = NULL, freq = NULL)]
                }
                
                to_rename <- which(grepl(pattern = "start", x = colnames(datatable)))
                colnames(datatable)[to_rename] <- "start"

                datatable[, `:=`(width = 2, seqnames = paste0("chr", chr))]
                return(datatable)
            }) 
            datatable <- rbindlist(datatable)
            self$true_breaks[[private$break_score]]$ref <- plyranges::as_granges(datatable)

            # expand breakpoint into k-mer
            self$true_breaks[[private$break_score]]$kmer <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref), 
                self$k-2
            )

            # discard any out-of-bounds ranges
            check.longrange <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref), 
                private$ranges$long_range-1
            )
            to.keep <- which(start(check.longrange) > 0)

            self$true_breaks[[private$break_score]]$kmer <- 
                self$true_breaks[[private$break_score]]$kmer[to.keep]
            self$true_breaks[[private$break_score]]$ref <- 
                self$true_breaks[[private$break_score]]$ref[to.keep]
            
            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description 
        #' Extracts control breakpoint locations.
        #' @return None.
        get_control_breaks = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Sampling control breakpoints"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # extract control regions
            control.folder <- paste0(private$bp_dir, "/control_coordinates/chr")

            # get breakpoints for each chromosome
            datatable <- lapply(private$which_chr, function(chr){
                datatable <- fread(
                    paste0(control.folder, chr, ".csv"), 
                    showProgress = FALSE
                )
                datatable[, seqnames := paste0("chr", chr)]
                return(datatable)
            })
            datatable <- rbindlist(datatable)
            datatable <- plyranges::as_granges(datatable)

            # sample rows proportional to width
            to.sample.max <- length(self$true_breaks[[private$break_score]]$ref)
            set.seed(self$seed)
            sample.rows <- sample(length(datatable), size = to.sample.max, 
                                  replace = TRUE, prob = width(datatable))
            sample.rows <- datatable[sample.rows]

            # sample within each chosen range
            sample.rows <- as.data.table(sample.rows)
            sample.rows[, position := sample(start:end, size = 1, replace = FALSE), 
                        by = 1:dim(sample.rows)[1]]
            sample.rows[, start := position]
            sample.rows[, end := position+1]
            sample.rows[, c("width", "strand", "position") := NULL]
            datatable <- plyranges::as_granges(sample.rows)

            # extract mid-range distance of true breakpoint locations
            true.bp.expanded <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref), 
                private$ranges$long_range
            )

            # exclude non-overlapping control breakpoints from expanded true breakpoint regions
            control.bp.expanded <- plyranges::stretch(
                plyranges::anchor_center(datatable), 
                private$ranges$long_range
            )

            control.bp.expanded <- plyranges::filter_by_non_overlaps(
                control.bp.expanded, 
                true.bp.expanded
            )            

            # shift expanded control breakpoints back into original 2-base breakpoint location
            datatable <- plyranges::shift_right(plyranges::stretch(
                plyranges::anchor_center(control.bp.expanded), 
                -private$ranges$long_range
            ), shift = 1L)            
            datatable <- arrange(datatable, seqnames, start)
            datatable <- unique(datatable)

            # expand control breakpoint into k-mer
            self$control_breaks[[private$break_score]]$ref <- datatable
            self$control_breaks[[private$break_score]]$kmer <- plyranges::stretch(
                plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref), 
                self$k-2
            )

            # discard any out-of-bounds ranges
            check.longrange <- plyranges::stretch(
                plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref), 
                private$ranges$long_range-1
            )

            to.keep <- sapply(paste0("chr", private$which_chr), function(chrs){
                return(match(filter(check.longrange, 
                    seqnames == chrs & start > 0 & end <= private$chrs_len[[chrs]]),
                    check.longrange
                ))
            }, USE.NAMES = FALSE)
            to.keep <- unlist(to.keep)

            self$control_breaks[[private$break_score]]$kmer <- 
                self$control_breaks[[private$break_score]]$kmer[to.keep]
            self$control_breaks[[private$break_score]]$ref <- 
                self$control_breaks[[private$break_score]]$ref[to.keep]

            # randomly sample subset of true breakpoints to 
            # achieve a balance of true and control breakpoints.
            # Optionally, can create an imbalance for crash test demonstrations.
            true.bp <- self$true_breaks[[private$break_score]]$ref
            to.sample.max <- length(
                self$control_breaks[[private$break_score]]$ref
            )*private$true_prop

            set.seed(self$seed)
            sample.rows <- sample(length(true.bp), size = to.sample.max, 
                                  replace = TRUE, prob = width(true.bp))
            sample.rows <- true.bp[sample.rows]
            true.bp <- arrange(sample.rows, seqnames, start)
            true.bp <- unique(true.bp)
            
            # expand control breakpoint into k-mer
            self$true_breaks[[private$break_score]]$ref <- true.bp
            self$true_breaks[[private$break_score]]$kmer <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref), 
                self$k-2
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description 
        #' Extracts true and control breakpoint scores.
        #' @param data GRanges of breakpoints. 
        #' @param ref.data GRanges of breakpoints expanded into k-mers.
        #' @param type Character vector of "breaks" or "controls".
        #' @return None.
        get_break_scores = function(type){
            # progress message
            t1 <- Sys.time()
            cur.msg <- ifelse(type == "breaks",
                              "Getting scores for true breaks",
                              "Getting scores for control breaks")
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            data <- switch(type,
                "breaks" = self$true_breaks[[private$break_score]]$kmer,
                "controls" = self$control_breaks[[private$break_score]]$kmer
            )
            ref.data <- switch(type,
                "breaks" = self$true_breaks[[private$break_score]]$ref,
                "controls" = self$control_breaks[[private$break_score]]$ref
            )

            # import query table of z-scores or probability ratios
            datatable.category <- private$extended_preparam_table$category
            datatable.category <- paste0(
                datatable.category, 
                "_", private$break_score
            )
            datatable <- private$extended_preparam_table[, -"category"]

            # get sequences of k-mer breakages        
            datatable.seq <- Biostrings::getSeq(private$ref, data)

            # match corresponding k-mer in query table
            query.kmers <- Biostrings::DNAStringSet(colnames(datatable))
            match.kmers <- match(datatable.seq, query.kmers)
            # drop remaining NAs as sequences contain "N"
            to.remove <- is.na(match.kmers)
            match.kmers <- match.kmers[!to.remove]
            datatable.seq <- datatable.seq[!to.remove]
            datatable.kmer <- data[!to.remove]
            datatable.ref <- ref.data[!to.remove]

            # extract columns based on k-mer matches
            if(private$scores_with_kmers){
                datatable <- datatable[, ..match.kmers]
                rownames(datatable) <- datatable.category
                datatable <- t(datatable)
                datatable <- data.table(kmer = rownames(datatable), scores = datatable)
                setnames(datatable, c("kmer", datatable.category))
            } else {
                datatable <- as.matrix(datatable[, ..match.kmers])
                colnames(datatable) <- NULL
                rownames(datatable) <- datatable.category
            }

            if(type == "breaks"){
                self$true_breaks[[private$break_score]]$ref <- datatable.ref
                self$true_breaks[[private$break_score]]$kmer <- datatable.kmer
                if(private$scores_with_kmers){
                    self$true_breaks[[private$break_score]]$scores <- datatable
                } else {
                    self$true_breaks[[private$break_score]]$scores <- t(datatable)
                }
            } else if(type == "controls"){
                self$control_breaks[[private$break_score]]$ref <- datatable.ref
                self$control_breaks[[private$break_score]]$kmer <- datatable.kmer
                if(private$scores_with_kmers){
                    self$control_breaks[[private$break_score]]$scores <- datatable
                } else {
                    self$control_breaks[[private$break_score]]$scores <- t(datatable)
                }
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description 
        #' Extracts the short, mid and long-range GRanges around the central breakpoint.
        #' @return None.
        expand_ranges = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Expanding GRanges into short, mid and long-ranges"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # true breakpoints
            self$true_breaks_expanded$short_range <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref), 
                private$ranges$short_range-1
            )

            self$true_breaks_expanded$mid_range <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref), 
                private$ranges$mid_range-1
            )

            self$true_breaks_expanded$long_range <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref), 
                private$ranges$long_range-1
            )

            # change end-position if out of bounds
            col_names <- colnames(as.data.table(self$true_breaks_expanded$long_range))
            if(private$regression & "Breaks" %in% col_names){
                self$true_breaks_expanded$long_range <- self$true_breaks_expanded$long_range %>% 
                    as_tibble() %>% 
                    dplyr::mutate(
                        seqnames = as.character(seqnames),
                        max_end = as.numeric(private$chrs_len[seqnames])
                    ) %>% 
                    dplyr::mutate(
                        end = ifelse(end <= max_end, end, max_end)
                    ) %>%
                    dplyr::select(seqnames, start, end, Breaks) %>% 
                    plyranges::as_granges()
            }

            # control breakpoints
            if(!private$regression & private$get_controls){
                self$control_breaks_expanded$short_range <- plyranges::stretch(
                    plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref), 
                    private$ranges$short_range-1
                )

                self$control_breaks_expanded$mid_range <- plyranges::stretch(
                    plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref), 
                    private$ranges$mid_range-1
                )

                self$control_breaks_expanded$long_range <- plyranges::stretch(
                    plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref), 
                    private$ranges$long_range-1
                )

                # change end-position if out of bounds
                col_names <- colnames(as.data.table(self$control_breaks_expanded$long_range))
                if(private$regression & "Breaks" %in% col_names){
                    self$control_breaks_expanded$long_range <- self$control_breaks_expanded$long_range %>% 
                        as_tibble() %>% 
                        dplyr::mutate(
                            seqnames = as.character(seqnames),
                            max_end = as.numeric(private$chrs_len[seqnames])
                        ) %>% 
                        dplyr::mutate(
                            end = ifelse(end <= max_end, end, max_end)
                        ) %>%
                        dplyr::select(seqnames, start, end, Breaks) %>% 
                        plyranges::as_granges()
                }
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        }
    )
)