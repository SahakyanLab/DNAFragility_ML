Breakpoints <- R6::R6Class(
    classname = "Breakpoints",
    public = list(
        #' @field k Numeric vector of k-mer size.
        k = 4,

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

        initialize = function(k, exp, seed, assembly){
            if(!missing(k)) self$k <- k
            if(!missing(exp)) private$bp_dir <- paste0("../data/experiments/", exp)
            if(!missing(seed)) self$seed <- seed
            if(!missing(assembly)) private$ref <- private$get_ref_seq(assembly)
        },

        #' @description
        #' Function extracts true and control breakpoints.
        #' @param break_score Character Vector of "all", "zscore" or "ratio".
        #' @param scores_with_kmers Boolean. If TRUE, scores return with k-mer col.
        #' @return None.
        get_breaks = function(break_score = "all",
                              scores_with_kmers = FALSE){
            private$scores_with_kmers <- scores_with_kmers

            # get range cutoffs
            private$get_ranges()

            run_functions = function(){
                msg <- paste0("Extracting for breakage score: ",
                              private$break_score)
                l <- paste0(rep(".", 70-nchar(msg)))
                cat(msg, l, "\n", sep = "")

                private$get_true_breaks()
                private$get_break_scores(type = "breaks")
                private$get_control_breaks()
                private$get_break_scores(type = "controls")
            }

            if(break_score == "all"){
                private$break_score_all <- TRUE
                for(break_score in c("zscore", "ratio")){
                    private$break_score <- break_score
                    run_functions()
                }
            } else if(break_score == "zscore" || break_score == "ratio"){
                private$break_score <- break_score
                private$break_score_all <- FALSE
                run_functions()
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

        #' @field break_score_all Boolean. If TRUE, extract both "zscore" and "ratios".
        break_score_all = FALSE,

        #' @field scores_with_kmers Boolean. If TRUE, breakage scores return with k-mer col.
        scores_with_kmers = FALSE,

        #' @field ref BSgenome reference genome assembly.
        ref = NULL,

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
        },

        #' @description 
        #' Extracts the short, mid and long-range upper limits around the central breakpoint.
        #' @return None.
        get_ranges = function(){
            datatable <- fread(
                paste0("../data/Ranges_cutoffs.csv"),
                select = c("Cluster", paste0("kmer_", self$k))
            )

            # convert values to ceiling
            datatable[, paste0("kmer_", self$k) := ceiling(datatable[, 2])]

            # extract each range
            private$ranges$short_range <- as.numeric(datatable["short.range", 2, on = "Cluster"])
            private$ranges$mid_range <- as.numeric(datatable["mid.range", 2, on = "Cluster"])
            private$ranges$long_range <- as.numeric(datatable["long.range", 2, on = "Cluster"])
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

            self$true_breaks_expanded$kmer_occupancy <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref), 
                4L
            )

            self$true_breaks_expanded$qm_left_kmer <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref %>% 
                                            dplyr::mutate(width = 1)),
                6L
            )

            self$true_breaks_expanded$qm_right_kmer <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref %>% 
                                            dplyr::mutate(width = 1) %>% 
                                            plyranges::shift_right(1)),
                6L
            )

            # control breakpoints
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

            self$control_breaks_expanded$kmer_occupancy <- plyranges::stretch(
                plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref), 
                4L
            )

            self$control_breaks_expanded$qm_left_kmer <- plyranges::stretch(
                plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref %>% 
                                            dplyr::mutate(width = 1)),
                6L
            )

            self$control_breaks_expanded$qm_right_kmer <- plyranges::stretch(
                plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref %>% 
                                            dplyr::mutate(width = 1) %>% 
                                            plyranges::shift_right(1)),
                6L
            )

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
            datatable <- lapply(1:22, function(chr){
                datatable <- fread(paste0(bp.folder, chr, ".csv"), 
                            showProgress = FALSE)
                setnames(datatable, "start")
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

            # extract mid-range distance of true breakpoint locations
            true.bp <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref), 
                min(private$ranges$mid_range*2, 200)
            )

            # true.bp.left <- plyranges::shift_left(
            #     plyranges::anchor_start(self$true_breaks[[private$break_score]]$ref), 
            #     private$ranges$long_range
            # )
            # true.bp.right <- plyranges::shift_right(
            #     plyranges::anchor_end(self$true_breaks[[private$break_score]]$ref), 
            #     private$ranges$long_range
            # )
            # true.bp <- c(true.bp.left, true.bp.right)
            # true.bp <- arrange(true.bp, seqnames, start)
            # true.bp <- plyranges::unanchor(true.bp)

            # extract control regions
            control.folder <- paste0(private$bp_dir, "/control_coordinates/chr")

            # get breakpoints for each chromosome
            datatable <- lapply(1:22, function(chr){
                datatable <- fread(
                    paste0(control.folder, chr, ".csv"), 
                    showProgress = FALSE
                )
                datatable[, seqnames := paste0("chr", chr)]
                return(datatable)
            })
            datatable <- rbindlist(datatable)
            datatable <- plyranges::as_granges(datatable)

            # find non-overlapping ranges to sample from 
            datatable <- plyranges::filter_by_non_overlaps(datatable, true.bp)
            to.sample.max <- length(true.bp)*1.5

            # sample rows proportional to width
            set.seed(self$seed)
            sample.rows <- sample(length(datatable), size = to.sample.max, 
                                  replace = TRUE, prob = width(datatable))
            sample.rows <- datatable[sample.rows]

            # sample uniformly within each chosen range
            sample.rows <- as.data.table(sample.rows)
            sample.rows[, position := sample(start:end, 
                                             size = 1, 
                                             replace = FALSE), 
                        by = 1:dim(sample.rows)[1]]
            sample.rows[, start := position]
            sample.rows[, end := position+1]
            sample.rows[, c("width", "strand", "position") := NULL]
            datatable <- plyranges::as_granges(sample.rows)
            datatable <- plyranges::filter_by_non_overlaps(datatable, true.bp)
            datatable <- arrange(datatable, seqnames, start)
            datatable <- unique(datatable)

            # expand breakpoint into k-mer
            self$control_breaks[[private$break_score]]$ref <- sort(unique(datatable))
            self$control_breaks[[private$break_score]]$kmer <- plyranges::stretch(
                plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref), 
                self$k-2
            )

            # discard any out-of-bounds ranges
            check.longrange <- plyranges::stretch(
                plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref), 
                private$ranges$long_range-1
            )
            to.keep <- which(start(check.longrange) > 0)
            self$control_breaks[[private$break_score]]$kmer <- 
                self$control_breaks[[private$break_score]]$kmer[to.keep]
            self$control_breaks[[private$break_score]]$ref <- 
                self$control_breaks[[private$break_score]]$ref[to.keep]

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
            datatable <- fread(
                paste0("../data/kmertone/QueryTable/QueryTable_",
                       "kmer-", self$k, "_",
                       private$break_score, ".csv"), 
                showProgress = FALSE
            )
            datatable.category <- datatable$category
            datatable.category <- paste0(datatable$category, "_", private$break_score)
            datatable <- datatable[, -"category"]

            # drop any infinite and NA values
            temp <- t(datatable)
            datatable <- t(temp[is.finite(matrixStats::rowSums2(temp, na.rm = TRUE)), ])
            datatable <- as.data.table(datatable)

            # get sequences of k-mer breakages        
            datatable.seq <- Biostrings::getSeq(private$ref, data)

            # match corresponding k-mer in query table
            query.kmers <- Biostrings::DNAStringSet(colnames(datatable))

            # forward k-mers
            match.kmers <- match(datatable.seq, query.kmers)
            # reverse k-mers
            any.na <- is.na(match.kmers)
            rev.kmers <- Biostrings::reverseComplement(datatable.seq[any.na])
            match.rev.kmers <- match(rev.kmers, query.kmers)
            # combine results
            match.kmers[any.na] <- match.rev.kmers
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
        }
    )
)