Breakpoints <- R6::R6Class(
    classname = "Breakpoints",
    public = list(
        #' @field label Character vector of the species label.
        label = NULL,

        #' @field k Numeric vector of k-mer size.
        k = 8,

        #' @field bp_dir Numeric vector to fix the seed.
        seed = 1234,

        #' @field df_bp Data.table of positions to extract.
        df_bp = NULL,

        #' @field true_breaks GRanges of true breakpoints.
        true_breaks = NULL,

        #' @field true_breaks_expanded GRanges list of break ranges.
        true_breaks_expanded = NULL,

        #' @field out Numeric vector. If -1, skips to the next chunk of the sequence for analysis. 
        out = 1,

        initialize = function(fasta_sequence, label, k, seed, break_score, df_bp){
            if(!missing(fasta_sequence)) private$ref <- fasta_sequence
            if(!missing(label)) self$label <- label
            if(!missing(k)) self$k <- k
            if(!missing(seed)) self$seed <- seed
            if(!missing(break_score)) private$break_score <- break_score
            if(!missing(df_bp)) self$df_bp <- df_bp
        },

        #' @description
        #' Function extracts true and control breakpoints.
        #' @param break_type Character vector of any of c("Biological", "Enzymatic", "High_frequency").  
        #'  If user knows the type of breakage but does not have any range effect values, the maximum
        #'  range effect cross similar types of breakages will be taken. If no type is given, user 
        #'  has to input 3 range effects into the ranges parameter.
        #' @param ranges Numeric vector of length 3, where position 1 is the short range effect, 
        #'  position 2 is the medium range effect, and position 3 is the long range effect.
        #' @return None.
        get_breaks = function(break_type = "Biological", ranges = NULL){
            # get range cutoffs
            private$get_ranges(break_type = break_type, ranges = ranges)

            # extend query table
            private$kmer_window <- self$k
            private$get_extended_querytable(kmer_size = private$kmer_window)

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
                private$get_break_scores(type = "breaks")
                if(self$out == -1) return(-1)
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
            
            if(self$out == -1) return(-1)
            private$expand_ranges()
        }
        
    ),
    private = list(
        #' @field ranges Numeric vector of upper limit of the short, mid and long ranges.
        ranges = NULL,

        #' @field break_score Character vector of "zscore" or "ratios".
        break_score = NULL,

        #' @field ref BSgenome reference genome assembly.
        ref = NULL,

        #' @field kmer_window Numeric vector of the window size for counting k-mers.
        kmer_window = NULL,

        #' @field kmer_ref Data.Table of forward and reverse complement k-mers.
        kmer_ref = NULL,

        #' @field kmer_list Character vector of all k-mers.
        kmer_list = NULL,

        #' @field extended_preparam_table Data Table of the extended query table.
        extended_preparam_table = NULL,

        #' @field datatable_seq GRanges object of expanded DNA sequences. 
        datatable_seq = NULL,

        #' @description 
        #' Extracts the short, mid and long-range upper limits around the central breakpoint.
        #' @param break_type Character vector of any of c("Biological", "Enzymatic", "High_frequency"). 
        #'  If user knows the type of breakage but does not have any range effect values, the maximum
        #'  range effect cross similar types of breakages will be taken. If no type is given, user 
        #'  has to input 3 range effects into the ranges parameter.
        #' @param ranges Numeric vector of length 3, where position 1 is the short range effect, 
        #'  position 2 is the medium range effect, and position 3 is the long range effect.
        #' @return None.
        get_ranges = function(break_type = "Biological", ranges = NULL){
            round_to_nearest_even <- function(x) round(x/2)*2

            if(!is.null(ranges)){
                ranges <- unname(ranges)
                all.ranges <- list(
                    short = min(ranges),
                    medium = median(ranges),
                    long = max(ranges)
                )
            } else {
                datatable <- fread(paste0(
                    "../data/range_effects/", 
                    "MaxValuesFromClustersByType.csv"
                ))

                # extract each range
                row.id <- which(datatable$break_type == break_type)
                all.ranges <- datatable[row.id, -"break_type"]
                all.ranges <- as.list(all.ranges)
            }

            # round to the nearest even number
            all.ranges <- lapply(all.ranges, round_to_nearest_even)

            private$ranges$short_range <- all.ranges$short
            private$ranges$mid_range <- all.ranges$medium
            private$ranges$long_range <- all.ranges$long
        },

        #' @description 
        #' Generate k-mer table for reference when counting k-mer frequencies.
        #' @return None.
        generate_kmer_table = function(){
            all.kmers <- do.call(
                data.table::CJ, 
                rep(list(c("A", "C", "G", "T")), private$kmer_window)
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
        #' @param kmer_size Numeric vector. Size of the k-mer. 
        #' @return Data table of the query table.
        get_extended_querytable = function(kmer_size, custom_cols = TRUE){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Retrieving table of breakage parameters"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            # cat(paste0(cur.msg, l))

            extend_querytable <- function(dat, bp.of.interest, kmer_size){
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
                    rep(list(c("A", "C", "G", "T")), kmer_size)
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

            check.file <- paste0(
                "../data/kmertone/QueryTable/Extended_QueryTable_",
                "kmer-", private$kmer_window, "_",
                private$break_score, ".csv"
            )
            if(file.exists(check.file)){
                private$extended_preparam_table <- fread(
                    check.file,
                    showProgress = FALSE
                )

                if(private$only_breaks){
                    tfbs.ind <- which(!grepl(pattern = "^TFBS_", private$extended_preparam_table$category))
                    private$extended_preparam_table <- private$extended_preparam_table[tfbs.ind,]

                    # epigenome marks
                    epigenome.ind <- which(!grepl(
                        pattern = "Epigenome|ATACseq|Chipseq|Dnaseseq|FAIREseq", 
                        private$extended_preparam_table$category
                    ))
                    private$extended_preparam_table <- private$extended_preparam_table[epigenome.ind,]
                }
            } else {
                private$generate_kmer_table()

                # import all values from pre-parameterised table
                preparam.table <- fread(
                    paste0("../data/kmertone/QueryTable/QueryTable_",
                        "kmer-", private$kmer_window, "_",
                        private$break_score, ".csv"), 
                    showProgress = FALSE
                )
                category_col <- preparam.table[, "category"]

                if(private$only_breaks){
                    cols_to_rm <- !grepl(
                        pattern = paste0(
                            "^TFBS_|^Epigenome|^ATACseq|",
                            "^Chipseq|^Dnaseseq|^FAIREseq"
                        ),
                        x = category_col$category,
                        ignore.case = TRUE
                    )
                    category_col <- category_col[cols_to_rm,]
                    preparam.table <- preparam.table[cols_to_rm,]
                }

                res <- lapply(preparam.table$category, function(b){
                    extend_querytable(
                        dat = preparam.table, 
                        bp.of.interest = b,
                        kmer_size = kmer_size
                    )
                })
                res <- do.call(rbind, res)
                private$extended_preparam_table <- cbind(category_col, res)

                fwrite(
                    private$extended_preparam_table,
                    file = paste0(
                        "../data/kmertone/QueryTable/Extended_QueryTable_",
                        "kmer-", private$kmer_window, "_",
                        private$break_score, ".csv"
                    ),
                    showProgress = FALSE
                )
            }

            total.time <- Sys.time() - t1
            # cat("DONE! --", signif(total.time[[1]], 2), 
            #     attr(total.time, "units"), "\n")
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

            datatable <- self$df_bp
            datatable[, `:=`(seqnames = self$label, width = 2)]
            setnames(datatable, c("start", "seqnames", "width"))
            datatable[, start.pos := start]

            self$true_breaks[[private$break_score]]$ref <- plyranges::as_granges(datatable)
            
            # expand breakpoint into k-mer
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
            datatable.seq <- BSgenome::getSeq(private$ref, data)

            # check if all NAs or not
            any_non_nas <- which(!grepl(pattern = "N", x = datatable.seq))
            if(length(any_non_nas) == 0){
                self$out <- -1
                return(self$out)
            }

            # match corresponding k-mer in query table
            query.kmers <- Biostrings::DNAStringSet(colnames(datatable))
            match.kmers <- match(datatable.seq, query.kmers)
            
            # drop remaining NAs
            to.remove <- is.na(match.kmers)
            match.kmers <- match.kmers[!to.remove]
            datatable.seq <- datatable.seq[!to.remove]
            datatable.kmer <- data[!to.remove]
            datatable.ref <- ref.data[!to.remove]

            # re-format data
            datatable <- as.matrix(datatable[, ..match.kmers])
            colnames(datatable) <- NULL
            rownames(datatable) <- datatable.category

            if(type == "breaks"){
                self$true_breaks[[private$break_score]]$ref <- datatable.ref
                self$true_breaks[[private$break_score]]$kmer <- datatable.kmer
                self$true_breaks[[private$break_score]]$scores <- t(datatable)
            } else if(type == "controls"){
                self$control_breaks[[private$break_score]]$ref <- datatable.ref
                self$control_breaks[[private$break_score]]$kmer <- datatable.kmer
                self$control_breaks[[private$break_score]]$scores <- t(datatable)
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
                private$ranges$short_range-2
            )

            self$true_breaks_expanded$mid_range <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref), 
                private$ranges$mid_range-2
            )

            private$datatable_seq$mid_range$breaks <- BSgenome::getSeq(
                private$ref, 
                self$true_breaks_expanded$mid_range
            )

            self$true_breaks_expanded$long_range <- plyranges::stretch(
                plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref), 
                private$ranges$long_range-2
            )

            private$datatable_seq$long_range$breaks <- BSgenome::getSeq(
                private$ref, 
                self$true_breaks_expanded$long_range
            )

            # control breakpoints
            if(!private$regression & private$get_controls){
                self$control_breaks_expanded$short_range <- plyranges::stretch(
                    plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref), 
                    private$ranges$short_range-2
                )

                self$control_breaks_expanded$mid_range <- plyranges::stretch(
                    plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref), 
                    private$ranges$mid_range-2
                )

                private$datatable_seq$mid_range$control <- BSgenome::getSeq(
                    private$ref, 
                    self$control_breaks_expanded$mid_range
                )

                self$control_breaks_expanded$long_range <- plyranges::stretch(
                    plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref), 
                    private$ranges$long_range-2
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

                private$datatable_seq$long_range$control <- BSgenome::getSeq(
                    private$ref, 
                    self$control_breaks_expanded$long_range
                )
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        }
    )
)