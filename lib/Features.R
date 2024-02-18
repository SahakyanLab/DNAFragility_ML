Features <- R6::R6Class(
    classname = "Features",
    inherit = Breakpoints,
    public = list(
        #' @field feature_matrix Data.Table of all features.
        feature_matrix = NULL,

        #' @field train_matrix Data.Table of all features for training set.
        train_matrix = NULL,

        #' @field test_matrix Data.Table of all features for test set.
        test_matrix = NULL,

        #' @field validation_matrix Data.Table of all features for validation set.
        validation_matrix = NULL,

        #' @field standardise_stats List of stats extracted from the training data.
        standardise_stats = NULL,

        #' @field custom_test_path character vector. Path to file of test file.
        custom_test_path = "",

        #' @field true_bp_table data.table of true positions broken.
        true_bp_table = NULL,

        initialize = function(fast_matrix, fasta_sequence, label, k, seed, break_score, df_bp,
                              regression, get_controls, num_cores, any_predictors, 
                              only_breaks){
            if(!missing(regression)) private$regression <- regression
            if(!missing(get_controls)) private$get_controls <- get_controls
            if(!missing(num_cores)) private$num_cores <- num_cores
            if(!missing(any_predictors)) private$any_predictors <- any_predictors
            if(!missing(only_breaks)) private$only_breaks <- only_breaks
            if(!missing(fast_matrix)) private$fast_matrix <- fast_matrix

            super$initialize(
                fasta_sequence = fasta_sequence,
                label = label,
                k = k, 
                seed = seed,
                break_score = break_score,
                df_bp = df_bp
            )
        },

        #' @description 
        #' Extract features for use in machine learning model.
        #' @param break_type Character vector of any of c("biological", "enzymatic", "high_frequency"). 
        #'  If user knows the type of breakage but does not have any range effect values, the maximum
        #'  range effect cross similar types of breakages will be taken. If no type is given, user 
        #'  has to input 3 range effects into the ranges parameter.
        #' @param ranges Numeric vector of length 3, where position 1 is the short range effect, 
        #'  position 2 is the medium range effect, and position 3 is the long range effect.
        #' @param FEAT_LONGRANGE_PREPARAMS Boolean. If True, gets summed k-mer pre-param z-scores.
        #' @param FEAT_PROCESS_LIKE_IMAGE Boolean. If True, will extract "pixels" per rolling k-mer position.
        #' @param FEAT_G4_REGEX Boolean. If TRUE, this feature gets extracted.
        #' @param g4_type Character vector. c("GPQS", "GQS", "APQS", "AQS").
        #' @param FEAT_GC_COUNT Boolean. If TRUE, this feature gets extracted.
        #' @param FEAT_KMER_COUNTS Boolean. If TRUE, this feature gets extracted.
        #' @param kmer_window Numeric vector. K-mer size to extract frequency counts.
        #' @param crash_test Boolean. If TRUE, extracts 3-mer count in same window as break scores.
        #' @param FEAT_VIENNA_RNA Boolean. If TRUE, this feature gets extracted.
        #' @param sliding_window Numeric vector. Local sub-sequence window size to run viennaRNA on.
        #' @param nuc_type Character vector. c("DNA", "RNA").
        #' @param FEAT_RNAfold.CALL Character vector. Root directory of the viennaRNA programme.
        #' @param maxloopsize Numeric vector of maximum loop size for G4 sequences.
        #' @param DNA_SHAPE Boolean. If TRUE, this feature gets extracted.
        #' @return None.
        get_features = function(break_type = "Biological", ranges = NULL,
                                FEAT_PROCESS_LIKE_IMAGE = FALSE,
                                FEAT_LONGRANGE_PREPARAMS = FALSE,
                                FEAT_G4_REGEX = TRUE, g4_type = "GPQS",
                                FEAT_GC_COUNT = TRUE,
                                FEAT_KMER_COUNTS = TRUE, crash_test = FALSE, kmer_window = 3,
                                FEAT_VIENNA_RNA = TRUE, sliding_window = NULL, nuc_type = "DNA",
                                RNAfold.CALL = "/Users/paddy/opt/anaconda3/bin/RNAfold",
                                maxloopsize = 12,
                                FEAT_DNA_SHAPE = TRUE){
            start.time <- Sys.time()

            # extracts true and control breakpoints
            self$get_breaks(break_type = break_type, ranges = ranges)
            if(self$out == -1) return(-1)
            col.ind <- NULL

            if(FEAT_PROCESS_LIKE_IMAGE){
                private$kmer_window <- 8
                private$generate_kmer_table()
                private$process_like_image(data = "breaks")
                true.breaks.mat <- private$sequence_image_encoding$true_breaks
                control.breaks.mat <- NULL

                if(!private$regression & private$get_controls){
                    private$process_like_image(data = "control")
                    control.breaks.mat <- private$sequence_image_encoding$control_breaks
                }
            } else {
                # private vars
                private$maxloopsize <- maxloopsize

                # generate true breakpoint table
                true.breaks.mat <- cbind(
                    self$true_breaks$zscore$scores,
                    self$true_breaks$ratio$scores
                )

                # Query table contains breaks, tfbs, epigenome marks and QM params. Split it up.
                # TFBS
                tfbs.ind <- which(grepl(pattern = "^TFBS_", colnames(true.breaks.mat)))
                private$column_names$TFBS <- colnames(true.breaks.mat)[tfbs.ind]

                # QM parameters
                qm.ind <- which(grepl(pattern = "^[A-Z]_ds.dEhof_ds_ss", colnames(true.breaks.mat)))
                private$column_names$QM_PARAMETERS <- colnames(true.breaks.mat)[qm.ind]

                # G4seq map
                g4.ind <- which(grepl(pattern = "^G4seq", colnames(true.breaks.mat)))
                private$column_names$G4MAP <- colnames(true.breaks.mat)[g4.ind]

                # epigenome marks
                epigenome.ind <- which(grepl(
                    pattern = "Epigenome|ATACseq|Chipseq|Dnaseseq|FAIREseq", 
                    colnames(true.breaks.mat)
                ))
                private$column_names$EPIGENOME_MARKS <- colnames(true.breaks.mat)[epigenome.ind]            

                # remaining ones are the breakage scores
                breaks.ind <- which(!grepl(
                    pattern = paste0(
                        "^TFBS_|^[A-Z]_ds.dEhof_ds_ss|^G4seq|",
                        "Epigenome|ATACseq|Chipseq|Dnaseseq|FAIREseq"
                    ), 
                    colnames(true.breaks.mat)
                ))
                private$column_names$BREAKS <- colnames(true.breaks.mat)[breaks.ind]

                # generate control breakpoint table
                control.breaks.mat <- cbind(
                    self$control_breaks$zscore$scores,
                    self$control_breaks$ratio$scores
                )
                
                if(FEAT_LONGRANGE_PREPARAMS){
                    private$kmer_window <- 6
                    private$generate_kmer_table()
                    private$get_longrange_preparams(data = "breaks")
                    if(!private$regression & private$get_controls){
                        private$get_longrange_preparams(data = "control")
                    }

                    true.breaks.mat <- cbind(
                        true.breaks.mat, 
                        private$longrange_preparams$true_breaks
                        # private$kmer_counts$true_breaks
                    )
                    if(!private$regression & private$get_controls){
                        control.breaks.mat <- cbind(
                            control.breaks.mat, 
                            private$longrange_preparams$control_breaks
                            # private$kmer_counts$control_breaks
                        ) 
                    }
                }

                if(FEAT_G4_REGEX){
                    private$find_g4_regex(data = "breaks", g4_type = g4_type)
                    if(!private$regression & private$get_controls){
                        private$find_g4_regex(data = "control", g4_type = g4_type)
                    }

                    col.ind <- c("G4seq_counts" = ncol(true.breaks.mat)+1)
                    private$column_names$G4_REGEX <- "G4seq_counts"
                    true.breaks.mat <- cbind(true.breaks.mat, private$g4_regex$true_breaks)

                    if(!private$regression & private$get_controls){
                        control.breaks.mat <- cbind(control.breaks.mat, private$g4_regex$control_breaks)
                    }
                }

                if(FEAT_GC_COUNT){
                    private$get_gc_count(data = "breaks")
                    if(!private$regression & private$get_controls){
                        private$get_gc_count(data = "control")
                    }

                    col.ind <- c(
                        col.ind,
                        "GC_content" = ncol(true.breaks.mat)+1
                    )
                    private$column_names$GC_CONTENT <- "GC_content"
                    private$column_names$GC_COUNT <- c("GC_content")

                    true.breaks.mat <- cbind(
                        true.breaks.mat,
                        private$gc_content$true_breaks
                        # private$singleton_content$true_breaks
                    )

                    if(!private$regression & private$get_controls){
                        control.breaks.mat <- cbind(
                            control.breaks.mat,
                            private$gc_content$control_breaks
                            # private$singleton_content$control_breaks
                        )
                    }
                }

                if(FEAT_KMER_COUNTS){
                    private$kmer_window <- kmer_window
                    private$generate_kmer_table()
                    private$get_kmer_counts(data = "breaks", crash_test = crash_test)
                    if(!private$regression & private$get_controls){
                        private$get_kmer_counts(data = "control", crash_test = crash_test)
                    }

                    true.breaks.mat <- cbind(true.breaks.mat, private$kmer_counts$true_breaks)
                    if(!private$regression & private$get_controls){
                        control.breaks.mat <- cbind(control.breaks.mat, private$kmer_counts$control_breaks)
                    }
                }

                if(FEAT_VIENNA_RNA){
                    private$sliding_window <- sliding_window
                    private$run_viennaRNA_fold(
                        data = "breaks", 
                        RNAfold.CALL = RNAfold.CALL, 
                        nuc_type = nuc_type
                    )
                    if(!private$regression & private$get_controls){
                        private$run_viennaRNA_fold(
                            data = "control", 
                            RNAfold.CALL = RNAfold.CALL, 
                            nuc_type = nuc_type
                        )
                    }

                    true.breaks.mat <- cbind(true.breaks.mat, private$viennaRNA$true_breaks)
                    if(!private$regression & private$get_controls){
                        control.breaks.mat <- cbind(control.breaks.mat, private$viennaRNA$control_breaks)
                    }
                }

                if(FEAT_DNA_SHAPE){
                    private$get_dnashape(data = "breaks")
                    if(!private$regression & private$get_controls) private$get_dnashape(data = "control")

                    true.breaks.mat <- cbind(true.breaks.mat, private$dna_shape$true_breaks)
                    if(!private$regression & private$get_controls){
                        control.breaks.mat <- cbind(control.breaks.mat, private$dna_shape$control_breaks)
                    }
                }
            }

            # progress message
            t1 <- Sys.time()
            cur.msg <- "Concatenating all features into one matrix"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            if(is.null(control.breaks.mat)){
                feature.mat <- true.breaks.mat
            } else {
                feature.mat <- rbind(true.breaks.mat, control.breaks.mat)
            }
            feature.mat <- as.data.table(feature.mat)
            if(length(col.ind) > 0) colnames(feature.mat)[col.ind] <- names(col.ind)

            if(private$regression){
                predictor <- mcols(self$true_breaks[[private$break_score]]$kmer)
                predictor <- as.data.table(predictor)
                predictor <- predictor[, 1]
                setnames(predictor, "predictor")
            } else {
                if(private$any_predictors){
                    if(private$get_controls){
                        predictor <- switch(private$break_score,
                            "zscore" = c(
                                rep("YES", dim(self$true_breaks$zscore$scores)[1]), 
                                rep("NO", dim(self$control_breaks$zscore$scores)[1])
                            ),
                            "ratio" = c(
                                rep("YES", dim(self$true_breaks$ratio$scores)[1]), 
                                rep("NO", dim(self$control_breaks$ratio$scores)[1])
                            )
                        )
                    } else {
                        predictor <- mcols(self$true_breaks[[private$break_score]]$kmer)
                        predictor <- as.data.table(predictor)
                        predictor <- predictor[, 1]
                        setnames(predictor, "predictor")
                        predictor <- predictor[, predictor := ifelse(predictor == 1, "YES", "NO")]$predictor
                    }
                    predictor <- as.factor(predictor)
                }
            }

            if(!private$any_predictors) predictor <- NULL
            feature.mat <- cbind(predictor, feature.mat)

            # need to save breakpoint coordinate used for this feature matrix
            if(private$any_predictors){
                if(!private$regression){
                    if(private$get_controls){
                        feat.mat.split <- split(feature.mat, predictor)
                        feat.mat.split.no <- feat.mat.split$NO
                        feat.mat.split.yes <- feat.mat.split$YES
                        feat.mat.split.yes.ind <- which(complete.cases(feat.mat.split.yes))
                        feat.mat.split.no.ind <- which(complete.cases(feat.mat.split.no))

                        ctrl.bp.to.save <- as.data.table(
                            self$control_breaks[[private$break_score]]$ref[feat.mat.split.no.ind]
                        )
                        ctrl.bp.to.save <- data.table(
                            seqnames = ctrl.bp.to.save$seqnames,
                            start = ctrl.bp.to.save$start,
                            width = 2,
                            Breaks = "NO"
                        )

                        true.bp.to.save <- as.data.table(
                            self$true_breaks[[private$break_score]]$ref[feat.mat.split.yes.ind]
                        )
                        true.bp.to.save <- data.table(
                            seqnames = true.bp.to.save$seqnames,
                            start = true.bp.to.save$start,
                            width = 2,
                            Breaks = "YES"
                        )
                    } else {
                        feat.mat.split.yes.ind <- which(feature.mat$predictor == "YES")
                        true.bp.to.save <- as.data.table(
                            self$true_breaks[[private$break_score]]$ref[feat.mat.split.yes.ind]
                        )
                        true.bp.to.save[, `:=`(end = NULL, width = 2, strand = NULL)]
                        true.bp.to.save[, Breaks := NULL]
                        true.bp.to.save[, Breaks := "YES"]
                        
                        feat.mat.split.no.ind <- which(feature.mat$predictor == "NO")
                        ctrl.bp.to.save <- as.data.table(
                            self$true_breaks[[private$break_score]]$ref[feat.mat.split.no.ind]
                        )
                        ctrl.bp.to.save[, `:=`(end = NULL, width = 2, strand = NULL)]
                        ctrl.bp.to.save[, Breaks := NULL]
                        ctrl.bp.to.save[, Breaks := "NO"]
                    }
                    
                    dir.create(
                        path = "../data/feature_matrix/control_break_locations/",
                        showWarnings = FALSE,
                        recursive = TRUE
                    )
                    fwrite(
                        ctrl.bp.to.save,
                        file = paste0(
                            "../data/feature_matrix/control_break_locations/",
                            private$exp, "_kmer-", self$k, "_seed-", self$seed, "_",
                            ifelse(private$break_score_all, "all", private$break_score), 
                            ".csv"
                        )
                    )

                    dir.create(
                        path = "../data/feature_matrix/true_break_locations/",
                        showWarnings = FALSE,
                        recursive = TRUE
                    )
                    fwrite(
                        true.bp.to.save,
                        file = paste0(
                            "../data/feature_matrix/true_break_locations/",
                            private$exp, "_kmer-", self$k, "_seed-", self$seed, "_",
                            ifelse(private$break_score_all, "all", private$break_score), 
                            ".csv"
                        )
                    )

                    self$true_bp_table <- rbind(true.bp.to.save, ctrl.bp.to.save)
                    setorder(self$true_bp_table, start)
                } else {
                    feat.mat.split.yes.ind <- which(complete.cases(feature.mat))
                
                    true.bp.to.save <- as.data.table(
                        self$true_breaks[[private$break_score]]$ref[feat.mat.split.yes.ind]
                    )
                    true.bp.to.save[, `:=`(end = NULL, width = 2, strand = NULL)] 
                    self$true_bp_table <- true.bp.to.save

                    dir.create(
                        path = "../data/feature_matrix/true_break_locations/",
                        showWarnings = FALSE,
                        recursive = TRUE
                    )
                    fwrite(
                        true.bp.to.save,
                        file = paste0(
                            "../data/feature_matrix/true_break_locations/",
                            private$exp, "_kmer-", self$k, "_seed-", self$seed, "_",
                            ifelse(private$break_score_all, "all", private$break_score), 
                            ".csv"
                        )
                    )
                }
            }
    
            if(private$any_predictors){
                if(private$regression | !private$get_controls){
                    if(private$regression){
                        self$feature_matrix <- feature.mat[feat.mat.split.yes.ind]
                    } else {
                        self$feature_matrix <- feature.mat
                    }
                } else {
                    feat.mat.split.yes <- feat.mat.split.yes[feat.mat.split.yes.ind]
                    feat.mat.split.no <- feat.mat.split.no[feat.mat.split.no.ind]
                    self$feature_matrix <- rbind(feat.mat.split.yes, feat.mat.split.no)
                }
            } else {
                self$feature_matrix <- feature.mat
                
                # start_pos <- mcols(self$true_breaks[[private$break_score]]$ref)$start.pos
                # ind <- match(start_pos, self$df_bp$start)
                # self$true_bp_table <- self$df_bp[ind,]

                self$true_bp_table <- as_tibble(self$true_breaks[[private$break_score]]$ref) %>% 
                    dplyr::select(start, seqnames, width, start.pos) %>% 
                    as.data.table()
                self$true_bp_table[, start.pos := NULL]
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")

            final.time <- Sys.time() - start.time
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cat("Final time taken:", final.time[[1]], 
                attr(final.time, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
        }
    ),
    private = list(
        #' @field fast_matrix Boolean. If TRUE, will use RcppArmadillo for matrix calculations.
        fast_matrix = FALSE,

        #' @field sliding_window Numeric vector of the window to slide across the sequence.
        sliding_window = NULL,

        #' @field maxloopsize Numeric vector of maximum loop size for G4 sequences.
        maxloopsize = NULL,

        #' @field sequence_image_encoding Matrix of rolling k-meric enrichment value 
        #'  akin to pixel values per position.
        sequence_image_encoding = NULL,

        #' @field longrange_preparams Matrix of summed k-meric enrichment z-scores
        #'  from the pre-parameterised tables scanned in a rolling 1-bp sliding window
        #'  within the long-range regions.
        longrange_preparams = NULL,

        #' @field only_breaks Boolean. If TRUE, will scan only using breakage data.
        only_breaks = TRUE,

        #' @field g4_regex List of G4 regex matrix matches.
        g4_regex = NULL,

        #' @field gc_content List of GC content.
        gc_content = NULL,

        #' @field singleton_content List of ATGC counts.
        singleton_content = NULL,

        #' @field kmer_counts List of k-mer frequency matrix.
        kmer_counts = NULL,

        #' @field viennaRNA Matrix Array of viennaRNA calculations.
        viennaRNA = NULL,

        #' @field dna_shape Matrix Array of DNAShapeR calculations.
        dna_shape = NULL,

        #' @field exp Character vector of experiment name.
        exp = NULL,

        #' @field column_names Character vector of all column names in feature matrix.
        column_names = NULL,

        #' @field backup_train_matrix Data.Table of the training feature matrix for backup.
        backup_train_matrix = NULL,
        
        #' @field backup_train_matrix Data.Table of the testing feature matrix for backup.
        backup_test_matrix = NULL,        

        #' @field exp_to_remove_in_table Character vector of the dataset to remove which the 
        #'  model is being trained on.
        exp_to_remove_in_table = "",

        #' @field regression Boolean. If True, will perform a regression if data is for regression.
        regression = FALSE,

        #' @field get_controls Boolean. If True, will extract both true and false breakpoints 
        #'  if available for the classification model.
        get_controls = TRUE,

        #' @field num_cores Numeric vector. Indicate the CPUs to use for parallel processing.
        num_cores = 4,

        #' @field any_predictors Boolean. If TRUE, will treat without any predictor column. 
        any_predictors = TRUE,

        #' @description 
        #' Scan along the long-range sequence and extract rolling
        #' k-meric enrichment values akin to pixel values per position.
        #' @param data Character vector of "breaks" or "control".
        #' @return None.
        process_like_image = function(data){
            # progress message
            t1 <- Sys.time()
            cur.msg <- ifelse(
                data == "breaks",
                "Encoding rolling k-mer position within LR true regions",
                "Encoding rolling k-mer position within control regions"
            )
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")

            datatable.seq <- switch(data,
                "breaks" = private$datatable_seq$long_range$breaks,
                "control" = private$datatable_seq$long_range$control
            )

            # normalise by experiment
            private$get_extended_querytable(kmer_size = private$kmer_window)
            category_col <- private$extended_preparam_table[, "category"]
            preparam.table <- apply(
                private$extended_preparam_table[, -"category"], 1, 
                scale, center = TRUE, scale = TRUE
            )
            preparam.table <- as.data.table(t(preparam.table))
            setnames(preparam.table, private$kmer_list)
            preparam.table <- cbind(category_col, preparam.table)
            preparam.vector <- as.numeric(colSums(preparam.table[, -"category"])) 

            compressed_image <- TRUE
            if(compressed_image){
                # extract kmers
                stops <- (width(datatable.seq[1])-private$kmer_window+1)
                
                with_progress({
                    # `stops` is the total number of steps/iterations
                    p <- progressor(steps = stops)
                    
                    norm_zscore_function <- function(
                        i, datatable_seq, kmer_window, 
                        kmer_list, preparam_vector, p
                        ){
                        # Signal progress at the beginning of each iteration
                        p()

                        kmers_extracted <- subseq(
                            datatable_seq,
                            start = i, 
                            end = i+kmer_window-1
                        )

                        # get norm z-scores
                        kmer_ind <- match(paste0(kmers_extracted), kmer_list)

                        # for neural network: one row per sequence
                        norm_zscore <- preparam_vector[kmer_ind]

                        # replace NAs with zero
                        norm_zscore[is.na(norm_zscore)] <- 0
                        
                        # # k-mer index map
                        # kmers_extracted <- get_kmer_indices(
                        #     private$kmer_list, 
                        #     paste0(kmers_extracted) # this is the bottleneck
                        # )

                        # # get norm z-scores
                        # norm_zscore <- preparam.vector[kmers_extracted]

                        # # replace NAs with zero
                        # norm_zscore[is.na(norm_zscore)] <- 0

                        return(norm_zscore)
                    }

                    norm_zscore <- furrr::future_map(
                        .x = 1:stops, 
                        .f = norm_zscore_function, 
                        datatable_seq = datatable.seq, 
                        kmer_window = private$kmer_window, 
                        kmer_list = private$kmer_list, 
                        preparam_vector = preparam.vector,
                        p = p  # pass the progressor object to the function
                    )
                    norm_zscore_matrix <- do.call(cbind, norm_zscore)
                })

                
                # # end position of each k-mer
                # stops <- (width(datatable.seq[1])-private$kmer_window+1)
                # norm_zscore <- furrr::future_map(
                #     .x = 1:stops, 
                #     .f = norm_zscore_function, 
                #     datatable_seq = datatable.seq, 
                #     kmer_window = private$kmer_window, 
                #     kmer_list = private$kmer_list, 
                #     preparam_vector = preparam.vector
                # )
                # norm_zscore_matrix <- do.call(cbind, norm_zscore)

                if(data == "breaks"){
                    private$sequence_image_encoding$true_breaks <- norm_zscore_matrix
                } else if(data == "control"){
                    private$sequence_image_encoding$control_breaks <- norm_zscore_matrix
                }
            } else {
                preparam.table.mat <- preparam.table[, -"category"]
                preparam.table.mat <- as.matrix(preparam.table.mat)
                # dummy column for sequences with N
                preparam.table.mat <- cbind(preparam.table.mat, 0)
                colnames(preparam.table.mat) <- paste0("V", 1:ncol(preparam.table.mat))

                dir.create(
                    path = "image_feathers",
                    showWarnings = FALSE
                )
                iters <- length(datatable.seq)
                seq_len <- width(datatable.seq[1])
                stops <- (seq_len-private$kmer_window+1)

                with_progress({
                    # `stops` is the total number of steps/iterations
                    p <- progressor(steps = stops)
                    
                    get_image_per_seq <- function(
                        seq_num, i, stop, seq_len, datatable_seq,
                        kmer_window, kmer_list, preparam_vector, last_col_n, p
                        ){
                        # Signal progress at the beginning of each iteration
                        p()

                        # Extract the specified sequence
                        sequence <- datatable_seq[[seq_num]]

                        # extract rolling k-mer 
                        kmers_extracted <- substring(
                            text = paste0(sequence),
                            first = 1:stops,
                            last = kmer_window:seq_len
                        )

                        # get norm z-scores
                        kmer_ind <- match(kmers_extracted, kmer_list)

                        # replace NAs with zero
                        seq_with_n <- which(grepl(
                            pattern = "N", 
                            x = kmers_extracted,
                            ignore.case = TRUE
                        ))
                        kmer_ind[seq_with_n] <- last_col_n

                        # for neural network: one matrix per sequence
                        norm_zscore_table <- preparam_vector[, kmer_ind]
                        colnames(norm_zscore_table) <- NULL

                        # save as feather
                        # temp <- as.data.table(norm_zscore_table)
                        # setnames(temp, new_colnames)
                        # write_feather(temp, paste0("./image_feathers/image_", seq_num, ".feather"))

                        # con <- file("binary_images/all_images.bin", "wb")
                        # writeBin(as.double(as.matrix(norm_zscore_table)), con, size = 8)  # 8 bytes for double
                        # close(con)

                        return(norm_zscore_table)
                    }

                    # suppressPackageStartupMessages(suppressWarnings(library(feather)))
                    # iters = 500
                    # t1 = Sys.time()
                    norm_zscore <- furrr::future_map(
                        .x = 1:iters, 
                        .f = get_image_per_seq, 
                        stop = stops, 
                        seq_len = seq_len,
                        datatable_seq = datatable.seq, 
                        kmer_window = private$kmer_window, 
                        kmer_list = private$kmer_list, 
                        preparam_vector = preparam.table.mat,
                        last_col_n = ncol(preparam.table.mat),
                        new_colnames = paste0("V", 1:stops),
                        p = p  # pass the progressor object to the function
                    )
                    # total.time = Sys.time() - t1
                    # print(total.time)

                    # (as.numeric(total.time)*60)/iters*length(datatable.seq)/60/60
                })
                                

            }

            cat(paste0(cur.msg, l))
            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },
        
        #' @description 
        #' Scan along the long-range sequence and extract rolling 
        #' k-meric enrichment values from pre-parameterised maps.
        #' @param data Character vector of "breaks" or "control".
        #' @return None.
        get_longrange_preparams = function(data){
            # progress message
            t1 <- Sys.time()
            cur.msg <- ifelse(
                data == "breaks",
                "Getting pre-parameterised k-mer vals within LR true regions",
                "Getting pre-parameterised k-mer vals within LR control regions"
            )
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            datatable.seq <- switch(data,
                "breaks" = private$datatable_seq$long_range$breaks,
                "control" = private$datatable_seq$long_range$control
            )

            # normalise by experiment
            private$get_extended_querytable(kmer_size = private$kmer_window)
            category_col <- private$extended_preparam_table[, "category"]
            preparam.table <- apply(
                private$extended_preparam_table[, -"category"], 1, 
                scale, center = TRUE, scale = TRUE
            )
            
            if(private$only_breaks){
                which_rows <- match(
                    gsub(
                        pattern = "_zscore$|_ratio$",
                        replacement = "",
                        x = private$column_names$BREAKS
                    ), 
                    category_col$category
                )

                preparam.table <- preparam.table[, which_rows]       
                category_col <- category_col[which_rows]
            }
            preparam.mat <- as.matrix(preparam.table)
            
            # require finite numeric values for matrix multiplications
            preparam.mat[is.na(preparam.mat)] <- 0
            preparam.mat[is.infinite(preparam.mat)] <- 0

            #' 1. get rolling k-mers in sliding window size of 1 bp
            # for some reason, the oligonucleotideFrequency function
            # causes a segmentation fault if I pass it too many 
            # sequences. To be safe, I do it in chunks
            datatable.seq.split <- split(
                datatable.seq, 
                ceiling(seq_along(datatable.seq)/500000)
            )

            # only keep first lexicologically occurring k-mer
            fwd.ind <- match(private$kmer_ref$kmer, private$kmer_list)
            rev.ind <- match(private$kmer_ref$rev.comp, private$kmer_list)
            first.lex.kmer <- 1:nrow(private$kmer_ref)
            preparam.mat <- preparam.mat[first.lex.kmer, ]

            kmer.counts <- lapply(1:length(datatable.seq.split), function(x){
                kmer.counts <- Biostrings::oligonucleotideFrequency(
                    x = datatable.seq.split[[x]], 
                    width = private$kmer_window,
                    step = 1
                )

                # count occurrence on minus strands
                kmer.counts <- 
                    kmer.counts[, fwd.ind, drop = FALSE] + 
                    kmer.counts[, rev.ind, drop = FALSE]

                # only keep first lexicologically occurring k-mer            
                kmer.counts <- kmer.counts[, first.lex.kmer, drop = FALSE]

                #'Match k-mers with the pre-parameterised table
                # The matrix multiplication is the # of times a given 
                # k-mer per sample came up, and multiplied by the 
                # corresponding k-mer in the preparam table per sample.
                # Thus, the final value is simply a summed k-meric enrichment z-score.
                # output format: samples (rows) vs. breakage source (cols)
                if(private$fast_matrix){
                    kmer.enrich.all <- eigenMatrixMultiply(kmer.counts, preparam.mat)
                } else {
                    kmer.enrich.all <- kmer.counts %*% preparam.mat
                }

                return(kmer.enrich.all)
            })

            # long-range breakage z-scores
            kmer.counts.all <- do.call(rbind, kmer.counts)
            colnames(kmer.counts.all) <- paste0(category_col$category, "_LR_sum")

            # filter for non-zero k-mer occurrences
            kmer.counts.all[is.na(kmer.counts.all)] <- 0
            kmer.counts.all[is.infinite(kmer.counts.all)] <- 0

            if(data == "breaks"){
                private$longrange_preparams$true_breaks <- kmer.counts.all
                private$column_names$LONGRANGE_PREPARAMS <- colnames(kmer.counts.all)
            } else if(data == "control"){
                private$longrange_preparams$control_breaks <- kmer.counts.all
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },        

        #' @description 
        #' Find G4 regular expression in supplied DNA sequence.
        #' @param data Character vector of "breaks" or "control".
        #' @param g4_type Character vector of "GPQS", "GQS", "APQS" or "AQS".
        #' @return None.
        find_g4_regex = function(data, g4_type){
            # progress message
            t1 <- Sys.time()
            cur.msg <- ifelse(
                data == "breaks",
                "Finding G4s with regex within true break regions",
                "Finding G4s with regex within control regions"
            )
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            datatable.seq <- switch(data,
                "breaks" = private$datatable_seq$mid_range$breaks,
                "control" = private$datatable_seq$mid_range$control
            )
            datatable.seq <- paste0(datatable.seq)

            if(g4_type == "GPQS" | g4_type == "GQS"){
                plus.strand.regex  <- paste0("([gG]{3,}[NATGCnatgc]{1,", 
                                             private$maxloopsize,"}){3,}[gG]{3,}")
                minus.strand.regex <- paste0("([cC]{3,}[NATGCnatgc]{1,", 
                                             private$maxloopsize,"}){3,}[cC]{3,}")
            } else if(g4_type == "APQS" | g4_type == "AQS"){
                plus.strand.regex  <- paste0("([aA]{3,}[NATGCnatgc]{1,", 
                                             private$maxloopsize,"}){3,}[aA]{3,}")
                minus.strand.regex <- paste0("([tT]{3,}[NATGCnatgc]{1,", 
                                             private$maxloopsize,"}){3,}[tT]{3,}")
            } else if(g4_type != "GPQS" & g4_type != "APQS" & 
                      g4_type != "GQS"  & g4_type != "AQS"  & g4_type != "other"){
                stop("Non-recognisable g4_type input.")
            }

            plus.strand <- stringr::str_count(
                string = datatable.seq, 
                pattern = plus.strand.regex
            )
            minus.strand <- stringr::str_count(
                string = datatable.seq, 
                pattern = minus.strand.regex
            )

            if(data == "breaks"){
                private$g4_regex$true_breaks <- plus.strand+minus.strand
            } else if(data == "control"){
                private$g4_regex$control_breaks <- plus.strand+minus.strand
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description 
        #' Calculate GC content.
        #' @param data Character vector of "breaks" or "control".
        #' @return None.
        get_gc_count = function(data){
            # progress message
            t1 <- Sys.time()
            cur.msg <- ifelse(
                data == "breaks",
                "Calculating GC content within true break regions",
                "Calculating GC content within control regions"
            )
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            datatable.seq <- switch(data,
                "breaks" = private$datatable_seq$long_range$breaks,
                "control" = private$datatable_seq$long_range$control
            )
            genome.len <- width(datatable.seq)

            # count base frequencies
            letter.counts <- Biostrings::letterFrequency(
                datatable.seq, letters = "ACGT", OR = 0
            )
            letter.counts.norm <- letter.counts/genome.len
            gc.content <- letter.counts.norm[, "G"]+letter.counts.norm[, "C"]
            g.minus.c <- letter.counts.norm[, "G"]-letter.counts.norm[, "C"]

            if(data == "breaks"){
                private$gc_content$true_breaks <- gc.content
                private$singleton_content$true_breaks <- letter.counts
                private$column_names$SINGLETON <- c("A", "C", "G", "T")
            } else if(data == "control"){
                private$gc_content$control_breaks <- gc.content
                private$singleton_content$control_breaks <- letter.counts
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description 
        #' Generate frequency table of k-mer occurrences.
        #' @param data Character vector of "breaks" or "control".
        #' @return None.
        get_kmer_counts = function(data, crash_test = FALSE){
            # progress message
            t1 <- Sys.time()
            cur.msg <- ifelse(
                data == "breaks",
                "Calculating k-mer counts within true break regions",
                "Calculating k-mer counts within control regions"
            )
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            if(crash_test){
                expand.breaks.into.kmers <- switch(data,
                    "breaks" = plyranges::stretch(
                        plyranges::anchor_center(self$true_breaks[[private$break_score]]$ref), 
                        self$k
                    ),
                    "control" = plyranges::stretch(
                        plyranges::anchor_center(self$control_breaks[[private$break_score]]$ref), 
                        self$k
                    )
                )
                datatable.seq <- Biostrings::getSeq(private$ref, expand.breaks.into.kmers)
            } else {
                datatable.seq <- switch(data,
                    "breaks" = private$datatable_seq$long_range$breaks,
                    "control" = private$datatable_seq$long_range$control
                )
            }

            kmer.counts <- Biostrings::oligonucleotideFrequency(
                x = datatable.seq,
                width = private$kmer_window,
                step = 1
            )
            all.kmers <- colnames(kmer.counts)
            fwd.ind <- match(private$kmer_ref$kmer, all.kmers)
            rev.ind <- match(private$kmer_ref$rev.comp, all.kmers)
            kmer.counts <- 
                kmer.counts[, fwd.ind, drop = FALSE] +
                kmer.counts[, rev.ind, drop = FALSE]

            if(data == "breaks"){
                private$kmer_counts$true_breaks <- kmer.counts
                private$column_names$KMER_COUNTS <- colnames(kmer.counts)
            } else {
                private$kmer_counts$control_breaks <- kmer.counts
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Sliding window of viennaRNA folding energies.
        #' @param data Character vector of "breaks" or "control".
        #' @param nuc_type Character vector of "DNA" or "RNA".
        #' @param RNAfold.CALL Character vector of root directory of the viennaRNA programme.
        #' @return None.
        run_viennaRNA_fold = function(data, nuc_type, RNAfold.CALL){
            # progress message
            t1 <- Sys.time()
            cur.msg <- ifelse(
                data == "breaks",
                "Predicting folding energies within true break regions",
                "Predicting folding energies within control regions"
            )
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            datatable.seq <- switch(data,
                "breaks" = private$datatable_seq$mid_range$breaks,
                "control" = private$datatable_seq$mid_range$control
            )

            #' @description
            #' RNA parameters are described in
            #' Mathews DH, Disney MD, Childs JL, Schroeder SJ, Zuker M, Turner DH. (2004).
            #' Incorporating chemical modification constraints into a dynamic programming algorithm
            #' for prediction of RNA secondary structure. Proc Natl Acad Sci USA 101(19):7287-92.
            #' Ts will be treated as Us, but not converted into Us in the output.
            #' @param nuc_type Character vector of "DNA" or "RNA". 
            #' @param RNAfold.CALL Character vector of root directory of the viennaRNA programme.
            #' @return Matrix Array of three features per expanded k-mer around the breakpoint.
            run_vrna <- function(nuc_type, RNAfold.CALL){
                run.command <- switch(nuc_type,
                    "RNA" = paste0(RNAfold.CALL, 
                                " --MEA -p0 -d2 --noLP --noPS < temp_seq.db > temp_rnafold.txt"),
                    "DNA" = paste0(RNAfold.CALL, " --noconv -P ../data/parameters/dna_mathews.par",
                                " --MEA -p0 -d2 --noLP --noPS < temp_seq.db > temp_rnafold.txt")
                )
                system(run.command)

                # read results and format it
                temp.out <- fread("temp_rnafold.txt", 
                                  fill = TRUE, 
                                  header = FALSE,
                                  showProgress = FALSE)

                # minimum free energy
                interval <- seq(3, dim(temp.out)[1], 5)
                viennaRNA.MFE <- gsub("^\\S* |[][()]", "", temp.out$V1[interval])
                viennaRNA.MFE <- as.numeric(viennaRNA.MFE)

                # free energy of ensemble
                interval <- seq(4, dim(temp.out)[1], 5)
                viennaRNA.FEE <- gsub("^.*=+ |[kcal/mol]", "", temp.out$V1[interval])
                viennaRNA.FEE <- as.numeric(viennaRNA.FEE)

                # frequency of mfe structure in ensemble
                interval <- seq(5, dim(temp.out)[1], 5)
                viennaRNA.FRQ <- gsub("frequency of mfe structure in ensemble|;", "", 
                                      temp.out$V1[interval])
                viennaRNA.FRQ <- as.numeric(viennaRNA.FRQ)

                # combine results
                viennaRNA.mat <- cbind(viennaRNA.MFE, viennaRNA.FEE, viennaRNA.FRQ)
                global.cols <<- c("viennaRNA.MFE", "viennaRNA.FEE", "viennaRNA.FRQ")
                colnames(viennaRNA.mat) <- global.cols
                return(viennaRNA.mat)
            }

            if(!is.null(private$sliding_window)){
                seq.len <- width(datatable.seq)[1]
                interval.nums <- ceiling(seq.len/private$sliding_window)-1

                start.ind <- 1
                end.ind <- start.ind+private$sliding_window
                for(ind in 1:interval.nums){
                    start.ind <- c(start.ind, end.ind[ind]+1)
                    if(ind < interval.nums){
                        if(seq.len <= start.ind[ind+1]+private$sliding_window){
                            end.ind <- c(end.ind, seq.len)
                            break
                        }
                        end.ind <- c(end.ind, start.ind[ind+1]+private$sliding_window)
                    } else {
                        end.ind <- c(end.ind, seq.len)
                    }
                }

                output <- pbapply::pblapply(1:length(start.ind), function(ind){
                    datatable.seq.temp <- Biostrings::DNAStringSet(stringr::str_sub(
                        string = datatable.seq, 
                        start = start.ind[ind], 
                        end = end.ind[ind]
                    ))
                    names(datatable.seq) <- paste0("str_", 1:length(datatable.seq))
                    Biostrings::writeXStringSet(datatable.seq.temp, filepath = "temp_seq.db")
                    v.mat <- run_vrna(nuc_type = nuc_type, RNAfold.CALL = RNAfold.CALL)
                    fwrite(as.data.table(v.mat), file = paste0("temp_viennaRNA_mat_", ind, ".csv"))
                    invisible(file.remove(c("temp_seq.db", "temp_rnafold.txt")))
                    return(NULL)
                })

                to.import <- list.files(path = "./", pattern = "temp_viennaRNA_mat_")
                viennaRNA.mat <- lapply(global.cols, function(cols){
                    out <- lapply(1:length(to.import), function(file){
                        datatable <- fread(to.import[file], select = cols)
                        colnames(datatable) <- NULL
                        datatable <- as.matrix(datatable)
                        return(datatable)
                    })
                    out <- do.call(cbind, out)
                    min.val <- matrixStats::rowMins(as.matrix(out), na.rm = TRUE)
                    max.val <- matrixStats::rowMaxs(as.matrix(out), na.rm = TRUE)
                    out <- matrix(data = c(min.val, max.val), ncol = 2)
                    colnames(out) <- c(paste0(cols, ".min"), paste0(cols, ".max"))
                    return(out)
                })
                viennaRNA.mat <- do.call(cbind, viennaRNA.mat)
            }

            # run viennaRNA on full DNA sequence
            names(datatable.seq) <- paste0("str_", 1:length(datatable.seq))
            Biostrings::writeXStringSet(datatable.seq, filepath = "temp_seq.db")
            v.mat <- run_vrna(nuc_type = nuc_type, RNAfold.CALL = RNAfold.CALL)
            colnames(v.mat) <- paste0(colnames(v.mat), ".full")

            if(is.null(private$sliding_window)){
                if(data == "breaks"){
                    private$viennaRNA$true_breaks <- v.mat
                } else if(data == "control"){
                    private$viennaRNA$control_breaks <- v.mat
                }
            } else {
                if(data == "breaks"){
                    private$viennaRNA$true_breaks <- cbind(v.mat, viennaRNA.mat)
                } else if(data == "control"){
                    private$viennaRNA$control_breaks <- cbind(v.mat, viennaRNA.mat)
                }
            }

            to.remove <- list.files(path = "./", pattern = "temp*")
            private$column_names$VIENNA_RNA <- colnames(private$viennaRNA$true_breaks)
            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
            invisible(file.remove(to.remove))
        },

        #' @description 
        #' Predicts DNA shape based on sliding window k-mer sequences.
        #' @param data Character vector of "breaks" or "control".
        #' @return None.
        get_dnashape = function(data){
            # progress message
            t1 <- Sys.time()
            cur.msg <- ifelse(
                data == "breaks",
                "Predicting DNA shape within true break regions",
                "Predicting DNA shape within control regions"
            )
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))            

            datatable.seq <- switch(data,
                "breaks" = self$true_breaks_expanded$mid_range,
                "control" = self$control_breaks_expanded$mid_range
            )

            # predict DNA shape
            suppressMessages(DNAshapeR::getFasta(
                datatable.seq, 
                private$ref, 
                width = as.numeric(width(datatable.seq)[1]),
                filename = "temp.fa"))
            pred <- suppressMessages(DNAshapeR::getShape("temp.fa"))

            # compute summary statistics
            shape.name <- c("HelT", "MGW", "ProT", "Roll")
            shape.pred <- lapply(1:length(shape.name), function(x){
                row.mean <- matrixStats::rowMeans2(pred[[x]], na.rm = TRUE)
                row.sd <- matrixStats::rowSds(pred[[x]], na.rm = TRUE)
                combined <- cbind(row.mean, row.sd)
                colnames(combined) <- c(
                    paste0(shape.name[x], ".mean"), 
                    paste0(shape.name[x], ".sd")
                )
                return(combined)
            })

            if(data == "breaks"){
                private$dna_shape$true_breaks <- do.call(cbind, shape.pred)
                private$column_names$DNA_SHAPE <- colnames(private$dna_shape$true_breaks)
            } else if(data == "control"){
                private$dna_shape$control_breaks <- do.call(cbind, shape.pred)
            }
            
            files.to.remove <- list.files(
                path = "./", 
                pattern = "temp", 
                full.names = TRUE
            )
            invisible(file.remove(files.to.remove))

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")            
        }
    )
)