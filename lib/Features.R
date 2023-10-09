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

        initialize = function(k, exp, seed, scores_with_kmers, true_prop, assembly, 
                              break_score, ranges, exp_to_remove_in_table, 
                              regression, get_controls, which_chr, num_cores){
            if(!missing(exp)) private$exp <- exp
            if(!missing(exp_to_remove_in_table)){
                private$exp_to_remove_in_table <- exp_to_remove_in_table
            }
            if(!missing(regression)) private$regression <- regression
            if(!missing(get_controls)) private$get_controls <- get_controls

            # if only want to extract features for specific chromosome
            if(!missing(which_chr)) private$which_chr <- which_chr
            if(!missing(num_cores)) private$num_cores <- num_cores

            super$initialize(
                k = k, 
                exp = exp, 
                seed = seed,
                scores_with_kmers = scores_with_kmers,
                true_prop = true_prop,
                assembly = assembly,
                break_score = break_score
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
        #' @param only_breaks Boolean. If TRUE, will scan only using breakage data.
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
        #' @param SAVE_OUTPUT Boolean. If TRUE, feature matrix gets backed-up as csv file.
        #' @return None.
        get_features = function(break_type = "biological", ranges = NULL,
                                FEAT_PROCESS_LIKE_IMAGE = FALSE,
                                FEAT_LONGRANGE_PREPARAMS = FALSE, only_breaks = TRUE,
                                FEAT_G4_REGEX = TRUE, g4_type = "GPQS",
                                FEAT_GC_COUNT = TRUE,
                                FEAT_KMER_COUNTS = TRUE, crash_test = FALSE, kmer_window = 3,
                                FEAT_VIENNA_RNA = TRUE, sliding_window = NULL, nuc_type = "DNA",
                                RNAfold.CALL = "/Users/paddy/opt/anaconda3/bin/RNAfold",
                                maxloopsize = 12,
                                FEAT_DNA_SHAPE = TRUE,
                                SAVE_OUTPUT = TRUE){
            start.time <- Sys.time()

            # extracts true and control breakpoints
            self$get_breaks(break_type = break_type, ranges = ranges)
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

                # remove the dataset that is present in the true breakpoint tables
                if(private$exp_to_remove_in_table != ""){
                    true.breaks.mat <- true.breaks.mat[,-which(grepl(
                        pattern = private$exp_to_remove_in_table, 
                        x = colnames(true.breaks.mat))
                    )]
                }

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
                # remove the dataset that is present in the control breakpoint tables
                if(private$exp_to_remove_in_table != ""){
                    control.breaks.mat <- control.breaks.mat[,-which(grepl(
                        pattern = private$exp_to_remove_in_table, 
                        x = colnames(control.breaks.mat))
                    )]                
                }

                if(FEAT_LONGRANGE_PREPARAMS){
                    private$kmer_window <- 8
                    private$generate_kmer_table()
                    private$get_longrange_preparams(data = "breaks", only_breaks = only_breaks)
                    if(!private$regression & private$get_controls){
                        private$get_longrange_preparams(data = "control", only_breaks = only_breaks)
                    }

                    true.breaks.mat <- cbind(
                        true.breaks.mat, 
                        private$longrange_preparams$true_breaks
                    )
                    if(!private$regression & private$get_controls){
                        control.breaks.mat <- cbind(
                            control.breaks.mat, 
                            private$longrange_preparams$control_breaks
                        ) 
                    }
                }

                if(FEAT_G4_REGEX){
                    private$find_g4_regex(data = "breaks", g4_type = g4_type)
                    if(!private$regression & private$get_controls){
                        private$find_g4_regex(data = "control", g4_type = g4_type)
                    }

                    col.ind <- c("g4seq.counts" = ncol(true.breaks.mat)+1)
                    private$column_names$G4_REGEX <- "g4seq.counts"
                    true.breaks.mat <- cbind(true.breaks.mat, private$g4_regex$true_breaks)

                    if(!private$regression & private$get_controls){
                        control.breaks.mat <- cbind(control.breaks.mat, private$g4_regex$control_breaks)
                    }
                }

                if(FEAT_GC_COUNT){
                    private$get_gc_count(data = "breaks")
                    if(!private$regression & private$get_controls) private$get_gc_count(data = "control")

                    col.ind <- c(
                        col.ind,
                        "gc.content" = ncol(true.breaks.mat)+1,
                        "gc.skew" = ncol(true.breaks.mat)+2
                    )
                    private$column_names$GC_CONTENT <- "gc.content"
                    private$column_names$GC_COUNT <- c("gc.content", "gc.skew")

                    true.breaks.mat <- cbind(
                        true.breaks.mat,
                        private$gc_content$true_breaks,
                        private$gc_skew$true_breaks,
                        private$singleton_content$true_breaks
                    )

                    if(!private$regression & private$get_controls){
                        control.breaks.mat <- cbind(
                            control.breaks.mat,
                            private$gc_content$control_breaks,
                            private$gc_skew$control_breaks,
                            private$singleton_content$control_breaks
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
                    predictor <- predictor[, predictor := ifelse(predictor, "YES", "NO")]$predictor
                }
                predictor <- as.factor(predictor)
            }
            feature.mat <- cbind(predictor, feature.mat)

            # need to save breakpoint coordinate used for this feature matrix
            if(!private$regression & private$get_controls){
                feat.mat.split <- split(feature.mat, predictor)
                feat.mat.split.no <- feat.mat.split$NO
                feat.mat.split.yes <- feat.mat.split$YES
                feat.mat.split.yes.ind <- which(complete.cases(feat.mat.split.yes))
                feat.mat.split.no.ind <- which(complete.cases(feat.mat.split.no))

                ctrl.bp.to.save <- as.data.table(
                    self$control_breaks[[private$break_score]]$ref[feat.mat.split.no.ind]
                )
                ctrl.bp.to.save[, `:=`(end = NULL, width = 2, strand = NULL)]
                
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
            } else {
                feat.mat.split.yes.ind <- which(complete.cases(feature.mat))
            }
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

            if(private$regression | !private$get_controls){
                self$feature_matrix <- feature.mat[feat.mat.split.yes.ind]
            } else {
                feat.mat.split.yes <- feat.mat.split.yes[feat.mat.split.yes.ind]
                feat.mat.split.no <- feat.mat.split.no[feat.mat.split.no.ind]
                self$feature_matrix <- rbind(feat.mat.split.yes, feat.mat.split.no)
            }

            if(is.null(private$backup_feature_matrix)){
                private$backup_feature_matrix <- self$feature_matrix
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")

            if(SAVE_OUTPUT){
                # progress message
                t1 <- Sys.time()
                cur.msg <- "Backing-up results as RData files"
                l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                cat(paste0(cur.msg, l))

                saveRDS(
                    private$column_names, 
                    file = paste0("../data/feature_matrix/", private$exp, 
                                  "_kmer-", self$k, "_seed-", self$seed, "_",
                                  ifelse(private$break_score_all == TRUE,
                                  "all", private$break_score), 
                                  "-features_COLNAMES.RData"),
                )

                write_parquet(
                    self$feature_matrix, 
                    file = paste0("../data/feature_matrix/", private$exp, 
                                  "_kmer-", self$k, "_seed-", self$seed, "_",
                                  ifelse(private$break_score_all == TRUE,
                                  "all", private$break_score), 
                                  "-features.parquet")
                )

                to.keep <- c("predictor", private$column_names$BREAKS)
                write_parquet(
                    self$feature_matrix[, ..to.keep], 
                    file = paste0("../data/feature_matrix/", private$exp, 
                                  "_kmer-", self$k, "_seed-", self$seed, "_", 
                                  ifelse(private$break_score_all == TRUE,
                                  "all", private$break_score), 
                                  "-breakscores.parquet")
                )

                total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")
            }

            final.time <- Sys.time() - start.time
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cat("Final time taken:", final.time[[1]], 
                attr(final.time, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
        },

        #' @description
        #' Import feature matrix from csv file if exists.
        #' @param custom_test_path character vector. Path to file of test file.
        #' @return None.
        get_features_from_file = function(custom_test_path = ""){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Importing feature matrix from RData file"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            file.to.import <- paste0(
                "../data/feature_matrix/", private$exp, 
                "_kmer-", self$k, "_seed-", self$seed, "_",
                ifelse(private$break_score_all, 
                       "all", private$break_score), 
                "-features.RData"
            )

            if(file.exists(file.to.import)){
                self$feature_matrix <- readRDS(file.to.import)
                private$column_names <- readRDS(
                    paste0("../data/feature_matrix/", private$exp, 
                           "_kmer-", self$k, "_seed-", self$seed, "_", 
                           ifelse(private$break_score_all, 
                                  "all", private$break_score), 
                           "-features_COLNAMES.RData")
                )
            } else {
                stop("File doesn't exist.")
            }

            if(is.null(private$backup_feature_matrix)){
                private$backup_feature_matrix <- self$feature_matrix
            }

            if(is.null(private$backup_train_matrix)){
                private$backup_train_matrix <- self$train_matrix
                private$backup_test_matrix <- self$test_matrix
            }

            if(custom_test_path != ""){
                self$custom_test_path <- custom_test_path
                private$backup_test_matrix <- self$test_matrix <- readRDS(custom_test_path)
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Select subset of columns from feature matrix. 
        #' @param cols Character Vector of columns to select:
        #' c("all", "only_breaks", "only_triplets", "pre_parameterised").
        #' @return None.
        select_columns = function(cols = "all"){
            if(is.null(private$backup_feature_matrix)){
                private$backup_feature_matrix <- self$feature_matrix
            }

            if(is.null(private$backup_train_matrix)){
                private$backup_train_matrix <- self$train_matrix
                private$backup_test_matrix <- self$test_matrix
            }

            to.filter <- switch(cols,
                "only_breaks" = private$column_names$BREAKS,
                "only_triplets" = private$column_names$KMER_COUNTS,
                "only_singleton" = private$column_names$SINGLETON,
                "only_gc_content" = private$column_names$GC_CONTENT,
                "only_gc" = private$column_names$GC_COUNT,
                "singleton_and_gc" = unlist(c(
                    private$column_names$SINGLETON,
                    private$column_names$GC_CONTENT),
                    use.names = FALSE),
                "triplets_and_gc" = unlist(c(
                    private$column_names$KMER_COUNTS,
                    private$column_names$GC_CONTENT),
                    use.names = FALSE),                    
                "pre_parameterised" = unlist(c(
                    private$column_names$BREAKS,
                    private$column_names$QM_PARAMETERS,
                    private$column_names$TFBS,
                    private$column_names$G4MAP,
                    private$column_names$LONGRANGE_PREPARAMS
                    # private$column_names$EPIGENOME_MARKS
                ), use.names = FALSE),
                "pre_parameterised_triplets_gc" = unlist(c(
                    private$column_names$BREAKS,
                    private$column_names$QM_PARAMETERS,
                    private$column_names$TFBS,
                    private$column_names$G4MAP,
                    private$column_names$LONGRANGE_PREPARAMS,
                    private$column_names$KMER_COUNTS,
                    private$column_names$GC_CONTENT,
                    private$column_names$SINGLETON
                ), use.names = FALSE),
                "breaks_with_longrange" = unlist(c(
                    private$column_names$BREAKS,
                    private$column_names$LONGRANGE_PREPARAMS),
                    use.names = FALSE),                    
                "all" = unlist(c(
                    private$column_names$BREAKS,
                    private$column_names$QM_PARAMETERS,
                    private$column_names$TFBS,
                    private$column_names$G4MAP,
                    private$column_names$EPIGENOME_MARKS,
                    private$column_names$LONGRANGE_PREPARAMS,
                    private$column_names$DNA_SHAPE,
                    private$column_names$KMER_COUNTS,
                    private$column_names$GC_CONTENT,
                    private$column_names$SINGLETON),
                    use.names = FALSE
                )
            )        

            if(!"predictor" %in% to.filter) to.filter <- c("predictor", to.filter)
            if(!length(to.filter) > 0) to.filter <- unlist(private$column_names, 
                                                           use.names = FALSE)
            
            self$feature_matrix <- private$backup_feature_matrix[, ..to.filter]  

            # only when another test data set is used
            if(self$custom_test_path != ""){
                private$backup_test_matrix <- self$test_matrix <- self$test_matrix[, ..to.filter]
            }                                                                
                                                        
            # self$train_matrix <- private$backup_train_matrix[, ..to.filter]
            # self$test_matrix <- private$backup_test_matrix[, ..to.filter]
        },

        #' @description 
        #' Downsample the size of the training set.
        #' @param proportion Numeric vector for the full data subset.
        #' @return None.
        downsample_training = function(proportion = 0.5){
            if(!nrow(self$feature_matrix) > 1){
                stop("Need to first extract the features!")
            }
            if(nrow(self$train_matrix) != nrow(private$backup_train_matrix)){
                # self$feature_matrix <- private$backup_feature_matrix
                self$train_matrix <- private$backup_train_matrix              
            }

            set.seed(self$seed)

            # Split dataset based on partition
            index_2 <- caret::createDataPartition(
                self$train_matrix$predictor, 
                p = proportion, 
                list = FALSE,
                times = 1
            )
            self$train_matrix <- self$train_matrix[index_2, ]
        },

        #' @description 
        #' Split the feature matrix into training and test sets.
        #' @param train.split Numeric vector for the training set proportion.
        #' @param val.split Numeric vector for the validation set proportion.
        #' @param split_for_validation Boolean. If TRUE, split full feature
        #' matrix into train, test and validation sets.
        #' @return None.
        train_test_split = function(train.split = 0.7, val.split = 0.6,
                                    split_for_validation = FALSE){
            if(!nrow(self$feature_matrix) > 1){
                stop("Need to first extract the features!")
            }

            # progress message
            t1 <- Sys.time()
            cur.msg <- "Generating train/test split of feature matrix"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            set.seed(self$seed)

            # Split dataset based on partition
            index_2 <- caret::createDataPartition(
                self$feature_matrix$predictor, 
                p = train.split, 
                list = FALSE,
                times = 1
            )
            self$train_matrix <- self$feature_matrix[index_2, ]

            if(split_for_validation){
                val.plus.test <- self$feature_matrix[-index_2, ]

                # Split remaining for test and validation sets
                index_2 <- caret::createDataPartition(
                    val.plus.test$predictor, 
                    p = val.split, 
                    list = FALSE,
                    times = 1
                )
                self$validation_matrix <- val.plus.test[index_2, ]
                self$test_matrix <- val.plus.test[-index_2, ]
            } else {
                if(self$custom_test_path == ""){
                    self$test_matrix <- self$feature_matrix[-index_2, ]
                }
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description 
        #' Standardise the training and test sets. 
        #' @param all Boolean to standardise all or not.
        #' @param triplets Boolean to standardise triplets or not.
        #' @return None.
        standardise = function(cols){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Standardising train and test sets"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            if(is.null(nrow(self$train_matrix))){
                stop("Need to split the feature matrix into train/test")
            }

            #' @description 
            #' Filter train/test feature matrix for standardisation.
            #' @param data Data.Table of train or test set.
            #' @param all Boolean to standardise all or not.
            #' @param triplets Boolean to standardise triplets or not.
            #' @return List of column order of unfiltered matrix, 
            #' columns to not standardise, and columns to standardise.
            filter_matrix <- function(data, cols){
                original.order <- colnames(data)
                to.filter <- switch(cols,
                    "only_breaks" = private$column_names$BREAKS,
                    "only_triplets" = private$column_names$KMER_COUNTS,
                    "only_singleton" = private$column_names$SINGLETON,
                    "only_gc_content" = private$column_names$GC_CONTENT,
                    "only_gc" = private$column_names$GC_COUNT,
                    "singleton_and_gc" = unlist(c(
                        private$column_names$SINGLETON,
                        private$column_names$GC_CONTENT),
                        use.names = FALSE),
                    "triplets_and_gc" = unlist(c(
                        private$column_names$KMER_COUNTS,
                        private$column_names$GC_CONTENT),
                        use.names = FALSE),                             
                    "pre_parameterised" = unlist(c(
                        private$column_names$BREAKS,
                        private$column_names$QM_PARAMETERS,
                        private$column_names$TFBS,
                        private$column_names$G4MAP,
                        private$column_names$LONGRANGE_PREPARAMS
                        # private$column_names$EPIGENOME_MARKS
                    ), use.names = FALSE),
                    "pre_parameterised_triplets_gc" = unlist(c(
                        private$column_names$BREAKS,
                        private$column_names$QM_PARAMETERS,
                        private$column_names$TFBS,
                        private$column_names$G4MAP,
                        private$column_names$LONGRANGE_PREPARAMS,
                        private$column_names$KMER_COUNTS,
                        private$column_names$GC_CONTENT,
                        private$column_names$SINGLETON
                    ), use.names = FALSE),
                    "breaks_with_longrange" = unlist(c(
                        private$column_names$BREAKS,
                        private$column_names$LONGRANGE_PREPARAMS),
                        use.names = FALSE),                        
                    "all" = unlist(c(
                        private$column_names$BREAKS,
                        private$column_names$QM_PARAMETERS,
                        private$column_names$TFBS,
                        private$column_names$G4MAP,
                        private$column_names$EPIGENOME_MARKS,
                        private$column_names$LONGRANGE_PREPARAMS,
                        private$column_names$DNA_SHAPE,
                        private$column_names$KMER_COUNTS,
                        private$column_names$GC_COUNT,
                        private$column_names$SINGLETON),
                        use.names = FALSE
                    )
                )
                to.filter <- original.order[!(original.order %in% to.filter)]
                if(!"predictor" %in% to.filter) to.filter <- c("predictor", to.filter)
                cols.to.norm <- data[, -..to.filter]
                return(list(original.order, to.filter, cols.to.norm))
            }

            # standardise training set
            filtered.train <- filter_matrix(data = self$train_matrix, cols = cols)
            to.remove <- filtered.train[[2]]
            output_norm <- scale(filtered.train[[3]], center = TRUE, scale = TRUE)
            output <- cbind(self$train_matrix[, ..to.remove], output_norm)
            setcolorder(output, filtered.train[[1]])

            # save standardised statistics
            self$train_matrix <- output
            self$standardise_stats$mu <- attr(output_norm, "scaled:center")
            self$standardise_stats$sigma <- attr(output_norm, "scaled:scale")

            # standardise test data based on train
            filtered.test <- filter_matrix(data = self$test_matrix, cols = cols)
            to.remove <- filtered.test[[2]]
            output_norm <- scale(
                filtered.test[[3]],
                center = self$standardise_stats$mu,
                scale = self$standardise_stats$sigma
            )
            output <- cbind(self$test_matrix[, ..to.remove], output_norm)
            setcolorder(output, filtered.test[[1]])
            self$test_matrix <- output

            # mu <- as.numeric(matrixStats::colMeans2(filtered.train[[3]], na.rm = TRUE))
            # sigma <- as.numeric(matrixStats::colSds(filtered.train[[3]], na.rm = TRUE))
            # output <- apply(filtered.train[[3]], 2, scale, center = TRUE, scale = TRUE)
            # output <- cbind(self$train_matrix[, ..to.remove], output)
            # setcolorder(output, filtered.train[[1]])

            # self$train_matrix <- output
            # self$standardise_stats$mu <- mu
            # self$standardise_stats$sigma <- sigma

            # # standardise test set based on test set stats
            # filtered.test <- filter_matrix(data = self$test_matrix, 
            #                                cols = cols)
            # to.remove <- filtered.test[[2]]
            # centering <- sweep(filtered.test[[3]], 2, self$standardise_stats$mu, 
            #                    "-", check.margin = FALSE)
            # scaling <- sweep(centering, 2, self$standardise_stats$sigma, "/", 
            #                  check.margin = FALSE)
            # output <- cbind(self$test_matrix[, ..to.remove], scaling)
            # setcolorder(output, filtered.test[[1]])
            # self$test_matrix <- output

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        }
    ),
    private = list(
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

        #' @field g4_regex List of G4 regex matrix matches.
        g4_regex = NULL,

        #' @field gc_content List of GC content.
        gc_content = NULL,

        #' @field gc_skew List of GC skew.
        gc_skew = NULL,

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

        #' @field backup_feature_matrix Data.Table of Feature Matrix for backup.
        backup_feature_matrix = NULL,

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
            cat(paste0(cur.msg, l))

            datatable.seq <- switch(data,
                "breaks" = Biostrings::getSeq(
                    private$ref, 
                    self$true_breaks_expanded$long_range
                ),
                "control" = Biostrings::getSeq(
                    private$ref, 
                    self$control_breaks_expanded$long_range
                )
            )

            # normalise by experiment
            category_col <- private$extended_preparam_table[, "category"]
            preparam.table <- apply(
                private$extended_preparam_table[, -"category"], 1, 
                scale, center = TRUE, scale = TRUE
            )
            preparam.table <- as.data.table(t(preparam.table))
            setnames(preparam.table, private$kmer_list)
            preparam.table <- cbind(category_col, preparam.table)
            preparam.vector <- as.numeric(colSums(preparam.table[, -"category"]))

            # end position of each k-mer
            stops <- (width(datatable.seq[1])-private$kmer_window+1)

            # extract kmers
            norm_zscore_function <- function(
                i, datatable_seq, kmer_window, kmer_list, preparam_vector
                ){
                kmers_extracted <- subseq(
                    datatable_seq,
                    start = i, 
                    end = i + kmer_window - 1
                )

                # get norm z-scores
                kmer_ind <- match(paste0(kmers_extracted), kmer_list)
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
                preparam_vector = preparam.vector
            )
            norm_zscore_matrix <- do.call(cbind, norm_zscore)

            # ##################################################################
            # out = dslabs::read_mnist(
            #     path = NULL,
            #     download = FALSE,
            #     destdir = tempdir(),
            #     url = "https://www2.harvardx.harvard.edu/courses/IDS_08_v2_03/",
            #     keep.files = TRUE
            # )

            # x_train = out$train$images
            # as_tibble(x_train)

            # y_labels = out$train$labels
            # as_tibble(y_labels)
            # ##################################################################

            if(data == "breaks"){
                private$sequence_image_encoding$true_breaks <- norm_zscore_matrix
            } else if(data == "control"){
                private$sequence_image_encoding$control_breaks <- norm_zscore_matrix
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },
        
        #' @description 
        #' Scan along the long-range sequence and extract rolling 
        #' k-meric enrichment values from pre-parameterised maps.
        #' @param data Character vector of "breaks" or "control".
        #' @param only_breaks Boolean. If TRUE, will scan only using breakage data.
        #' @return None.
        get_longrange_preparams = function(data, only_breaks = TRUE){
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
                "breaks" = Biostrings::getSeq(
                    private$ref, 
                    self$true_breaks_expanded$long_range
                ),
                "control" = Biostrings::getSeq(
                    private$ref, 
                    self$control_breaks_expanded$long_range
                )
            )

            # import all values from pre-parameterised table
            preparam.table <- fread(
                paste0("../data/kmertone/QueryTable/QueryTable_",
                    "kmer-", private$kmer_window, "_",
                    private$break_score, ".csv"), 
                showProgress = FALSE
            )

            if(only_breaks){
                preparam.table <- preparam.table[match(
                    gsub(
                        pattern = "_zscore$|_ratio$",
                        replacement = "",
                        x = private$column_names$BREAKS
                    ), 
                    preparam.table$category
                )]             
            }
            preparam.mat <- as.matrix(preparam.table[, -"category"])
            # require finite numeric values for matrix multiplications
            preparam.mat[is.na(preparam.mat)] <- 0
            preparam.mat[is.infinite(preparam.mat)] <- 0
            preparam.mat.transpose <- t(preparam.mat)

            #' 1. get rolling k-mers in sliding window size of 1 bp
            # for some reason, the oligonucleotideFrequency function
            # causes a segmentation fault if I pass it too many 
            # sequences. To be safe, I do it in chunks
            datatable.seq.split <- split(
                datatable.seq, 
                ceiling(seq_along(datatable.seq)/10000)
            )

            fwd.ind <- match(private$kmer_ref$kmer, private$kmer_list)
            rev.ind <- match(private$kmer_ref$rev.comp, private$kmer_list)
            first.lex.kmer <- 1:nrow(private$kmer_ref)
            
            kmer.counts <- pbapply::pblapply(1:length(datatable.seq.split), function(x){
                kmer.counts <- Biostrings::oligonucleotideFrequency(
                    x = datatable.seq.split[[x]], 
                    width = private$kmer_window,
                    step = 1
                )

                # count occurrence on minus strands
                kmer.counts <- kmer.counts[, fwd.ind]+kmer.counts[, rev.ind]

                # only keep first lexicologically occurring k-mer            
                kmer.counts <- kmer.counts[, first.lex.kmer]

                #' 2. match k-mers with the pre-parameterised table
                # The matrix multiplication is the # of times a given 
                # k-mer per sample came up, and multiplied by the 
                # corresponding k-mer in the preparam table per sample.
                # Thus, the final value is simply a summed k-meric enrichment z-score.
                # output format: samples (rows) vs. breakage source (cols)
                # kmer.enrich.all <- kmer.counts %*% preparam.mat.transpose
                kmer.enrich.all <- eigenMatrixMultiply(kmer.counts, preparam.mat.transpose)
                return(kmer.enrich.all)
            })
            kmer.counts.all <- do.call(rbind, kmer.counts)
            colnames(kmer.counts.all) <- paste0(preparam.table$category, "_LR_sum")

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
            # cat("DONE! --", signif(total.time[[1]], 2), 
            #     attr(total.time, "units"), "\n")
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
                "breaks" = Biostrings::getSeq(
                    private$ref, 
                    self$true_breaks_expanded$mid_range
                ),
                "control" = Biostrings::getSeq(
                    private$ref, 
                    self$control_breaks_expanded$mid_range
                )
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
        #' Calculate GC content and GC skew.
        #' @param data Character vector of "breaks" or "control".
        #' @return None.
        get_gc_count = function(data){
            # progress message
            t1 <- Sys.time()
            cur.msg <- ifelse(
                data == "breaks",
                "Calculating GC content/skew within true break regions",
                "Calculating GC content/skew within control regions"
            )
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            datatable.seq <- switch(data,
                "breaks" = Biostrings::getSeq(private$ref, self$true_breaks_expanded$long_range),
                "control" = Biostrings::getSeq(private$ref, self$control_breaks_expanded$long_range)
            )
            genome.len <- width(datatable.seq)

            # count base frequencies
            letter.counts <- Biostrings::letterFrequency(datatable.seq, letters = "ACGT", OR = 0)
            letter.counts.norm <- letter.counts/genome.len
            gc.content <- letter.counts.norm[, "G"]+letter.counts.norm[, "C"]
            g.minus.c <- letter.counts.norm[, "G"]-letter.counts.norm[, "C"]
            gc.skew <- mapply('/', g.minus.c, gc.content)

            if(data == "breaks"){
                private$gc_content$true_breaks <- gc.content
                private$gc_skew$true_breaks <- gc.skew
                private$singleton_content$true_breaks <- letter.counts
                private$column_names$SINGLETON <- c("A", "C", "G", "T")
            } else if(data == "control"){
                private$gc_content$control_breaks <- gc.content
                private$gc_skew$control_breaks <- gc.skew
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
                    "breaks" = Biostrings::getSeq(private$ref, self$true_breaks_expanded$long_range),
                    "control" = Biostrings::getSeq(private$ref, self$control_breaks_expanded$long_range)
                )
            }

            kmer.counts <- Biostrings::oligonucleotideFrequency(datatable.seq, private$kmer_window)
            all.kmers <- colnames(kmer.counts) 
            fwd.ind <- match(private$kmer_ref$kmer, all.kmers)
            rev.ind <- match(private$kmer_ref$rev.comp, all.kmers)
            kmer.counts <- kmer.counts[, fwd.ind]+kmer.counts[, rev.ind]

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
                "breaks" = Biostrings::getSeq(private$ref, self$true_breaks_expanded$mid_range),
                "control" = Biostrings::getSeq(private$ref, self$control_breaks_expanded$mid_range)
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