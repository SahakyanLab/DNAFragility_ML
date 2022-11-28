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

        initialize = function(k, exp, seed, assembly, 
                              break_score, scores_with_kmers){
            if(!missing(exp)) private$exp <- exp
            super$initialize(
                k = k, 
                exp = exp, 
                seed = seed,
                assembly = assembly
            )

            # extracts true and control breakpoints
            self$get_breaks(
                break_score = break_score,
                scores_with_kmers = scores_with_kmers
            )
        },

        #' @description 
        #' Extract features for use in machine learning model.
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
        #' @param DNA_SHAPE Boolean. If TRUE, this feature gets extracted.
        #' @param FEAT_TFBS_EUCLIDEAN_DISTANCE Boolean. If TRUE, this feature gets extracted.
        #' @param FEAT_OCCUPANCY_SCORES Boolean. If TRUE, this feature gets extracted.
        #' @param FEAT_QM_PARAMETERS Boolean. If TRUE, this feature gets extracted.
        #' @param SAVE_OUTPUT Boolean. If TRUE, feature matrix gets backed-up as csv file.
        #' @return None.
        get_features = function(FEAT_G4_REGEX = TRUE, g4_type = "GPQS",
                                FEAT_GC_COUNT = TRUE,
                                FEAT_KMER_COUNTS = TRUE, crash_test = FALSE, kmer_window = 3,
                                FEAT_VIENNA_RNA = TRUE, sliding_window = NULL, nuc_type = "DNA",
                                RNAfold.CALL = "/Users/paddy/opt/anaconda3/bin/RNAfold",
                                FEAT_DNA_SHAPE = TRUE,
                                FEAT_TFBS_EUCLIDEAN_DISTANCE = FALSE,
                                FEAT_OCCUPANCY_SCORES = TRUE,
                                FEAT_QM_PARAMETERS = TRUE,
                                SAVE_OUTPUT = TRUE){
            start.time <- Sys.time()

            true.breaks.mat <- cbind(
                self$true_breaks$zscore$scores,
                self$true_breaks$ratio$scores
            )                                        
            private$column_names$BREAKS <- colnames(true.breaks.mat)
            control.breaks.mat <- cbind(
                self$control_breaks$zscore$scores,
                self$control_breaks$ratio$scores
            )
            col.ind <- NULL

            if(FEAT_G4_REGEX){
                private$find_g4_regex(data = "breaks", g4_type = g4_type)
                private$find_g4_regex(data = "control", g4_type = g4_type)

                col.ind <- c("g4seq.counts" = ncol(true.breaks.mat)+1)
                private$column_names$G4_REGEX <- "g4seq.counts"
                true.breaks.mat <- cbind(true.breaks.mat, private$g4_regex$true_breaks)
                control.breaks.mat <- cbind(control.breaks.mat, private$g4_regex$control_breaks)
            }

            if(FEAT_GC_COUNT){
                private$get_gc_count(data = "breaks")
                private$get_gc_count(data = "control")

                col.ind <- c(
                    col.ind,
                    "gc.content" = ncol(true.breaks.mat)+1,
                    "gc.skew" = ncol(true.breaks.mat)+2
                )
                private$column_names$GC_COUNT <- c("gc.content", "gc.skew")

                true.breaks.mat <- cbind(
                    true.breaks.mat,
                    private$gc_content$true_breaks,
                    private$gc_skew$true_breaks,
                    private$singleton_content$true_breaks
                )

                control.breaks.mat <- cbind(
                    control.breaks.mat,
                    private$gc_content$control_breaks,
                    private$gc_skew$control_breaks,
                    private$singleton_content$control_breaks
                )
            }

            if(FEAT_KMER_COUNTS){
                private$kmer_window <- kmer_window
                private$generate_kmer_table()
                private$get_kmer_counts(data = "breaks", crash_test = crash_test)
                private$get_kmer_counts(data = "control", crash_test = crash_test)

                true.breaks.mat <- cbind(true.breaks.mat, private$kmer_counts$true_breaks)
                control.breaks.mat <- cbind(control.breaks.mat, private$kmer_counts$control_breaks)
            }

            if(FEAT_VIENNA_RNA){
                private$sliding_window <- sliding_window
                private$run_viennaRNA_fold(data = "breaks", RNAfold.CALL = RNAfold.CALL, 
                                           nuc_type = nuc_type)
                private$run_viennaRNA_fold(data = "control", RNAfold.CALL = RNAfold.CALL, 
                                           nuc_type = nuc_type)

                true.breaks.mat <- cbind(true.breaks.mat, private$viennaRNA$true_breaks)
                control.breaks.mat <- cbind(control.breaks.mat, private$viennaRNA$control_breaks)
            }

            if(FEAT_DNA_SHAPE){
                private$get_dnashape(data = "breaks")
                private$get_dnashape(data = "control")

                true.breaks.mat <- cbind(true.breaks.mat, private$dna_shape$true_breaks)
                control.breaks.mat <- cbind(control.breaks.mat, private$dna_shape$control_breaks)
            }

            if(FEAT_TFBS_EUCLIDEAN_DISTANCE){
                private$get_tfbs_euclidean_distance(data = "breaks")
                private$get_tfbs_euclidean_distance(data = "control")

                true.breaks.mat <- cbind(
                    true.breaks.mat, 
                    private$tfbs_euclidean_distance$true_breaks
                )
                control.breaks.mat <- cbind(
                    control.breaks.mat, 
                    private$tfbs_euclidean_distance$control_breaks
                )
            }

            if(FEAT_OCCUPANCY_SCORES){
                find.upper.lim <- c(2, 4, 6, self$k)
                upper.limit <- max(find.upper.lim[find.upper.lim <= 6])

                private$get_occupancy_scores(data = "breaks", k = upper.limit)
                private$get_occupancy_scores(data = "control", k = upper.limit)

                true.breaks.mat <- cbind(
                    true.breaks.mat, 
                    private$occupancy_scores$true_breaks
                )
                control.breaks.mat <- cbind(
                    control.breaks.mat, 
                    private$occupancy_scores$control_breaks
                )
            }

            if(FEAT_QM_PARAMETERS){
                private$generate_kmer_table(k = 7)
                private$get_qm_parameters(data = "breaks")
                private$get_qm_parameters(data = "control")

                true.breaks.mat <- cbind(
                    true.breaks.mat, 
                    private$qm_parameters$true_breaks
                )
                control.breaks.mat <- cbind(
                    control.breaks.mat, 
                    private$qm_parameters$control_breaks
                )
            }

            # progress message
            t1 <- Sys.time()
            cur.msg <- "Concatenating all features into one matrix"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            feature.mat <- rbind(true.breaks.mat, control.breaks.mat)
            feature.mat <- as.data.table(feature.mat)
            if(length(col.ind) > 0) colnames(feature.mat)[col.ind] <- names(col.ind)
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
            predictor <- as.factor(predictor)
            feature.mat <- cbind(predictor, feature.mat)
            self$feature_matrix <- feature.mat[complete.cases(feature.mat)]
            if(is.null(private$backup_feature_matrix)){
                private$backup_feature_matrix <- self$feature_matrix
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")

            if(SAVE_OUTPUT){
                # progress message
                t1 <- Sys.time()
                cur.msg <- "Backing-up results as csv files"
                l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                cat(paste0(cur.msg, l))

                saveRDS(
                    private$column_names, 
                    file = paste0("../data/feature_matrix/", private$exp, 
                                  "_kmer-", self$k, "_", 
                                  ifelse(private$break_score_all == TRUE,
                                  "all", private$break_score), 
                                  "-features_COLNAMES.RData"),
                )

                fwrite(
                    x = self$feature_matrix, 
                    file = paste0("../data/feature_matrix/", private$exp, 
                                  "_kmer-", self$k, "_", 
                                  ifelse(private$break_score_all == TRUE,
                                  "all", private$break_score), 
                                  "-features.csv"),
                    showProgress = FALSE
                )

                to.keep <- c("predictor", private$column_names$BREAKS)
                fwrite(
                    x = self$feature_matrix[, ..to.keep], 
                    file = paste0("../data/feature_matrix/", private$exp, 
                                  "_kmer-", self$k, "_", 
                                  ifelse(private$break_score_all == TRUE,
                                  "all", private$break_score), 
                                  "-breakscores.csv"),
                    showProgress = FALSE
                )

                total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")
            }

            final.time <- Sys.time() - start.time
            cat("Final time taken:", t[[1]], attr(t, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
        },

        #' @description
        #' Import feature matrix from csv file if exists.
        #' @return None.
        get_features_from_csv = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Importing feature matrix from csv"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            file.to.import <- paste0(
                "../data/feature_matrix/", private$exp, 
                "_kmer-", self$k, "_", 
                ifelse(private$break_score_all, 
                       "all", private$break_score), 
                "-features.csv"
            )

            if(file.exists(file.to.import)){
                self$feature_matrix <- fread(
                    file.to.import,
                    showProgress = FALSE
                )

                private$column_names <- readRDS(
                    paste0("../data/feature_matrix/", private$exp, 
                           "_kmer-", self$k, "_", 
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

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Select subset of columns from feature matrix. 
        #' @param cols Character Vector of columns to select:
        #' c("all", "only_breaks", "only_triplets", "complex").
        #' @return None.
        select_columns = function(cols = "all"){
            if(is.null(private$backup_feature_matrix)){
                private$backup_feature_matrix <- self$feature_matrix
            }
            
            to.filter <- switch(cols,
                "only_breaks" = private$column_names$BREAKS,
                "only_triplets" = private$column_names$KMER_COUNTS,
                "only_singleton" = private$column_names$SINGLETON,
                "singleton_and_gc" = unlist(c(
                    private$column_names$SINGLETON,
                    private$column_names$GC_COUNT),
                    use.names = FALSE),
                "complex" = unlist(c(
                    private$column_names$BREAKS,
                    private$column_names$VIENNA_RNA,
                    private$column_names$DNA_SHAPE,
                    private$column_names$OCCUPANCY_SCORES,
                    private$column_names$QM_PARAMETERS),
                    use.names = FALSE),
                "all" = unlist(private$column_names, use.names = FALSE)
            )

            if(!"predictor" %in% to.filter) to.filter <- c("predictor", to.filter)
            if(!length(to.filter) > 0) to.filter <- unlist(private$column_names, 
                                                           use.names = FALSE)
            self$feature_matrix <- private$backup_feature_matrix[, ..to.filter]
        },

        #' @description 
        #' Split the feature matrix into training and test sets.
        #' @param proportion Numeric vector for the full data subset.
        #' @return None.
        setup_crash_test = function(proportion = 0.5){
            if(!nrow(self$feature_matrix) > 1){
                stop("Need to first extract the features!")
            }
            if(nrow(self$feature_matrix) != nrow(private$backup_feature_matrix)){
                self$feature_matrix <- private$backup_feature_matrix
            }

            set.seed(self$seed)

            # Split dataset based on partition
            index_2 <- caret::createDataPartition(
                self$feature_matrix$predictor, 
                p = proportion, 
                list = FALSE,
                times = 1
            )
            self$feature_matrix <- self$feature_matrix[index_2, ]
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
                self$test_matrix <- self$feature_matrix[-index_2, ]
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
                    "complex" = unlist(c(
                        private$column_names$BREAKS,
                        private$column_names$VIENNA_RNA,
                        private$column_names$DNA_SHAPE,
                        private$column_names$OCCUPANCY_SCORES,
                        private$column_names$QM_PARAMETERS),
                        use.names = FALSE),
                    "all" = unlist(private$column_names, use.names = FALSE)
                )
                to.filter <- original.order[!(original.order %in% to.filter)]
                if(!"predictor" %in% to.filter) to.filter <- c("predictor", to.filter)
                cols.to.norm <- data[, -..to.filter]
                return(list(original.order, to.filter, cols.to.norm))
            }

            # standardise training set
            filtered.train <- filter_matrix(data = self$train_matrix, 
                                            cols = cols)
            to.remove <- filtered.train[[2]]
            mu <- as.numeric(colMeans(filtered.train[[3]], na.rm = TRUE))
            sigma <- as.numeric(apply(filtered.train[[3]], 2, sd, na.rm = TRUE))
            output <- apply(filtered.train[[3]], 2, scale, center = TRUE, scale = TRUE)
            output <- cbind(self$train_matrix[, ..to.remove], output)
            setcolorder(output, filtered.train[[1]])

            self$train_matrix <- output
            self$standardise_stats$mu <- mu
            self$standardise_stats$sigma <- sigma

            # standardise test set based on test set stats
            filtered.test <- filter_matrix(data = self$test_matrix, 
                                           cols = cols)
            to.remove <- filtered.test[[2]]
            centering <- sweep(filtered.test[[3]], 2, self$standardise_stats$mu, 
                               "-", check.margin = FALSE)
            scaling <- sweep(centering, 2, self$standardise_stats$sigma, "/", 
                             check.margin = FALSE)
            output <- cbind(self$test_matrix[, ..to.remove], scaling)
            setcolorder(output, filtered.test[[1]])
            self$test_matrix <- output

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        }
    ),
    private = list(
        #' @field kmer_window Numeric vector of the window size for counting k-mers.
        kmer_window = NULL,

        #' @field kmer_ref Data.Table of forward and reverse complement k-mers.
        kmer_ref = NULL,

        #' @field kmer_ref Data.Table of forward and reverse complement heptamers.
        heptamer_ref = NULL,

        #' @field sliding_window Numeric vector of the window to slide across the sequence.
        sliding_window = NULL,

        #' @field maxloopsize Numeric vector of maximum loop size for G4 sequences.
        maxloopsize = 12,

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

        #' @field tfbs_euclidean_distance Data.Table of normalised
        #' euclidean distance calculations.
        tfbs_euclidean_distance = NULL,

        #' @field occupancy_scores Data.Table of k-mer enrichment and depletion scores
        #' extracted from protein occupancy profiles.
        occupancy_scores = NULL,

        #' @field qm_parameters Data.Table of average heptamer quantum mechanical parameters.
        qm_parameters = NULL,

        #' @field exp Character vector of experiment name.
        exp = NULL,

        #' @field column_names Character vector of all column names in feature matrix.
        column_names = NULL,

        #' @field backup_feature_matrix Data.Table of Feature Matrix for backup.
        backup_feature_matrix = NULL,

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
                "breaks" = Biostrings::getSeq(private$ref, self$true_breaks_expanded$mid_range),
                "control" = Biostrings::getSeq(private$ref, self$control_breaks_expanded$mid_range)
            )

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
        #' Generate k-mer table for reference when counting k-mer frequencies.
        #' @return None.
        generate_kmer_table = function(k){
            if(missing(k)) k <- private$kmer_window
            all.kmers <- do.call(data.table::CJ, 
                                 rep(list(c("A", "C", "G", "T")), k))
            all.kmers <- all.kmers[, do.call(paste0, .SD)]
            kmer_ref <- data.table(
                'kmer' = all.kmers,
                'rev.comp' = as.character(
                    Biostrings::reverseComplement(Biostrings::DNAStringSet(all.kmers))
                ))
            kmer_ref[, `:=`(cond = ifelse(seq(1:nrow(.SD)) < match(kmer, rev.comp), 
            TRUE, ifelse(kmer == rev.comp, TRUE, FALSE)))]
            kmer_ref <- kmer_ref[cond == TRUE, .(kmer, rev.comp)]

            if(k == 7){
                private$heptamer_ref <- kmer_ref[, "kmer"]
            } else {
                private$kmer_ref <- kmer_ref
            }
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
                    "breaks" = Biostrings::getSeq(private$ref, self$true_breaks_expanded$mid_range),
                    "control" = Biostrings::getSeq(private$ref, self$control_breaks_expanded$mid_range)
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
            cur.msg <- ifelse(
                data == "breaks",
                "Predicting DNA shape within true break regions",
                "Predicting DNA shape within control regions"
            )
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")

            datatable.seq <- switch(data,
                "breaks" = self$true_breaks_expanded$mid_range,
                "control" = self$control_breaks_expanded$mid_range
            )

            # predict DNA shape
            DNAshapeR::getFasta(datatable.seq, private$ref, 
                                width = as.numeric(width(datatable.seq)[1]), 
                                filename = "temp.fa")
            pred <- DNAshapeR::getShape("temp.fa")

            # compute summary statistics
            shape.name <- c("HelT", "MGW", "ProT", "Roll")
            shape.pred <- lapply(1:length(shape.name), function(x){
                row.mean <- matrixStats::rowMeans2(pred[[x]], na.rm = TRUE)
                row.sd <- matrixStats::rowSds(pred[[x]], na.rm = TRUE)
                combined <- cbind(row.mean, row.sd)
                colnames(combined) <- c(paste0(shape.name[x], ".mean"), 
                                        paste0(shape.name[x], ".sd"))
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
        },

        #' @description 
        #' Extracts the normalised euclidean distances between
        #' the PWMs from core homo sapiens TFBS and PWMs of 4^k-mers.
        #' @param data Character vector of "breaks" or "control".
        #' @return None.
        get_tfbs_euclidean_distance = function(data){
            # progress message
            t1 <- Sys.time()
            cur.msg <- ifelse(
                data == "breaks",
                "Extracting euc.dist of TFBS/k-mers within true break regions",
                "Extracting euc.dist of TFBS/k-mers within control regions"
            )
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            datatable.seq <- switch(data,
                "breaks" = Biostrings::getSeq(private$ref, self$true_breaks[[private$break_score]]$kmer),
                "control" = Biostrings::getSeq(private$ref, self$control_breaks[[private$break_score]]$kmer)
            )

            # import tfbs euclidean distance file
            tfbs.datatable <- fread(
                file = paste0("../data/tfbs/PCs_QueryTable-kmer_", self$k,".csv"),
                # file = paste0("../data/tfbs/QueryTable-kmer_", self$k,".csv")
                showProgress = FALSE
            )

            # match corresponding k-mer in query table
            tfbs.kmers <- Biostrings::DNAStringSet(tfbs.datatable$kmer)

            # forward k-mers
            match.kmers <- match(datatable.seq, tfbs.kmers)
            # reverse k-mers
            any.na <- is.na(match.kmers)
            rev.kmers <- reverseComplement(datatable.seq[any.na])
            match.rev.kmers <- match(rev.kmers, tfbs.kmers)
            # combine results
            match.kmers[any.na] <- match.rev.kmers

            # extract rows based on k-mer matches
            datatable <- tfbs.datatable[match.kmers, -"kmer"]

            if(data == "breaks"){
                private$tfbs_euclidean_distance$true_breaks <- datatable
                private$column_names$TFBS_EUCLIDEAN_DISTANCE <- colnames(datatable)
            } else if(data == "control"){
                private$tfbs_euclidean_distance$control_breaks <- datatable
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description 
        #' Extract the enrichment and depletion scores of each k-mer
        #' for each protein occupancy profile.
        #' @param data Character vector of "breaks" or "control".
        #' @return None.
        get_occupancy_scores = function(data, k){
            # progress message
            t1 <- Sys.time()
            cur.msg <- ifelse(
                data == "breaks",
                "Extracting occupancy k-mer scores within true break regions",
                "Extracting occupancy k-mer scores within control regions"
            )
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            datatable.seq <- switch(data,
                "breaks" = Biostrings::getSeq(
                    private$ref, 
                    self$true_breaks_expanded$kmer_occupancy
                ),
                "control" = Biostrings::getSeq(
                    private$ref, 
                    self$control_breaks_expanded$kmer_occupancy
                )
            )

            # import tfbs euclidean distance file
            datatable <- fread(
                file = paste0("../data/kmertone/QueryTable/QueryTable_occupancy_",
                              "kmer-", 6, "_",
                              private$break_score, ".csv"), 
                showProgress = FALSE
            )
            datatable.category <- datatable$category
            datatable.category <- paste0(datatable$category, "_", private$break_score)
            datatable <- datatable[, -"category"]
            
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

            # # extract columns based on k-mer matches
            datatable <- as.matrix(datatable[, ..match.kmers])
            colnames(datatable) <- NULL
            rownames(datatable) <- datatable.category
            datatable <- as.data.table(t(datatable))

            if(data == "breaks"){
                private$occupancy_scores$true_breaks <- datatable
                private$column_names$OCCUPANCY_SCORES <- colnames(datatable)
            } else if(data == "control"){
                private$occupancy_scores$control_breaks <- datatable
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description 
        #' Extract the QM pre-parameterised k-mers.
        #' @param data Character vector of "breaks" or "control".
        #' @return None.
        get_qm_parameters = function(data){
            # progress message
            t1 <- Sys.time()
            cur.msg <- ifelse(
                data == "breaks",
                "Extracting QM params for 7-mers within true break regions",
                "Extracting QM params for 7-mers within control regions"
            )
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            left.seq <- switch(data, 
                "breaks" = Biostrings::getSeq(private$ref, self$true_breaks_expanded$qm_left_kmer),
                "control" = Biostrings::getSeq(private$ref, self$control_breaks_expanded$qm_left_kmer)
            )
            right.seq <- switch(data, 
                "breaks" = Biostrings::getSeq(private$ref, self$true_breaks_expanded$qm_right_kmer),
                "control" = Biostrings::getSeq(private$ref, self$control_breaks_expanded$qm_right_kmer)
            )

            # import tfbs euclidean distance file
            qm.table <- fread(
                file = "../data/QM_parameters/compiled_data/denergy.txt",
                showProgress = FALSE
            )
            # only keep diff in heat of formation, and ionisation potential features
            to.keep <- colnames(qm.table)[grepl(
                pattern = "seq|dEhof|dIP", 
                x = colnames(qm.table)
            )]
            qm.table <- qm.table[, ..to.keep]

            #' @description 
            #' Process forward and reverse k-mers
            #' @param qm.datatable Data.Table of quantum mechanical parameters.
            #' @param which.seq GRanges object of the heptamers to process.
            #' @return None.
            process.kmers <- function(qm.datatable, which.seq){
                # only keep first occurring kmer in lexicological order 
                qm.datatable <- qm.datatable[match(private$heptamer_ref$kmer, qm.datatable$seq)]
                query.kmers <- Biostrings::DNAStringSet(qm.datatable$seq)

                # forward k-mers
                match.kmers <- match(which.seq, query.kmers)
                # reverse k-mers
                any.na <- is.na(match.kmers)
                rev.kmers <- Biostrings::reverseComplement(which.seq[any.na])
                match.rev.kmers <- match(rev.kmers, query.kmers)
                # combine results
                match.kmers[any.na] <- match.rev.kmers

                # extract columns based on k-mer matches
                datatable <- qm.datatable[match.kmers, -"seq"]
                return(as.matrix(datatable))
            }

            left.seq <- process.kmers(qm.datatable = qm.table, which.seq = left.seq)
            right.seq <- process.kmers(qm.datatable = qm.table, which.seq = right.seq)
            avg.seq <- (left.seq+right.seq)/2
            avg.seq <- as.data.table(avg.seq)

            if(data == "breaks"){
                private$qm_parameters$true_breaks <- avg.seq
                private$column_names$QM_PARAMETERS <- colnames(avg.seq)
            } else if(data == "control"){
                private$qm_parameters$control_breaks <- avg.seq
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        }
    )
)