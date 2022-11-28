KmerTable <- R6::R6Class(
    classname = "KmerTable",
    public = list(
        #' @field k Numeric vector of k-mer size.
        k = 4,

        #' @field score Character vector of "zscore" or "ratios".
        break_score = "zscore",

        #' @field statistic Character vector of summation technique
        #' for a group of k-mer scores per breakage type.
        statistic = "mean",

        #' @field exp Character vector of experiment name.
        exp = NULL,

        #' @field kmer_list Character vector of k-mers.
        kmer_list = NULL,

        #' @field kmer_ref Data.table of first occurring 
        #' kmer in lexicological order.
        kmer_ref = NULL,

        #' @field kmer_ref_occupancy Data.table of first occurring 
        #' kmer in lexicological order for protein occupancies.
        kmer_ref_occupancy = NULL,

        #' @field leftout_table Data.table of breakage scores
        #' for the experiment left out.
        leftout_table = NULL,

        initialize = function(k, break_score, statistic, 
                              exp, all_exp, group_exp){
            if(!missing(k)) self$k <- k
            if(!missing(break_score)) self$break_score <- break_score
            if(!missing(statistic)) self$statistic <- statistic
            if(!missing(exp)) self$exp <- exp
            if(!missing(group_exp)) private$group_exp <- group_exp

            # generate k-mer table
            private$generate_table()

            # get org file
            private$get_org(all_exp = all_exp)
        },

        #' @description 
        #' Generates k-mer probability ratios or z-score query table.
        #' @return None.
        generate_querytable = function(){
            # progress message
            t1 <- Sys.time()
            cur.msg <- paste0("Generating ", self$break_score, " ", 
                              self$k, "-mer table")
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # Query_Table for k-mer breakage propensities
            df <- lapply(private$categories, function(category){
                return(private$get_data(data = private$org, 
                                        category = category,
                                        occupancy = FALSE))
            })
            df <- rbindlist(df)
            df <- as_tibble(df) %>% 
                tidyr::spread(key = "kmer", value = "value")
            fwrite(
                x = df, 
                file = paste0("../data/kmertone/QueryTable/",
                                "QueryTable_kmer-", 
                                self$k, "_", self$break_score, ".csv")
            )

            df <- lapply(private$category_leftout, function(category){
                return(private$get_data(data = private$org_leftout, 
                                        category = category,
                                        occupancy = FALSE))
            })
            df <- rbindlist(df)
            df <- as_tibble(df) %>% 
                tidyr::spread(key = "kmer", value = "value")
            fwrite(
                x = df, 
                file = paste0("../data/kmertone/QueryTable/",
                                "QueryTable_Leftout_kmer-", 
                                self$k, "_", self$break_score, ".csv")
            )

            # Query_Table for protein occupancies
            df <- lapply(private$categories_occupancy, function(category){
                return(private$get_data(data = private$org_occupancy, 
                                        category = category,
                                        occupancy = TRUE))
            })
            df <- rbindlist(df)
            df <- as_tibble(df) %>% 
                tidyr::spread(key = "kmer", value = "value")
            fwrite(
                x = df, 
                file = paste0("../data/kmertone/QueryTable/",
                                "QueryTable_occupancy_kmer-", 
                                self$k, "_", self$break_score, ".csv")
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description 
        #' Get the breakage scores for the experiment that was left out.
        #' @return None.
        get_leftout_table = function(){
            datatable <- fread(
                file = paste0("../data/kmertone/QueryTable/",
                              "QueryTable_Leftout_kmer-", 
                              self$k, "_", self$break_score, ".csv")
            )
            datatable <- cbind(
                self$kmer_ref, 
                as.data.table(t(datatable)[-1])
            )
            setnames(datatable, c("kmer", self$break_score))
            self$leftout_table <- datatable
        }
    ),
    private = list(
        #' @field categories Character vector of all breakage types.
        categories = NULL,

        #' @field categories_occupancy Character vector of all occupancy types.
        categories_occupancy = NULL,

        #' @field org Data.table of the full org_file.csv
        org = NULL,

        #' @field org_occupancy Data.table of the full org_file.csv occupancies.
        org_occupancy = NULL,

        #' @field category_leftout Character vector of the breakage type left out.
        category_leftout = NULL,

        #' @field org_leftout Data.Table of the experiment type left out.
        org_leftout = NULL,

        #' @field group_exp Boolean. If TRUE, groups breakage exp by similar types.
        group_exp = NULL,

        #' @description 
        #' A utility function to generate k-mers.
        #' Then, only keeps first occurring k-mer in lexicological order.
        #' @return None.
        generate_table = function(){
            k.mers <- do.call(data.table::CJ, 
                              rep(list(c("A", "C", "G", "T")), self$k))
            self$kmer_list <- k.mers[, do.call(paste0, .SD)]
            
            rev.comp <- as.character(
                Biostrings::reverseComplement(Biostrings::DNAStringSet(self$kmer_list))
            )
            kmer_ref <- data.table('kmer' = self$kmer_list, 'rev.comp' = rev.comp)
            kmer_ref[, cond := ifelse(seq(1:nrow(.SD)) < match(kmer, rev.comp), 
            TRUE, ifelse(kmer == rev.comp, TRUE, FALSE))]
            self$kmer_ref <- kmer_ref[cond == TRUE, .(kmer)]
        },
        
        #' @description 
        #' A utility function to get subset of org.file's breakage type.
        #' @return None.
        get_org = function(all_exp){
            org.file <- fread("../data/org_file.csv", showProgress = FALSE)
            org.file[, bp.exp := paste0(Fragmentation_type, "/", 
                                        Experiment_folder)]

            # discard rows for use in machine learning model
            to.discard <- which(grepl(
                pattern = self$exp, 
                x = org.file$Fragmentation_type, 
                ignore.case = TRUE
            ))
            if(length(to.discard) > 0){
                private$org_leftout <- org.file[to.discard, ]
                if(private$group_exp){
                    private$category_leftout <- private$org_leftout$Category_main
                } else {
                    private$category_leftout <- private$org_leftout$bp.exp
                }
                if(!all_exp) org.file <- org.file[-to.discard, ]
            }

            # for breakpoints
            private$org <- org.file[DSB_Map == TRUE, ]
            if(private$group_exp){
                categories <- unique(private$org$Category_main)
            } else {
                categories <- private$org$bp.exp
            }
            private$categories <- sort(categories[which(categories != "")])

            # for occupancy profiles
            private$org_occupancy <- org.file[DSB_Map == FALSE, ]
            if(private$group_exp){
                categories_occupancy <- unique(private$org_occupancy$Category_main)
            } else {
                categories_occupancy <- private$org_occupancy$bp.exp
            }
            private$categories_occupancy <- sort(
                categories_occupancy[which(categories_occupancy != "")]
            )
        },

        #' @description
        #' A utility function to calculate summarised probability ratios 
        #' or z-scores per group of breakage source.
        #' @param category Character vector of a single breakage source.
        #' @return A data.table of k-mer scores per breakage category.
        get_data = function(data, category, occupancy = FALSE){
            if(private$group_exp){
                org.filtered <- data[Category_main == category, bp.exp]
            } else {
                org.filtered <- data[bp.exp == category, bp.exp]
            }

            # mean/median of same breakage type
            out <- lapply(1:length(org.filtered), function(x){
                return(private$summarise_score(file = org.filtered[x],
                                               occupancy = occupancy))
            })

            out <- switch(self$statistic,
                "mean" = {
                    sapply(1:length(out[[1]]), function(x){
                        mean(sapply(out, `[[`, x), na.rm = TRUE)
                    })
                },
                "median" = {
                    sapply(1:length(out[[1]]), function(x){
                        median(sapply(out, `[[`, x), na.rm = TRUE)
                    })
                }
            )

            datatable <- data.table(
                kmer = self$kmer_ref$kmer,
                value = out,
                category = as.factor(category)
            )
            return(datatable)
        },

        #' @description 
        #' A utility function to calculate k-mer probability ratios or z-scores.
        #' @param file Character vector of the specified breakage type.
        #' @return A data.table of scores.
        summarise_score = function(file, occupancy){
            files <- list.files(
                path = ifelse(occupancy,
                              paste0("../data/occupancy_profiles/", 
                                     file, "/kmertone/results"),
                              paste0("../data/kmertone/", file)),
                pattern = paste0("score_", self$k),
                full.names = TRUE
            )

            data.set <- fread(
                file = files,
                sep = ",", header = TRUE, showProgress = FALSE,
                select = c("case", "control", "z")
            )

            data.set[, kmer := self$kmer_list]
            data.set <- data.set[!which(is.na(data.set$z)), ]
            data.set <- data.set[!which(is.infinite(data.set$z)), ]

            # only keep first occurring kmer in lexicological order 
            data.set <- data.set[match(self$kmer_ref$kmer, data.set$kmer)]

            output <- switch(self$break_score,
                "ratio" = {
                    norm.case <- data.set$case/sum(data.set$case, 
                                                   na.rm = TRUE)
                    norm.control <- data.set$control/sum(data.set$control, 
                                                         na.rm = TRUE)
                    ratio <- norm.case/norm.control 
                },
                "zscore" = data.set$z
            )
            return(output)
        }
    )
)