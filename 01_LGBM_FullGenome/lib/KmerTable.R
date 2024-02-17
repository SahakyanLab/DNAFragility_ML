KmerTable <- R6::R6Class(
    classname = "KmerTable",
    public = list(
        #' @field k Numeric vector of k-mer size.
        k = 8,

        #' @field score Character vector of "zscore" or "ratios".
        break_score = "zscore",

        #' @field statistic Character vector of summation technique
        #' for a group of k-mer scores per breakage type.
        statistic = "mean",

        #' @field exp Character vector of experiment name.
        exp = "",

        #' @field kmer_list Character vector of k-mers.
        kmer_list = NULL,

        #' @field kmer_ref Data.table of first occurring 
        #' kmer in lexicological order.
        kmer_ref = NULL,

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
                return(private$get_data(
                    data = private$org, 
                    category = category
                ))
            })
            df <- rbindlist(df)
            df <- as_tibble(df) %>% 
                tidyr::pivot_wider(
                    names_from = "kmer",
                    values_from = "value"
                )
            fwrite(
                x = df, 
                file = paste0(
                    "../data/kmertone/QueryTable/QueryTable_kmer-", 
                    self$k, "_", self$break_score, ".csv"
                )
            )

            if(length(private$category_leftout) > 0){
                df <- lapply(private$category_leftout, function(category){
                    return(private$get_data(
                        data = private$org_leftout, 
                        category = category
                    ))
                })
                df <- rbindlist(df)
                df <- as_tibble(df) %>% 
                    tidyr::pivot_wider(
                        names_from = "kmer",
                        values_from = "value"
                    )
                fwrite(
                    x = df, 
                    file = paste0(
                        "../data/kmertone/QueryTable/QueryTable_Leftout_kmer-", 
                        self$k, "_", self$break_score, ".csv"
                    )
                )
            }

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

        #' @field org Data.table of the full org_file.csv
        org = NULL,

        #' @field category_leftout Character vector of the breakage type left out.
        category_leftout = NULL,

        #' @field org_leftout Data.Table of the experiment type left out.
        org_leftout = NULL,

        #' @field group_exp Boolean. If TRUE, groups breakage exp by similar types.
        group_exp = TRUE,

        #' @description 
        #' A utility function to generate k-mers.
        #' Then, only keeps first occurring k-mer in lexicological order.
        #' @return None.
        generate_table = function(){            
            k.mers <- do.call(data.table::CJ, rep(list(c("A", "C", "G", "T")), self$k))
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
            org.file[, bp.exp := paste0(`Fragmentation type`, "/", `Experiment folder`)]

            # discard rows for use in machine learning model
            to.discard <- NULL
            if(self$exp != ""){
                to.discard <- which(grepl(
                    pattern = self$exp, 
                    x = org.file$Category_Main, 
                    ignore.case = TRUE
                ))
            }

            if(length(to.discard) > 0){
                private$org_leftout <- org.file[to.discard, ]
                if(private$group_exp){
                    private$category_leftout <- private$org_leftout$Category_Main
                } else {
                    private$category_leftout <- private$org_leftout$bp.exp
                }
                if(!all_exp) org.file <- org.file[-to.discard, ]
            }

            # for breakpoints
            private$org <- org.file[`DSB Map` == TRUE, ]
            if(private$group_exp){
                categories <- unique(private$org$Category_Main)
            } else {
                categories <- private$org$bp.exp
            }
            private$categories <- sort(categories[which(categories != "")])
        },

        #' @description
        #' A utility function to calculate summarised probability ratios 
        #' or z-scores per group of breakage source.
        #' @param category Character vector of a single breakage source.
        #' @return A data.table of k-mer scores per breakage category.
        get_data = function(data, category){
            if(private$group_exp){
                org.filtered <- data[Category_Main == category, bp.exp]
            } else {
                org.filtered <- data[bp.exp == category, bp.exp]
            }

            # mean/median of same breakage type
            out <- lapply(1:length(org.filtered), function(x){
                return(private$summarise_score(file = org.filtered[x]))
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
        summarise_score = function(file){
            files <- list.files(
                path = paste0("../data/kmertone/", file),
                pattern = paste0("score_", self$k),
                full.names = TRUE
            )

            data.set <- fread(
                file = files,
                sep = ",", header = TRUE, showProgress = FALSE,
                select = c("case", "control", "z")
            )
            data.set[, kmer := self$kmer_list]

            #' @description
            #'  Impute missing values by 5 nearest neighbour of k-mer similarity.
            #' @param df Data.table of the query table.
            #' @param column Character vector of the column to impute.
            #' @param to_impute Numeric vector of the values to impute. 
            #' @return None.
            impute_missing_vals <- function(df, column, to_impute){
                # generate rolling k-mers
                which_kmers <- df$kmer[to_impute]

                for(ind in 1:length(which_kmers)){
                    kmer <- which_kmers[ind]
                    dat <- df[, ..column][[1]]

                    # replace with 5 nearest neighbour result
                    mat <- cbind(df$kmer, rep(kmer, nrow(df)))

                    # global (NW) alignment
                    levdist_mat <- LevenshteinLoop(mat)[,1]

                    # first five nearest neighbours of finite values 
                    # and excluding the first result as it's the query kmer.
                    k <- 5
                    nn_ind <- order(levdist_mat)[2:length(levdist_mat)]
                    pot_candidates <- dat[nn_ind]
                    target_val <- mean(pot_candidates[is.finite(pot_candidates)][1:k])
                    dat[to_impute[ind]] <- target_val

                    if(column == "z") df[, z := dat]
                    if(column == "ratio") df[, ratio := dat]
                }
            }

            # if any infinite or NAs, impute using smaller k-mers
            to_impute <- which(!is.finite(data.set$z))
            if(length(to_impute) > 0){
                impute_missing_vals(
                    df = data.set, 
                    column = "z", 
                    to_impute = to_impute
                )    
            }

            # get probability ratios
            data.set[, ratio := {
                norm.case <- data.set$case/sum(data.set$case, na.rm = TRUE)
                norm.control <- data.set$control/sum(data.set$control, na.rm = TRUE)
                ratio <- norm.case/norm.control 
            }]

            to_impute <- which(!is.finite(data.set$ratio))
            if(length(to_impute) > 0){
                impute_missing_vals(
                    df = data.set, 
                    column = "ratio", 
                    to_impute = to_impute
                )    
            }

            # only keep first occurring kmer in lexicological order 
            data.set <- data.set[match(self$kmer_ref$kmer, data.set$kmer)]

            output <- switch(action,
                "ratio" = data.set$ratio,
                "zscore" = data.set$z
            )
            return(output)
        }
    )
)