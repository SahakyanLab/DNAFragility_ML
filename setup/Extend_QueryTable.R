Breakpoints <- R6::R6Class(
    classname = "Breakpoints",
    public = list(
        initialize = function(break_score){
            if(!missing(break_score)) private$break_score <- break_score
        },

        #' @description
        #'  Extend table to have all k-mers not just lexicographically occurring ones.
        #' @return Data table of the query table.
        get_extended_tables = function(k = c(2,4,6,8)){
            for(kmer in k){
                private$kmer_window <- kmer
                private$get_extended_querytable(kmer_size = private$kmer_window)
            }
        }        
    ),
    private = list(
        #' @field break_score Character vector of "zscore" or "ratios".
        break_score = NULL,

        #' @field kmer_window Numeric vector of the window size for counting k-mers.
        kmer_window = NULL,

        #' @field kmer_ref Data.Table of forward and reverse complement k-mers.
        kmer_ref = NULL,

        #' @field kmer_list Character vector of all k-mers.
        kmer_list = NULL,

        #' @field extended_preparam_table Data Table of the extended query table.
        extended_preparam_table = NULL,

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
        get_extended_querytable = function(kmer_size){
            # progress message
            t1 <- Sys.time()
            cur.msg <- paste0("Retrieving table of breakage parameters for ", kmer_size, "-mer")
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

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

            private$generate_kmer_table()

            # import all values from pre-parameterised table
            preparam.table <- fread(
                paste0("./DNAfrAIlib/QueryTable_",
                    "kmer-", private$kmer_window, "_",
                    private$break_score, ".csv"), 
                showProgress = FALSE
            )
            category_col <- preparam.table[, "category"]
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

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        }
    )
)

suppressPackageStartupMessages(suppressWarnings(library(data.table)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(gtools)))
data.table::setDTthreads(threads = 1)

Breaks <- Breakpoints$new(break_score = "zscore")
Breaks$get_extended_tables(k = c(6,8))