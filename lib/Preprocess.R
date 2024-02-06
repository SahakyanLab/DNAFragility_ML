Preprocess <- R6::R6Class(
    classname = "Preprocess",
    public = list(
        #' @field fasta_sequence DNAString instance of the species sequence.
        fasta_sequence = NULL,

        #' @field label Character vector of the species label.
        label = NULL,

        #' @field k Numeric vector of k-mer size.
        k = 8,

        #' @field break_type Character vector of any of c("Biological", "Enzymatic", "High_frequency"). 
        #'  If user knows the type of breakage but does not have any range effect values, the maximum
        #'  range effect cross similar types of breakages will be taken. If no type is given, user 
        #'  has to input 3 range effects into the ranges parameter.
        break_type = "Biological",

        #' @field df_bp Data.table of positions to extract.
        df_bp = NULL,

        #' @field Character vector. Type of structural variant. 
        SV_type = "",

        initialize = function(fasta_sequence, label, k, break_type, ranges, SV_type){
            if(!missing(fasta_sequence)) self$fasta_sequence <- fasta_sequence
            if(!missing(label)) self$label <- label
            if(!missing(k)) self$k <- k
            if(!missing(break_type)) self$break_type <- break_type
            if(!missing(SV_type)) self$SV_type <- SV_type

            # extract 3 ranges of sequence influences
            private$get_ranges(ranges = ranges)
        },

        #' @description
        #'  Pad sequences and extract single nucleotide positions for feature extractions.
        #' @return None.
        preprocess_sequences = function(){
            if(self$SV_type != "SNP"){
                private$pad_fasta_sequence()
            }

            private$extract_positions()
        }
    ),
    private = list(
        #' @field ranges Numeric vector of upper limit of the short, mid and long ranges.
        ranges = NULL,

        #' @description 
        #' Extracts the short, mid and long-range upper limits around the central breakpoint.
        #' @param ranges Numeric vector of length 3, where position 1 is the short range effect, 
        #'  position 2 is the medium range effect, and position 3 is the long range effect.
        #' @return None.
        get_ranges = function(ranges){
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
                row.id <- which(datatable$break_type == self$break_type)
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
        #'  Pad the fasta sequence on the left and right to not have out-of-bounds errors.
        #' @return None.
        pad_fasta_sequence = function(){
            pad_seq <- rep("N", (private$ranges$long_range/2))

            self$fasta_sequence <- c(
                Biostrings::DNAString(paste(pad_seq, collapse = "")), 
                self$fasta_sequence,
                Biostrings::DNAString(paste(pad_seq, collapse = ""))
            )
        },

        #' @description
        #'  Extract each position except for the padded sequence regions.
        #' @return None.
        extract_positions = function(){
            if(self$SV_type == "SNP"){
                # self$df_bp <- data.table(start.pos = ceiling(width(self$fasta_sequence) / 2))
                self$df_bp <- data.table(start.pos = ceiling(nchar(self$fasta_sequence) / 2))
            } else {
                start_pos <- (private$ranges$long_range/2)+(self$k/2)
                end_pos <- length(self$fasta_sequence)-(private$ranges$long_range/2)-(self$k/2)
                self$df_bp <- data.table(start.pos = start_pos:end_pos)

                self$fasta_sequence <- Biostrings::DNAStringSet(self$fasta_sequence)
            }
            names(self$fasta_sequence) <- self$label
        }
    )
)