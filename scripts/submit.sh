#!/bin/bash
# Rscript Cor.R "$(pwd)/" "all" 8 "all"
Rscript Cor.R "$(pwd)/" "all" 6 "all"
Rscript Cor.R "$(pwd)/" "all" 4 "all"

Rscript Cor.R "$(pwd)/" "zscore" 8 "all"
Rscript Cor.R "$(pwd)/" "zscore" 6 "all"
Rscript Cor.R "$(pwd)/" "zscore" 4 "all"

Rscript Cor.R "$(pwd)/" "ratio" 8 "all"
Rscript Cor.R "$(pwd)/" "ratio" 6 "all"
Rscript Cor.R "$(pwd)/" "ratio" 4 "all"