#!/usr/bin/bash

pwd="$(pwd)/"
homer_path="/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/09_HOMER/lib"

# run homer motif discoveries
Rscript Process.R $pwd "K562_DMSO_DSBs" $homer_path
Rscript Process.R $pwd "K562_Top2_mediated_DSBs" $homer_path