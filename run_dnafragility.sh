#!/usr/bin/bash

pwd="$(pwd)/"
RNAfold_path="/home/imm/hert6114/ViennaRNA-2.6.4/src/bin/RNAfold"

# 1. run COSMIC analysis
bash ./COSMIC/submit.sh

# 2. run ML demonstration studies
bash ./00_ML_proof_of_concept/submit.sh $RNAfold_path

# 3. run generalised fragility model studies
bash ./01_LGBM_FullGenome/submit.sh $RNAfold_path

# 4. check if fwd and rev comp sequence predictions are identical
bash ./08_Forward_vs_RevComp/submit.sh $RNAfold_path

# 5. run virus fragility
bash ./02_Virus_fragility/submit.sh $RNAfold_path

# 6. run delta fragility
bash ./05_DeltaFragility/submit.sh $RNAfold_path