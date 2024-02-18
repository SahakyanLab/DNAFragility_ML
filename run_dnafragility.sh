#!/usr/bin/bash

RNAfold_path="/home/imm/hert6114/ViennaRNA-2.6.4/src/bin/RNAfold"

# 1. run COSMIC analysis
cd ./COSMIC/scripts
bash submit.sh
cd ../../

# # 2. run ML demonstration studies
# cd ./00_ML_proof_of_concept/
# bash submit.sh $RNAfold_path
# cd ../../

# # 3. run generalised fragility model studies
# cd ./01_LGBM_FullGenome/
# bash submit.sh $RNAfold_path
# cd ../../

# # 4. check if fwd and rev comp sequence predictions are identical
# cd ./08_Forward_vs_RevComp/
# bash submit.sh $RNAfold_path
# cd ../../

# # 5. run virus fragility
# cd ./02_Virus_fragility/
# bash submit.sh $RNAfold_path
# cd ../../

# # 6. run delta fragility
# cd ./05_DeltaFragility/
# bash submit.sh $RNAfold_path
# cd ../../