#!/bin/bash

# extract feature matrix of the full human genome
for batch in {1..8}
do
    batch_start=$(( (batch - 1) * 5 + 1 ))
    batch_end=$(( batch * 5 ))

    touch "run_batch_$batch.sh"
    echo "#!/bin/bash" > "run_batch_$batch.sh"

    for sub_batch in $(seq $batch_start $batch_end)
    do
        echo "Rscript 01_Get_Features_SV.R $sub_batch" >> "run_batch_$batch.sh"
    done
done