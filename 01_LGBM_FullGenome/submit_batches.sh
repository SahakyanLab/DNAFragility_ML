#!/bin/bash

# extract feature matrix of the full human genome
for chr in {1..22}
do
    touch "run_batch_$chr.sh"
    echo "#!/bin/bash" > "run_batch_$chr.sh"
    echo "Rscript 01_Extract_Features_for_Model.R FALSE $chr" >> "run_batch_$chr.sh"
done