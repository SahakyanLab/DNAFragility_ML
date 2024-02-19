# DNAFragility

Development and application of the generalised DNA fragility prediction engine based only on sequence-context.

## Setup

Clone the project:

```
git clone https://github.com/SahakyanLab/DNAFragility_ML.git
```

Please follow the instructions below on how to acquire the public datasets, setup the directory stucture, and software necessary to run all the studies from the publication.  At the end of this `README` file, you can find two separate bash script commands that runs the majority of the setup and runs the calculations sequentially. 

## 1. Software requirements
The resource-demanding computations were performed on a single NVIDIA RTX A6000 GPU with 40GB RAM. The developed workflows and analyses employed the [R programming language 4.3.2](https://www.r-project.org/) and [Python 3.9.12](https://www.python.org/).

Please run the below script to install the latest versions of the R and Python packages necessary to perform the calculations and analyses. 

```bash
bash ./setup/install_packages.sh
```

Please also download and install the below software.

### Edlib
* Please clone the repo from [this link](https://github.com/Martinsos/edlib) (Edlib >= 1.2.7). Place the [edlib.h](https://github.com/Martinsos/edlib/tree/master/edlib/include) and [edlib.cpp](https://github.com/Martinsos/edlib/tree/master/edlib/src) into [lib/edlib/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/lib/edlib) folder.

### Secondary structure folding parameter file
* Please download the [DNA parameter file](https://github.com/ViennaRNA/ViennaRNA/blob/master/misc/dna_mathews2004.par) and place it into [data/parameters](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/data/parameters) folder.

### `ggpattern`

Please note, if you are using Ubuntu, you may have trouble installing the [ggpattern](https://github.com/trevorld/ggpattern) R package. However, the below steps has worked for us. 

1. sudo apt-get install libmagick++-dev
2. sudo apt install libgdal-dev
3. sudo apt-get install -y libudunits2-dev
4. install.packages("units")
5. install.packages("sf")
6. install.packages("gridpattern")
7. install.packages("ggpattern")

## 2. Public files to download
### Cancer-associated DNA strand breaks
We retrieved all the somatic mutation data of both the non-coding and coding regions associated with cancer from the [Catalogue of Somatic Mutations in Cancer (COSMIC) database](https://cancer.sanger.ac.uk/cosmic/), including Non-Coding Variants, Cancer Gene Census, and Breakpoints (structural variants) datasets obtained from release v98, May 2023.

These versions can be downloaded from the following links:
* [Non-Coding Variants](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v98/noncodingvariantstsv) as "Cosmic_NonCodingVariants_PROCESSED.tsv"
* [Coding Variants](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v98/mutantcensus) as "Cosmic_MutantCensus_PROCESSED.tsv"
* [Cancer Gene Census](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v98/cancergenecensus) as "Cosmic_CancerGeneCensus_v98_GRCh38.tsv"
* [Breakpoints](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v98/breakpoints) as "Cosmic_Breakpoints_v98_GRCh38.tsv".
* [Classification](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v98/classification) as "Cosmic_Classification_v98_GRCh38.tsv"

Unpack and extract the relevant files. Place the contents into [COSMIC/data/COSMIC/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/COSMIC/data/COSMIC) folder. Please note, we renamed the above first two files with the "PROCESSED" suffix, as the files were very large due to the SNPs, hence, we removed them. We suggest you do this too, unless you have sufficient memory to load and process them all.

### ClinVar SVs and SNPs
We obtained SVs and SNPs from the variant_summary.txt.gz file downloaded from the [ClinVar database](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz) accessed on December 6<sup>th</sup> in 2023 that had a clinically associated pathogenic or benign label. Please note, this file gets updated weekly.

Unpack and extract the relevant files. Place the contents into [04_ClinVar/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/04_ClinVar) folder. 

### Annotation of genomic and genic features on the human genome

* [Gene annotation](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz)
* [Centromere and Pericentromere](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.0.bed)
* [RepeatMasker](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed)
* [Telomere](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_telomere.bed)
* [Housekeeping genes from the HRT Atlas](https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv) as "Human_Housekeeping_Genes.csv"
* [Gene-centric chromosomal fragile sites from HumCFS](https://webs.iiitd.edu.in/raghava/humcfs/fragile_site_bed.zip). Unpack the individual files into the [COSMIC/data/annotations/chromosome_bed](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/COSMIC/data/COSMIC/data/annotations/chromosome_bed) folder.
* [Sites of G4 structures](https://github.com/SahakyanLab/G4Damage/blob/master/raw_data/PQSdata.txt.gz) as "G4_PQS.txt". Then, convert this file to "G4_PQS.csv".

To download CpG islands and Isochores from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables/), please select the following:

* [CpG Islands.](https://genome.ucsc.edu/cgi-bin/hgTables/) clade: Mammal, genome: Human, assembly: Jan 2022 (T2T CHM13v2.0/hs1), group: All Tracks, track: CpG Islands, table: hub_3671779_cpgIslandExtUnmasked, output format: BED - browser extensible data, output filename: output_CpG_Islands.csv.
* [Isochores.](https://genome.ucsc.edu/cgi-bin/hgTables/) clade: Mammal, genome: Human, assembly: May 2004 (NCBI35/hg17), group: All Tracks, track: Isochores, table: ct_Isochores_9145, output format: BED - browser extensible data, output filename: iso_hg17.bb.

Unpack and extract the relevant files from above. Place the contents into [COSMIC/data/annotations/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/COSMIC/data/annotations) folder. 

* [Cancer driver genes from COSMIC relased v98, May 2023](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v98/cancergenecensus). Unpack and extract the relevant files. Place the contents into [COSMIC/data/COSMIC/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/COSMIC/data/COSMIC) folder. 

### Chromothripsis breakpoint events

We obtained the chromothripsis breakpoint cases from ChromothripsisDB. Please download the dataset from Download -> Full Dataset -> [Chromothripsis case data](http://cailab.labshare.cn/ChromothripsisDB/download/)

Unpack and extract the relevant files from above. Place the contents into [03_Chromothripsis/data](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/03_Chromothripsis/data) folder.

### Transcription Factor data

We retrieved 247 core-validated vertebrate transcript factor binding sites (TFBS) from the [JASPAR 2024 database](https://jaspar.elixir.no/downloads/).

* Download the [vertebrate dataset from here](https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_redundant_pfms_jaspar.txt)
* Download the [bed files from here](https://jaspar.elixir.no/download/data/2024/bed.tar.gz) as "jaspar_beds"

Unpack and extract the relevant files from above. Place the contents into [data/TFBS/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/data/JASPAR/) folder. 

## 3. Liftover files

We processed all datasets in the reference genome version used as per the deposition. When doing comparative analysis, we lifted the genomic coordinates over to the latest T2T genome assembly. 

* [hg17 to hg38](https://hgdownload.cse.ucsc.edu/goldenpath/hg17/liftOver/#:~:text=hg17ToHg38.over.chain.gz)
* [hg19 to hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/#:~:text=hg19ToHg38.over.chain.gz)
* [hg38 to hg19](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz)
* [hg38 to T2T](https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38-chm13v2.over.chain.gz)

Unpack and extract the relevant files. Place the contents into [COSMIC/data/liftover/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/COSMIC/data/liftover) folder. 

## 4. Reference sequences

We processed all datasets in the reference genome version used as per the deposition. For Kmertone, the individual fasta files were needed. This GitHub repo is dependent on the results of [DNAFragility_dev](https://github.com/SahakyanLab/DNAFragility_dev), where the reference genomes are downloaded already.

## 5. `DNAfrAIlib` feature library

The genomic sequence-based octameric features can be downloaded from the [DNAfrAIlib](https://github.com/SahakyanLab/DNAfrAIlib) repo. The quantum mechanical hexameric parameters can be downloaded from [DNAkmerQM](https://github.com/SahakyanLab/DNAkmerQM/tree/master/6-mer).

This has been automatically setup if you run the below bash script.

```bash
bash get_feature_lib.sh
```

## 6. Notes on `00_ML_proof_of_concept` folder

To run the `00_ML_proof_of_concept` work, you need to have two datasets downloaded and processed following the method from [DNAFragility_dev](https://github.com/SahakyanLab/DNAFragility_dev). The demonstration used in the paper is based on [this study](https://doi.org/10.1016/j.molcel.2019.05.015) with data deposited on [the GEO database](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121742). We specifically used [DMSO-treated, endogenous DNA fragility in K562 cells](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3444988). You can also run it on the [etoposide-treated DNA fragility in K562 cells enriched at topoisomerase II sites](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3444989).

For any ML task, you require the genomic sequence range of influence for each of the short-, medium-, and long-range effects. Depending on the dataset used, some datasets had to be pre-processed to handle 5'-3' DNA strand breaks. Hence, running the full [DNAFragility_dev](https://github.com/SahakyanLab/DNAFragility_dev) study beforehand is strongly advised.

Alternatively, if you wish to skip the [DNAFragility_dev](https://github.com/SahakyanLab/DNAFragility_dev) process, and just want to process these DNA strand breaks for the present study, please run the below bash script.

```bash
bash get_MLdemo_datasets.sh
```

## 7. Other notes

* All cpp files are interfaced *via* the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) library in R with `omp.h` when possible. Please ensure you have this installed.
* `RcppArmadillo.h` and `RcppEigen.h` are necessary for the feature extraction process. Please ensure you have this installed. By default, will not use it in case you have not installed it.
* Various model predictions have been deposited if the compressed file size was within the GitHub file size limit. If you wish to view and/or use them, please `gunzip` the files.
* While this repo can run as a standalone study, the results are dependent on [DNAFragility_dev](https://github.com/SahakyanLab/DNAFragility_dev) and when possible, we have deposited the necessary dependent files. 
* When you run the `run_dnafragility.sh` bash script, you will need to include the path to the viennaRNA `RNAFold` programme as the first argument. Some operating systems allow you to interface it directly *via* `RNAfold`, others require the literal path to the programme.

## 8.Run all setup files

If you wish to run all setups, including all the aforementioned bash scripts, please run the below bash script.

```bash
bash run_all_setup_files.sh
```

## 9. Run the full DNAFragility_ML study

Please note that many of the calculations were computationally intensive, particularly the `01_LGBM_FullGenome` and `05_DeltaFragility` folders. Most things were run in parallel in smaller batches. However, if you submit the below bash script, it runs all scripts sequentially. This can **take several months** to complete. 
Most tasks take up several tens to hundreds of GBs worth of RAM. The entire study requires between 2-4 TB of hard drive space. 

You may need to monitor your memory usage, memory cache, and swap to ensure calculations run smoothly.

**Arguments**

* `Rnafold_path` path to the RNAfold function for secondary structure prediction.
* `fast_matrix` If TRUE, will use fast RcppArmadillo matrix calculations. Default FALSE.

```bash
bash run_dnafragility.sh $RNAfold_path $fast_matrix
```