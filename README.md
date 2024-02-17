# DNAFragility

Development and application of the generalised DNA fragility prediction engine based only on sequence-context.

## 1. Software requirements
The resource-demanding computations were performed on a single NVIDIA RTX A6000 GPU with 40GB RAM. The developed workflows and analyses employed the [R programming language 4.3.2](https://www.r-project.org/) and [Python 3.9.12](https://www.python.org/).

Please run the below script to install the latest versions of the R and Python packages necessary to perform the calculations and analyses. 

```bash
bash ./setup/install_packages.sh
```

Please also download and install the below other software.

### Edlib
* Please clone the repo from [this link](https://github.com/Martinsos/edlib) (Edlib >= 1.2.7). Place the [edlib.h](https://github.com/Martinsos/edlib/tree/master/edlib/include) and [edlib.cpp](https://github.com/Martinsos/edlib/tree/master/edlib/src) into [lib/edlib/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/lib/edlib) and [01_LGBM_FullGenome/lib/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/01_LGBM_FullGenome/lib/edlib)

### phmap.hpp via gtl
* Please clone the repo from [this link](https://github.com/greg7mdp/gtl). Place the contents of gtl into [lib/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/lib) and [01_LGBM_FullGenome/lib/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/lib)

## 2. Public files to download
### Cancer-associated DNA strand breaks
We retrieved all the somatic mutation data of both the non-coding and coding regions associated with cancer from the [Catalogue of Somatic Mutations in Cancer (COSMIC) database](https://cancer.sanger.ac.uk/cosmic/), including Non-Coding Variants, Cancer Gene Census, and Breakpoints (structural variants) datasets obtained from release v98, May 2023.

These versions can be downloaded from the following links:
* [Non-Coding Variants](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v98/noncodingvariantstsv) as "Cosmic_NonCodingVariants_PROCESSED.tsv"
* [Cancer Gene Census](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v98/cancergenecensus) as "Cosmic_MutantCensus_PROCESSED.tsv"
* [Breakpoints](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v98/breakpoints) as "Cosmic_Breakpoints_v98_GRCh38.tsv".
* [Classification](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v98/classification) as "Cosmic_Classification_v98_GRCh38.tsv"

Unpack and extract the relevant files. Place the contents into [COSMIC/data/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/COSMIC/data/COSMIC) folder. Please note, we renamed the above first two files with the "PROCESSED" suffix, as the files were very large due to the SNPs, hence, we removed them. We suggest you do this too, unless you have sufficient memory to load and process them all.

### ClinVar SVs and SNPs
We obtained SVs and SNPs from the variant_summary.txt.gz file downloaded from the [ClinVar database](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz) accessed on December 6<sup>th</sup> in 2023 that had a clinically associated pathogenic or benign label. Please note, this file gets updated weekly.

Unpack and extract the relevant files. Place the contents into [04_ClinVar/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/04_ClinVar) folder. 

### Annotation of genomic and genic features on the human genome

* [Gene annotation](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz)
* [Centromere and Pericentromere](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.0.bed)
* [RepeatMasker](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed)
* [Telomere](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_telomere.bed)
* [Housekeeping genes from the HRT Atlas](https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv) as "Human_Housekeeping_Genes.csv"
* [Chromosomal fragile sites from HumCFS](https://webs.iiitd.edu.in/raghava/humcfs/fragile_site_bed.zip). Unpack the individual files into the [COSMIC/data/annotations/chromosome_bed](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/COSMIC/data/COSMIC/data/annotations/chromosome_bed) folder.
* [Sites of G4 structures](https://github.com/SahakyanLab/G4Damage/blob/master/raw_data/PQSdata.txt.gz) as "G4_PQS.txt". Then, convert this file to "G4_PQS.csv".

To download CpG islands and Isochores from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables/), please select the following:

* [CpG Islands.](https://genome.ucsc.edu/cgi-bin/hgTables/) clade: Mammal, genome: Human, assembly: Jan 2022 (T2T CHM13v2.0/hs1), group: All Tracks, track: CpG Islands, table: hub_3671779_cpgIslandExtUnmasked, output format: BED - browser extensible data, output filename: output_CpG_Islands.csv.
* [Isochores.](https://genome.ucsc.edu/cgi-bin/hgTables/) clade: Mammal, genome: Human, assembly: May 2004 (NCBI35/hg17), group: All Tracks, track: Isochores, table: ct_Isochores_9145, output format: BED - browser extensible data, output filename: iso_hg17.bb.

Unpack and extract the relevant files. Place the contents into [COSMIC/data/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/COSMIC/data/annotations) folder. 

* [Cancer driver genes from COSMIC relased v98, May 2023](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v98/cancergenecensus)

Unpack and extract the relevant files. Place the contents into [COSMIC/data/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/COSMIC/data/COSMIC) folder. 

## 3. Liftover files

We processed all datasets in the reference genome version used as per the deposition. When doing comparative analysis, we lifted the genomic coordinates over to the latest T2T genome assembly. 

* [hg17 to hg38](https://hgdownload.cse.ucsc.edu/goldenpath/hg17/liftOver/#:~:text=hg17ToHg38.over.chain.gz)
* [hg19 to hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/#:~:text=hg19ToHg38.over.chain.gz)
* [hg38 to T2T](https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38-chm13v2.over.chain.gz)

Unpack and extract the relevant files. Place the contents into [COSMIC/data/liftover/](https://github.com/SahakyanLab/DNAFragility_ML/tree/master/COSMIC/data/liftover) folder. 

## 4. Reference sequences

We processed all datasets in the reference genome version used as per the deposition. For Kmertone, the individual fasta files were needed. This GitHub repo is dependent on the results of [DNAFragility_dev](https://github.com/SahakyanLab/DNAFragility_dev), where the reference genomes are downloaded already.