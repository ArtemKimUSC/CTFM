# CT-FM - Cell Type Fine-Mapping
CT-FM is a statistical genetics method to identify causal cell types underlying diseases and complex traits. We leverage cell type specific candidate cis-regulatory elements (cCREs) together with GWAS summary statistics while accounting for co-regulation of cis-regulatory effects across diverse cell types.

Contacts:<br /> 
Artem Kim artemkim@usc.edu<br />
Steven Gazal gazal@usc.edu<br />

# Installation

## Pre-requisites

Make sure to have anaconda/conda and mamba installed in your environment <br />
https://www.anaconda.com/download <br />
https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html

## 1. Clone the repository
`git clone https://github.com/ArtemKimUSC/CTFM/`


## 2. Install conda environments
CT-FM uses two independent conda environments. <br /> The first one - `ldsc` - is used for pre-processing of GWAS summary statistic files and for running Stratified LD score regression (S-LDSC) method. The second one - `ctfm` - is used to run SuSiE fine-mapping on S-LDSC results.
`cd CTFM`
`mamba env create -f install/ctfm.yml`<br />
`mamba env create -f install/ldsc.yml` <br />

## 3. Test conda environments 
Once created, try to activate the two conda environments and load the susieR library and make sure there are no errors.

-Testing ldsc <br />

`conda activate ldsc` <br />
`python ~/bin/ldsc/ldsc.py -h` <br />
`conda deactivate` <br />

-Testing ctfm <br />

`conda activate ctfm`<br />
`R`<br />
`library(susieR)`<br />
`q()`<br />

## 4. Download reference files
To run CT-FM, you need to download the 1000 Genomes reference files and default S-LDSC annotations. 
**Warning : the download size is big (~20 GB)**

`bash install/download_reference.sh`

The script will download all necessary files and put them in the `CTFM/data/` folder


# Running the default CT-FM pipeline

The following workflow uses the CT-FM pipeline with 927 cCRE annotations of diverse human cell types, as it was performed in the CT-FM paper : **arxiv link**

## 1. Process GWAS summary statistics
The goal here is to take your GWAS sumstats and convert it to a hg19 S-LDSC friendly format using the script 1_create_sumstats.pl
The code below assumes that it is launched from the root directory and that the unprocessed GWAS sumstats is gzipped and already located in the sumstats_in/ folder.

`conda activate ldsc` <br />

`perl scripts/1_create_sumstats.pl --filein sumstats_in/YOUR_SUMSTATS \   # precise the name of your GWAS sumstats omitting the .gz extension` <br />
 `--fileout sumstats_ready/YOUR_NEW_SUMSTATS \ # output sumstats, put it in the sumstats_ready/ directory for downstream analyses`<br />
 `--chr  \`<br />
 `--pos  \`<br />
 `--beta  \`<br />
 `--se  \`<br />
 `--A1  \`<br />
 `--A2  \`<br />
 `--hg  \`<br />
 ` --N`<br />
`conda deactivate`

For this script, we need to precise:<br />
-the input GWAS summary statistics file (the one you put in the sumstats_in folder): `--filein` argument <br />
-the output processed GWAS summary statistics file (put it in the sumstats_ready folder): `--fileout` argument<br />
-the columns (**column numbers**) corresponding to the chromosome (`--chr`), position (`--pos`), beta (`--beta`) and standard error (`--se`) estimates, alleles A1 (`--A1`) and A2 (`--A2`)<br />
-the genome version of the GWAS summary statistics (`--hg 38` or `--hg 19`) and the samplesize column number (`--N`)<br />


### Additional arguments: <br />

The script was also designed to take as input other GWAS sumstats information, for example: <br /> 

Z-scores: add `--Z` argument (column number) and remove `--beta` and `--se` arguments <br />
OR values: add `--OR` argument (column number) and remove `--beta` and `--se` arguments <br />
p-values: add `--pval` argument (column number) and remove `--beta` and `--se` arguments. Note: this option will convert P values into unsigned Z scores, which allows heritability partitioning analyses but not genetic correlation analyses. If OR is available, you can precise both `--p-val` and `--OR` arguments to estimate signed effects. <br />
total sample size of GWAS: add `--myN` argument (total samplesize number) and remove `--N` argument <br />
total sample size of cases and controls: add `--N1` (number of cases) and `--N2` (number of controls) and remove `--N` / `--myN` arguments <br />

## 2. Run S-LDSC with 927 cCRE annotations

We will use the `2_launch_SLDSC_default.sh` script which will launch 5 analyses (running time ~10-15 hours for each) to run S-LDSC on your GWAS sumstats and on ~1000 different CTS annotations. <br />
The main code for S-LDSC analysis is stored in `scripts/sldsc.sh`<br />
The script to launch is `scripts/2_launch_SLDSC.sh` => If you work on a cluster with SLURM-like job submissions, I recommend you modify this script to launch 5 analyses as 5 separate jobs.<br />


`conda activate ldsc`<br />
`bash scripts/2_launch_SLDSC.sh $YOUR_NEW_SUMSTATS_NAME  # precise the name of your sumstats file omitting the "sumstats.gz" part`<br />
`conda deactivate`<br />

The output will be stored in `CTFM/out/SLDSC/YOUR_NEW_SUMSTATS/`<br />


## 3. Launch fine-mapping for your S-LDSC results 
We will use the R script `3_launch_susie.R` to perform SuSiE fine-mapping for the analyzed GWAS (`$YOUR_NEW_SUMSTATS` argument) and 927 S-LDSC annotations.<br />


`conda activate ctfm`<br />
`Rscript scripts/3_launch_susie.R $directory $YOUR_NEW_SUMSTATS #precise the work directory in which CT-FM was downloaded and the name of your sumstats file`<br />
`conda deactivate`<br />


=> The results are composed of 3 files containing the initial S-LDSC results (`_CTFM_sldsc.txt`), CT-FM PIP values (`_CTFM_pips.txt`) and Credible sets (`_CTFM_CS.txt`) -> stored in `out/susie/`






