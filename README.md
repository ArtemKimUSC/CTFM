# CT-FM - Cell Type Fine-Mapping
CT-FM is a statistical genetics method to identify causal cell types underlying diseases and complex traits. We leverage cell type specific candidate cis-regulatory elements (cCREs) together with GWAS summary statistics while accounting for co-regulation of cis-regulatory effects across diverse cell types.

Contacts:<br /> 
Artem Kim artemkim@usc.edu<br />
Steven Gazal gazal@usc.edu<br />

# Installation

## Pre-requisites <br />

Make sure to have anaconda/conda and mamba installed in your environment <br />
https://www.anaconda.com/download <br />
https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html

## 1. Clone the repository <br />


```bash
git clone https://github.com/ArtemKimUSC/CTFM/
```


## 2. Install conda environment

```bash
cd CTFM
mamba env create -f install/ctfm.yml
```

## 3. Test conda environment
Once created, try to activate the conda environment, test LDSC build and the susieR library to make sure there are no errors.

```bash
conda activate ctfm
python ~/bin/ldsc/ldsc.py -h
R
library(susieR)
q()
conda deactivate
```

## 4. Download reference files for CT-FM
To run CT-FM, you need to download the 1000 Genomes reference files and default S-LDSC annotations.<br />
**Warning : the download size is big (~20 GB)**

```bash
bash install/download_reference.sh
```

The script will download all necessary files and put them in the `CTFM/data/` folder

## 5. (optional) Download reference files for CT-FM-SNP

You will need the bed files of the initial annotations to run CT-FM-SNP

```bash
bash install/download_ctfmsnp_reference.sh
```
**Warning : the download size is big (~850 Mb archive, ~17 GB uncompressed)**



# Running the default CT-FM pipeline

The following workflow uses the CT-FM pipeline with 927 cCRE annotations of diverse human cell types, as it was performed in the CT-FM paper : **arxiv link**

## 1. Process GWAS summary statistics
The goal here is to take your GWAS sumstats and convert it to a hg19 S-LDSC friendly format using the script `1_create_sumstats.pl`
The code below assumes that it is launched from the root directory and that the unprocessed GWAS sumstats is gzipped and already located in the sumstats_in/ folder.
```bash
conda activate ctfm

perl scripts/1_create_sumstats.pl --filein sumstats_in/YOUR_SUMSTATS \   # precise the name of your GWAS sumstats omitting the .gz extension
 --fileout sumstats_ready/YOUR_NEW_SUMSTATS \ # output sumstats, put it in the sumstats_ready/ directory for downstream analyses
 --chr  \
 --pos  \
 --beta  \
 --se  \
 --A1  \
 --A2  \
 --hg  \
  --N

conda deactivate
```

For this script, we need to precise:<br />
-the input GWAS summary statistics file (the one you put in the sumstats_in folder): `--filein` argument <br />
-the output processed GWAS summary statistics file (put it in the sumstats_ready folder): `--fileout` argument<br />
-the columns (**column numbers**) corresponding to the chromosome (`--chr`), position (`--pos`), beta (`--beta`) and standard error (`--se`) estimates, alleles A1 (`--A1`) and A2 (`--A2`)<br />
-the genome version of the GWAS summary statistics (`--hg 38` or `--hg 19`) and the samplesize column number (`--N`)<br />


### Additional arguments: <br />

The script was also designed to take as input other GWAS sumstats information, for example: <br /> 

**Z-scores**: add `--Z` argument (column number) and remove `--beta` and `--se` arguments <br />
**OR values**: add `--OR` argument (column number) and remove `--beta` and `--se` arguments <br />
**p-values**: add `--pval` argument (column number) and remove `--beta` and `--se` arguments. Note: this option will convert P values into unsigned Z scores, which allows heritability partitioning analyses but not genetic correlation analyses. If OR is available, you can precise both `--p-val` and `--OR` arguments to estimate signed effects. <br />
**total sample size of GWAS**: add `--myN` argument (total samplesize number) and remove `--N` argument <br />
**total sample size of cases and controls**: add `--N1` (number of cases) and `--N2` (number of controls) and remove `--N` / `--myN` arguments <br />

## 2. Run S-LDSC with 927 cCRE annotations

We will use the `2_launch_SLDSC_default.sh` script which will launch 5 analyses (running time ~10-15 hours for each) to run S-LDSC on your GWAS sumstats with 927 different cCRE annotations. <br />
The main code for S-LDSC analysis is stored in `scripts/sldsc.sh`<br />
The script to launch is `scripts/2_launch_SLDSC.sh`<br /> => **If you work on a cluster with SLURM-like job submissions, I recommend you modify this script to launch 5 analyses as 5 separate jobs.** <br />

```bash
conda activate ctfm
bash scripts/2_launch_SLDSC.sh $YOUR_NEW_SUMSTATS_NAME  # precise the name of your sumstats file omitting the "sumstats.gz" part
conda deactivate
```


The output will be stored in `CTFM/out/SLDSC/YOUR_NEW_SUMSTATS/`<br />


## 3. Launch fine-mapping for your S-LDSC results 
We will use the R script `3_launch_susie.R` to perform SuSiE fine-mapping for the analyzed GWAS (`$YOUR_NEW_SUMSTATS` argument) and 927 S-LDSC annotations.<br />

```bash
conda activate ctfm
Rscript scripts/3_launch_susie.R $directory $YOUR_NEW_SUMSTATS #precise the work directory in which CT-FM was downloaded and the name of your sumstats file
conda deactivate
```

=> The results are composed of 3 files containing the initial S-LDSC results (`_CTFM_sldsc.txt`), CT-FM PIP values (`_CTFM_pips.txt`) and Credible sets (`_CTFM_CS.txt`) -> stored in `out/susie/`

# Running the default CT-FM-SNP pipeline

The following workflow uses the CT-FM-SNP pipeline with 927 cCRE annotations of diverse human cell types, as it was performed in our paper : **arxiv link** <br />

Make sure you have downloaded the necessary bed files (step 5 of **Installation**)


## 1. Get overlapping annotations 

We will use the script `4_CTFMSNP_default_annots.sh` to retreive, for each SNP, the overlapping annotations. The script takes as input a list of SNPs to analyze in hg19 bed format, as shown in the example below:<br />

`chr1	100818727	100818728	rs17420882`

Importantly, the 4th column needs to correspond to the name of the SNP (will be used for output). You also need to precise the output directory:

```bash
conda activate ctfm
bash scripts/4_CTFMSNP_default_annots.sh $SNP_LIST.bed $OUTPUT_DIRECTORY
conda deactivate
```
The script will generate, for each SNP, a file `$SNP.out` containing overlapping S-LDSC annotations.


## 2. Run CT-FM-SNP on S-LDSC results 
For each SNP, the script `5_CTFMSNP_default_susie.R` runs fine-mapping using S-LDSC results of overlapping annotations for the trait of interest.

```bash
conda activate ctfm
Rscript scripts/5_CTFMSNP_default_susie.R $WORKING_DIRECTORY $TRAIT $SNP.out
conda deactivate
```

The script requires 3 arguments:
`$WORKING_DIRECTORY` - the directory CT-FM was downloaded into <br />
`$TRAIT` - the trait of interest, must have the same name as the directory with S-LDSC results in `out/SLDSC/` <br />
`$SNP.out` - the SNP to analyze, genereated in step 1. precise the full path to the file <br />


The results will be stored in `out/ctfmsnp/$TRAIT/$SNP_pips.txt` and `out/ctfmsnp/$TRAIT/$SNP_CS.txt` - PIP values for each overlapping annotation and the inferred credible sets.

# Running CT-FM with custom annotations

In development.




