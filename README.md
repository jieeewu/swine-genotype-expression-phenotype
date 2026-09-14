# Swine Genotype Expression Phenotype (SGEP)

## Overview

This repository contains the scripts and workflow used in our study to integrate genotype, gene expression, and phenotype data in swine, enabling downstream association analyses between genetic variation, transcriptomic regulation, and phenotypic traits.

The pipeline includes data preprocessing, expression quantification, and genotype–expression–phenotype association analyses, as described in the accompanying manuscript.

## Pipeline Overview
The following diagram illustrates the overall pipeline of the project, including data preprocessing and core processing logic.

![Pipeline Overview](Pipeline/Pipeline.png)



## Overview of statistical analyses

- Data preprocessing and quality control

  - Genotype data

  - RNA-seq data

- Genotype–expression association analyses
      - eQTL mapping and functional characterization analysis
     - Tissue sharing patterns

-  Expression-mediated genetic effects

  -  Mediation analysis
  -  Pathway analysis

- Genotype–phenotype association analyses
  - Genome-Wide Association Study (GWAS)
  - Enrichment analysis of eQTLs within GWAS QTLs

  - Transcriptome-wide association study (TWAS)

- Causal inference and colocalization analyses
  - Mendelian randomization
  - Colocalization analysis
    

## Individual-level PrediXcan Analysis Using Pig TWAS Prediction Models

We provide the necessary shell scripts and corresponding data to facilitate the use of pre-built pig PrediXcan prediction models for transcriptome-wide association studies (TWAS).


### Model Description

These prediction models were trained using pig eQTL data based on the [PredictDB](https://github.com/hakyimlab/PredictDB-Tutorial) framework.

- Species: Sus scrofa
- Genome assembly: Sscrofa11.1
- Model type: Elastic Net
- Tissues:
  - Muscle (Mu)
  - Liver (Li)
  - Abdominal adipose (AF)
  - Backfat (BF)
    

### Prerequisites

Ensure you have Python 3.7 or higher. Set up a Python environment following the requirements of [MetaXcan](https://github.com/hakyimlab/MetaXcan):

numpy
scipy
pandas
statsmodels
patsy
h5py
cyvcf2
bgen_reader

You can also download the `PrediXcanAssociation.py` and `Predict.py` scripts from the [Association_study folder](https://github.com/jieeewu/swine-genotype-expression-phenotype/tree/main/05_Association_study) or [MetaXcan](https://github.com/hakyimlab/MetaXcan).

### Input Files

1. **Predictive Models**: Download the corresponding tissue-specific `*.db` files in the [Tissue_DB](Tissue_DB/) folder (`${tissue}_models.db`). These `*.db` files are obtained using the [PredictDB tutorial](https://github.com/hakyimlab/PredictDB-Tutorial). Each database contains SNP weights for predicting genetically regulated gene expression in a specific tissue.
2. **Genotype Files**: Phased genotype files in VCF format, generated from WGS or imputed genotype data (imputation server：[SWIM](http://106.13.12.181:9088/#/home)). Variant IDs must match the identifiers used in the prediction models. The SNP ID format should be `chrom_position_ref` (e.g., `1_502855_C`; Sscrofa11.1).
3. **Phenotype Files**: `phenotype_name.list` and `${phenos_file}`. The phenotype file should contain FID, IID, and phenotype columns (e.g., pheno1, pheno2, etc.)

### Output Files

1. `${tissue}_prediction_output.txt`: Individual-level genetically predicted gene expression matrix.
2. 2. `${tissue}_prediction_summary_output.txt`: Prediction performance information including number of SNPs used for each gene.
3. `${tissue}_${pheno}_association.txt`: Association between genetically predicted expression and phenotype.

### Running the Shell Script
Individual level PrediXcan: [introduction, tutorials and manual](https://github.com/hakyimlab/MetaXcan/wiki/Individual-level-PrediXcan:-introduction,-tutorials-and-manual)

```bash
# Mu: Muscle
# Li: Liver
# AF: Abdominal adipose
# BF: Backfat

module load GCCcore/8.2.0 Python/3.7.2
mkdir -p ${Predixcan}/Predict_output ${Predixcan}/PrediXcanAssociation 

for tissue in Mu Li AF BF 
do
printf "Predict expression\n\n"
python3 $METAXCAN/Predict.py \
--model_db_path $db/${tissue}_models.db \
--vcf_genotypes ${tissue}_phased_for_Predixcan.vcf.gz \
--vcf_mode genotyped \
--prediction_output ${Predixcan}/Predict_output/${tissue}_prediction_output.txt \
--prediction_summary_output ${Predixcan}/Predict_output/${tissue}_prediction_summary_output.txt \
--verbosity 9 \
--throw

printf "association\n\n"
mkdir -p ${Predixcan}/PrediXcanAssociation/${tissue}

for pheno in `cat phenotype_name.list`
do
python3 $METAXCAN/PrediXcanAssociation.py \
--expression_file ${Predixcan}/Predict_output/${tissue}_prediction_output.txt \
--input_phenos_file ${phenos_file} \
--input_phenos_column $pheno \
--output ${Predixcan}/PrediXcanAssociation/${tissue}/${tissue}_${pheno}_association.txt \
--verbosity 9 \
--throw
done
done


```

## Summary-level TWAS Analysis
S-PrediXcan uses GWAS summary statistics to perform TWAS.

### Summary-level S-PrediXcan

Required files:

- Tissue-specific model database (*.db)
- LD covariance file (*.txt.gz)
- GWAS summary statistics


Please refer to the MetaXcan [S-PrediXcan documentation](https://github.com/hakyimlab/MetaXcan/wiki/S-PrediXcan-Command-Line-Tutorial).

### LD Covariance Files (for S-PrediXcan)
LD covariance files are required only for S-PrediXcan analysis.
The covariance files were calculated using genotype data from the pig eQTL reference population used for model training.

For each tissue-specific model:

- Muscle:
  - Mu_models.db
  - Muscle_cov.txt.gz

- Liver:
  - Li_models.db
  - Liver_cov.txt.gz

These covariance files are required for summary-level TWAS analysis using [S-PrediXcan](https://github.com/hakyimlab/MetaXcan/wiki/S-PrediXcan-Command-Line-Tutorial).

### GWAS Summary Statistics Format

For S-PrediXcan analysis, users should provide GWAS summary statistics containing SNP identifiers, alleles, effect sizes, standard errors, and p-values.

Example:
SNP A1 A2 beta se p
1_502855_C C T 0.0231 0.0102 0.0235
1_503001_A A G -0.0154 0.0087 0.0764
1_503245_G G A 0.0312 0.0121 0.0101

Column description:

- SNP: Variant identifier matching the prediction model (e.g., chrom_position_ref)
- A1: Effect allele
- A2: Non-effect allele
- beta: GWAS effect size for A1 allele
- se: Standard error of beta
- p: Association P-value

### Running S-PrediXcan
  
```bash
#!/bin/bash

GWAS=pig_trait_GWAS.txt

for tissue in Mu Li AF BF
do

echo "Running S-PrediXcan for ${tissue}"

python SPrediXcan.py \
--model_db_path ${tissue}_models.db \
--covariance ${tissue}_cov.txt.gz \
--gwas_file ${GWAS} \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column beta \
--pvalue_column p \
--output ${tissue}_TWAS.txt

done
```

## Code Availability

The source code is available at:
https://github.com/jieeewu/swine-genotype-expression-phenotype

A versioned release of the code used in this study has been archived in Zenodo
and assigned a DOI: https://doi.org/10.5281/zenodo.22275711

## Citation

If you use this pipeline or code in your work, please cite our manuscript and the archived code:

- **Publication:**.
Jie Wu*, Ming Yang*, Zebin Zhang*, Enqin Zheng*, Zhanwei Zhuang, Shenping Zhou, Cineng Xu, Yibin Qiu, Donglin Ruan, Jianping Quan, Rongrong Ding, ..., Wen Huang#, Jie Yang#, Zhenfang Wu#. *Integrative Systems Genetics Analysis Advances Elucidation of the Genetic Basis of Complex Traits in Pigs*. *Nature Communications* (Accepted).

------

## Contact

For questions or issues regarding the pipeline, please open an issue in this repository or contact the corresponding author.

![Visitors](https://api.visitorbadge.io/api/visitors?path=https%3A%2F%2Fgithub.com%2Fjieeewu%2Fswine-genotype-expression-phenotype&label=VISITOR&countColor=%23263759)


