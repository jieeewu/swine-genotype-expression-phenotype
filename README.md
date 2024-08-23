# swine-genotype-expression-phenotype

## TWAS Predictive Model Usage Guide

We provide the necessary shell scripts and corresponding data to facilitate the direct use of pre-built TWAS predictive models.

### Note: PrediXcan is Suitable for Individual-Level Data

### Prerequisites

Ensure you have Python 3.5 or higher. Set up a virtual environment according to the [MetaXcan](https://github.com/hakyimlab/MetaXcan) requirements:

- numpy (>=1.11.1)
- scipy (>=0.18.1)
- pandas (>=0.18.1)
- patsy (>=0.5.0)
- statsmodels (>=0.8.0)
- h5py (>=2.7.1)
- bgen_reader (>=3.0.3)
- cyvcf2 (>=0.8.0)

### Input Files

1. **Database Files**: Download the corresponding tissue-specific `*.db` files from our site as per the [PredictDB tutorial](https://github.com/hakyimlab/PredictDB-Tutorial).
2. **Genotype Files**: Phased genotype files in VCF format, which can be either WGS data or imputed genotype data (imputation server：[SWIM](http://106.13.12.181:9088/#/home)).
3. **Phenotype Files**: Files containing phenotype data for association analysis, with columns for FID, IID, pheno1, pheno2, etc.

### Output Files

1. `${tissue}_prediction_output.txt`
2. `${tissue}_prediction_summary_output.txt`
3. `${tissue}_${pheno}_association.txt`

### Running the Shell Script

```bash
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
--input_phenos_file ${phenos} \
--input_phenos_column $pheno \
--output ${Predixcan}/PrediXcanAssociation/${tissue}/${tissue}_${pheno}_association.txt \
--verbosity 9 \
--throw
done
done

