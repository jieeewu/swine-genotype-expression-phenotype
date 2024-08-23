# swine-genotype-expression-phenotype

## TWAS Prediction Model Usage Guide

To facilitate the use of pre-built TWAS prediction models, we provide the necessary shell code along with the corresponding data files.

### Requirements

- **Python 3.5 or higher**: Ensure you have Python 3.5 or a higher version installed.
- **Virtual Environment**: Set up a virtual environment according to the [MetaXcan](https://github.com/hakyimlab/MetaXcan).

### Input Files

1. **Model Database Files**: Download the `*.db` files for the corresponding tissues from our website.
2. **Genotype VCF Files**: Phased genotype data in VCF format.
3. **Phenotype Files**: The phenotype file should contain columns in the following order: `FID`, `IID`, `pheno1`, `pheno2`, ...

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

