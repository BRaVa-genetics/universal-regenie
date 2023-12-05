# universal-regenie

## Usage

### 01_step1_fitNULLGLMM.sh

```
bash 01_step1_fitNULLGLMM.sh \
    --genotypePlink "ukb_cal_allChrs.bed" \
    --phenoFile in/pheno_list \
    --phenoCol $pheno \
    --covarColList $covariates \
    --categCovarColList "${categorical_covariates}" \
    --sampleIDs in/sample_ids \
    --snp_ids in/snp_ids \
    --sampleIDCol "IID" \
    --outputPrefix ${output_prefix} \
    --isSingularity false \
    --traitType $trait_type
```

- refer to https://rgcgithub.github.io/regenie/options/#input_1 for more input details
- output will be a .list file and a .loco file.

### 02_step2_fitALTGLMM.sh

```
bash 02_step2_SPAtests_variant_and_gene.sh \
    --chr $chrom \
    --bsize $bsize \
    --testType $test_type \
    --phenoType $pheno_type \
    --plink "in/exome" \
    --locoFile in/loco_file \
    --groupFile in/group_file \
    --pheno $pheno \
    --phenoFile in/pheno_file \
    --covariates $covariates \
    --annotations $annotations \
    --outputPrefix $output_prefix \
    --isSingularity false \
    --subSampleFile in/subsample_file

```

- output results will be in a .regenie file.


