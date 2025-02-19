#!/bin/bash

source ./run_container.sh

POSITIONAL_ARGS=()

SINGULARITY=false
SAMPLEIDCOL="IID"
OUT="out"
TRAITTYPE=""
PLINK=""
PHENOFILE=""
PHENOCOL=""
COVARCOLLIST=""
CATEGCOVARCOLLIST=""
WD=$(pwd)

while [[ $# -gt 0 ]]; do
  case $1 in
    -o|--outputPrefix)
      OUT="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--isSingularity)
      shift # past argument
      shift # past value
      ;;
    -t|--traitType)
      TRAITTYPE="$2"
      if ! ( [[ ${TRAITTYPE} == "quantitative" ]] || [[ ${TRAITTYPE} == "binary" ]] ); then
        echo "Trait type is not in {quantitative,binary}"
        exit 1
      fi
      shift # past argument
      shift # past value
      ;;
    -p|--genotypePlink)
      GENOTYPE_PLINK="$2"
      shift # past argument
      shift # past value
      ;;
    --phenoFile)
      PHENOFILE="$2"
      shift # past argument
      shift # past value
      ;;
    --phenoCol)
      PHENOCOL="$2"
      shift # past argument
      shift # past value
      ;;
    -c|--covarColList)
      COVARCOLLIST="$2"
      shift # past argument
      shift # past value
      ;;
    --categCovarColList)
      CATEGCOVARCOLLIST="$2"
      shift # past argument
      shift # past value
      ;;
    --sampleIDs)
      SAMPLEIDS="$2" 
      shift
      shift
      ;; 
    --snp_ids)
      snp_ids="$2" 
      shift
      shift
      ;; 
    -i|--sampleIDCol)
      SAMPLEIDCOL="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      echo "usage: 01_step1_fitNULLGLMM.sh
  required:
    -t,--traitType: type of the trait {quantitative,binary}.
    --genotypePlink: plink filename prefix of bim/bed/fam files. These must be present in the working directory at ./in/plink_for_vr_bed/
    --sparseGRM: filename of the sparseGRM .mtx file. This must be present in the working directory at ./in/sparse_grm/
    --sparseGRMID: filename of the sparseGRM ID file. This must be present in the working directory at ./in/sparse_grm/
    --phenoFile: filename of the phenotype file. This must be present in the working directory at ./in/pheno_files/
    --phenoCol: the column names of the phenotype to be analysed in the file specified in --phenoFile.
  optional:
    -o,--outputPrefix:  output prefix of the SAIGE step 1 output.
    -s,--isSingularity (default: false): is singularity available? If not, it is assumed that docker is available.
    -c,--covarColList: comma separated column names (e.g. age,pc1,pc2) of continuous covariates to include as fixed effects in the file specified in --phenoFile.
    --categCovarColList: comma separated column names of categorical variables to include as fixed effects in the file specified in --phenoFile.
    --sampleIDCol (default: IID): column containing the sample IDs in the phenotype file, which must match the sample IDs in the plink files.
      "
      shift # past argument
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

# Checks
if [[ ${TRAITTYPE} == "" ]]; then
  echo "traitType not set"
  exit 1
fi

if [[ ${SAMPLEIDS} != "" ]]; then
  SAMPLEIDS=${HOME}/$SAMPLEIDS
fi

if [[ ${PHENOFILE} == "" ]]; then
  echo "phenoFile not set"
  exit 1
fi

if [[ ${PHENOCOL} == "" ]]; then
  echo "phenoCol not set"
  exit 1
fi

if [[ $OUT = "out" ]]; then
  echo "Warning: outputPrefix not set, setting outputPrefix to ${PHENOCOL}. Check that this will not overwrite existing files."
  OUT="${PHENOCOL}"
fi

if [[ $COVARCOLLIST = "" ]]; then
  echo "Warning: no continuous fixed effeect covariates included."
fi

if [[ $CATEGCOVARCOLLIST = "" ]]; then
  echo "Warning: no categorical fixed effeect covariates included."
fi

echo "OUT               = ${OUT}"
echo "SINGULARITY       = ${SINGULARITY}"
echo "TRAITTYPE         = ${TRAITTYPE}"
echo "PLINK             = ${PLINK_WES}.{bim/bed/fam}"
echo "PHENOFILE         = ${PHENOFILE}"
echo "PHENOCOL          = ${PHENOCOL}"
echo "COVARCOLLIST      = ${COVARCOLLIST}"
echo "CATEGCOVARCOLLIST = ${CATEGCOVARCOLLIST}"
echo "SAMPLEIDS         = ${SAMPLEIDS}"
echo "SAMPLEIDCOL       = ${SAMPLEIDCOL}"

if [[ "$PHENOCOL" =~ .*"-".* || "$PHENOCOL" =~ .*",".* || "$PHENOCOL" =~ .*"=".* ]]; then
  echo "Phenotype name cannot contain \"-\" or \",\" or \"=\""
  exit 1
fi

# For debugging
set -exo pipefail

## Set up directories
WD=$( pwd )

# Get number of threads
n_threads=$(( $(nproc --all) - 1 ))

# Get inverse-normalize flag if trait_type=="quantitative"
if [[ ${TRAITTYPE} == "quantitative" ]]; then
  echo "Quantitative trait passed to regenie, perform IRNT"
  INVNORMALISE=TRUE
else
  echo "Binary trait passed to regenie"
  INVNORMALISE=FALSE
fi

bed="${HOME}/${GENOTYPE_PLINK}"
bed="${bed%.*}"

awk '{print $1, $1}' ${SAMPLEIDS} > sampleids
sed -i '1i FID IID' sampleids

awk -v colnames="$COVARCOLLIST" '
BEGIN{ FS=OFS="\t" }
FNR==1{
  split(colnames, cols, ",");
  for(i in cols){
    for(j=1; j<=NF; j++){
      if($j==cols[i]){
        col[i]=j;
      }
    }
  }
  gsub(",", OFS, colnames); print "FID IID", colnames
}
{
  if (FNR > 1) {
    printf "%s %s", $1, $1
    for(i in col){
      if ($(col[i]) == "") {
        printf " NA"
      } else {
        printf " %s", $(col[i])
      }
    }
    print ""
  }
}
' ${HOME}/${PHENOFILE} > covariates.tsv

awk -v colnames="$PHENOCOL" '
BEGIN{ FS=OFS="\t" }
FNR==1{
  split(colnames, cols, ",");
  for(i in cols){
    for(j=1; j<=NF; j++){
      if($j==cols[i]){
        col[i]=j;
      }
    }
  }
  gsub(",", OFS, colnames); print "FID IID", colnames
}
{
  if (FNR > 1) {
    printf "%s %s", $1, $1
    for(i in col){
      if ($(col[i]) == "") {
        printf " NA"
      } else {
        printf " %s", $(col[i])
      }
    }
    print ""
  }
}
' ${HOME}/${PHENOFILE} > phenotypes.tsv

head phenotypes.tsv
head covariates.tsv
head sampleids

# Assuming trait_type is already set in your environment
if [ "$trait_type" == "quantitative" ]; then
    trait_flag="--qt"
else
    trait_flag="--bt"
fi

cmd="""regenie \
          --step 1 \
          --bed $bed \
          --extract "${HOME}/${snp_ids}" \
          --phenoFile ${HOME}/phenotypes.tsv \
          --phenoCol ""${PHENOCOL}"" \
          --covarFile ${HOME}/covariates.tsv \
          --covarColList ""${COVARCOLLIST}"" \
          --catCovarList=""${CATEGCOVARCOLLIST}"" \
          --keep ${HOME}/sampleids \
          --threads=${n_threads} \
          --bsize 1000 \
          $trait_flag \
          --out ${HOME}/$OUT
"""

run_container

ls -alt