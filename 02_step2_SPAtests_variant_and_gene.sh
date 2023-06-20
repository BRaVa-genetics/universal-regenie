#!/bin/bash

source ./run_container.sh

POSITIONAL_ARGS=()

SINGULARITY=false
OUT="out"
TESTTYPE=""
PLINK=""
VCF=""
MODELFILE=""
GROUPFILE=""

saige_version="1.1.8"

while [[ $# -gt 0 ]]; do
  case $1 in
    -o|--outputPrefix)
      OUT="$2"
      shift # past argument
      shift # past value
      ;;
    --chr)
      CHR="$2"
      shift # past argument
      shift # past value
      ;;
    --testType)
      TESTTYPE="$2"
      if ! ( [[ ${TESTTYPE} = "variant" ]] || [[ ${TESTTYPE} = "group" ]] ); then
        echo "Test type is not in {variant,group}"
        exit 1
      fi
      shift # past argument
      shift # past value
      ;;
    -s|--isSingularity)
      shift # past argument
      shift # past value
      ;;
    -p|--plink)
      PLINK="$2"
      shift # past argument
      shift # past value
      ;;
    --vcf)
      VCF="$2"
      shift # past argument
      shift # past value
      ;;
    -m|--modelFile)
      MODELFILE="$2"
      shift # past argument
      shift # past value
      ;;
    -g|--groupFile)
      GROUPFILE="$2"
      shift # past argument
      shift # past value
      ;;
    --annotations)
      ANNOTATIONS="$2"
      shift # past argument
      shift # past value
      ;;
    --subSampleFile)
      SUBSAMPLES="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      echo "usage: 02_step2_SPAtests_variant_and_gene.sh
  required:
    --testType: type of test {variant,group}.
    -p,--plink: plink filename prefix of bim/bed/fam files. These must be present in the working directory at ./in/plink_for_vr_bed/
    --vcf vcf exome file. If a plink exome file is not available then this vcf file will be used. These must be present in the working directory at ./in/vcf/
    --modelFile: filename of the model file output from step 1. This must be in relation to the working directory.
	--chr: chromosome to test.
  optional:
    -o,--outputPrefix:  output prefix of the SAIGE step 2 output.
    -s,--isSingularity (default: false): is singularity available? If not, it is assumed that docker is available.
    -g,--groupFile: required if group test is selected. Filename of the annotation file used for group tests. This must be in relation to the working directory.
    --annotations: required if group test is selected. comma seperated list of annotations to test found in groupfile.
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
if [[ ${TESTTYPE} == "" ]]; then
  echo "Test type not set"
  exit 1
fi

if [[ ${PLINK} == "" ]] && [[ ${VCF} == "" ]]; then
  echo "plink files plink.{bim,bed,fam} and vcf not set"
  exit 1
fi

if [[ ${MODELFILE} == "" ]]; then
  echo "model file not set"
  exit 1
fi

if [[ $GROUPFILE == "" ]] && [[ ${TESTTYPE} == "group" ]]; then
  echo "attempting to run group tests without an annotation file"
  exit 1
fi

if [[ $ANNOTATIONS == "" ]] && [[ ${TESTTYPE} == "group" ]]; then
  echo "attempting to run group tests without selected annotations"
  exit 1
fi

if [[ $SUBSAMPLES != "" ]]; then
  SUBSAMPLES="${HOME}/${SUBSAMPLES}"
fi

if [[ $OUT = "out" ]]; then
  echo "Warning: outputPrefix not set, setting outputPrefix to 'out'. Check that this will not overwrite existing files."
fi

echo "OUT               = ${OUT}"
echo "TESTTYPE          = ${TESTTYPE}"
echo "SINGULARITY       = ${SINGULARITY}"
echo "PLINK             = ${PLINK}.{bim/bed/fam}"
echo "MODELFILE         = ${MODELFILE}"
echo "GROUPFILE         = ${GROUPFILE}"
echo "ANNOTATIONS"      = ${ANNOTATIONS}

# For debugging
set -exo pipefail

## Set up directories
WD=$( pwd )

# Get number of threads
n_threads=$(( $(nproc --all) - 1 ))

## Set up directories
WD=$( pwd )


if [[ "$TESTTYPE" = "variant" ]]; then
  echo "variant testing"
  min_mac=0.5
  GROUPFILE=""
else
  echo "gene testing"
  min_mac=20
  GROUPFILE="${HOME}/${GROUPFILE}"
fi

if [[ ${PLINK} != "" ]]; then
  PLINK="${HOME}/${PLINK}"
  BED=${PLINK}".bed"
  BIM=${PLINK}".bim"
  FAM=${PLINK}".fam"
  VCF=""
elif [[ ${VCF} != "" ]]; then 
  BED=""
  BIM=""
  FAM="" 
  VCF="${VCF}"
else
  echo "No plink or vcf found!"
  exit 1
fi

awk '{print $1, $1}' ${SAMPLEIDS} > sampleids
sed -i '1i FID IID' sampleids

IFS=',' read -r -a array <<< "$COVARCOLLIST"

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
}
{
  printf "%s %s", $1, $2
  for(i in col){
    printf " %s", $col[i]
  }
  print ""
}
' ${HOME}/${PHENOFILE} > covariates.tsv

IFS=',' read -r -a array <<< "$PHENOCOL"

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
}
{
  printf "%s %s", $1, $2
  for(i in col){
    printf " %s", $col[i]
  }
  print ""
}
' ${HOME}/${PHENOFILE} > phenotypes.tsv

cmd="regenie \
  --step 2 \
  --bed $PLINK \
  --phenoFile ${HOME}/phenotypes.tsv \
  --phenoCol ""${PHENOCOL}"" \
  --covarFile ${HOME}/covariates.tsv \
  --covarColList ""${COVARCOLLIST}"" \
  --catCovarList=""${CATEGCOVARCOLLIST}"" \
  --keep ${HOME}/sampleids \
  --firth --approx --pThresh 0.1 \
  --bsize 400 \
  --pred ${HOME}/${MODELFILE} \
  --threads 4 \
  --out ${HOME}/${OUT}
"

echo "Running variant based tests for all variants in with MAC > 20"

run_container
