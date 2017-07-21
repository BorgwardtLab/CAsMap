#!/bin/bash
# Uncomment both for selective regression tests, e.g. only exact non-WY method
#STEM_WY=""; STEM_NONWY="exact"

ROOT=$(cd "$(dirname $0)/.."; pwd) # OS X and BSD compatible absolute path
OUT_DIR="${ROOT}/test/output/plink"
PRE_RM_OUT_DIR=false
ECHO_CMD=false
EXPECTED_OUT_DIR="${ROOT}/test/expected_output/plink"



####
# FAIS

INPUT_DIR="${ROOT}/data/fais/test_plink_big"
DATA_OPT="-plink ${INPUT_DIR}/sample_data"
# Uncomment to test if ETH format gives same results
#DATA_OPT="-eth ${INPUT_DIR}/data.txt ${INPUT_DIR}/label.txt"
source "${ROOT}/test/_regression_test.sh"



####
# FastCMH

STEM_WY=""; STEM_NONWY="fastcmh"

INPUT_DIR="${ROOT}/data/fastcmh/test_plink"
DATA_OPT="-plink ${INPUT_DIR}/sample_data ${INPUT_DIR}/sample_covariate.txt"
source "${ROOT}/test/_regression_test.sh"
