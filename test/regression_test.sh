#!/bin/bash

ROOT=$(cd "$(dirname $0)/.."; pwd) # OS X and BSD compatible absolute path
OUT_DIR="${ROOT}/test/output"
PRE_RM_OUT_DIR=false
ECHO_CMD=false
EXPECTED_OUT_DIR="${ROOT}/test/expected_output"



####
# FAIS

# Uncomment both for selective regression tests
#STEM_WY="chi"; STEM_NONWY="chi" # only Chi methods
#STEM_WY="" # only non-WY methods

INPUT_DIR="${ROOT}/data/fais/test"
DATA_OPT="-eth ${INPUT_DIR}/data.txt ${INPUT_DIR}/label.txt"
source "${ROOT}/test/_regression_test.sh"



####
# FastCMH

STEM_WY=""; STEM_NONWY="fastcmh"

#INPUT_DIR="${ROOT}/data/fastcmh/r_demo"
# inconsistent covariates number error
#DATA_OPT="-eth ${INPUT_DIR}/data.txt ${INPUT_DIR}/label.txt ${INPUT_DIR}/cov.dat"
# one covariate output == FAIS Chi output
#DATA_OPT="-eth ${INPUT_DIR}/data.txt ${INPUT_DIR}/label.txt ${INPUT_DIR}/cov0.txt"

INPUT_DIR="${ROOT}/data/fastcmh/test"
DATA_OPT="-eth ${INPUT_DIR}/data.txt ${INPUT_DIR}/label.txt ${INPUT_DIR}/cov.txt"
# Uncomment to test on a larger data set, that finds no significant intervals
#DATA_OPT="-eth ${INPUT_DIR}/data2.txt ${INPUT_DIR}/label.txt ${INPUT_DIR}/cov.txt"
source "${ROOT}/test/_regression_test.sh"



####
# FACS

STEM_WY=""; STEM_NONWY=""; STEM_ISET="facs"

INPUT_DIR="${ROOT}/data/facs/test"
DATA_OPT="-eth ${INPUT_DIR}/tictactoe_data.txt ${INPUT_DIR}/tictactoe_labels.txt ${INPUT_DIR}/tictactoe_covariates.txt"
# Uncomment to test on a larger data set
# EXPECTED_OUT_DIR="${EXPECTED_OUT_DIR}/big"
# OUT_DIR="${OUT_DIR}/big"
# INPUT_DIR="${ROOT}/data/facs/test_eth_big"
# DATA_OPT="-eth ${INPUT_DIR}/X_file.txt ${INPUT_DIR}/Y_file.txt ${INPUT_DIR}/C_file.txt"
# ALPHA=0.1
source "${ROOT}/test/_regression_test.sh"
