#!/bin/bash

# set -x

if [ -z ${OUT_DIR} ]; then
    echo "Error: set env var OUT_DIR with a path to a directory where files will be output (note: if dir exists it will be removed first)."
    exit 1
fi
if [ -z ${EXPECTED_OUT_DIR} ]; then
    echo "Error: set env var EXPECTED_OUT_DIR with a path to a directory where expected output files are."
    exit 1
fi
if [[ ${DATA_OPT} != "-eth "* ]] && [[ ${DATA_OPT} != "-plink "* ]]; then
    echo "Error: set env var DATA_OPT with a -eth or -plink option and input data paths."
    exit 1
fi
if [ -z ${ROOT} ]; then
    echo "Error: set env var ROOT with a path to a main repo directory."
    exit 1
fi


if [ -z ${PRE_RM_OUT_DIR} ]; then
    PRE_RM_OUT_DIR=true
fi
if ${PRE_RM_OUT_DIR}; then
    if [ -d "${OUT_DIR}/" ]; then
        rm -rf "${OUT_DIR}/"
    fi
fi
mkdir -p "${OUT_DIR}"


# Note: ${..+x} is necessary since STEM_* might be an empty string
if [ -z ${STEM_NONWY+x} ]; then
    STEM_NONWY="exact chi"
fi
if [ -z ${STEM_WY+x} ]; then
    STEM_WY="exact chi"
fi
if [ -z ${STEM_ISET+x} ]; then
    STEM_ISET=""
fi

if [ -z ${LMAX} ]; then
    LMAX=0
fi
if [ -z ${ALPHA} ]; then
    ALPHA=0.05
fi
if [ -z ${NPERM} ]; then
    NPERM=50
fi
if [ -z ${SEED} ]; then
    SEED=100
fi


if [ -z ${ECHO_CMD} ]; then
    ECHO_CMD=false
fi


diff_output() {
    local file_pref=$1
    for file_path in $(ls "${EXPECTED_OUT_DIR}/${file_pref}"*); do
        file=$(basename ${file_path})
        if [ -f "${OUT_DIR}/${file}" ]; then
            echo "diff \"${file}\""
            diff "${OUT_DIR}/${file}" "${EXPECTED_OUT_DIR}/${file}"
        else
            echo "diff fail: missing output file \"${file}\""
        fi
    done
}

eval_error() {
    local log_path=$1
    cat $log_path
    exit 1
}

for stem in ${STEM_NONWY}; do
    out_pref="output_${stem}";
    log_file="last_run_${stem}.log"
    log_path="${OUT_DIR}/${log_file}"
    echo run significant_interval_search_${stem}
    cmd="\"${ROOT}/executables/significant_interval_search_${stem}\" ${DATA_OPT} ${ALPHA} ${LMAX} \"${OUT_DIR}/${out_pref}\" > \"${log_path}\""
    ${ECHO_CMD} && echo ${cmd}
    eval ${cmd} || eval_error ${log_path}
    diff_output $out_pref
    diff_output $log_file
done
for stem in ${STEM_WY}; do
    out_pref="output_wy_${stem}";
    log_file="last_run_wy_${stem}.log"
    log_path="${OUT_DIR}/${log_file}"
    echo run significant_interval_search_wy_${stem};
    cmd="\"${ROOT}/executables/significant_interval_search_wy_${stem}\" -seed ${SEED} ${DATA_OPT} ${ALPHA} ${NPERM} ${LMAX} \"${OUT_DIR}/${out_pref}\" > \"${log_path}\""
    ${ECHO_CMD} && echo ${cmd}
    eval ${cmd} || eval_error ${log_path}
    diff_output $out_pref
    diff_output $log_file
done
for stem in ${STEM_ISET}; do
    out_pref="output_iset_${stem}";
    log_file="last_run_iset_${stem}.log"
    log_path="${OUT_DIR}/${log_file}"
    echo run significant_itemset_search_${stem};
    cmd="\"${ROOT}/executables/significant_itemset_search_${stem}\" ${DATA_OPT} ${ALPHA} ${LMAX} \"${OUT_DIR}/${out_pref}\" > \"${log_path}\""
    ${ECHO_CMD} && echo ${cmd}
    eval ${cmd} || eval_error ${log_path}
    diff_output $out_pref
    diff_output $log_file
done
