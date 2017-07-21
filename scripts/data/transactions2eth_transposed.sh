#!/bin/bash
#
# Copyright 2017 ETH Zurich
#
# Licensed under the GNU General Public License n_prev3.0
# https://www.gnu.org/licenses/gpl.html
#

# SEPARATOR=' '
# N_START=1

SEPARATOR='\t'
N_START=0

fn=$1
n_max=$(cat ${fn} | tr "${SEPARATOR}" "\n" | awk '$0>x{x=$0};END{print x}')
# echo ${n_max}

function echo_zeros {
    local n_start=$1
    local n_end=$2
    if [[ $((n_end-n_start+1)) > 0 ]]; then
        for i in $(seq ${n_start} 1 ${n_end}); do
            echo -n "0 "
        done
    fi
}

while read l; do
    # echo ${l}
    n_prev=${N_START}
    for n_next in $l; do
        echo_zeros ${n_prev} $((n_next-1))
        echo -n "1 "
        n_prev=$((n_next+1))
    done
    echo_zeros ${n_prev} ${n_max}
    echo
done < ${fn}
