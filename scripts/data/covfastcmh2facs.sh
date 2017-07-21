#!/bin/bash
#
# Copyright 2017 ETH Zurich
#
# Licensed under the GNU General Public License v3.0
# https://www.gnu.org/licenses/gpl.html
#

v=0
for n in $(cat $1); do
    for i in $(seq $n); do
        echo $v
    done
    v=$((v+1))
done
