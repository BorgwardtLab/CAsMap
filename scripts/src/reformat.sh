#!/bin/bash
#
# Copyright 2017 ETH Zurich
#
# Licensed under the GNU General Public License v3.0
# https://www.gnu.org/licenses/gpl.html
#

clang-format -i -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 80}" $*
echo reformated files: $*
