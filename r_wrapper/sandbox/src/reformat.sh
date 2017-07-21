#!/bin/bash

clang-format -i -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 80}" $*
echo reformated files: $*
