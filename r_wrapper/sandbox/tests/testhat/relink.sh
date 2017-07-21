#!/bin/bash

for l in $(find . -type l); do
    fn=$(basename ${l});
    ln -sf $(find ../../../test/plink -name ${fn}) ${fn}
done
