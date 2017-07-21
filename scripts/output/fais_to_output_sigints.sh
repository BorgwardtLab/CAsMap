#!/bin/bash
FAIS=$1
if [ -z "${FAIS}" ]; then
    echo "Error: no fais file given"
    exit 1
fi

FAIS_EXPECTED_HEADER="l,tau,a,x,P-value"
FAIS_HEADER=$(head -n 1 ${FAIS})
# check header
if [ "${FAIS_HEADER}" != "${FAIS_EXPECTED_HEADER}" ]; then
    echo "Error: expected header \"${FAIS_EXPECTED_HEADER}\"; are you sure it's a fais file?"
    exit 1
fi

echo "start,end,a,x,P-value"
for row in $(tail -n +2 ${FAIS}); do
	length=$(echo ${row} | cut -f 1 -d ",")
	start=$(echo ${row} | cut -f 2 -d ",")
	rest=$(echo ${row} | cut -f 3,4,5 -d ",")
	echo ${start},$((start+length-1)),${rest}
done