#!/bin/bash

if [ $# -eq 0 ]; then
    echo "ERROR! usage: dnadiff_workflow.sh bestMash_hits_db.tab out_suffix"
    exit 1
fi

while read col1 col2 rem
do
    echo "dnadiff ${col2} ${col1} -p ${col1%%.fa*}_${2}" 
    dnadiff ${col2} ${col1} -p ${col1%%.fa*}_${2}
done < ${1}
