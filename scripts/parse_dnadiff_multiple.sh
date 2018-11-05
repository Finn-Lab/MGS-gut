#!/bin/bash

if [ $# -eq 0 ]; then
    echo "usage: parse_dnadiff_workflow.sh bin_folder/ _DB [suffix/prefix] use_prefix"
    exit 1
fi

for x in ${1}/*${2}.report
do 
    echo ${x}
    parse_dnadiff.py ${x} ${1}/bestMash_hits${2}.tab ${2} ${3} ${4} > ${x%%.report}_parsed.tab
done
