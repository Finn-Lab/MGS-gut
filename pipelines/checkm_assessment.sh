#!/bin/bash

if [ $# -eq 0 ]; then
    echo ""
    echo "Runs the checkm workflow (tree, tree_qa, lineage_set, analyze and qa)"
    echo "Might need to 'cd' into directory if input path is too long"
    echo ""
    echo "Example usage:"
    echo "checkm_quality.sh bin_folder/ fa prefix"
    exit 1
fi

path=${1}
ext=${2}
prefix=${3}

if [ -d "${path}/tmp" ]
then
    rm -rf "${path}/tmp"
fi

mkdir "${path}/tmp"

# clean output files before running
rm -rf "${path}/checkm_output"
rm -f "${path}/checkm.log"
rm -f "${path}/bins_"*
rm -f "${path}/marker_file"

# checkm tree
bsub -M 85000 -n 8 -o "${path}/checkm.log" -J "checkm_tree_${path}" "checkm tree -t 8 -x ${ext} ${path} ${path}/checkm_output --tmpdir ${path}/tmp"

# checkm tree_qa
bsub -o "${path}/checkm.log" -M 5000 -J checkm_tree_qa_${path} -w "ended(checkm_tree_${path})" "checkm tree_qa ${path}/checkm_output --tab_table -f ${path}/bins_taxonomy.tab --tmpdir ${path}/tmp"

# checkm lineage_set
bsub -M 5000 -o "${path}/checkm.log" -J checkm_lineage_set_${path} -w "ended(checkm_tree_qa_${path})" "checkm lineage_set ${path}/checkm_output/ ${path}/marker_file --tmpdir ${path}/tmp"

# checkm analyze
bsub -M 50000 -n 8 -o "${path}/checkm.log" -J checkm_analyze_${path} -w "ended(checkm_lineage_set_${path})" "checkm analyze -t 8 -x ${ext} ${path}/marker_file ${path} ${path}/checkm_output --tmpdir ${path}/tmp"

# checkm qa
bsub -M 10000 -o "${path}/checkm.log" -J checkm_qa_${path} -w "ended(checkm_analyze_${path})" "checkm qa ${path}/marker_file ${path}/checkm_output --tab_table -f ${path}/bins_qa.tab --tmpdir ${path}/tmp; rm -rf ${path}/tmp; parse_checkm.py ${path}/bins_qa.tab ${path}/bins_taxonomy.tab ${prefix} > ${path}/checkm_results.tab"
