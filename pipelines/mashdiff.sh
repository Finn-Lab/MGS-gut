#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

This script runs Mash between a set of bins and a reference database. \
It then uses the best matching reference to perform whole-genome alignment with dnadiff.

OPTIONS:
   -i      Absolute path of input directory containing bins (run ID must be immediately before 'metaspades')
   -r      Reference file to use for mash (.msh)
   -s      Database suffix to add
   -p      Prefix for output files ("use_prefix" if using the metaspades file structure) 

EXAMPLE:
./script.sh -i bin_folder/ -r hgr_files.msh -s HGR -p use_prefix
EOF
}

bins=
db_file=
db=
prefix=

while getopts “i:r:s:p:” OPTION
do
     case ${OPTION} in
         i)
             bins=${OPTARG}
             ;;
         r)
             db_file=${OPTARG}
             ;;
         s)
             db=${OPTARG}
             ;;
         p)
             prefix=${OPTARG}
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

if [[ -z ${bins} ]] || [[ -z ${db_file} ]] || [[ -z ${db} ]]
then
     echo "ERROR : Please supply required arguments"
     usage
     exit 1
fi

bsub -M 5000 -n 8 -J "mash_dist_${bins}" -o ${bins}/mash_dist_${db}.tab "mash dist -p 8 ${db_file} ${bins}/*.fa"

bsub -o /dev/null -J "bestMash_${bins}" -w "ended(mash_dist_${bins})" "bestMash.py ${bins}/mash_dist_${db}.tab > ${bins}/bestMash_hits_${db}.tab"

bsub -o ${bins}/dnadiff_${db}.log -M 20000 -w "ended(bestMash_${bins})" "dnadiff_multiple.sh ${bins}/bestMash_hits_${db}.tab ${db}; parse_dnadiff_multiple.sh ${bins} _${db} suffix ${prefix}"
