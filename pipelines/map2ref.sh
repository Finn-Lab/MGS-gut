#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Map reads against a custom reference database"
    echo ""
    echo "Notes:" 
    echo "- run 'bwa index' on reference FASTA beforehand"
    echo "- headers in FASTA file must include genome name (use rename_multifasta_prefix.py script)"
    echo "" 
    echo "usage: script.sh input_1.fastq(gz) input_2.fastq(gz) ref.fasta outprefix"
    echo ""
    exit 1
fi

# variables
ref=${3}
refix=$(basename ${ref%%.fa*})
reads=${1}
reads2=${2}
outprefix=${4}
readname=$(basename ${outprefix})

# initial mapping
bsub -n 8 -M 40000 -oo ${outprefix}_${refix}_mapping.log -J bwa_${readname}_${refix} \
"bwa mem -t 8 ${ref} ${reads} ${reads2} | samtools view -@ 7 -uS - -o ${outprefix}_${refix}_unsorted.bam"

# sort bam file
bsub -n 8 -M 40000 -o ${outprefix}_${refix}_mapping.log  -J sort_${readname}_${refix} -w "ended(bwa_${readname}_${refix})" \
"samtools sort -@ 7 ${outprefix}_${refix}_unsorted.bam -o ${outprefix}_${refix}_sorted.bam"

# extract unique counts
bsub -n 8 -M 5000 -o ${outprefix}_${refix}_mapping.log  -J unique_${readname}_${refix} -w "ended(sort_${readname}_${refix})" \
"samtools view -@ 7 -q 1 -f 2 -u ${outprefix}_${refix}_sorted.bam -o ${outprefix}_${refix}_unique_sorted.bam"
bsub -n 8 -o ${outprefix}_${refix}_mapping.log  -J index-unique_${readname}_${refix} -w "ended(unique_${readname}_${refix})" \
"samtools index -@ 7 ${outprefix}_${refix}_unique_sorted.bam"
bsub -M 10000 -o ${outprefix}_${refix}_mapping.log  -J depth-unique_${readname}_${refix} -w "ended(index-unique_${readname}_${refix})" \
"samtools idxstats ${outprefix}_${refix}_unique_sorted.bam > ${outprefix}_${refix}_unique_depth.tab; samtools depth ${outprefix}_${refix}_unique_sorted.bam > ${outprefix}_${refix}_unique_depth-pos.tab; parse_bwa-depth.py ${outprefix}_${refix}_unique_depth.tab ${outprefix}_${refix}_unique_depth-pos.tab > ${outprefix}_${refix}_unique.tab; rm -rf ${outprefix}_${refix}_unique_sorted.ba* ${outprefix}_${refix}_unique_depth*"

# extract total counts
bsub -n 8 -o ${outprefix}_${refix}_mapping.log  -J index_${readname}_${refix} -w "ended(sort_${readname}_${refix})" \
"samtools index -@ 7 ${outprefix}_${refix}.bam"
bsub -M 10000 -o ${outprefix}_${refix}_mapping.log  -J depth_${readname}_${refix} -w "ended(index_${readname}_${refix})" \
"samtools idxstats ${outprefix}_${refix}_sorted.bam > ${outprefix}_${refix}_depth.tab; samtools depth ${outprefix}_${refix}_sorted.bam > ${outprefix}_${refix}_depth-pos.tab; parse_bwa-depth.py ${outprefix}_${refix}_depth.tab ${outprefix}_${refix}_depth-pos.tab > ${outprefix}_${refix}_total.tab; rm -rf ${outprefix}_${refix}_unsorted.bam ${outprefix}_${refix}_sorted.bam.bai ${outprefix}_${refix}_depth*"
