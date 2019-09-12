Analysis of Metagenomic Species (MGS)
=====================================

Scripts used for characterizing metagenome-assembled genomes (MAGs) used in the following publication:

A Almeida, AL Mitchell, M Boland, SC Forster, GB Gloor, A Tarkowska, TD Lawley and RD Finn (2019) [A new genomic blueprint of the human gut microbiota](https://www.nature.com/articles/s41586-019-0965-1). <i>Nature</i> <b>568</b>, 499–504

Associated data can also be found in: ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/umgs_analyses/

## mashdiff.sh

Compare a set of genomes against a reference database. Selects best representative for whole-genome alignment.

<b>Requirements:</b>
* Mash (tested v2.0)
* MUMmer (tested v3.23)
* Python 2.7

Coded for running within LSF cluster environments. 

<b>Usage:</b> 
```
mashdiff.sh -i genome_folder/ -r reference.msh -s db_name -p output_prefix
```
Arguments:  
`-i` folder containing the genomes to analyse in FASTA format  with `.fa` extension  
`-r` reference file `.msh` generated with `mash sketch`  
`-s` user-defined name for the database (e.g. refseq)  
`-p` user-defined prefix to label the query genomes in the output (e.g. gut)  

<b>Notes:</b>
- The `scripts/` directory needs to be part of your `$PATH` system variable
- The output for each query genome is a `*dbname_parsed.tab` file containing the dnadiff and Mash results. Column headers in the resultant file are: 
```
Query name / Reference name / Ref length / % Ref covered / Query length / % Query aligned / ANI / Mash distance
```

## checkm_assessment.sh

Runs the CheckM `lineage_wf` workflow with the recommended `tree_qa` step for more detailed taxonomic assignment.

<b>Requirements:</b>
* CheckM (tested v1.0.7-1.0.10)
* Python 2.7

Coded for running within LSF cluster environments. 

<b>Usage:</b>
```
checkm_assessment.sh genome_folder/ fa output_prefix
```
Positional arguments:  
1: folder containing the genomes to analyse in FASTA format  
2: extension of the FASTA files to be analysed in the `genome_folder/`  
3: user-defined prefix to label the query genomes in the output (e.g. gut)  

<b>Notes:</b>
- The `scripts/` directory needs to be part of your `$PATH` system variable
- Output is a `checkm_parsed.tab` file with the taxonomy results from `tree_qa` combined with the quality scores determined with `lineage_wf`

## map2ref.sh

Map metagenomic reads against a genome database using BWA.

<b>Requirements:</b>
* BWA (tested v0.7.16a-r1181)
* Samtools (tested v1.5)
* Python 2.7

Coded for running within LSF cluster environments. 

<b>Usage:</b>
```
map2ref.sh input_1.fastq(gz) input_2.fastq(gz) ref-db.fasta out_prefix
```
Positional arguments:  
1 and 2: Forward and reverse reads of metagenome to query (fastq or fastq.gz)  
2: Genome database indexed with `bwa index`, where FASTA headers are in the following structure: `>genome-name_1...`  
3: user-defined output prefix to save output files (e.g. results/metagenome1) 

<b>Notes:</b>
- The `scripts/` directory needs to be part of your `$PATH` system variable
- Results will be stored in the `out_prefix_ref-db_total.tab` and `out_prefix_ref-db_unique.tab` files. The former contains the counts/coverage/depth/variation for all reads mapped per genome, while the latter only takes into account the uniquely mapped reads.

## Other analysis and plotting scripts

<b>scripts/</b>
* `concat2tree.py`: Concatenate protein sequence alignments and build tree with either RAxML or FastTree.
* `count_taxa.py`: Take a taxonomy tabular file and count the number of genomes per taxon.
* `rename_multifasta_prefix.py`: Rename a multi-fasta file based on a user-defined prefix.
* ... remaining scripts are part of the `mashdiff.sh`, `checkm_assessment.sh` and `map2ref.sh` pipelines.

<b>R/</b>
* `antismash_extfig8.R` : Plot counts of biosynthetic gene clusters (BGCs) calculated with antiSMASH.
* `bwa_geo-prevalence_fig4a.R` : Characterize geographic distribution of MGS.
* `bwa_pan-metagenome_fig4d.R` : Build accumulation curve of the number of species as a function of samples.
* `bwa_prev_fig2b_extfig7c.R` : Determine overall prevalence of MGS.
* `bwa_thresholds_extfig7.R` : Define thresholds for species presence/absence based on BWA results.
* `funcs_phy-assoc_fig5b.R` : Find functions (GO slim terms) differentially abundant between two sets of species.
* `gprop_pca_fig5a.R` : Perform Principal Component Analysis (PCA) using Genome Properties.
* `kegg-cats_extfig9b.R` : Calculate proportion of KEGG functions differentially abundant.
* `mags-cluster_extfig5.R` : Cluster MAGs based on Mash distances.
* `mags-quality_extfig2.R` : Plot CheckM quality scores.
* `mashdiff-counts_fig1b.R` : Evaluate mashdiff results in terms of total reference matches.
* `mashdiff-hist_scatter_fig1a.R` : Scatterplots and histograms to visualize mashdiff results.
* `mgs-quality_extfig6.R` : Evaluate MGS quality scores.
* `phylo-diversity_fig3b.R` : Calculate phylogenetic diversity from a newick tree file.
* `read-class_fig4c.R` : Analyse sourmash classification results.
* `taxa_counts_fig2a.R` : Stacked plots of taxa proportions.
* `virfinder_analysis.R` : Detect viral contigs in a FASTA file using VirFinder.
