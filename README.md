Analysis of Metagenomic Species (MGS)
=====================================

Scripts used for characterizing metagenome-assembled genomes (MAGs).

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

## Other analysis and plotting scripts

<b>scripts/</b>
* `concat2tree.py`: Concatenate protein sequence alignments and build tree with either RAxML or FastTree.
* `count_taxa.py`: Take a taxonomy.tab file and count the number of genomes per taxon.
* ... remaining scripts are part of the `mashdiff.sh` and `checkm_assessment.sh` pipelines.

<b>R/</b>
* `aai-cogs_extfig4c.R` : Assess distribution of average amino acid identities between marker genes.
* `antismash_extfig7.R` : Plot counts of biosynthetic gene clusters (BGCs) calculated with antiSMASH.
* `bwa_geo-prevalence_fig4a.R` : Characterize geographic distribution of MGS.
* `bwa_pan-metagenome_fig4d.R` : Build accumulation curve of the number of species as a function of samples.
* `bwa_prev_fig2b_extfig6c.R` : Determine overall prevalence of MGS.
* `bwa_thresholds_extfig6.R` : Define thresholds for species presence/absence based on BWA results.
* `funcs_phy-assoc_fig5b.R` : Find functions (GO slim terms) differentially abundant between two sets of species.
* `gprop_pca_fig5a.R` : Perform Principal Component Analysis (PCA) using Genome Properties.
* `kegg-cats_extfig8b.R` : Calculate proportion of KEGG functions differentially abundant.
* `mags-cluster_extfig4.R` : Cluster MAGs based on Mash distances.
* `mags-quality_extfig2.R` : Plot CheckM quality scores.
* `mashdiff-counts_fig1b.R` : Evaluate mashdiff results.
* `mashdiff-hist_scatter_fig1a.R` : Scatterplots and histograms to visualize mashdiff results.
* `mgs-quality_extfig5.R` : Evaluate MGS quality scores.
* `phylo-diversity_fig3b.R` : Calculate phylogenetic diversity from a newick tree file.
* `read-class_fig4c.R` : Analyse sourmash classification results.
* `taxa_counts_fig2a.R` : Stacked plots of taxa proportions.
* `virfinder_analysis.R` : Detect viral contigs in a FASTA file using VirFinder.
