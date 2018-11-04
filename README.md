Analysing Metagenomic Species (MGS) from the human gut
============================================

Scripts for characterizing metagenome-assembled genomes (MAGs).

<b>Pipelines:</b>
* mashdiff.sh
* checkm_assessment.sh

Pipelines were written to be run in LSF cluster environments (bsub launched jobs), so will need to be edited if the goal is to run them locally.

## mashdiff.sh

Compare genomes against a reference database. Selects best representative for whole-genome alignment.

<b>Requirements:</b>
* Mash (tested v2.0)
* MUMmer (tested v3.23)
* Python 2.7

<b>Usage:</b> 
> mashdiff.sh -i genome_folder/ -r reference.msh -s db_name -p output_prefix

<b>Notes:</b>
- A database file with the target genomes needs to be first generated with <i>mash sketch</i>
- The output for each genome is a \*db_parsed.tab file containing the dnadiff and Mash results. Column headers in the resultant file are: 
> Query name / Reference name / Ref length / % Ref covered / Query length / % Query aligned / ANI / Mash distance

## checkm_assessment.sh

Runs the CheckM lineage_wf with the recommended qa_tree step.

<b>Analysis and plotting scripts:</b>
* aai-cogs_extfig4c.R
* antismash_extfig7.R
* bwa_geo-prevalence_fig4a.R
* bwa_pan-metagenome_fig4c.R
* bwa_prev_fig2b_extfig6c.R
* bwa_thresholds_extfig6.R
* funcs_phy-assoc_fig5b.R
* gprop_pca_fig5a.R
* kegg-cats_extfig8b.R
* mags-cluster_extfig4.R
* mags-quality_extfig2.R
* mashdiff-counts_fig1b.R
* mashdiff-hist_scatter_fig1a.R
* mgs-quality_extfig5.R
* phylo-diversity_fig3b.R
* read-class_fig4b.R
* taxa_counts_fig2a.R
