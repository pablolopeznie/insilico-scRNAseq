# insilico-scRNAseq
---------------------------
scRNA-seq analysis pipeline
---------------------------

@author: Pablo López Nieto
@date: 30-12-2025
@email_contact: pablolopeznto@gmail.com

-----------
DESCRIPTION
-----------

This directory contains all the required files and scripts to reproduce 
a complete analysis carried out from GSE271406 single-cell dataset in R language. 

Make sure you configure your R session working directory to this path.

--------
PIPELINE
--------

1.- Open download_createobjects_scdata.R and run the whole script.
	-> It will download the count matrices, feature and barcode info
	from GEO database, for all 13 samples.
	-> Datasets will be read in 10X format and stored into individual
	Seurat objects, which will be returned in individual lists for
	healthy individuals (1 list of 3 objects) and patients (1 list of 10
	elements).

	WARNING: Downloaded raw files will remain in you directory unless you
	explicitly erase them (this pipeline is space-agnostic).

2.- Run CITEseq_patients_QC.Rmd
	-> It will perform QC and preprocessing for all 10 patient samples,
	returning a list with qc-passed seurat objects (patients_CITE_seu_filt.rda)

3.- Run CITEseq_healthydonors_QC+analysis.Rmd
	IMPORTANT DISCLAIMER: This script does not only perform QC and preprocessing
	for healthy donor data, it is able to run a complete CITE-seq analysis for
	this 3-sampled dataset. This was used as a first approximation to the study the
	preexisting cellular composition. 
	Nevertheless, if you are only interested in the integrated analysis (patients 
	+ healthy donors), run all the steps necessary to generate the preprocessed Seurat 
	object list (healthy_CITE_seurat_preprocessed.rda, chunk 15).

4.- Explore integrated_CITEseq_analysis.Rmd and enjoy the integrated analysis! 
You will manage with:

	- Batch effect correction (multisample integration)
	- Cell clustering
	- Combined RNA and ADT expression analysis
	- Annotation of cells with a reference dataset
	- Celltype abundance comparisons
	- Pseudobulk analysis of individual populations between conditions
	- Gene Set Enrichment Analysis between celltypes
	
	and many many plots! :)

*Note*:
Since the integration processes are computationally intensive given the big amount of cells in the analysis, there exists an alternative script (server_scrnaseq.R) which is devoted to be run in a cluster or computing server. It will create a file with the integrated Seurat object (hp_pat_integrated.Rda), which can be directly loaded into the corresponding Rmd to keep on with downstream analysis. In case you have trouble with running the integration directly from the source chunks in the Rmd file, use this script.
