# scRNA-seq Analysis Pipeline

## Installation

### Dependencies
- NextFlow >= 23.10.1
- Apptainer(singularity); tested with version 1.2.5-1.el7

### Script and container download
#### Script
Request has to be asked beforehand, or can be run on the server under FOLDERNAME

```
git clone https://github.com/parmejohn/scRNA-pipeline-UofM.git
```

#### Apptainer (Singularity)
Image can be downloaded from [Sylabs](https://cloud.sylabs.io/library/parmejohn/uofm/scrnaseq_singularity) or pulled using 

```
apptainer pull --arch amd64 library://parmejohn/uofm/scrnaseq_singularity:latest
```
Apptainer and singularity are the same, but apptainer is the new name. Also to note, when using apptainer, library may not be set up which could lead to pull errors if using the command above.


## Usage
Run the script from the downloaded directory

```
nextflow run scRNA_pipeline.nf \
	--indir <input folder to sorted cellranger counts> \
	--outdir <output folder> \
	--species <species name> \
	--bind <bind file path locations> \
	--reference_seurat [path(s) to labelled reference seurat object] \
	--beginning_cluster [inferred earliest celltype] \
	--clusters_optimal [optimal number of clusters] \
	--resolution [a chosen resolution] \
	--run_escape [True/False] \
	-with-apptainer path_to_image/scrnaseq_singularity_tmp.sif
```

### Arguments
- indir: Input directory with folders of condition names, which have cellranger count results in them for each sample (REQUIRED)
	- DEFAULT: ''
- outdir: Where to place analysis/ folder (REQUIRED)
	- DEFAULT: ''
- species: Species name; Must lowercase with no spaces eg. musmusculus homosapiens (REQUIRED)
	- DEFAULT: ''
- bind: Paths where input or output data is moving. This is because apptainer(singularity) without root permission needs explicit file paths. Multiple paths can be defined with a ',' to seperate them. (REQUIRED)
	- DEFAULT: ''
- reference_seurat: Path (or paths seperated by ',') to a labelled reference seurat object. The object(s) will be used to automatically label your own dataset. Note, you should still manually check your clusters to see if the identification was what was expected
	- DEFAULT: 'none'
- beginning_cluster: If you know what celltype is would be the earliest in a differentiation pathway.
	- DEFAULT: ''
- clusters_optimal: Number of optimal clusters/dimensionality for machine learning methods. If set to 0, will be calculated automatically through NbClust() Calinski-Harabasz (ch) index.
	- DEFAULT: 0
- resolution: Desired resolution, increasing this will increase the number of clusters and vice versa. WILL IMPLEMENT AUTOMATICALLY THROUGH CLUSTTREE LATER
	- DEFAULT: 1
- run_escape: Performs escape analysis to provide single-cell GSEA results. IMPORTANT: This analysis takes a substantial amount of time and memory, so ensure that you provide enough of both before setting to 'True'. With Nextflow capabilities as well, it is fine to run the pipeline once without escape and then with, since the pipeline will cache the rest of the data.
	- Default: False
- with-apptainer: path to downloaded image (REQUIRED)
	- When this option is not set, you will use your native environment, which may be missing packages, dependecies, or have version mismatches.

### Outputs
Analysis folder:

```
analysis/
├── data
│   ├── optimal_clusters.txt
│   ├── qc
│   │   ├── CONDITION1
│   │   │   ├── [sample_1]_soupx
│   │   │   │   ├── barcodes.tsv
│   │   │   │   ├── genes.tsv
│   │   │   │   └── matrix.mtx
│   │   │   └── ...
│   │   └── ...
│   ├── sce_slingshot.rds
│   ├── sc_integrated_milo_traj.rds
│   ├── se_filtered_list.rds
│   ├── se_filtered_singlets_list.rds
│   ├── se_integrated_auto_label.rds
│   ├── se_integrated_dimred.rds
│   ├── se_integrated_escape_norm.rds
│   ├── se_integrated_escape.rds
│   ├── se_integrated.rds
│   ├── se_list_raw.rds
│   ├── se_markers_presto_integrated.txt
│   └── ti
│       ├── ti_gene_clusters_slingPseudotime_\*.txt
│       ├── ...
└── plots
    ├── conserved_marker_unlabelled.pdf
    ├── da
    │   ├── milo_DA_DE_heatmap_\*.pdf
    │   ├── milo_DA_fc_distribution.pdf
    │   ├── milo_DA_umap.pdf
    │   ├── milo_pval_distribution.pdf
    │   └── milo_volcano_plot.pdf
    ├── deseq2
    │   ├── deseq2_cluster_[CLUSTER]_[CONDITION1]_vs_[CONDITION2].pdf
    │   ├── ...
    ├── gsea
    │   ├── comparative
    │   │   ├── gsea_cluster_[CLUSTER]_[CONDITION1]_vs_[CONDITION2].pdf
    │   │   ├── ...
    │	└── escape
    │       ├── escape_heatmap_top5.pdf
    │       ├── [CLUSTER]
    │       │   ├── GEYSER_PLOT_[1_path].pdf
    │	    │   ├──...
    │       ├── ...
    ├── integrated_elbow_plot.pdf
    ├── integrated_umap_grouped.pdf
    ├── integrated_umap_labelled.pdf
    ├── integrated_umap_split.pdf
    ├── integrated_umap_unlabelled.pdf
    ├── qc
    │   ├── [sample1]_soupx_nGenes_nUMI.pdf
    │   ├── [sample1]_soupx_percent_mt.pdf
    │   ├── ...
    ├── reference_marker_mapping_heatmap.pdf
    ├── ti
    │   ├── ti_de_slingPseudotime_*.pdf
    │   └── ti_start_smooth.pdf OR ti_no_start_not_smooth.pdf
    └── top3_markers_expr_heatmap.pdf
```

#### Data descriptions:
- optimal_clusters.txt: Contains the number of optimal clusters used for machine learning algorithms throughout the pipline
- qc/: Ambient corrected Cellranger counts using SoupX
- sce_slingshot.rds: SingleCellExperiment R object containing pseudotime calculations for trajectory inference
- sc_integrated_milo_traj.rds: Integrated Seurat object with miloR calculated neighborhoods
- se_filtered_list.rds: Basic QC'd list of Seurat objects; removal of low quality cells from MAD cutoffs of nGenes, nUMI, and percentage of mitochondrial reads
- se_filtered_singlets_list.rds: List of Seurat objects with removed doublets from ScDblFinder
- se_integrated_auto_label.rds: Integrated Seurat object from automatic labelling from reference Seurat object(s)
- se_integrated_dimred.rds: Integrated Seurat object with Seurat dimensional reduction techniques saved (FindClusters, RunPCA, RunUMAP)
- se_integrated_escape_norm.rds: Integrated Seurat object with UCell enrichment scores calculated for each cell, which is normalized in respect to the number of genes
- se_integrated_escape.rds: Integrated Seurat object with UCell enrichment scores calculated for each cell
- se_integrated.rds: Seurat object that has all samples over all conditions integrated in order to be comparable downstream
- se_list_raw.rds: Seurat object list of corrected for ambient RNA
- se_markers_presto_integrated.txt: Results from FindAllMarkers (presto implementation), where differential expression (1 cluster against the rest) is calculated for each gene.
- ti/:
	- ti_gene_clusters_slingPseudotime_*.txt: Gene names for each cluster in the differential expressed genes according to pseudotime

#### Plot descriptions:
- conserved_marker_unlabelled.pdf: Conserved markers between conditions; aids in identifying cluseters manually
- da/: Differential abundance analysis through the miloR package. This will test which condition is more prevelant in the different clusters
	- milo_DA_DE_heatmap_\*.pdf: Differentially expressed genes that meet a differential abundance cutoff between conditions
	- milo_DA_fc_distribution.pdf: Beeswarm plot showing fold-change distribution
	- milo_DA_umap.pdf: Single-cell clustered UMAP from Seurat, which is now labelled with the differential abundance across conditions. Each vertice is the amount of cells in a given neighborhood, and each edge is the number of cells shared between the neighbourhoods.
	- milo_pval_distribution.pdf: Uncorrected P-value distribution, should have a right skew (anti-conservative) distribution
	- milo_volcano_plot.pdf: Visualizes the the SpatialFDR values to see if any neighbourhoods make the cutoff
- deseq2/: Finds differentially expressed genes between conditions for each cluster; calculated with pseudobulk
- gsea/: Gene-set enrichment analyses
	- comparative/: DESeq2 (pseudobulk) GSEA results between conditions for each cluster
	- escape/: UCell GSEA results calculated on clusters of cells one by one
		- escape_heatmap_top5.pdf: Aggregated enrichment scores for each cluster, for the top five GO pathways
		- GEYSER_PLOT_[path].pdf: Shows individual GO paths with individual cells for each cluster. Central dot = median. Thick/thin lines = 66/95% interval summaries
- integrated_elbow_plot.pdf: Rank PCs based on variance percentage, view the elbow for choosing an optimal clustering number
- integrated_umap_grouped.pdf: UMAP of integrated samples separated by conditions. Check if integrated properly and for batch effects
- integrated_umap_labelled.pdf: UMAP of integrated samples with automatic labelling from reference Seurat object(s)
- integrated_umap_split.pdf: UMAP of integrated samples separated by samples
- integrated_umap_unlabelled.pdf: UMAP of integrated samples
- qc/: Preliminary QC graphs to remove low quality cells through MAD values
	- [sample1]_soupx_nGenes_nUMI.pdf: Filtered by number of genes and number of UMIs. High amount = potential multiplets; Low amount = Lysed or ambient RNA cells
	- [sample1]_soupx_percent_mt.pdf: Filtered by the percentage of mitochodrial reads in the cell. High amount = Lysed cell
- reference_marker_mapping_heatmap.pdf: If provided reference Seurat object(s), a prediction score heatmap will be made to visualize the cell labelling to the unnamed clusters
- ti/:
	- ti_de_slingPseudotime_*.pdf = Differentially expressed genes for across pseudotime values calculated by slingshot
	- ti_start_smooth.pdf: Trajectory inference with smooth principal curves; produced when a starting cluster/celltype is explicitly mentioned
	- ti_no_start_not_smooth.pdf: Trajectory inference, useful for trying to figure out if cells are differentiating from one cluster to another; direction is not known though
- top3_markers_expr_heatmap.pdf: Top 3 genes for each cluster from FindAllMarkers call

### Example
From Cisplatin-treated and non-treated mice data.

Counts folder:
```
counts
├── CTRL
│   ├── pilot_study_C1 -> /home/projects/CIO/yard/run_cellranger_count/pilot_study_C1
│   └── pilot_study_C2 -> /home/projects/CIO/yard/run_cellranger_count/pilot_study_C2
└── TREAT
    ├── pilot_study_T1 -> /home/projects/CIO/yard/run_cellranger_count/pilot_study_T1
    └── pilot_study_T2 -> /home/projects/CIO/yard/run_cellranger_count/pilot_study_T2
```
Above I am using symlinks to save space (the symlink path also needs to be included in the binded paths).Ensure that REAL paths are used, an error will occur if attempting to use a symlink path name since apptainer/singularity does not run with root access and needs to be explicitly told where the files can be found.

```
nextflow run scRNA_pipeline.nf \
	--indir /home/projects/sc_pipelines/counts/ \
	--outdir /home/projects/sc_pipelines/test_run_nf_1 \
	--species musmusculus \
	--bind /home/projects/,/home/projects/CIO/yard/run_cellranger_count \
	--reference_seurat /home/projects/sc_pipelines/analysis/data/se_michalski.rds \
	--beginning_cluster Osteoblasts \ 
	-with-apptainer /home/phamj7@med.umanitoba.ca/bin/scrnaseq_singularity.sif
```
