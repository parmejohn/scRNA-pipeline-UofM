# scRNA-seq Analysis Pipeline

## Installation

### Dependencies
- NextFlow >= 23.10.1
- Apptainer(singularity); tested with version 1.2.5-1.el7

### Script

```
git clone https://github.com/parmejohn/scRNA-pipeline-UofM.git
```

## Usage
Run the script from the downloaded directory

```
nextflow run scRNA_pipeline.nf \
	--indir <input folder to sorted cellranger counts> \
	--outdir <output folder> \
	--species <species name> \
	--bind <bind file path locations> \
	--reference_seurat [path(s) to labelled reference seurat object] \
	--sc_atac [true/false] \
	--replicates [true/false] \
	--run_sling [true/false] \
	--main_time [true/false] \
	--beginning_cluster [inferred earliest celltype] \
	--clusters_optimal [optimal number of clusters] \
	--resolution [0-9*] \
	--run_escape [true/false] \
	--pathways [PATHWAY1,PATHWAY2,...] \
	--co_conditions [X1,X2,...] \
	--reduced_dim [integrated.cca | integrated.mnn | harmony] \
	--test_data [0-9*] \
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
- sc_atac: Will perform joint integration of scRNA-seq and scATAC-seq data. Also will perform a series of scATAC analyses. Must have multiome 10X data available in the indir.
	- DEFAULT: false
- replicates: Signifies if analyses has replicates available per main condition. If false, will do the bare minimum analyses.
	- DEFAULT: true
- main_time: If the main condition being tested for is time, this parameter will specify it for the downstream analyses.
	- DEFAULT: false
- run_escape: Performs escape analysis to provide single-cell GSEA results. IMPORTANT: This analysis takes a substantial amount of time and memory, so ensure that you provide enough of both before setting to 'true'. With Nextflow capabilities as well, it is fine to run the pipeline once without escape and then with  '-resume', since the pipeline will cache the rest of the data.
	- DEFAULT: false
- pathways: List of desired pathways/phrases separated by ',' to search on the Gene onotology: biological processes for MsigDB for the escape analysis. Setting this option will make it so the top 5 pathways will not be printed during escape. Used if '--run_escape true'
	- DEFAULT: 'none'
- run_sling: Performs slingshot trajectory inference analysis
	- DEFAULT: false
- beginning_cluster: If you know what celltype is would be the earliest in a differentiation pathway. Used if '--run_slingshot true'
	- DEFAULT: ''
- clusters_optimal: Number of optimal clusters/dimensionality for machine learning methods. If set to 0, will be calculated automatically through NbClust() Calinski-Harabasz (ch) index.
	- DEFAULT: 0
- resolution: Desired resolution, increasing this will increase the number of clusters and vice versa.
	- DEFAULT: 1
- co_conditions: List of co-conditions for each sample separated by ','. Your sample name should contain an underscore for each condition, including the main one being tested for. For instance, SAMPLE1_MUT_DAY50_treated, could have the arguement of '--co_condition time,treatment'
	- DEFAULT: 'none'
- reduced_dim: Batch correction method to carry out while integrating samples. Possible arguements are "integrated.cca (CCA), integrated.mnn (FastMNN), or harmony"
	- DEFAULT: 'integrated.cca'
- test_data: Number of cells to downsample for, if 0 the pipeline will take all cells into account
	- DEFAULT: 0

### Outputs
##### QC
<details>
<summary>Click to expand</summary>
<br>

- Ambient RNA removal with soupX
- Removing low quality cells using MAD thresholds based off ...
	- Percentage of mitochondrial reads
	- Number of UMIs
 	- Number of genes
- Doublet removal with scDblFinder
- Files
	- Data
		- *_soupx
		- se_list_raw.rds
		- se_filtered_list.rds
	 	- se_filtered_singlets_list.rds
		- se_filtered_doublets_list.rds
	 - Plots
		- *_nGenes_nUMI.pdf
	 	- *_percent_mt.pdf
		- *_nucleosome_signal
	 	- *_tss.pdf
	  	- *_ncount_atac.pdf
</details>

#### Seurat processing and dimensional reduction
<details>
<summary>Click to expand</summary>
<br>

- Normalize, find variable genes, and scale the data
- Perform integration with desired batch correction tool
	-  Joint integration of scRNA-seq and scATAC-seq data if present
- Perform dimensional reduction and clustering
- Files
	- Data
		- se_integrated.rds
		- se_integrated_dimred.rds
		- optimal_clusters.txt
	- Plots
 		- integrated_umap_grouped.pdf
		- integrated_umap_split.pdf
		- integrated_umap_unlabelled.pdf
		- percent_cells_group_unlabelled.pdf
</details>

##### Cell marker identification
<details>
<summary>Click to expand</summary>
<br>

- Find differentially expressed cell markers per cluster
- Performs reference marker mapping if reference seurat was provided
- Files
	- Data
		- se_markers_presto_integrated.txt
		- se_integrated_auto_label.rds
	- Plots
 		- top3_markers_expr_heatmap.pdf
		- conserved_marker_unlabelled.pdf
		- integrated_umap_labelled.pdf
		- percent_cells_group_labelled.pdf
</details>

#### Pseudo-bulk analysis
<details>
<summary>Click to expand</summary>
<br>

- Averages gene expression for each cluster and performs differential gene expression with DESeq2
- Using the log2FC values as rank, performs GSEA using the fgsea package
- Will perform pairwise comparisons, but matching the like co-condition with like
	- For example, given sample1_CTRL_1hr, sample2_CTRL_2hr, sample1_TREAT_1hr, sample2_TREAT_2hr the pairs will be sample1_CTRL_1hr vs sample1_TREAT_1hr and sample2_CTRL_2hr vs sample2_TREAT_2hr
- Files
	- Data
		- deseq2_cluster_CONDITION1_vs_CONDITION2.txt
	- Plots
 		- deseq2_cluster_CONDITION1_vs_CONDITION2.pdf
		- gsea_cluster_CONDITION1_vs_CONDITION2.pdf
</details>

#### Trajectory inference
##### Slingshot
<details>
<summary>Click to expand</summary>
<br>
	
- Files
	- Data
		- ti_de_between_group.txt
		- ti_DEGs_qval_full_lineage_[0-9].txt
		- sce_slingshot.rds
 	- Plots
  		- ti_no_start_not_smooth.pdf
    		- ti_start_smooth.pdf
      		- ti_deg_between_group.pdf
		- ti_de_lineage[0-9].pdf
</details>

##### Tempora
<details>
<summary>Click to expand</summary>
<br>

- Files
	- Data
  		- se_integrated_tempora_seurat_v3.rds
	- Plots
 		- tempora_screeplot_CONDITION.pdf
		- tempora_inferred_lineages_CONDITION.pdf
</details>

##### Psupertime
<details>
<summary>Click to expand</summary>
<br>

- Files
	- Data
	- Plots
</details>

#### Escape
<details>
<summary>Click to expand</summary>
<br>

- Files
	- Data
	- Plots
</details>

#### CellChat
<details>
<summary>Click to expand</summary>
<br>

- Files
	- Data
	- Plots
</details>

#### miloR
<details>
<summary>Click to expand</summary>
<br>

- Files
	- Data
	- Plots
</details>

#### scATAC
<details>
<summary>Click to expand</summary>
<br>

- Files
	- Data
	- Plots
</details>

#### Data descriptions:
- optimal_clusters.txt: Contains the number of optimal clusters used for machine learning algorithms throughout the pipline
- qc/: Ambient corrected Cellranger counts using SoupX
- sce_slingshot.rds: SingleCellExperiment R object containing pseudotime calculations for trajectory inference
- sc_integrated_milo_traj.rds: Integrated Seurat object with miloR calculated neighborhoods
- se_filtered_list.rds: Basic QC'd list of Seurat objects; removal of low quality cells from MAD cutoffs of nGenes, nUMI, and percentage of mitochondrial reads
- se_filtered_singlets_list.rds: List of Seurat objects with removed doublets from ScDblFinder
- se_integrated_auto_label.rds: Integrated Seurat object from automatic labelling from reference Seurat object(s)
- se_integrated_dimred.rds: Integrated Seurat object with Seurat dimensional reduction techniques saved (FindClusters, RunPCA, RunUMAP)
- se_integrated_escape_norm.rds: Integrated Seurat object with UCell enrichment scores calculated for each cell, with normalized in respect to the number of genes and unormalized values
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
Cisplatin-treated and non-treated mice data.

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
Above I am using symlinks to save space (the symlink path also needs to be included in the binded paths). Ensure that REAL paths are used, an error will occur if attempting to use a symlink path name since apptainer/singularity does not run with root access and needs to be explicitly told where the files can be found.

```
nextflow run scRNA_pipeline.nf \
	--indir /home/projects/sc_pipelines/counts/ \
	--outdir /home/projects/sc_pipelines/test_run_nf_1 \
	--species musmusculus \
	--bind /home/projects/,/home/projects/CIO/yard/run_cellranger_count \
	--reference_seurat /home/projects/sc_pipelines/analysis/data/se_michalski.rds \
	--run_sling true \
	--beginning_cluster Osteoblasts
```

## References
### Bioinformatics Tools
- batchelor: Haghverdi, L., Lun, A. T. L., Morgan, M. D. & Marioni, J. C. Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. Nat Biotechnol 36, 421–427 (2018).
- Comparing batch corrections: Luecken, M. D. et al. Benchmarking atlas-level data integration in single-cell genomics. Nat Methods 19, 41–50 (2022).
- DESeq2: Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550 (2014).
- Escape: Borcherding, N. et al. Mapping the immune environment in clear cell renal carcinoma by single-cell genomics. Commun Biol 4, 1–11 (2021).
- fgsea: Korotkevich, G. et al. Fast gene set enrichment analysis. 060012 Preprint at https://doi.org/10.1101/060012 (2021).
- harmony: Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16, 1289–1296 (2019).
- miloR: Dann, E., Henderson, N. C., Teichmann, S. A., Morgan, M. D. & Marioni, J. C. Differential abundance testing on single-cell data using k-nearest neighbor graphs. Nat Biotechnol 40, 245–253 (2022).
- presto: Korsunsky, I., Nathan, A., Millard, N. & Raychaudhuri, S. Presto scales Wilcoxon and auROC analyses to millions of observations. 653253 Preprint at https://doi.org/10.1101/653253 (2019).
- scDblFinder: Germain, P.-L., Lun, A., Meixide, C. G., Macnair, W. & Robinson, M. D. Doublet identification in single-cell sequencing data using scDblFinder. Preprint at https://doi.org/10.12688/f1000research.73600.2 (2022).
- Seurat: Hao, Y. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat Biotechnol 42, 293–304 (2024).
- Slingshot: Street, K. et al. Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. BMC Genomics 19, 477 (2018).
- SoupX: Young, M. D. & Behjati, S. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. GigaScience 9, giaa151 (2020).
- Tempora: Tran, T. N. & Bader, G. D. Tempora: Cell trajectory inference using time-series single-cell RNA sequencing data. PLOS Computational Biology 16, e1008205 (2020).
- tradeSeq: Van den Berge, K. et al. Trajectory-based differential expression analysis for single-cell sequencing data. Nat Commun 11, 1201 (2020).

