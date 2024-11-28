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
	--main_time [true/false] \
	--beginning_cluster [inferred earliest celltype] \
	--clusters_optimal [optimal number of clusters] \
	--resolution [0-9*] \
	--run_escape [true/false] \
	--pathways [PATHWAY1,PATHWAY2,...] \
	--co_conditions [X1,X2,...] \
	--reduced_dim [integrated.cca | integrated.mnn | harmony] \
	--test_data [0-9*] \
	--run_da [true/false] \
	--run_trajectory_inference [true/false] \
	--run_cellchat [true/false]
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
- run_trajectory_inference: Performs trajectory inference analyses
	- DEFAULT: false
- run_da: Performs miloR differential abundance analysis
	- DEFAULT: false
- run_cellchat: Performs CellChat cell-cell communication analysis
	- DEFAULT: false

### Outputs
##### QC
<details>
<summary>Click to expand</summary>
<br>

- Ambient RNA removal with soupX
- Removing low quality cells using median average deviation (MAD) thresholds based off ...
	- Percentage of mitochondrial reads
	- Number of UMIs
 	- Number of genes
- Doublet removal with scDblFinder
- Files
	- Data
		- *_soupx: SoupX directory output after removal of ambient RNA
		- se_list_raw.rds: Seurat object list of corrected for ambient RNA
		- se_filtered_list.rds: Basic QC'd list of Seurat objects; removal of low quality cells from MAD cutoffs of nGenes, nUMI, and percentage of mitochondrial reads
	 	- se_filtered_singlets_list.rds: List of Seurat objects with removed doublets from ScDblFinder
		- se_filtered_doublets_list.rds: List of Seurat objects with doublets from ScDblFinder
	 - Plots
		- *_nGenes_nUMI.pdf: Filtered by number of genes and number of UMIs. High amount = potential multiplets; Low amount = Lysed or ambient RNA cells
	 	- *_percent_mt.pdf: Filtered by the percentage of mitochodrial reads in the cell. High amount = Lysed cell
		- *_nucleosome_signal: Filters the ratio of mono-nucleosomal to nucleosome-free fragments, also known as a signal-to-noise ratio. Expect to see general decrease in peaks
	 	- *_tss.pdf: Filtered by the transcription start site (TSS) enrichment score. The score of the framents are centered at the TSS. Expect higher enrichment of fragments around the TSS
	  	- *_ncount_atac.pdf: Few reads should be excluded because of low sequencing depth, high levesl can represent doublets, nuclei clumps, etc.
</details>

#### Seurat processing and dimensional reduction
<details>
<summary>Click to expand</summary>
<br>

- Normalize, find variable genes, and scale the data
- Perform integration with desired batch correction tool
	-  Joint integration of scRNA-seq and scATAC-seq data if present
- Perform dimensional reduction and clustering
	- Automatically determines the optimal number of k to use for clustering by NBClust with the euclidean distance option
- Does not perform SCTransform at the moment
- Files
	- Data
		- se_integrated.rds: Seurat object that has all samples over all conditions integrated in order to be comparable downstream
		- se_integrated_dimred.rds: Integrated Seurat object with Seurat dimensional reduction techniques saved (FindClusters, RunPCA, RunUMAP)
		- optimal_clusters.txt: Contains the number of optimal clusters used for machine learning algorithms throughout the pipline
	- Plots
 		- integrated_umap_grouped.pdf: UMAP of integrated samples separated by conditions. Check if integrated properly and for batch effects
		- integrated_umap_split.pdf: UMAP of integrated samples separated by samples
		- integrated_umap_unlabelled.pdf: UMAP of integrated samples
		- percent_cells_group_unlabelled.pdf: Cell cluster proportions for each main condition; Unlabelled clusters
</details>

##### Cell marker identification
<details>
<summary>Click to expand</summary>
<br>

- Find differentially expressed cell markers per cluster
- Performs reference marker mapping if reference seurat was provided
- Files
	- Data
		- se_markers_presto_integrated.txt: Results from FindAllMarkers (presto implementation), where differential expression (1 cluster against the rest) is calculated for each gene.
		- se_integrated_auto_label.rds: Integrated Seurat object from automatic labelling from reference Seurat object(s)
	- Plots
 		- top3_markers_expr_heatmap.pdf: Top 3 genes for each cluster from FindAllMarkers call
		- conserved_marker_unlabelled.pdf: Conserved markers between conditions; aids in identifying cluseters manually
		- integrated_umap_labelled.pdf: UMAP of integrated samples with automatic labelling from reference Seurat object(s)
		- percent_cells_group_labelled.pdf: Cell cluster proportions for each main condition; Labelled clusters
  		- reference_marker_mapping_heatmap.pdf: If provided reference Seurat object(s), a prediction score heatmap will be made to visualize the cell labelling to the unnamed clusters
</details>

#### Pseudo-bulk analysis
<details>
<summary>Click to expand</summary>
<br>

- Averages gene expression for each cluster and performs differential gene expression with DESeq2
- Using the log2FC values as rank, performs GSEA using the fgsea package
	- This method uses ALL genes, not an over-representation analyses method; DEGs would affect the pathways, but will have to filter the GSEA for them
- Will perform pairwise comparisons, but matching the like co-condition with like
	- For example, given sample1_CTRL_1hr, sample2_CTRL_2hr, sample1_TREAT_1hr, sample2_TREAT_2hr the pairs will be sample1_CTRL_1hr vs sample1_TREAT_1hr and sample2_CTRL_2hr vs sample2_TREAT_2hr
- Files
	- Data
		- deseq2_cluster_CONDITION-X_vs_CONDITION-Y.txt: Tab-delimited file with DESeq2 results, unfiltered
	- Plots
 		- deseq2_cluster_CONDITION-X_vs_CONDITION-Y.pdf: Finds differentially expressed genes between conditions for each cluster; calculated with pseudobulk
		- gsea_cluster_CONDITION-X_vs_CONDITION-Y.pdf: DESeq2 (pseudobulk) GSEA results between conditions for each cluster
</details>

#### Trajectory inference
##### Slingshot
<details>
<summary>Click to expand</summary>
<br>
	
- Uses a minimum spanning tree (MST) on a given UMAP to form the initial cluster connections
- Depending on whether there is a starting celltype specified (--beginning_cluster NAME), it will run in 2 different ways
	- No start: No pseudotimes are calculated, since the direction is unknown
	- Start: The MST will be smoothed using principal curves and then each cell will have a pseudotime calculated based on their orthogonal distance from the curve
- Files
	- Data
		- ti_de_between_GROUPING.txt: List of differentially expressed genes between the given condition. Sorted in the same order as the corresponding heatmap. [Tutorial](https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html)
		- ti_DEGs_qval_full_lineage_[0-9].txt: List of differentially expressed genes according to pseudotime in a predicted lineage. 
		- sce_slingshot.rds: SingleCellExperiment R object containing pseudotime calculations for trajectory inference
 	- Plots
  		- ti_no_start_not_smooth.pdf: Trajectory inference, useful for trying to figure out if cells are differentiating from one cluster to another; direction is not known though
    		- ti_start_smooth.pdf: Trajectory inference with smooth principal curves; produced when a starting cluster/celltype is explicitly mentioned
      		- ti_deg_between_GROUPING.pdf: Trajectory inference of differentially expressed genes between the given condition. Rows are genes, columns are cells, and the columns can also be split by each predicted lineage. |Log2FC| >= 2 and padj < 0.05
		- ti_de_lineage[0-9].pdf: Differentially expressed genes for across pseudotime values calculated by slingshot. Only the 50 genes are listed max in the heatmap. |Log2FC| >= 2 and padj < 0.05
</details>

##### Tempora (Deprecated for now)
<details>
<summary>Click to expand</summary>
<br>
	
- Uses time-series scRNA-seq
- Uses pathway analyses and time variable to order cells
	- Uses mutual information rank and ARACNE to calculate the ranks in all cluster pairs and assigns the order of edges made
	- ARACNE in this case infers direct regulatory relationships between transcriptional regulator proteins and target genes
- Files
	- Data
  		- se_integrated_tempora_seurat_v3.rds: Tempora object; Seurat3 converted from Seurat5
	- Plots
 		- tempora_screeplot_GROUPING.pdf: Used to visualize and determine the optimal number of clusters in tempora
		- tempora_inferred_lineages_GROUPING.pdf: Displays the predicted trajectory of cell types from early to late cell types; Please use literature to confirm results
</details>

##### Psupertime
<details>
<summary>Click to expand</summary>
<br>

- Uses time-series scRNA-seq
- Identifies ordering coefficients for each genes, where the value is equal to the contribution the gene has to cell ordering (based on their sequential expression)
	- This coefficient is multipled with the gene expression matrix to get the pseudotime values for each cell
- Files
	- Data
 		- psuper_top_20_genes_CONDITION-X_CLUSTER.txt: List of top 20 genes with the largest absolute coefficient values
	- Plots
		- psuper_boxplot_compare_dist_CLUSTER.pdf: Showing distribution of psupertimes for a given cluster label for each timepoint. T-test performed for significance test 
		- psuper_density_pseudotime_CONDITION-X_CLUSTER.pdf: Density plots which is sorted by psupertime
  		- psuper_gene_coefficients_CONDITION-X_CLUSTER.pdf: Top gene coefficient scores in the given condition
    		- psuper_gene_coefficients_CONDITION-X_CLUSTER_CONDITION-Y_genes.pdf: Looking at the top genes in CONDITION-X with the coefficients seen in CONDITION-Y
      		- psuper_top_20_genes_over_pseudotime_CONDITION-X_CLUSTER.pdf = CONDITION-X gene expression of the top rated coefficients
        	- psuper_top_20_genes_over_pseudotime_CONDITION-X_CLUSTER_CONDITION-Y_genes.pdf: CONDITION-Y gene expression of the CONDITION-X top rated coefficients genes
         	- psuper_train_results_CONDITION-X_CLUSTER.pdf: Machine learning statistics from the psupertime functions
</details>

#### Escape
<details>
<summary>Click to expand</summary>
<br>

- Utilizes the heterogeneity of scRNA-seq data to perform GSEA
- Multiple methods of calculating enrichment scores, but UCell is used here for a heuristic
	- UCell = Based on Mann Whitney U Statistic
 		- Score calculated relatively to gene expression in individual cells
		- Creates a ranked list of genes for EACH cell, so the process will take a while
- As of the moment there are only visualizations, but [differential enrichment](https://www.borch.dev/uploads/screpertoire/articles/running_escape) can be calculated 
- Files
	- Data
 		- se_integrated_escape_norm.rds: Integrated Seurat object with UCell enrichment scores calculated for each cell, with normalized in respect to the number of genes and unormalized values
	- Plots
 		- escape_heatmap_PATHWAY-PATTERN.pdf: Aggregated single-cell GSEA results for pathways containing a phrase from '--pathways' if it was used in the intial command
   		- escape_heatmap_top5.pdf: Aggregated enrichment scores for each cluster, for the top five GO pathways
   		- PATHWAY_geyser.pdf: Shows individual GO paths with individual cells for each cluster. Central dot = median. Thick/thin lines = 66/95% interval summaries
</details>

#### miloR
<details>
<summary>Click to expand</summary>
<br>

- Detects changes in compositions between multicondition scRNA-seq datasets. [Tutorial](https://www.bioconductor.org/packages/devel/bioc/vignettes/miloR/inst/doc/milo_gastrulation.html) with more in-depth info for the analyses
	- https://marionilab.github.io/miloR/articles/milo_gastrulation.html
- Simplified Steps:
	1. k-nearest neighbours (KNN) graph is built from the batch correction integration method that was chosen (default: integrated_cca from Seurat)
	2. Define what the a neighbourhood will be for each cell (index)
 		- This will form edges connecting the index to surrounding cells
   		- Want average neighbourhood size over 5 x N_samples (This needs manual checks)
     		- Refined by calculating the median profile for a given neighbourhood by selecting the vertex (cell) with the highest number of triangle in the neighbourhood
       			- This may reduce multiple neighbourshoods to prevent oversampling
     	3. Count the amount of cells from each sample are in each neighbourhood
      	4. Define the experimental design (this will be between the main condition and then co-conditions as well)
       	5. Testing for DA using a negative binomial generalized linear model implementation in edgeR between neighbourhoods
    		- Contains Fold-change and FDR for each neightbourhood to determine significant differential abundance
       		- SpatialFDR = Corrected pvals for multiple testing from the overlap between neighbourhoods
       	 		- Weighted version of BH method in the KNN graph
- Files
	- Data
		- sc_integrated_milo_traj.rds: Integrated Seurat object with miloR calculated neighborhoods
		- da_CLUSTER_markers_avg_logfc_GROUPING.txt: FindNhoodGroupMarkers results; 1 vs all differential gene expresssion for each neighbourhood group. Added 'avg_logFC' to use in GSEA
		- da_CLUSTER_markers_by_neighborhood_GROUPING_expr_matrix.txt: Expression matrix of each significant neighbourhood
		- da_CLUSTER_markers_by_neighborhood_GROUPING.txt: Gene list from DEG heatmap. Cluster[0-10] is shown for downstream filtering
		- da_diff_test_GROUPING.txt: Differential abundance testing from testNhoods. Determins if there is a significant differential abundace between the GROUPING
		- milo_gsea_cluster_CLUSTER_GROUPING.txt: List of GSEA results from fgsea. Unfiltered
	- Plots
 		- milo_DA_DE_heatmap_CLUSTER_GROUPING.pdf: Differentially expressed genes that meet a differential abundance cutoff between conditions; logFC used to determine if a gene is differentially expressed is from using da_CLUSTER_markers_avg_logfc_GROUPING 'avg_logFC' value
   		- milo_DA_fc_distribution_GROUPING.pdf: Beeswarm plot showing fold-change distribution
     		- milo_DA_umap_GROUPING.pdf: Single-cell clustered UMAP from Seurat, where each point is now gathered into neighbourhoods, edges show shared cells between neighrbourhoods, and colour is the logFC value from testNhoods
       		- milo_gsea_cluster_CLUSTER_GROUPING.pdf: GSEA results that uses the da_CLUSTER_markers_avg_logfc_GROUPING 'avg_logFC' values as a ranked list.
         	- milo_pval_distribution_GROUPING.pdf: Uncorrected P-value distribution, should have a right skew (anti-conservative) distribution
          	- milo_volcano_plot_GROUPING.pdf: Visualizes the the SpatialFDR values to see if any neighbourhoods make the cutoff
</details>

#### CellChat
<details>
<summary>Click to expand</summary>
<br>

- Finding over-expressed ligand-receptor pairs and pathways using a communication probability
	- https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
	- Communication probability is built off gene expression, the ligand-receptor interactions from the curated CellChat database, and law of mass action
		- Genes are chosen if they are overexpressed ligands and/or receptors
- Significant cell-cell communication is calculated by performing a permutation test on the communication probability
- For cellchat_compare_*heatmap.pdf, the top bar plot shows total signal strength of a cell group (summarizes all signalling pathways for that cell group). The right bar plot will show the total signaling strength for a pathway (summarizes all cell groups for that signaling pathway)
- Files
	- Data
		 - cellchat_merged.rds: Merged CellChat object; used in all downstream analyses
		 - cellchat_object_list.rds: List of CellChat objects; length = # of samples
	- Plots
		 - cellchat_compare_all_signal_heatmap.pdf: Aggregated heatmap of outgoing and incoming signals;
		 - cellchat_compare_incoming_signal_heatmap.pdf: Heatmap of incoming signals
		 - cellchat_compare_outgoing_signal_heatmap.pdf: Heatmap of outgoing signals
		 - cellchat_differential_interaction_circle.pdf: Red colored edges show increased signaling in the second dataset (Sorted it lexicographical order). Size of the edge represents the differential number of interactions/strength
		 - cellchat_differential_interaction_heatmap.pdf: Red show increased signaling in the second dataset (Sorted it lexicographical order). Top bar plot = Incoming signal for the given cell group; Right bar plot = Outgoing signal for the given cell group
		 - cellchat_information_flow_compare.pdf: Information flow = sum of the communication probability amoung all pairs of cell groups in the pathway, or known as total weight of the network. Highlighted pathways display the significant pathways according to a paired Wilcoxon test.
		 - cellchat_interaction_summary_bar.pdf: Compares total number of interactions and strength of the cell-cell communication networks
		 - cellchat_num_interactions_circle.pdf: Shows the number of interactions between the different cell groups
		 - cellchat_population_send_receive.pdf: Identifying cell populations with significant changes in sending/receiving signals
		 - commun_prob/
			 - cellchat_CLUSTER_expression.pdf: Comparing communication probabilities from different ligand-receptor pairs. Each column is cell-cell specific, and the different colors represent the condition. If the title of the plot does not show "Unfiltered" that means the communication probability difference between two cell-cell for the two different conditions is ~50%
		 - signaling_pathways/
			 - PATHWAY
				 - cellchat_PATHWAY_expression.pdf: Shows the different expression of genes involved in the given pathway, if found in the dataset
				 - cellchat_PATHWAY_signal_path.pdf: Shows which cells are involved in the given PATHWAY. The inner thinner bars reprsent the target cell that receive the signal, and the size is proportional to the signal strength received
</details>

#### scATAC
<details>
<summary>Click to expand</summary>
<br>

- If the data contains corresponding scATAC-seq data for the same set of scRNA-seq cells, the data will be jointly analyzed
- QC: Added in filtering low-quality cells based on the total number of fragments in peaks, TSS enrichment score, and nucleosome signal
- Changed steps to dimensional reduction and clustering
	- Dimensional reduction
		- Term frequency-inverse document frequency (TF-IDF) normalization = Normalizes cell sequencing depth across the peaks
	 	- Feature selection = Chooses top percentage of peaks or remove peaks that are in less than n cells
		- Similar to PCA in scRNA-seq, perform singular value decomposition (SVD) on the normalized gene matrix
		- When performed in this order, it is called the latent semantic indexing (LSI)
	- Calculate the weighted nearest neighbours graph to combine the RNA and ATAC-seq information to be used in downstream clustering
- Files
	- Data
		 - CONDITION-X_vs_CONDITION-Y
			 - volcano_data
				 - dap_CELLTYPE_CONDITION-X_vs_CONDITION-Y.txt: Tab-delimited file with DESeq2 results, unfiltered
	- Plots
		 - CONDITION-X_vs_CONDITION-Y
			 - closest_gene_plots: Closest gene to one of the top differentially expressed peaks found from the DESeq2 results. Default is the top and bot 5 peaks. adj_pval < 0.05 and |avg_log2FC| >= 1
				 - scatac_closest_genes_dap_coverage_CELLTYPE_GENE_up.pdf: Displays the genome track with the corresponding RNA-seq, ATAC-seq, and possible ATAC-seq peak links
				 - scatac_closest_genes_dap_gex_CELLTYPE_GROUPING.pdf: Violin plots of the top 5 associated genes
				 - scatac_closest_genes_dap_gsea_CELLTYPE.pdf: After linking each peak to their closest gene, perform GSEA. Ranked list by log2FC
			 - motif_plots
				 - scatac_motif_CELLTYPE_GROUPING.pdf: Top (default: 6) over-represented motifs in the peak set
			 - volcano_plots
				 - scatac_volcano_CELLTYPE_CONDITION-X_vs_CONDITION-Y.pdf: Differentiall accessible peaks. adj_pval < 0.05 and |avg_log2FC| >= 1
</details>

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

