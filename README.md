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

#### Apptainer(Singularity)
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
	-with-apptainer path_to_image/scrnaseq_singularity_tmp.sif
```

### Arguements
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
- with-apptainer: path to downloaded image (REQUIRED)

### Outputs
Analysis folder:

```

```



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
Above I am using symlinks to save space (the symlink path also needs to be included in the binded paths)

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
