
// Required inputs; will be checked in the R scripts
params.indir = ''
params.outdir = ''
params.species= ''
params.bind = ''
params.plot_format = 'pdf'

//Optional inputs
params.reference_seurat = 'none'
params.beginning_cluster = 'none'
params.co_conditions = 'none'
params.reduced_dim = 'integrated.cca'
params.pathways = 'none'
params.main_time = false
params.replicates = true

//choosing which analyses to run
params.sc_atac = false
params.run_trajectory_inference = false
params.run_da = false
params.run_escape = false
params.run_cellchat = false
params.run_neuroestimator = false

//customizing parameters for different analyses
params.clusters_optimal = 0
params.resolution = 1
params.test_data = 0
params.mitochondrial_percent_cutoff = 0

include {AMBIENTRNAREMOVAL} from './modules/ambientremoval.nf'
include {FILTERLOWQUALDOUBLETS} from './modules/filterlowqualdoublets.nf'
include {INTEGRATEDATA} from './modules/integratedata.nf'
include {DIMENSIONALREDUCTION} from './modules/dimensionalreduction.nf'
include {IDENTIFYMARKERS} from './modules/identifymarkers.nf'
include {COMPARATIVEANALYSIS} from './modules/comparativeanalysis.nf'
include {TRAJECTORYINFERENCE} from './modules/trajectoryinference.nf'
include {ESCAPEANALYSIS} from './modules/escapeanalysis.nf'
include {DAANALYSIS} from './modules/daanalysis.nf'
include {TEMPORAANALYSIS} from './modules/temporaanalysis.nf'
include {PSUPERTIME} from './modules/psupertime.nf'
include {CELLCHAT} from './modules/cellchat.nf'
include {ATACANALYSES} from './modules/atacanalyses.nf'
include {NEUROESTIMATOR} from './modules/neuroestimator.nf'
include {NEUROESTIMATORPLOT} from './modules/neuroestimator.nf'
include {SUMMARYREPORT} from './modules/summaryreport.nf'

process INITIALIZEFOLDERS {
    input:
    path params.outdir
    
    output:
    //path '${params.outdir}/*'
    //outdata_ch = Channel.value('${params.outdir}/analysis/data/')
    //plotdata_ch = Channel.value('${params.outdir}/analysis/plots/')

    """
    #!/bin/bash

    if [! -d ${params.outdir}/analysis/]; then
        mkdir ${params.outdir}/analysis/
        mkdir ${params.outdir}/analysis/data/

        mkdir ${params.outdir}/analysis/plots
    fi
    """ 
}

process CONVERTSEURAT {
    debug true
    cache 'deep'
	containerOptions "--bind $params.bind"	    

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
    )
    
    input:
    path integrated
    
    output:
    path "se_integrated_v3.rds", emit: se_integrated_obj
	path "se_integrated_v3_counts_matrix.rds", emit: se_integrated_cnts

	script:
    """
	 ${projectDir}/src/ConvertSeurat5.R $integrated
    """ 
}

workflow {
    INITIALIZEFOLDERS(params.outdir)
    
    processed_co_condition = params.co_conditions.replaceAll(',', ' ')
    ambient_ch = AMBIENTRNAREMOVAL(params.indir, params.test_data) // array of conditions
    
    if (params.sc_atac){
  
      FILTERLOWQUALDOUBLETS(ambient_ch, params.species, processed_co_condition, "yes", params.indir, params.mitochondrial_percent_cutoff, params.plot_format)
      filtered_ch = FILTERLOWQUALDOUBLETS.out.se_filtered_singlets_list
      
    } else {

      FILTERLOWQUALDOUBLETS(ambient_ch, params.species, processed_co_condition, "no", params.indir, params.mitochondrial_percent_cutoff, params.plot_format)
      filtered_ch = FILTERLOWQUALDOUBLETS.out.se_filtered_singlets_list
    }
    
    // integrated seurat object
    integrated_ch = INTEGRATEDATA(filtered_ch, params.reduced_dim, processed_co_condition)

    // seurat dim red
    DIMENSIONALREDUCTION(integrated_ch, params.clusters_optimal, params.resolution, params.reduced_dim, params.plot_format)
    dimred_ch = DIMENSIONALREDUCTION.out.se_integrated_dimred

  	new_opt_clust = DIMENSIONALREDUCTION.out.clusters_optimal_n.splitText().map{it -> it.trim()}
  	println new_opt_clust

    // identify cell markers
    if (params.reference_seurat != 'none'){
        IDENTIFYMARKERS(dimred_ch, new_opt_clust, params.reference_seurat, params.plot_format)
	      identified_report = IDENTIFYMARKERS.out.report
        identified_ch = IDENTIFYMARKERS.out.se_integrated_auto_label
    } else {
        IDENTIFYMARKERS(dimred_ch, new_opt_clust, params.reference_seurat, params.plot_format)
        identified_ch = dimred_ch
	      identified_report = "no reference"
    }

	if (params.run_trajectory_inference){
		TRAJECTORYINFERENCE(identified_ch, params.beginning_cluster, params.plot_format)
		sling_ch = TRAJECTORYINFERENCE.out.slingshot_out_pdf
	} else {
		sling_ch = "SKIPPING SLINGSHOT ANALYSIS, RERUN WITH '--params.run_trajectory_inference true' if you desired those results"
	}
	
  m = params.co_conditions ==~ '.*time.*'
  assert m instanceof Boolean
  if (m && params.run_trajectory_inference){
    //TEMPORAANALYSIS(identified_ch, "no",  params.species, params.plot_format)
    //tempora_ch = TEMPORAANALYSIS.out.report
	tempora_ch = "SKIPPING TEMPORA, no time condition set, or --params.run_trajectory_inference false"    

    PSUPERTIME(identified_ch, "no", params.plot_format)
    psupertime_ch = PSUPERTIME.out.report
    
  } else if (params.main_time && params.run_trajectory_inference) {
    //TEMPORAANALYSIS(identified_ch, "yes", params.species, params.plot_format)
    //tempora_ch = TEMPORAANALYSIS.out.report
	tempora_ch = "SKIPPING TEMPORA, no time condition set, or --params.run_trajectory_inference false"    

    PSUPERTIME(identified_ch, "yes", params.plot_format)
    psupertime_ch = PSUPERTIME.out.report
    
  } else {
    tempora_ch = "SKIPPING TEMPORA, no time condition set, or --params.run_trajectory_inference false"
    psupertime_ch = "SKIPPING PSUPERTIME, no time condition set, or --params.run_trajectory_inference false"
  }
  
  if (params.replicates){
    COMPARATIVEANALYSIS(identified_ch, params.species, params.plot_format)
	  comparative_ch = COMPARATIVEANALYSIS.out.report
	  
	  if (params.run_da) {
      DAANALYSIS(identified_ch, new_opt_clust, params.reduced_dim, params.species, params.plot_format)
  	  daanalysis_ch = DAANALYSIS.out.report
	  } else {
	    daanalysis_ch = "--run_da false"
	  }
	  
	  if (params.run_cellchat) {
      CELLCHAT(identified_ch, params.species, params.plot_format)
  	  cellchat_ch = CELLCHAT.out.report
	  } else {
	    cellchat_ch = "--run_cellchat false"
	  }
	  
  } else {
    comparative_ch = "NO REPLICATES"
    daanalysis_ch = "NO REPLICATES"
    cellchat_ch = "NO REPLICATES"
  }
  
	if (params.run_escape){
	  processed_pathways = params.pathways.replaceAll(',', ' ')
		println "WARNING: This may take a while"
		ESCAPEANALYSIS(identified_ch, params.species, processed_pathways, params.plot_format)
		escape_ch = ESCAPEANALYSIS.out.se_integrated_escape
	} else {
		escape_ch = "SKIPPING ESCAPE ANALYSIS, RERUN WITH '--run_escape true' if you desired those results"
	}

	if (params.sc_atac){
		ATACANALYSES(identified_ch, params.species, params.plot_format)
		atac_ch = ATACANALYSES.out.report
	} else {
		atac_ch = "no atac-seq info provided"
	}

  CONVERTSEURAT(identified_ch)
	seurat3_ch = CONVERTSEURAT.out.se_integrated_cnts
	seurat3_obj_ch = CONVERTSEURAT.out.se_integrated_obj
  
  if (params.run_neuroestimator){
    neuroestimator_ch = NEUROESTIMATOR(seurat3_ch, params.species)
    NEUROESTIMATORPLOT(seurat3_obj_ch, neuroestimator_ch, params.plot_format)
	neuroestimatorplot_ch = NEUROESTIMATORPLOT.out.plot
  } else {
    neuroestimatorplot_ch = "run_neuroestimator false"
  }

    SUMMARYREPORT(
        comparative_ch,
        identified_report,
        sling_ch,
        daanalysis_ch,
        escape_ch,
        tempora_ch,
        psupertime_ch,
        cellchat_ch,
        atac_ch,
        neuroestimatorplot_ch,
        new_opt_clust
        )
}
