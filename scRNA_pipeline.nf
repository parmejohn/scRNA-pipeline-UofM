
// Required inputs; will be checked in the R scripts
params.indir = ''
params.outdir = ''
params.species= ''
params.bind = ''

//Optional inputs
params.reference_seurat = 'none'
params.clusters_optimal = 0
params.resolution = 1
params.beginning_cluster = 'none'
params.run_escape = false
params.run_sling = false
params.test_data = 0
params.co_conditions = 'none'
params.reduced_dim = 'integrated.cca'
params.pathways = 'none'

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

workflow {
    INITIALIZEFOLDERS(params.outdir)
    ambient_ch = AMBIENTRNAREMOVAL(params.indir, params.test_data) // array of conditions
    // ambient_ch.view() // outputs just the qc folder, for downstream qc, all samples need to be loaded as a list
    
    processed_co_condition = params.co_conditions.replaceAll(',', ' ')
    FILTERLOWQUALDOUBLETS(ambient_ch, params.species, processed_co_condition)
    filtered_ch = FILTERLOWQUALDOUBLETS.out.se_filtered_singlets_list
    
    // integrated seurat object
    integrated_ch = INTEGRATEDATA(filtered_ch, params.reduced_dim, processed_co_condition)

    // seurat dim red
    DIMENSIONALREDUCTION(integrated_ch, params.clusters_optimal, params.resolution, params.reduced_dim)
    dimred_ch = DIMENSIONALREDUCTION.out.se_integrated_dimred

	new_opt_clust = DIMENSIONALREDUCTION.out.clusters_optimal_n.splitText().map{it -> it.trim()}
	println new_opt_clust

    // identify cell markers
    if (params.reference_seurat != 'none'){
        IDENTIFYMARKERS(dimred_ch, new_opt_clust, params.reference_seurat)
        identified_ch = IDENTIFYMARKERS.out.se_integrated_auto_label
    } else {
        IDENTIFYMARKERS(dimred_ch, new_opt_clust, params.reference_seurat)
        identified_ch = dimred_ch
    }

	if (params.run_sling){
		TRAJECTORYINFERENCE(identified_ch, params.beginning_cluster)
		sling_ch = TRAJECTORYINFERENCE.out.slingshot_out_pdf
	} else {
		sling_ch = "SKIPPING SLINGSHOT ANALYSIS, RERUN WITH '--run_sling true' if you desired those results"
	}
	println sling_ch
	// TRAJECTORYINFERENCE(identified_ch, params.beginning_cluster)
	
  m = params.co_conditions ==~ '.*time.*'
  assert m instanceof Boolean
  if (m){
    TEMPORAANALYSIS(identified_ch)
    tempora_ch = TEMPORAANALYSIS.out.report
  } else {
    tempora_ch = "SKIPPING TEMPORA, no time co-condition set"
  }

    COMPARATIVEANALYSIS(identified_ch, params.species)

	if (params.run_escape){
	  processed_pathways = params.pathways.replaceAll(',', ' ')
		println "WARNING: This may take a while"
		ESCAPEANALYSIS(identified_ch, params.species, processed_pathways)
		escape_ch = ESCAPEANALYSIS.out.se_integrated_escape
	} else {
		escape_ch = "SKIPPING ESCAPE ANALYSIS, RERUN WITH '--run_escape true' if you desired those results"
	}
	println escape_ch

    DAANALYSIS(identified_ch, new_opt_clust, params.reduced_dim)

    SUMMARYREPORT(
        COMPARATIVEANALYSIS.out.report,
        sling_ch,
        DAANALYSIS.out.report,
        escape_ch,
        tempora_ch,
        "${params.outdir}/analysis/",
        new_opt_clust
        )
}
