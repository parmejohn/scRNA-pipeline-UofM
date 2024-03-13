
// Required inputs; will be checked in the R scripts
params.indir = ''
params.outdir = ''
params.species= ''
params.bind = ''

//Optional inputs
params.reference_seurat = ''
params.clusters_optimal = 0
params.resolution = 1
params.beginning_cluster = ''

//Global Variables
//outdata_ch = Channel.value()
//plotdata_ch = Channel.value()

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

process AMBIENTRNAREMOVAL {
    //debug true
    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true
    )

    containerOptions "--bind $params.bind"

    input:
    val params.indir

    output:
    path '*'
    
    script:
       """
        ${projectDir}/src/AmbientRNARemovalMain.R --i ${params.indir}
       """
}

process FILTERLOWQUALDOUBLETS {
    debug true

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/qc",
        mode: 'copy',
        overwrite: true,
        pattern: "*.pdf"
    )

    containerOptions "--bind $params.bind"

    input:
    path ambient_rmv

    output:
    path "*.pdf"
    path "se_filtered_singlets_list.rds", emit: se_filtered_singlets_list
    path "se_filtered_list.rds"
    path "se_list_raw.rds"
    
    script:
       """
        ${projectDir}/src/FilterLowQualityAndDoublets.R --i $ambient_rmv
       """
}

process INTEGRATEDATA {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true
    )

    containerOptions "--bind $params.bind"

    input:
    val rds

    output:
    path 'se_integrated.rds'
    
    script:
       """
        ${projectDir}/src/IntegrateSamplesMain.R --i $rds
       """
}

workflow {
    INITIALIZEFOLDERS(params.outdir)
    ambient_ch = AMBIENTRNAREMOVAL(params.indir) // array of conditions
    // ambient_ch.view() // outputs just the qc folder, for downstream qc, all samples need to be loaded as a list
    

    FILTERLOWQUALDOUBLETS(ambient_ch)
    filtered_ch = FILTERLOWQUALDOUBLETS.out.se_filtered_singlets_list
    //filtered_ch = Channel.fromPath( "${params.outdir}/analysis/data/se_filtered_singlets_list.rds", checkIfExists: true)
    
    // integrated seurat object
    integrated_ch = INTEGRATEDATA(filtered_ch)
    //integrated_ch.view()

}