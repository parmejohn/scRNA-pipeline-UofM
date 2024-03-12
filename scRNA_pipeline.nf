
// Required inputs; will be checked in the R scripts
params.indir = ''
params.outdir = ''
params.species= ''

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
    debug true
    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true
    )

    containerOptions '--bind /home/projects/,/home/projects/CIO/yard/run_cellranger_count'

    input:
    val params.indir

    output:
    path '*'
    
    script:
       """
        ${projectDir}/src/AmbientRNARemovalMain.R --i ${params.indir}
       """
}

workflow {
    INITIALIZEFOLDERS(params.outdir)
    ambient_ch = AMBIENTRNAREMOVAL(params.indir) // array of conditions
    ambient_ch.view()

    
}