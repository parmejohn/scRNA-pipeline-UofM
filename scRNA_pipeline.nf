
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
    if [! -d ${params.outdir}/analysis/]; then
        mkdir ${params.outdir}/analysis/
        mkdir ${params.outdir}/analysis/data/
        mkdir ${params.outdir}/analysis/data/qc

        mkdir ${params.outdir}/analysis/plots
    fi
    """ 
}

process AMBIENTRNAREMOVAL {
    input:
    val params.indir

    output:
    path '${params.outdir}/analysis/qc/*'
    
    script:
       """

       """
}

workflow {
    INITIALIZEFOLDERS(params.outdir)

}