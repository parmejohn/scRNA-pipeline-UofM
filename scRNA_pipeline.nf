
// Required inputs; will be checked in the R scripts
params.indir = ''
params.outdir = ''
params.species= ''
params.bind = ''

//Optional inputs
params.reference_seurat = 'none'
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

process DIMENSIONALREDUCTION {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
    )
    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.txt"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.pdf"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val clusters_optimal
    val resolution

    output:
    path "*.pdf"
    path "se_integrated_dimred.rds", emit: se_integrated_dimred
    path "optimal_clusters.txt", emit: clusters_optimal_n
    
    script:
       """
        ${projectDir}/src/SeuratDimRedMain.R --i $integrated -clusters_optimal $clusters_optimal -resolution $resolution
       """
}

process IDENTIFYMARKERS {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
    )
    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.txt"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.pdf"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val clusters_optimal
    path reference_seurat

    output:
    path "*.pdf"
    path "se_markers_presto_integrated.txt"
    path "se_integrated_auto_label.rds", emit: se_integrated_auto_label, optional: true
    
    script:
       """
        ${projectDir}/src/IdentifyCellMarkersMain.R --i $integrated -clusters_optimal $clusters_optimal -reference_seurat $reference_seurat
       """
}

process TRAJECTORYINFERENCE {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
    )
    publishDir (
        path: "$params.outdir/analysis/data/ti",
        mode: 'copy',
        overwrite: true,
        pattern: "*.txt"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/ti",
        mode: 'copy',
        overwrite: true,
        pattern: "*.pdf"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val beginning_cluster

    output:
    path "*.pdf"
    path "*.txt"
    path "*.rds"
    
    script:
       """
        ${projectDir}/src/TrajectoryInferenceMain.R --i $integrated -beginning_cluster $beginning_cluster
       """
}

process COMPARATIVEANALYSIS {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/plots/gsea/comparative",
        mode: 'copy',
        overwrite: true,
        pattern: "gsea_cluster*.pdf"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/deseq2",
        mode: 'copy',
        overwrite: true,
        pattern: "deseq2*.pdf"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val species

    output:
    path "*.pdf"
    //path "*.txt" // should include DEGs txt file?
    
    script:
       """
        ${projectDir}/src/ComparativeAnalysisMain.R --i $integrated --s $species
       """
}

process ESCAPEANALYSIS {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/gsea/",
        mode: 'copy',
        overwrite: true,
        pattern: "*escape/"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val species

    output:
    path "*.rds"
    path "*escape/"
    
    script:
       """
        ${projectDir}/src/EscapeMain.R --i $integrated --s $species
       """
}

process DAANALYSIS {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/da",
        mode: 'copy',
        overwrite: true,
        pattern: "*.pdf"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val clusters_optimal

    output:
    path "*.rds"
    path "*.pdf"
    
    script:
       """
        ${projectDir}/src/DifferentialAbundanceMain.R --i $integrated -clusters_optimal $clusters_optimal
       """
}

workflow {
    // TODO: change all val to path after -> dont know why some of these are working
    INITIALIZEFOLDERS(params.outdir)
    ambient_ch = AMBIENTRNAREMOVAL(params.indir) // array of conditions
    // ambient_ch.view() // outputs just the qc folder, for downstream qc, all samples need to be loaded as a list
    

    FILTERLOWQUALDOUBLETS(ambient_ch)
    filtered_ch = FILTERLOWQUALDOUBLETS.out.se_filtered_singlets_list
    //filtered_ch = Channel.fromPath( "${params.outdir}/analysis/data/se_filtered_singlets_list.rds", checkIfExists: true)
    
    // integrated seurat object
    integrated_ch = INTEGRATEDATA(filtered_ch)
    //integrated_ch.view()

    // seurat dim red
    DIMENSIONALREDUCTION(integrated_ch, params.clusters_optimal, params.resolution)
    dimred_ch = DIMENSIONALREDUCTION.out.se_integrated_dimred
    //DIMENSIONALREDUCTION.out.clusters_optimal_n.view()
    dimred_ch.view()
    opt_clust_file = file("${params.outdir}/analysis/data/optimal_clusters.txt", checkIfExists: true)
    new_opt_clust = opt_clust_file.text

    // identify cell markers
    identified_ch = Channel.of()
    println params.reference_seurat
    if (params.reference_seurat != 'none'){
        println "true"
        println params.reference_seurat
        IDENTIFYMARKERS(dimred_ch, new_opt_clust, params.reference_seurat)
        identified_ch = IDENTIFYMARKERS.out.se_integrated_auto_label
    } else {
        println "false"
        println params.reference_seurat
        IDENTIFYMARKERS(dimred_ch, new_opt_clust, params.reference_seurat)
        identified_ch = dimred_ch
    }
    identified_ch.view() // confirm the output
    //TRAJECTORYINFERENCE(identified_ch, params.beginning_cluster)


    COMPARATIVEANALYSIS(identified_ch, params.species)
    //ESCAPEANALYSIS(identified_ch, params.species)
    DAANALYSIS(identified_ch, new_opt_clust)

}