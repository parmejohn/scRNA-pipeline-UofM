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
    val reference_seurat

    output:
    path "*.pdf"
    path "se_markers_presto_integrated.txt"
    path "se_integrated_auto_label.rds", emit: se_integrated_auto_label, optional: true
    
    script:
       """
        ${projectDir}/src/IdentifyCellMarkersMain.R --i $integrated -clusters_optimal $clusters_optimal -reference_seurat $reference_seurat
       """
}