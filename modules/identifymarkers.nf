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
    publishDir (
        path: "$params.outdir/analysis/plots/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.svg"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val clusters_optimal
    val reference_seurat
    val plots_format

    output:
    path "*.pdf"
    path "*.svg"
    path "se_markers_presto_integrated.txt"
    path "se_integrated_auto_label.rds", emit: se_integrated_auto_label, optional: true
    val true, emit: report
    
    script:
       """
        ${projectDir}/src/IdentifyCellMarkersMain.R --i $integrated -clusters_optimal $clusters_optimal -reference_seurat $reference_seurat -plots_format $plots_format
       """
}