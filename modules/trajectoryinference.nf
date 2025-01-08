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
    publishDir (
        path: "$params.outdir/analysis/plots/ti",
        mode: 'copy',
        overwrite: true,
        pattern: "*.jpeg"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val beginning_cluster
    val plots_format

    output:
    path "*.pdf", emit: slingshot_out_pdf
    path "*.txt", optional: true
    path "*.rds", optional: true
    path "*.jpeg"

    // val true, emit: report
    
    script:
       """
        ${projectDir}/src/TrajectoryInferenceMain.R --i $integrated -beginning_cluster $beginning_cluster -plots_format $plots_format
       """
}