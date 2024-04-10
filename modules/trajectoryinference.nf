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
    path "*.pdf", emit: slingshot_out_pdf
    path "*.txt", optional: true
    path "*.rds", optional: true
    // val true, emit: report
    
    script:
       """
        ${projectDir}/src/TrajectoryInferenceMain.R --i $integrated -beginning_cluster $beginning_cluster
       """
}