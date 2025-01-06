process TEMPORAANALYSIS {
    debug true
    //cache 'deep'
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
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
        pattern: "*.svg"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val main_time
    val species
    val plots_format

    output:
    path "*.pdf"
    path "*.rds"
    val true, emit: report
    path "*.svg"

    script:
       """
        ${projectDir}/src/TemporaMain.R --i $integrated -main_time $main_time -species $species -plots_format $plots_format
       """
}