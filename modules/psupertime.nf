process PSUPERTIME {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/plots/ti/",
        mode: 'copy',
        overwrite: true,
        pattern: "*plots"
    )
    publishDir (
        path: "$params.outdir/analysis/data/ti/",
        mode: 'copy',
        overwrite: true,
        pattern: "*data"
    )
    
    containerOptions "--bind $params.bind"

    input:
    path integrated
    val main_time
    val plots_format

    output:
    path "psupertime_plots"
    path "psupertime_data"
    val true, emit: report

    
    script:
       """
        ${projectDir}/src/PsupertimeMain.R --i $integrated -main_time $main_time -plots_format $plots_format
       """
}