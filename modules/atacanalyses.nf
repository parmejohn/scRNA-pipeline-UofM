process ATACANALYSES {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*data"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/",
        mode: 'copy',
        overwrite: true,
        pattern: "*plots"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val species
    val plots_format

    output:
    path "*plots"
    path "*data"
    val true, emit: report
    
    script:
       """
        ${projectDir}/src/AtacAnalyses.R --i $integrated --s $species -plots_format $plots_format
       """
}
