process SUMMARYREPORT{
        publishDir (
        path: "$params.outdir/analysis/",
        mode: 'copy',
        overwrite: true
    )

    containerOptions "--bind $params.bind"

    input:
    val comparative_signal
    val identify_signal
    val trajectory_signal
    val da_signal
    val escape_signal
    val tempora_signal
    val psupertime_signal
    val cellchat_signal
    val atac_signal
    path analysis_dir
    val opt_clusters
    path output

    output:
    path "*.html"

    script:
        """
        cp ${projectDir}/src/SummaryReport.qmd $analysis_dir/SummaryReport.qmd
        quarto render $analysis_dir/SummaryReport.qmd -P data:$analysis_dir --to html
        rm $analysis_dir/SummaryReport.qmd
        """
}
