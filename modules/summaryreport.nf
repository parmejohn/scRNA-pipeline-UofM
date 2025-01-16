process SUMMARYREPORT{

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
    val opt_clusters

    script:
        """
        quarto render ${projectDir}/src/SummaryReport.qmd -P data:${params.outdir}/analysis/ --to html --output-dir ${params.outdir}/analysis/
        """
}
