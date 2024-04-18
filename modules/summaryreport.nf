process SUMMARYREPORT{
        publishDir (
        path: "$params.outdir/analysis/",
        mode: 'copy',
        overwrite: true
    )

    containerOptions "--bind $params.bind"

    input:
    val comparative_signal
    val trajectory_signal
    val da_signal
    val escape_signal
    val tempora_signal
    val analysis_dir
    val opt_clusters

    output:
    path "*.html"

    script:
        """
        #!/usr/local/bin/Rscript

        rmarkdown::render("${projectDir}/src/SummaryReport.R", params = list(data = "$analysis_dir", opt_c = "$opt_clusters"), output_dir = ".")
        """

}
