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

    output:
    path '*.html'

    script:
        """
        #!/usr/local/bin/Rscript

        rmarkdown::render("${projectDIr}/src/SummaryReport.R")
        """

}