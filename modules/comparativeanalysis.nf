process COMPARATIVEANALYSIS {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/plots/gsea/comparative",
        mode: 'copy',
        overwrite: true,
        pattern: "gsea_cluster*.pdf"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/deseq2",
        mode: 'copy',
        overwrite: true,
        pattern: "deseq2*.pdf"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/gsea/comparative",
        mode: 'copy',
        overwrite: true,
        pattern: "gsea_cluster*.svg"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/deseq2",
        mode: 'copy',
        overwrite: true,
        pattern: "deseq2*.svg"
    )
	  publishDir (
        path: "$params.outdir/analysis/data/deseq2",
        mode: 'copy',
        overwrite: true,
        pattern: "*.txt"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val species
    val plots_format

    output:
    path "*.pdf"
    path "*.txt" // DEGs txt files
    val true, emit: report
    path "*.svg"

    script:
       """
        ${projectDir}/src/ComparativeAnalysisMain.R --i $integrated --s $species -plots_format $plots_format
       """
}