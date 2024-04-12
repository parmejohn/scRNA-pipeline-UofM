process ESCAPEANALYSIS {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/gsea/",
        mode: 'copy',
        overwrite: true,
        pattern: "esc*"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val species

    output:
    path "*.rds", emit: se_integrated_escape
    path "escape/"
    // val true, emit: report
    
    script:
       """
        ${projectDir}/src/EscapeMain.R --i $integrated --s $species
       """
}
