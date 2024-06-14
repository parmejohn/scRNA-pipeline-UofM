process CELLCHAT {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/plots/",
        mode: 'copy',
        overwrite: true,
        pattern: "*plots"
    )
    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*data"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val species

    output:
    path "*plots"
    path "*data"
    val true, emit: report

    script:
       """
        ${projectDir}/src/CellChatMain.R --i $integrated --s $species
       """
}