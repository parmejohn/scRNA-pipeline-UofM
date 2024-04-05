process AMBIENTRNAREMOVAL {
    //debug true
    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true
    )

    containerOptions "--bind $params.bind"

    input:
    val indir
    val test_data

    output:
    path '*'
    
    script:
       """
        ${projectDir}/src/AmbientRNARemovalMain.R --i $indir -test_data $test_data
       """
}
