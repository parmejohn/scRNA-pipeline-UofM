process INTEGRATEDATA {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true
    )

    containerOptions "--bind $params.bind"

    input:
    val rds

    output:
    path 'se_integrated.rds'
    
    script:
       """
        ${projectDir}/src/IntegrateSamplesMain.R --i $rds
       """
}