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
    val reduced_dim

    output:
    path 'se_integrated.rds'
    
    script:
       """
        ${projectDir}/src/IntegrateSamplesMain.R --i $rds -reduced_dim $reduced_dim
       """
}