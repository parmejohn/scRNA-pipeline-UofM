process DAANALYSIS {
    debug true
    cache 'deep'
//	cache false

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/da",
        mode: 'copy',
        overwrite: true,
        pattern: "*.pdf"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val clusters_optimal

    output:
    path "*.rds"
    path "*.pdf"
    
    script:
       """
        ${projectDir}/src/DifferentialAbundanceMain.R --i $integrated -clusters_optimal $clusters_optimal
       """
}