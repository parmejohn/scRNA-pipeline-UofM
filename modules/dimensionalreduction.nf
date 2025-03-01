process DIMENSIONALREDUCTION {
    debug true
    cache 'deep'

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
    )
    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.txt"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.pdf"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.jpeg"
    )
    containerOptions "--bind $params.bind"

    input:
    path integrated
    val clusters_optimal
    val resolution
    val reduced_dim
    val plots_format

    output:
    path "*.pdf"
    path "se_integrated_dimred.rds", emit: se_integrated_dimred
    path "optimal_clusters.txt", emit: clusters_optimal_n
	  val true, emit: report
	  path "*.jpeg"
    
    script:
       """
        ${projectDir}/src/SeuratDimRedMain.R --i $integrated -clusters_optimal $clusters_optimal -resolution $resolution -reduced_dim $reduced_dim -plots_format $plots_format
       """
}
