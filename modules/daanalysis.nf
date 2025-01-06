process DAANALYSIS {
    debug true
    cache 'deep'
	//cache false
	

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
    )
    publishDir (
        path: "$params.outdir/analysis/data/da/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.txt"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/da",
        mode: 'copy',
        overwrite: true,
        pattern: "*.pdf"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/da",
        mode: 'copy',
        overwrite: true,
        pattern: "*.svg"
    )

    containerOptions "--bind $params.bind"

    input:
    path integrated
    val clusters_optimal
    val reduced_dim
    val species
    val plots_format

    output:
    path "*.rds"
    path "*.pdf"
    path "*.txt"
    val true, emit: report
    path "*.svg"

    
    script:
       """
        ${projectDir}/src/DifferentialAbundanceMain.R --i $integrated -clusters_optimal $clusters_optimal -reduced_dim $reduced_dim -species $species -plots_format $plots_format
       """
}
