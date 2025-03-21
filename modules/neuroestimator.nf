process NEUROESTIMATOR{

    containerOptions "--bind $params.bind"
    cache 'lenient'
    debug true
    label 'neuroestimator'
        
    publishDir (
      path: "$params.outdir/analysis/data/",
      mode: 'copy',
      overwrite: true,
      pattern:'*'
    )
    
    input:
    path integrated_ch
    val species
    
    output:
    path "neuroestimator_results.txt"
    
    script:
    """
	  ${projectDir}/src/Neuroestimator.R $integrated_ch $species
    """
}

process NEUROESTIMATORPLOT{

    containerOptions "--bind $params.bind"
    cache 'lenient'
    debug true
    label 'neuroestimator_plot'
    
    publishDir (
      path: "$params.outdir/analysis/plots/neuroestimator/",
      mode: 'copy',
      overwrite: true,
      pattern:'*'
    )
    
    input:
    path integrated_ch
    path neuroestimator_results
	val plots_format

    output:
    path "*.pdf", emit: plot
    path "*.jpeg"
    
    script:
    """
	  ${projectDir}/src/NeuroestimatorPlot.R $integrated_ch $neuroestimator_results $plots_format
    """
}
