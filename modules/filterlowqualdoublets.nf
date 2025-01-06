process FILTERLOWQUALDOUBLETS {
    debug true

    publishDir (
        path: "$params.outdir/analysis/data/",
        mode: 'copy',
        overwrite: true,
        pattern: "*.rds"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/qc",
        mode: 'copy',
        overwrite: true,
        pattern: "*.pdf"
    )
    publishDir (
        path: "$params.outdir/analysis/plots/qc",
        mode: 'copy',
        overwrite: true,
        pattern: "*.svg"
    )

    containerOptions "--bind $params.bind"

    input:
    path ambient_rmv
    val species
    val conditions
    val atac
    val og
    val mitochondrial_percent_cutoff
    val plots_format

    output:
    path "*.pdf"
    path "se_filtered_singlets_list.rds", emit: se_filtered_singlets_list
    path "se_filtered_doublets_list.rds"
    path "se_filtered_list.rds"
    path "se_list_raw.rds"
    path "*.svg"
    
    script:
       """
        ${projectDir}/src/FilterLowQualityAndDoublets.R --i $ambient_rmv -species $species -coconditions $conditions -atac $atac -original_files $og -mito_cutoff $mitochondrial_percent_cutoff -plots_format $plots_format
       """
}
