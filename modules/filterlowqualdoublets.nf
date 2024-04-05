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

    containerOptions "--bind $params.bind"

    input:
    path ambient_rmv
    val species

    output:
    path "*.pdf"
    path "se_filtered_singlets_list.rds", emit: se_filtered_singlets_list
    path "se_filtered_list.rds"
    path "se_list_raw.rds"
    
    script:
       """
        ${projectDir}/src/FilterLowQualityAndDoublets.R --i $ambient_rmv -species $species
       """
}