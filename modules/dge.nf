/*
 * Different Gene Expression Analysis Processes
 */

/*
 * Script written by NASA, parameterized by dataset
 */

process DGE_BY_DESEQ2 {
  conda 'envs/r_deseq2.yml'
  publishDir "${params.publishDirPath}/dge"

  input:
    tuple path(ISA_zip), path(organisms_csv), path(Rsem_gene_counts)
  output:
    tuple path("Unnormalized_Counts.csv"), path("Normalized_Counts.csv")
  script:
    """
    # create working directory
    mkdir workDir

    # create rsem counts directory and move them into it
    mkdir workDir/RSEM_GENE_COUNTS
    mv $Rsem_gene_counts workDir/RSEM_GENE_COUNTS

    # create metadata dir and move metadata in
    mkdir workDir/metaDir
    mv $ISA_zip $organisms_csv workDir/metaDir

    # cd into workdir and create output directories
    cd workDir
    mkdir norm_counts_output
    mkdir dge_output

    # run the script with R
    GLDS-104-norm_DGE_analysis.R
    """
}
