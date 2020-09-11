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
    # create rsem counts directory and move them into it
    mkdir RSEM_GENE_COUNTS
    mv $Rsem_gene_counts RSEM_GENE_COUNTS

    # create metadata dir and move metadata to script expected places
    mkdir metaDir
    mv $ISA_zip metaDir

    # create output directories
    mkdir norm_counts_output
    mkdir dge_output

    # run the script with R
    GLDS-104-norm_DGE_analysis.R
    """
}
