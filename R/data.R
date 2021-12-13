#' @title scRNA Gene Expression Example Data
#' @description scRNA data with labels for cancer types (from the UCI ML repository).
#'     Each row is a sample, and each column is a gene. 
#' @format A data frame with 801 rows and 2000 variables:
#' \describe{
#'   \item{gene_xxxx}{gene expression data for gene_xxxx for each sample}
#'}
#' @source \url{https://archive.ics.uci.edu/ml/datasets/gene+expression+cancer+RNA-Seq#}
"tumor_reduced"



#' @title scRNA Gene Expression Example Labels
#' @description Labels for cancer types that correspond to the rows of tumor_reduced.
#' @format A 801-length vector
#' \describe{
#'   \item{label}{one of BRCA, KIRC, COAD, LUAD, PRAD}
#'}
#' @source \url{https://archive.ics.uci.edu/ml/datasets/gene+expression+cancer+RNA-Seq#}
"TC"