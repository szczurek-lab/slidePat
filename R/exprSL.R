#'Formatting the correlation matrix into a table with 4 columns containing :
#'Column 1 : row names (variable 1 for the correlation test)
#'Column 2 : clumn names (variable 2 for the correlation test)
#'Column 3 : the correlation coefficients
#'Column 4 : the p-values of the correlation
#'
#'@param cormatrix matrix of the correlation coefficients
#'@param pmatrix  matrix of the correlation p-values
#'@importFrom tibble rownames_to_column
#'@importFrom tidyr gather
#'@importFrom dplyr left_join
#'
#'@export

flatCorMatrix <- function(cormatrix, pmatrix) {
   cor_r <- tibble::rownames_to_column(as.data.frame(cormatrix), var = "row")
   cor_r <- tidyr::gather(cor_r, "column", "cor", -1)
   cor_p <- tibble::rownames_to_column(as.data.frame(pmatrix), var = "row")
   cor_p <- tidyr::gather(cor_p, "column", "p", -1)
   cor_p_matrix <- dplyr::left_join(cor_r, cor_p, by = c("row", "column"))
}

#'Get the correlation matrix of gene expression for pair of genes
#'
#'@param gene.pair pairs of genes to be tested
#'@param expression data.frame; gene expression in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@importFrom Hmisc rcorr
#'@importFrom dplyr filter arrange
#'@importFrom rlang .data
#'@export

corrExprGenePair <- function(gene.pair, expression){
    gene_name1 <- as.vector(unlist(unique(gene.pair[1])))
    gene_name2 <- as.vector(unlist(unique(gene.pair[2])))
    gene_list <- append(gene_name1, gene_name2)
    data <- expression %>%
        dplyr::select(any_of(gene_list))
    corr <- Hmisc::rcorr(as.matrix(data, type = c("spearman")))
    corr_results<- flatCorMatrix(corr$r, corr$P)
    colnames(corr_results) <- c("geneA", "geneB","corr","pvalue")
    corr_results <- corr_results %>% 
        dplyr::filter(.data$geneA %in% gene_name1 & .data$geneB %in% gene_name2) %>% 
        dplyr::arrange(-corr) 
    return(corr_results)
}
