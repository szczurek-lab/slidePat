#'Formatting the correlation matrix into a table with 4 columns containing :
#'Column 1 : row names (variable 1 for the correlation test)
#'Column 2 : clumn names (variable 2 for the correlation test)
#'Column 3 : the correlation coefficients
#'Column 4 : the p-values of the correlation
#'
#'@param cormatrix matrix of the correlation coefficients
#'@param pmatrix  matrix of the correlation p-values
#'@importFrom tibble rownames_to_column
#'@importFrom tidyr pivot_longer
#'@importFrom dplyr left_join
#'
#'@export

flatCorMatrix <- function(cormatrix, pmatrix) {
   cor_r <- cormatrix %>% 
       as.data.frame() %>%
       tibble::rownames_to_column(var = "row") %>%
       tidyr::pivot_longer(!row, names_to = "column", values_to = "cor")
   cor_p <- pmatrix %>%
       as.data.frame() %>%
       tibble::rownames_to_column(var = "row") %>%
       tidyr::pivot_longer(!row, names_to = "column", values_to = "p")
   cor_matrix <- dplyr::left_join(cor_r, cor_p)
   return(cor_matrix)
}

#'Get the correlation matrix of gene expression for pair of genes
#'
#'@param gene.pair data.frame; pairs of genes to be tested. First column - gene A, second column - gene B
#'@param expression data.frame; gene expression in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@importFrom Hmisc rcorr
#'@importFrom dplyr filter select
#'@export

corrExprGenePair <- function(gene.pair, expression){
    geneA <- as.vector(unlist(unique(gene.pair[1])))
    geneB <- as.vector(unlist(unique(gene.pair[2])))
    gene_list <- base::append(geneA, geneB) %>% 
        unique()
    
    if(!(geneA %in% gene_list) || (geneB %in% gene_list)) {
        return(data.frame("geneA" = geneA, "geneB" = geneB, "corr" = NA, "pvalue" = NA))
    }
    
    data <- expression %>%
        dplyr::select(any_of(gene_list))
    
    corr <- Hmisc::rcorr(as.matrix(data, type = c("spearman")))
    corr_results<- flatCorMatrix(corr$r, corr$P)
    colnames(corr_results) <- c("geneA", "geneB","corr","pvalue")
    corr_results <- corr_results %>% 
        dplyr::filter(geneA %in% geneA & geneB %in% geneB) 
    return(corr_results)
}
