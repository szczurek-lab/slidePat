#'Prepare pair of genes to be tested
#'
#'@param geneA character; list of gene A
#'@param geneB character; list of gene B
#'@importFrom dplyr mutate_if rename filter
#'@return data frame
#'@export

genePairCombination <- function(geneA, geneB){
    gene_pairs <- expand.grid(geneA, geneB) %>% 
        dplyr::mutate_if(is.factor, as.character) %>%
        dplyr::rename(geneA = Var1, 
                      geneB = Var2) %>% 
        filter(!geneA == geneB) #remove pairs with equal genes in the pair
    return(gene_pairs)
}

