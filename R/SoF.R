#'Calculate Wicoxon test for 2 genes
#'
#'@param gene.pair data.frame; pairs of genes to be tested. First column - gene A, second column - gene B
#'@param alt data.frame; gene alterations in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@param expression data.frame; gene expression in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@importFrom stats wilcox.test
#'@importFrom dplyr select
#'@export

runWilcoxonFor2Genes <- function(gene.pair, alt, expression){
    gene_name1 <- as.vector(unlist(gene.pair[1]))
    gene_name2 <- as.vector(unlist(gene.pair[2]))
    #print(paste("runWilcoxonFor2Genes", gene_name1, ",", gene_name2, sep = " "))
    alt_gene1 <- alt %>% 
        dplyr::select(all_of(gene_name1))
    exp_gene2 <- expression %>% 
        dplyr::select(all_of(gene_name2))
    exp_gene2_alt_gene1_1 <- exp_gene2[alt_gene1 == 1]
    exp_gene2_alt_gene1_0 <- exp_gene2[alt_gene1 == 0]
    res <- stats::wilcox.test(exp_gene2_alt_gene1_0, exp_gene2_alt_gene1_1, alternative = c("less"))
    results <- cbind(gene.pair, res$p.value) %>% as.data.frame()
    names(results) <- c("geneA", "geneB", "p.value")
    return(results)
}

#'Survival of the fittest (SoF)
#'
#'@param gene.pair data.frame; pairs of genes to be tested. First column - gene A, second column - gene B
#'@param alt data.frame; gene alterations in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@param expression data.frame; gene expression in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@param cores number of cores to be used; default 1
#'@importFrom dplyr select filter rename arrange
#'@importFrom stats p.adjust
#'@importFrom parallel mclapply
#'@importFrom data.table rbindlist
#'@importFrom tibble column_to_rownames
#'@export

SoF <- function(gene.pair, alt, expression, cores = 1){
   
    gene_name1 <- as.vector(unlist(unique(gene.pair[1])))
    gene_name2 <- as.vector(unlist(unique(gene.pair[2])))
    alt <- alt %>% 
        dplyr::select(any_of(gene_name1)) %>% 
        dplyr::filter(rownames(alt) %in% rownames(expression))
    expression <- expression %>% 
        dplyr::select(any_of(gene_name2)) %>% 
        dplyr::filter(rownames(expression) %in% rownames(alt))
    expression <- expression[order(rownames(expression)), , drop=FALSE]
    alt <- alt[order(rownames(alt)), , drop=FALSE]
    
    gene.pairs <- split(gene.pair, row.names(gene.pair))
    pvalues2Genes <- parallel::mclapply(gene.pairs, mc.cores = cores, function(x) {runWilcoxonFor2Genes(x, alt, expression)})
    output <- data.table::rbindlist(pvalues2Genes) %>% as.data.frame()
    p.vals.test.bh <- stats::p.adjust(output$p.value, method = "BH") %>% as.data.frame()
    p.values <- merge(output, p.vals.test.bh, by = "row.names") %>%
        column_to_rownames(var = "Row.names") %>%
        dplyr::rename("p.value.bh" = ".") %>%
        dplyr::arrange(p.value.bh, p.value)

    return(p.values)
}
