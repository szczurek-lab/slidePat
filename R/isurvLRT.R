#'Calculate genotypes for pair of genes
#'
#'@param gene1 character; gene A
#'@param gene2 character; gene B
#'@param data data.frame; gene alteration with survival status
#'@param expression data.frame; gene expression in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@param thr expression threshold
#'@importFrom dplyr select
#'@export

icalculatePairGenotype <- function(gene1, gene2, data, expression, thr){
    # if (!gene1%in%colnames(alt)){ print(paste( "iSurvLRT: Gene1 ", gene1, " not in alteration data, returning NA", sep=""))
    #     return(NA)}
    # if (!gene2%in%colnames(expression)){ print(paste( "iSurvLRT: Gene2 ", gene2, " not in expression data, returning NA", sep=""))
    #     return(NA)}
    gene1 <- data %>% 
        dplyr::select(one_of(gene1))
    gene2 <- expression %>% 
        dplyr::select(one_of(gene2))
    gene2 <- ifelse(gene2 < thr, 1, 0)
    genotype_data <- cbind(gene1, gene2)
    genotype_data$genotype <- paste(genotype_data[ ,1], genotype_data[ ,2])
    return(genotype_data)
}

#' Run isurvLRT for gene pairs
#'
#'@param gene.pair data.frame; pairs of genes to be tested. First column - gene A, second column - gene B
#'@param data data.frame; gene alteration with survival status
#'@param expression data.frame; gene expression in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@param thr expression threshold
#'@importFrom dplyr select
#'@importFrom rlang .data
#'@return result table
#'@export

irunFor2Genes <- function(gene.pair, data, expression, thr){
    gene_name1 <- as.vector(unlist(unique(gene.pair[1])))
    gene_name2 <- as.vector(unlist(unique(gene.pair[2])))
    pairGenotype <- icalculatePairGenotype(gene_name1, gene_name2, data, expression, thr)
    #    if (is.na(pairGenotype)){ return(rep( NA, 3))}
    gene1 <- pairGenotype[,1] %>% data.frame()
    gene2 <- pairGenotype[,2] %>% data.frame()
    alive <- data %>% dplyr::select(alive)
    glast_follow_up <- data %>% 
        dplyr::select(.data$GlastFollowUp)
    #print(paste("runFor2Genes",gene_name1,",",gene_name2,sep =" ", thr))
    results <- checkLinearityPairwise(gene1, gene2, alive, glast_follow_up)
    as.matrix(results)
}

#'Calcultate threshold
#'
#'@param gene.pair data.frame; pairs of genes to be tested. First column - gene A, second column - gene B
#'@param expression data.frame; gene expression in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@return data frame
#'@importFrom dplyr select
#'@importFrom stats quantile
#'@export

threshold <- function(gene.pair, expression){
    gene_name2 <- as.vector(unlist(gene.pair[2]))
    if (!gene_name2 %in% colnames(expression)){
        # print(paste( "iSurvLRT: Gene2 ", gene.name2, " not in expression data, returning NA", sep=""))
        return(as.numeric(rep(NA, 10)))}
    gene2 <- expression %>% 
        dplyr::select(one_of(gene_name2))
    thr <- stats::quantile(gene2[,1], probs = seq(0.05, 0.5, 0.05))
}

#'Test the pair of genes
#'
#'@param gene.pair data.frame; pairs of genes to be tested. First column - gene A, second column - gene B
#'@param data data data.frame; gene alteration with survival status
#'@param expression data.frame; gene expression in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@importFrom dplyr select
#'@importFrom rlang .data
#'@return data frame
#'@export

isurvLRTPair <- function(gene.pair, data, expression) {
    # gene_name1 <- as.vector(unlist(unique(gene.pair[1])))
    # gene_name2 <- as.vector(unlist(unique(gene.pair[2])))
    # data <- data %>% dplyr::select(any_of(gene_name1), .data$alive, .data$GlastFollowUp)
    # expression <- expression %>% 
    #     dplyr::select(any_of(gene_name2))
    thrs <- threshold(gene.pair, expression)
    output <- t(rep(NA, 6)) %>% as.data.frame()
    names(output) <- c('statistic', 'pvalue', 'SLflag', 'CLflag', 'effect', 'thr_min')
    if(!is.na(thrs)){
        results <- t(sapply(thrs, function(thr){irunFor2Genes(gene.pair, data, expression, thr)}))
        colnames(results) <- c('statistic', 'pvalue', 'SLflag', 'CLflag', 'effect')
        if (!all(is.na( results[, "pvalue"] ))) {
            ind_min <- which.min(results[,"pvalue"])
            results_min <- results[ind_min,]
            thr_min <- thrs[ind_min]
            output <- as.data.frame(cbind(t(results_min), thr_min))
        }
        else{
            output <- output
        }
    }
    isurvLRT_results <- cbind(gene.pair, output) %>% data.frame()
    return(isurvLRT_results)
}

#'Test the pairs of genes
#'
#'@param gene.pair data.frame; pairs of genes to be tested. First column - gene A, second column - gene B
#'@param alt data.frame; gene alterations in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@param survival data.frame; survival status 1 - alive, 0 - dead
#'@param expression data.frame; gene expression in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@param cores number of cores to be used; default 1
#'@importFrom dplyr select filter arrange rename_with
#'@importFrom rlang .data
#'@importFrom parallel mclapply
#'@return data frame
#'@export

isurvLRT <- function(gene.pair, alt, survival, expression, cores = 1) {
    gene_name1 <- as.vector(unlist(unique(gene.pair[1])))
    gene_name2 <- as.vector(unlist(unique(gene.pair[2])))
    alt <- alt %>% 
        dplyr::select(any_of(gene_name1)) %>% 
        dplyr::filter(rownames(alt) %in% rownames(expression))
    data <- inputDataSurvLRT(alt, survival)
    expression <- expression %>% 
        dplyr::select(any_of(gene_name2)) %>% 
        dplyr::filter(rownames(expression) %in% rownames(data))
    expression <- expression[order(rownames(expression)), , drop=FALSE]
    data <- data[order(rownames(data)), , drop=FALSE]
    gene.pair <- gene.pair %>%
        dplyr::rename_with(~ c('geneA', 'geneB'), 1:2) %>%
        dplyr::filter(geneA %in% colnames(data) & geneB %in% colnames(expression))
    gene.pairs <- split(gene.pair, row.names(gene.pair))
    results <- parallel::mclapply(gene.pairs, mc.cores = cores, function(x) {isurvLRTPair(x, data, expression)})
    output <- as.data.frame(data.table::rbindlist(results)) %>%
        dplyr::arrange(-SLflag, pvalue)
}
