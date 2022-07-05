#'Calculate likelihood ratio test
#'
#'@param S count of patients in group with a specific mutational and survival status
#'@param H precalculated sums of logs of survival function values for patients within groups with specific mutations
#'@importFrom stats pchisq
#'@export

LRT <- function(S, H) {
    if(any(S < 5) || any(H == 0) || any(is.nan(H))) {
        
        return(c(NA, NA, NA))
    }
    TS <- numeric(3)
    TS[1] <- S[1] - S[4]
    TS[2] <- S[2] + S[4]
    TS[3] <- S[3] + S[4]
    delta11 <- -(2*H[2]*H[3]*TS[2] - H[1]*H[4]*TS[2] + H[1]*H[4]*TS[3] + H[2]*H[3]*TS[1] - sqrt((2*H[2]*H[3]*TS[2] - H[1]*H[4]*TS[2] + H[1]*H[4]*TS[3] + H[2]*H[3]*TS[1])^2-4*(H[2]*H[2]*H[3] - H[1]*H[2]*H[4])*(H[3]*TS[2]*TS[2] + H[3]*TS[1]*TS[2]))) / (2*(H[2]*H[2]*H[3] - H[1]*H[2]*H[4]))
    delta12 <- -(2*H[2]*H[3]*TS[2] - H[1]*H[4]*TS[2] + H[1]*H[4]*TS[3] + H[2]*H[3]*TS[1] + sqrt((2*H[2]*H[3]*TS[2] - H[1]*H[4]*TS[2] + H[1]*H[4]*TS[3] + H[2]*H[3]*TS[1])^2-4*(H[2]*H[2]*H[3] - H[1]*H[2]*H[4])*(H[3]*TS[2]*TS[2] + H[3]*TS[1]*TS[2]))) / (2*(H[2]*H[2]*H[3] - H[1]*H[2]*H[4]))
    delta21 <- (TS[2] - TS[3] + H[2] * delta11) / H[3]
    delta22 <- (TS[2] - TS[3] + H[2] * delta12) / H[3]
    delta01 <- -H[4] * delta11 * delta21 / (TS[2] + H[2] * delta11)
    delta02 <- -H[4] * delta12 * delta22 / (TS[2] + H[2] * delta12)
    delta31 <- delta11 * delta21 / delta01
    delta32 <- delta12 * delta22 / delta02
    ml0 <- -S[1] / H[1]
    ml1 <- -S[2] / H[2]
    ml2 <- -S[3] / H[3]
    ml3 <- -S[4] / H[4]
    if (any(c(delta01, delta11, delta21) < 0)) {
        T1 <- -Inf
    }
    else {
        T1 <- sum(S[1:4]*(log(c(delta01, delta11, delta21, delta31)) - log(c(ml0, ml1, ml2, ml3)))) + sum(H[1:4]*c(delta01 - ml0, delta11 - ml1, delta21 - ml2, delta31 - ml3))
    }
    if (any(c(delta02, delta12, delta22) < 0)) {
        T2 <- -Inf
    }
    else {
        T2 <- sum(S[1:4]*(log(c(delta02, delta12, delta22, delta32)) - log(c(ml0, ml1, ml2, ml3)))) + sum(H[1:4]*c(delta02 - ml0, delta12 - ml1, delta22 - ml2, delta32 - ml3))
    }
    if ((any(c(delta01, delta11, delta21) < 0)) & (any(c(delta02, delta12, delta22) < 0)) ){
        return( return(c(NA, NA, NA)))
    }
    c(-2*max(T1,T2), (1 - stats::pchisq(-2*max(T1,T2), df = 1)), 2*as.numeric(ml0*ml3 < ml1*ml2) - 1)
}

#'Test for pairwaise epistasis type
#'
#'@param gene1 character; gene A alterations
#'@param gene2 character; gene B alterations
#'@param alive status to indicate alive or dead; 1 - alive, 0 - dead
#'@param GlastFollowUp the probablity of surviving more than the corresponding number of months
#'@export

checkLinearityPairwise <- function(gene1, gene2, alive, GlastFollowUp) {
    S <- numeric(4)
    H <- numeric(4)
    case0 <- (gene1 == 0) & (gene2 == 0)
    case1 <- (gene1 == 0) & (gene2 == 1)
    case2 <- (gene1 == 1) & (gene2 == 0)
    case3 <- (gene1 == 1) & (gene2 == 1)
    S[1] <- sum(case0 & alive == 0)
    S[2] <- sum(case1 & alive == 0)
    S[3] <- sum(case2 & alive == 0)
    S[4] <- sum(case3 & alive == 0)
    H[1] <- sum(log(GlastFollowUp[case0]))
    H[2] <- sum(log(GlastFollowUp[case1]))
    H[3] <- sum(log(GlastFollowUp[case2]))
    H[4] <- sum(log(GlastFollowUp[case3]))
    output <- LRT(S, H)
    names(output) <- c('statistic', 'pvalue', 'SLflag')
    output
}

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
#'@param gene.pair pair of genes to be tested
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
#'@param gene.pair pair of genes to be tested
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
#'@param gene.pair pair of genes to be tested
#'@param data data data.frame; gene alteration with survival status
#'@param expression data.frame; gene expression in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@importFrom dplyr select
#'@importFrom rlang .data
#'@return data frame
#'@export

isurvLRTPair <- function(gene.pair, data, expression) {
    gene_name1 <- as.vector(unlist(unique(gene.pair[1])))
    gene_name2 <- as.vector(unlist(unique(gene.pair[2])))
    data <- data %>% dplyr::select(all_of(gene_name1), .data$alive, .data$GlastFollowUp)
    expression <- expression %>% 
        dplyr::select(all_of(gene_name2))
    thrs <- threshold(gene.pair, expression)
    output <- t(rep(NA, 4)) %>% as.data.frame()
    names(output) <- c('statistic', 'pvalue', 'SLflag', 'thr_min')
    if(!is.na(thrs)){
        results <- t(sapply(thrs, function(thr){irunFor2Genes(gene.pair, data, expression, thr)}))
        colnames(results) <- c('statistic', 'pvalue', 'SLflag')
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
#'@param gene.pair data.frame; pair of genes to be tested
#'@param alt data.frame; gene alterations in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@param survival data.frame; survival status 1 - alive, 0 - dead
#'@param expression data.frame; gene expression in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@param cores number of cores to be used; default 1
#'@import dplyr
#'@importFrom rlang .data
#'@importFrom parallel mclapply
#'@return data frame
#'@export

isurvLRT <- function(gene.pair, alt, survival, expression, cores = 1) {
    gene_name1 <- as.vector(unlist(unique(gene.pair[1])))
    gene_name2 <- as.vector(unlist(unique(gene.pair[2])))
    alt <- alt %>% dplyr::select_if(names(alt) %in% gene_name1)
    expression <- expression %>% 
        dplyr::select_if(names(expression) %in% gene_name2)
    data <- inputDataSurvLRT(alt, survival)
    expression <- expression %>% 
        dplyr::filter(rownames(expression) %in% rownames(data))
    data <- data %>% 
        dplyr::filter(rownames(data) %in% rownames(expression))
    expression <- expression[ order(rownames(expression)), , drop=FALSE]
    data <- data[ order(rownames(data)), , drop=FALSE]
    gene.pairs <- split(gene.pair, row.names(gene.pair))
    results <- parallel::mclapply(gene.pairs, mc.cores = cores, function(x) {isurvLRTPair(x, data, expression)})
    output <- as.data.frame(data.table::rbindlist(results))
    output <- output[order(-output$SLflag == 1, output$pvalue),]
}
