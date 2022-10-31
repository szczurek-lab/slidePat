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
#'@param alive status to indicate alive or dead 1 - alive, 0 - dead
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

    ml0 <- -S[1] / H[1]
    ml1 <- -S[2] / H[2]
    ml2 <- -S[3] / H[3]
    ml3 <- -S[4] / H[4]
    
    effect <- log(ml0) + log(ml3) - log(ml1) - log(ml2)
    CL_flag <- ifelse(ml0 > ml3, 1, -1) #CL_flag - clinicaly relevant 
    output <- c(output, CL_flag, effect)
    names(output) <- c('statistic', 'pvalue', 'SLflag', 'CLflag', 'effect')
    return(output)
}

#' Run SurvLRT for gene pairs
#'
#'@param gene.pair data.frame; pairs of genes to be tested. First column - gene A, second column - gene B
#'@param data data data.frame; gene alteration with survival status
#'@importFrom dplyr select
#'@return result table
#'@export

runFor2Genes <- function(gene.pair, data){
    gene_name1 <- as.vector(unlist(unique(gene.pair[1])))
    gene_name2 <- as.vector(unlist(unique(gene.pair[2])))
    #  print(paste("runFor2Genes",gene.name1,",",gene.name2,sep =" "))
    gene1 <- data %>% 
        dplyr::select(any_of(gene_name1))
    gene2 <- data %>% 
        dplyr::select(any_of(gene_name2))
    alive <- data %>% 
        dplyr::select(alive)
    glast_follow_up <- data %>% 
        dplyr::select(GlastFollowUp)
    results <- checkLinearityPairwise(gene1, gene2, alive, glast_follow_up) %>%
        t() %>%
        as.data.frame()
    results <- base::cbind(gene.pair, results)
    return(results)
}

#'Prepare alteration data with survival for survLRT and isurvLRT test
#'
#'@param alt data.frame gene alteration
#'@param survival data.frame input data with survival information
#'@importFrom tibble column_to_rownames
#'@return data frame
#'@export

inputDataSurvLRT <- function(alt, survival) {
    data <- base::merge(alt, survival, by = "row.names") %>%
        tibble::column_to_rownames(var = "Row.names")
}

#'Test the pair of genes
#'
#'@param gene.pair data.frame; pairs of genes to be tested. First column - gene A, second column - gene B
#'@param alt data.frame; gene alterations in tumor samples. The columns correspond to genes and the rows to tumor samples.
#'@param survival data.frame; survival status 1 - alive, 0 - dead
#'@param cores number of cores to be used; default 1
#'@importFrom dplyr rename_with filter arrange
#'@importFrom parallel mclapply
#'@importFrom data.table rbindlist
#'@return data frame
#'@export

survLRT <- function(gene.pair, alt, survival, cores = 1) {
    
    gene.pair <- gene.pair %>%
        dplyr::rename_with(~ c('geneA', 'geneB'), 1:2) %>%
        dplyr::filter(geneA %in% colnames(alt) & geneB %in% colnames(alt))
    data <- inputDataSurvLRT(alt, survival)
    gene.pairs <- base::split(gene.pair, row.names(gene.pair))
    results <- parallel::mclapply(gene.pairs, mc.cores = cores, function(x) {runFor2Genes(x, data)})
    output <- as.data.frame(data.table::rbindlist(results)) %>%
        dplyr::arrange(-SLflag, pvalue)
}
