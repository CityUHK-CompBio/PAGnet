# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

###################################################################################
###################Master Regulator Analysis based on PAGnet#######################
#load("data/Pavirnet_basicdata.Rda")




pagnet.mra <- function(rnet=PAGnet, tflist=tf, signature=t3ss, pValueCutoff=0.05, pAdjustMethod=NULL, minRegulonSize=5){
  universe <- unique(c(as.matrix(rnet[,1]),as.matrix(rnet[,2])))
  universe.size <- length(unique(c(as.matrix(rnet[,1]),as.matrix(rnet[,2]))))

  signature.size <- length(as.character(as.matrix(signature)))
  signature.innet <- intersect(as.character(as.matrix(signature)),universe)
  signature.innet.size <- length(signature.innet)

  mra_results<-c()

  for(i in tflist){
    itf <- i
    regulon <- rnet[which(rnet[,1] == as.character(as.matrix(itf))),]
    regulon.target <- as.character(as.matrix(regulon[,2]))
    regulon.size <- length(regulon.target)

    observed.signature <- intersect(regulon.target,signature.innet)
    observed.signature.size <- length(observed.signature)
    pval <- round(phyper(observed.signature.size - 1, m = signature.innet.size, n = universe.size - signature.innet.size, k = regulon.size, lower.tail=F ),digits = 4)
    mra_results <- rbind(mra_results,data.frame(itf,universe.size,regulon.size,signature.innet.size,observed.signature.size,pval))
  }
  colnames(mra_results) <- c("TF","network.size", "regulon.size", "signature.size","observed.signature.size", "Pvalue")
  #mra_results <- mra_results[which(mra_results[,5]  != 0),]
  if(!is.null(pAdjustMethod)){
    pvals.adj <- round(p.adjust(mra_results$Pvalue, method=pAdjustMethod),4)
    mra_results <- cbind(mra_results,pvals.adj)
    colnames(mra_results) <- c("TF","network.size", "regulon.size", "signature.size","observed.signature.size", "Pvalue","adjust.Pvalue")
    mra_results[which(mra_results[,7] == 0),7] <- c("< 1e-4")
  }
  mra_results <- data.frame(mra_results[order(mra_results[,6]),])

  if(!is.null(pValueCutoff)){
    mra_results <- mra_results[which(mra_results$Pvalue < pValueCutoff),]
  }
  mra_results[which(mra_results[,6] == 0),6] <- c("< 1e-4")

  return(mra_results)

}
