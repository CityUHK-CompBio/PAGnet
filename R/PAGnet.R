
###################################################################################
###################Master Regulator Analysis based on PAGnet#######################
#load("data/PAGnet.rda")

pagnet.mra <- function(rnet=PAGnet, tflist=tf, signature=qs, pValueCutoff=0.05, pAdjustMethod=NULL, minRegulonSize=5){
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
  mra_results[which(mra_results$Pvalue == 0),"Pvalue"] <- c("< 1e-4")

  return(mra_results)

}


