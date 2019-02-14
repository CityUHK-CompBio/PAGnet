\name{pagnet.mra}
\alias{pagnet.mra}
\title{Master regulator analysis based on PAGnet}
\usage{
pagnet.mra <- function(rnet=PAGnet, tflist=tf, signature=t3ss, pValueCutoff=0.05, pAdjustMethod=NULL)
}
\description{
  This function perform master regulator analysis (MRA) over a list of regulons from PAGnet or user uploaded regulator network.
}
\arguments{
  \item{rnet}{a regulatory network in Pseudomonas aeruginosa genomic netork.}
   \item{tflist}{a list of transcription factor in regulatory network.}
   \item{signature}{interested gene sets.}
   \item{pValueCutoff}{a single numeric value specifying the cutoff for p-values considered significant.}
   \item{pAdjustMethod}{a single character value specifying the p-calue adjustment method to be used.}

}

\examples{
  data(Pavirnet_basicdata)

  #Use PAGnet as regulatory network

  PAGnet <- network_inf[,c(1,5)]
  tf <- tf

  #Select T3SS related genes as signatures

  t3ss <- signature_example

  #run MRA
  pagnet.mra(rnet=PAGnet,tflist=tf,signature = t3ss, pValueCutoff = 0.05,pAdjustMethod="BH")

}
