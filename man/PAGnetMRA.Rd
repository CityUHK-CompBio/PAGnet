\name{pagnet.mra}
\alias{pagnet.mra}
\title{PAGnet: Pseudomonas aeruginosa genomic integrated regulatory network}
\usage{
  pagnet.mra(rnet=PAGnet, tflist=tf, signature=qs, pValueCutoff=0.05, pAdjustMethod=NULL)
}
\description{
  This function perform master regulator analysis (MRA) over a list of regulons from the default PAGnet or user uploaded regulator network for identification of key transcription factors mediating a biological process or pathway.
}
\arguments{
  \item{rnet}{a regulatory network in Pseudomonas aeruginosa genomic netork.The user can choose to use the default PAVIRnet or to upload their own regulatory network in a predefined format}
   \item{tflist}{a list of transcription factor in regulatory network.}
   \item{signature}{interested gene sets.}
   \item{pValueCutoff}{a single numeric value specifying the cutoff for p-values considered significant.}
   \item{pAdjustMethod}{a single character value specifying the p-calue adjustment method to be used (see 'p.adjust' fro details). }

}

\examples{
  data(PAGnet)

  #Use PAGnet as regulatory network
  #Select quorum sensing (QS) related genes as signatures
  #run MRA

  MRA_results <- pagnet.mra(rnet=PAGnet,tflist=tf,signature = qs, pValueCutoff = 0.05, pAdjustMethod="BH")

  #Having completed master regulator analysis, a table will be returned
  head(MRA_results)
}
