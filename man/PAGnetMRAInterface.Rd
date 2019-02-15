\name{pagnet.mra.interface}
\alias{pagnet.mra.interface}
\title{PAGnet: an integrated regulatory network in Pseudomonas aeruginosa}
\usage{
  pagnet.mra.interface()
}
\description{
 master regulator analysis for identification of key transcription factors mediating a biological process or pathway. First, the user can choose to use the default PAVIRnet or to upload their own regulatory network in a predefined format. Second, the user needs to specify a gene signature associated with a biological function or pathway of interest, either by selecting a gene set from public databases or uploading a user-customized gene list. In the current version, the platform provides gene sets in Gene Ontology (GO) and KEGG databases obtained from Pseudomonas Genome DB. Having completed master regulator analysis, a table will be returned with information about each transcription factor’s corresponding gene ID, gene name, number of target genes, total number of hits (all signature genes in the network), observed hits (signature genes in the TF’s regulon), and a p-value calculated based on a hypergeometric test. The table is sorted according to the statistical significance indicated by the p-values, and the top significant TFs can be prioritized as master regulators.
}

\examples{
  library(PAGnet)

  #To start interface, call pagnet.mra.interface
  pagnet.mra.interface()
}
