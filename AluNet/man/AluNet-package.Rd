\name{AluNet-package}
\alias{AluNet-package}
\alias{AluNet}
\docType{package}
\title{
  Identifies candidate Alu elements that modulate 3D chromatin structure 
  by combining Hi-C, ATAC-seq data with a database of Alu elements. It achieves this
  with a custom weighteed adjacnecy matrix calculation and scratch-implemented 
  Leiden clustering algorithm.
}
\description{
  Weights are determined by DNA contact probabilities, openness of chromatin
  to enforce that Alu elements are able to be transcribed, and the maximum
  number of expected Alu pairs in the same family between two DNA fragments.
  This assumes Alus in the same family are complementary and transcribable (i.e. 
  have not accumulated mutations that prevent transcription or non-complementarity).
}
\details{
  The most important functions are runLeiden which run the weighted adjacnecy matrix through the
  Leiden algorithm as defined by Traag et al. (2019). There are also functions for
  cleaning standard Fit-Hi-C data, ATAC-seq peak call tables, and Dfam database files.
}
\author{
Brandt Bessell, Xuyuan Zhang, and Sicheng Qian

Maintainer: Your Name <bbessell@umich.edu>
}
\references{
  Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. Scientific reports, 9(1), 1-12.
}
\keyword{ package }
