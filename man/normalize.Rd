\name{normalize}
\alias{normalize}
\title{Min-max normalization}
\description{Min-max normalization for data matrices and data frames
}

\usage{
normalize(x,type="all")
}

\arguments{
\item{x}{a data frame or data matrix}
\item{type}{type of normalizaion. "row" normalization row by row, "col" normalization column by colum, "all" normalization for the entire data frame/matrix (default)}
}

\examples{
  mtx<-matrix(rnorm(20),5,4)
  n_mtx<-normalize(mtx) # Fully normalized matrix
  r_mtx<-normalize(mtx,type="row") # Normalize row by row
  c_mtx<-normalize(mtx,type="col") # Normalize col by col
  print(n_mtx) # Print fully normalized matrix
}

\value{
Returns a normalized data.frame/matrix.
}

\author{
Zsolt T. Kosztyan, University of Pannonia

e-mail: kosztyan.zsolt@gtk.uni-pannon.hu
}

\keyword{matrix}
