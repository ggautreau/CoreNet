\encoding{utf-8}
\docType{methods}
\name{getPartition, corenet-method}
\alias{getPartition}
\title{getPartition}
\arguments{
  \item{X}{a corenet object}
  \item{nk}{an integer value corresponding to number of chosen classes. Default
  is 2.}
}

\value{
  A list describing a  partition: \code{nk} the number of clusters, \code{cluster} a numeric vector indicating to which class belongs each node, \code{Cik} a classification matrix with \code{nk} columns and \code{n} rows where \code{n} is the number of nodes, \code{mean} a matrix whose rows contain the centers of the clusters (\code{nk} rows with length \code{p} where \code{p} is the number of organisms).
}
\description{Method for accessing a partition within an object of class corenet.
}

\author{Christophe Ambroise \email{cambroise@genopole.cnrs.fr} and Catherine Matias.
}



\keyword{clustering}
\keyword{methods}