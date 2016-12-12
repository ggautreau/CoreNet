\name{randError}
\alias{randError}
\title{Rand Index between two partitions.
}
\description{
Function computing the adjusted rand Index between two partitions. 
}
\usage{
randError(x, y)
}
\arguments{
  \item{x}{First partition, vector of integers describing the class of each observations.
}
  \item{y}{Second partition vector describing the same observations as \code{x}.
}

}
\value{
Return a number less than 1 (1 meaning total agreement between partitions).
}
\references{
Lawrence Hubert and Phipps Arabie, \emph{Comparing partitions}, Journal of Classification, 2(1): 193-218, 1985.
}
\author{
Christophe Ambroise
}
\keyword{clustering}
