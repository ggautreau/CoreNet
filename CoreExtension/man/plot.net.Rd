\encoding{utf-8}
\docType{methods}
\name{plot,corenet-method}
\alias{plot.corenet}
\title{Plot method for a corenet object}

\arguments{
  \item{x}{output of a nem method of the
  \pkg{CoreExtension} package (\code{\link{nem.corenet}}).
  Must be of  class \code{\link{corenet}}.}

  \item{y}{used for S4 compatibility.}

  \item{type}{character vector to indicate the type of plot: \code{mean} (image of the resulting partition) or \code{BICp} (penalized pseudo-likelihood with BIC penalty values).
 Default is \code{BICp}.}
}
\description{
  Produces a plot of the solutions of a \code{corenet} object clustered via \code{nem}.
}
\examples{
data(pseudomonas)
result<-nem(pseudomonas,nkMin=1,nkMax=5)
plot(result)
plot(result,type="mean")   
}

\references{
D. Stanford, \emph{Fast Automatic Unsupervised Image Segmentation and Curve Detection in Spatial Point Processes}, PhD thesis, Dept. of Statistics, Univ. of Washington, Seattle, 1999.

Mô Van Dang, \emph{Classification de données spatiales : modèles probabilistes et critères de partitionnement},  Ph.D. thesis, Université de Technologie de Compiègne, 1998. \url{http://www.sudoc.fr/126401616}
}


