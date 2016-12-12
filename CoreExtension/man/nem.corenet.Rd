\encoding{utf-8}
\docType{methods}
\name{nem, corenet-method}
\alias{nem}
\title{Spatial Classification by Neighborhood EM algorithm of corenet object}
\arguments{
  \item{X}{a corenet object}
  \item{nk}{an integer value corresponding to the number of classes. Default
  is 2.}
\item{nkMin}{an integer value corresponding to minimal number of classes. Default value is 2.}
\item{nkMax}{an integer value corresponding to maximal number of classes. Default value is 2.}
\item{family}{a character string indicating which method is to be used. Possible values are "bernoulli" (default parameter), "normal".}
\item{algo}{a character string indicating which method is to be used. Possible values are "ncem" (neighborhood classifying em), "nem"  (neighborhood em) and "gem" (generalized em). Default value is "nem".}
  \item{beta}{a real positive number, that is used as scale parameter for importance of neighborhood  information. Default value is 1. Using beta=0 corresponds to relying on a multivariate Bernoulli mixture model. }
  \item{bmod}{a character string indicating the type of choice for scale 
  parameter beta. Possible values are "fix" (beta is fixed by the user)
  or "optim" (beta is optimised relying on an approximation of the likelihood).
 Default value is "fix".}
\item{iters}{an integer value, maximum number of iterations
  allowed. Should be at least 100. Default value is 2000.}
\item{nbInit}{an integer value, number of different starting points. Default value is 5.}
  \item{\ldots}{\ldots}

}

\value{
  An object of class \code{\link{corenet}} with modified \code{partitionList} and \code{partitionOptions} attributes. The \code{partitionList} is a  list of partitions. Each partition is described by the following items: \code{nk} the number of clusters, \code{cluster} a numeric vector indicating to which class belongs each node, \code{Cik} a classification matrix with \code{nk} columns and \code{n} rows where \code{n} is the number of nodes, \code{mean} a matrix whose rows contain the centers of the clusters (\code{nk} rows with length \code{p} where \code{p} is the number of organisms).
\code{partitionOptions} is the  list of options used for obtaining the \code{partitionList}.  
}
\description{
 Method for Function \code{nem} performs  clustering of the nodes (=reactions or genes) of a graph (metabolic graph or gene adjacencies graph) common to a set of organisms, according to a Markov random field model. The dependency graph between reactions is the metabolic graph where two reactions are neighbors if they share a metabolite. The dependency graph between genes links  any two genes whose positions are adjacent in the organism's genome.
 Each node (reaction or gene) is moreover described by a binary vector of presence/absence in each considered organism.
}

\references{
C. Ambroise, \emph{Approche probabiliste en classification automatique et contraintes de voisinage}, Ph.D. thesis, Université de Technologie de Compiègne, 1996.

C. Ambroise, M. Dang, and G. Govaert, \emph{Clustering of
    spatial data by the EM algorithm} in Geostatistics for Environmental
  Applications, A. Soares, J. Gomez-Hernandez, and R. Froidevaux, Eds.,
  pp. 493-504. Kluwer Academic Publisher, 1997.

Mô Van Dang, \emph{Classification de données spatiales : modèles probabilistes et critères de partitionnement},  Ph.D. thesis, Université de Technologie de Compiègne, 1998. \url{http://www.sudoc.fr/126401616}
}


\author{Christophe Ambroise, \email{cambroise@genopole.cnrs.fr} and Catherine Matias. Relying on preliminary code written by Philippe Huppé for \code{R} package \code{Manor}. Interfaces a \code{C} source code (\code{nem}) entirely written by Mô Van Dang. 
}


\keyword{clustering}
\keyword{methods}