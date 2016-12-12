\encoding{utf-8}
\name{corenet-class}
\docType{class}
\alias{corenet}
\alias{initialise,corenet-method}

\title{Class \code{"corenet"}}
\description{
This class of objects encodes a set of metabolic reactions or genes in different organisms.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("corenet", dir_files, union_graph)}. 
}
\section{Slots}{
  \describe{
    \item{\code{dir_files}:}{Object of class \code{"character"}. A directory containing a set of .sif files (Cytoscape format). Each .sif file should correspond to an organism and contain the list of metabolic reactions or genes WE SHOULD MODIFY THAT within this organism. Default value is ".".
 }
    \item{\code{union_graph}:}{Object of class \code{"character"}. A .sif file (Cytoscape format) that corresponds to the metabolic union graph of a set of organisms. If not given, the function constructs it by listing all the
reactions appearing in the set of organisms. }
    \item{\code{nodes}:}{Object of class \code{"character"}. 
    	A vector that contains the names  of 
each node (reaction or gene).}
    \item{\code{organisms}:}{Object of class \code{"character"}. 
    	A vector that contains the names  of 
each sif file describing an organism.}
    \item{\code{data}:}{Object of class \code{"matrix"}.  A binary matrix with number of rows equal to the total number of (unique) nodes and number of columns equal to the number of organisms (number of .sif files). Each entry indicates whether the reaction is present (1) or absent (0) in the organism.}
    \item{\code{neighbors}:}{Object of class \code{"numeric"}. 
    	A vector that concatenates information on
each node from entries \code{nodes}, numbered
as they appear in vector \code{nodes}. For each node, one 
concatenates its number, the number of neighbors it has and the numbers
of these neighbors. Example: if node number 3 has 4 neighbors that
are reactions 6, 7, 10 and 12, then one concatenates
(3,4,6,7,10,12). This is further concatenated for all nodes. }
    \item{\code{partitionList}:}{Object of class \code{"list"}. A list of partitions. Each partition is described by the following items: \code{nk} the number of clusters, \code{cluster} a numeric vector indicating to which class belongs each node, \code{Cik} a classification matrix with \code{nk} columns and \code{n} rows where \code{n} is the number of nodes, \code{mean} a matrix whose rows contains the centers of the clusters (\code{nk} rows with length \code{p} where \code{p} is the number of organisms).}
    \item{\code{partitionOptions}:}{Object of class \code{"list"}. Options of the \code{nem Method}
used for obtaining the \code{partitionList}. }
 }
}



\section{Methods}{This class comes with
  \describe{
    \item{nem}{a Markov Random Field model based clustering method (see \code{\link{nem.corenet}})}
    \item{plot}{a specific plotting method documented (see \code{\link{plot.corenet}})}
    \item{getPartition}{a method for accessing a given partition of an object of class \code{corenet}
(see \code{\link{getPartition.corenet}})}}
	 }
\author{
Christophe Ambroise and Catherine Matias

}
\note{
More details about Cytoscape may be found at their web site \url{http://www.cytoscape.org/}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("corenet")
## If there is a union_graph USE it
# myNet<-new("corenet",union_graph="union.sif",dir_files=".")

## If there is no union_graph it will be created with all sif files in dir_files
# myNet<-new("corenet",dir_files=".")

}
\seealso{
  See also \code{\link{plot.corenet}}, \code{\link{nem.corenet}}.
}

\keyword{classes}
