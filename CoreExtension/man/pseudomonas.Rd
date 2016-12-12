\name{pseudomonas}
\alias{pseudomonas}
\docType{data}
\title{Pseudomonas data set.
}
\description{
The data set contains data related to 13 species from Pseudomonas genus. The list of metabolic reactions from each species has been filtered out (cofactors and current compounds were filtered out). 
}
\usage{data(pseudomonas)}
\format{
The dataset contains an object called \code{pseudomonas} of class \code{corenet} with attributes
\itemize{
\item{reactions} is a vector containing the name of each reaction	
\item{organisms} is a vector containing the name of each sif file used for the organism description	
\item{data} is a
 binary matrix with 2009 rows corresponding to the total number of (unique) metabolic reactions for this set of organisms and 13 columns corresponding to the total number of organisms. Each entry indicates whether the reaction is present (1) or absent (0) in the corresponding organism.

\item{neighbors}  a numeric vector that concatenates information on
each reaction. Reactions are numbered
according to their order in \code{data} attribute. For each reaction, we
concatenate its number, the number of neighbors it has and the numbers
of these neighbors. Example: if reaction number 3 has 4 neighbors that
are reactions 6, 7, 10 and 12, then one concatenates
(3,4,6,7,10,12). This is further concatenated for all reactions. 
}

Complete list of species is \emph{Pseudomonas aeruginosa LESB58, Pseudomonas aeruginosa PAO1, Pseudomonas aeruginosa UCBPP-PA14, Pseudomonas entomophila L48, 
Pseudomonas fluorescens Pf-5, 
Pseudomonas fluorescens Pf0-1, 
Pseudomonas fluorescens SBW25, 
Pseudomonas putida F1, Pseudomonas putida GB-1, Pseudomonas putida KT2440, Pseudomonas syringae pv. phaseolicola 1448A, 
Pseudomonas syringae pv. syringae B728a, 
Pseudomonas syringae pv. tomato DC3000}.
}

\source{
Original .sif files were uploaded from \url{http://www.cytoscape.org/}.
}

\examples{
data(pseudomonas)
show(pseudomonas)
str(pseudomonas)
}
\keyword{datasets}
