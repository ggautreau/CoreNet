\encoding{utf-8}
\name{coreextension-package}
\alias{coreextension-package}
\alias{CoreExtension}
\docType{package}
\title{
Extended core characterization of a set of objects (metabolic reactions or genes) within a set of species.}
\description{
\code{CoreExtension} performs probabilistic model based clusterings of objects (= metabolic reactions or genes) within a set of organisms. Two underlying model may be used: 1) a simple Bernoulli mixture model using only information about presence/absence of an object within a species ; 2) a hidden Markov random field model that combines two types of information: the former presence/absence information plus a dependency structure on the set of objects.
The dependency graph between reactions is the metabolic graph where two reactions are neighbors if they share a metabolite. The dependency graph between genes links any two genes whose positions are adjacent in the organism's genome.
}

\details{
\tabular{ll}{
Package: \tab CoreExtension\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-06-28\cr
License: \tab GPL-2 \cr
}
}
\author{
Christophe Ambroise <cambroise@genopole.cnrs.fr>, Yolande Diaz <ydiaz@genopole.cnrs.fr> and Catherine Matias <catherine.matias@math.cnrs.fr>, 
relying on initial code by Mô Van Dang.

Maintainer: Christophe Ambroise <cambroise@genopole.cnrs.fr>
}
\references{
C. Ambroise,  Y. Diaz, M.  Fernandez, C. Gautier,  C.C. Klein,
  V. Lacroix, C. Matias, M.-F. Sagot, (2013). A MODIFIER !!!
 Manuscript.
}
\keyword{neighborhood EM, Bernoulli mixture}

\examples{
data(pseudomonas)

# pseudomonas is corenet object with data corresponding to a binary matrix 
# indicating presence/absence of a metabolic reaction in each organism
str(pseudomonas)

# Clustering of reactions in two groups without neighborhood information (beta=0)
res_without_nei <- nem(pseudomonas,bmod="fix",beta=0)
str(res_without_nei)
cluster_without_nei <- res_without_nei@partitionList[[1]]$cluster


# Running NEM clustering procedure with strong beta value (beta=1)
res_with_nei_bstrong <- nem(pseudomonas,beta=1)
cluster_with_nei_bstrong <- res_with_nei_bstrong@partitionList[[1]]$cluster

# Comparing clusters
##  confusion table     
table(cluster_without_nei,cluster_with_nei_bstrong)
## rand index 
randError(cluster_without_nei,cluster_with_nei_bstrong)
# Estimated parameters
res_with_nei_bstrong@partitionList[[1]]$mean

# Running NEM procedure with optimised beta value   
res_with_nei_boptim <- nem(pseudomonas,bmod="optim")
#Optimised value for beta
res_with_nei_boptim@partitionList[[1]]$beta

# Running NEM for finding the best partition using BICp
result<-nem(pseudomonas,nkMin=1,nkMax=5)
plot(result)
plot(result,type="mean")   

}
