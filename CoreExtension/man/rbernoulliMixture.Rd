\name{rbernoulliMixture}
\alias{rbernoulliMixture}
\title{
Simulation of a binary vector according to an r-states Bernoulli mixture model.
}
\description{
The function \code{rbernoulliMixture} simulates a binary vector according to an r-states Bernoulli mixture model.
}

\usage{
rbernoulliMixture(pi, theta, n)
}
\arguments{
  \item{pi}{a vector containing r groups proportions. Must sum
    to one.
  }
  \item{theta}{
a matrix of Bernoulli probabilities within each group (rows) and for
each coordinate (column). The number of rows must be equal to length of
pi and number of columns is length of vector to be simulated.  
}
  \item{n}{
an integer value, the size of data to be simulated.
}
}

\value{
  A list of four objects.
\itemize{
  \item{x} is generated dataset. A matrix with n rows and same number of
  columns as theta. 
\item{z} is a vector of length n, indicating the group of each observation.
\item{pi} is a  vector of initial group proportions.
\item{theta} is a matrix of initial parameters.
}
}

\author{
Christophe Ambroise.
}
\keyword{Bernoulli}
\keyword{Mixture}
