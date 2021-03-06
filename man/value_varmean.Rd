% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM.R
\name{value_varmean}
\alias{value_varmean}
\title{The Variance-Mean value function used in the EM algorithm}
\usage{
value_varmean(update_par = NULL, i = 1, para, Rr, Rf,
  M = NROW(Rr[[1]]), Tn = length(Rr), W, discount = 1/1.01,
  lambda = 1/2, returnW = FALSE, detail = FALSE)
}
\arguments{
\item{update_par}{Matrix. Each row contains the parameters at a time point.
The number of column are the sum of the number of parameter in beta function,
the number of parameters in the consumption function, and the number of weights of the assets minus 1
(the sum of the weights is one so only need to estimate N-1 weights).
The positions in the vector should also be in that order.}

\item{i}{Integer. The position of the current time point}

\item{Rr}{List of matrix. The rertuns of the assets.
Each list is for one time point and each column in the matrix is one assets.}

\item{Rf}{Scalar. The risk free rate. e.g. 1.01}

\item{M}{Integer. The number of instances for one assets at a time}

\item{Tn}{Integer. The length of time horizon. Not including 0.}

\item{W}{Matrix. Each column is one time point and each row is one observation.}

\item{discount}{Scalar. The rate of discount. e.g. 1/1.01}

\item{lambda}{The scalar that make the calculation easiler.}

\item{returnW}{If TRUE, return the wealth. If False, return the value from the value function.}
}
\description{
The Variance-Mean value function used in the EM algorithm
}
\author{
Yangzhuoran Yang
}
