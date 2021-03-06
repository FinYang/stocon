% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weights.R
\name{weights_lm}
\alias{weights_lm}
\title{Finding weights that minimise variance of return using linear regression}
\usage{
weights_lm(Rt)
}
\arguments{
\item{Rt}{Lists of return. e.g. output from \code{\link{sim_simple}}}
}
\value{
List of weights for each period
}
\description{
Finding weights that minimise variance of return using linear regression as described in
Fan, Zhang & Yu (2012)
}
\references{
Fan, J., Zhang, J., & Yu, K. (2012). Vast Portfolio Selection With Gross-Exposure Constraints. Journal of the American Statistical Association, 107(498), 592–606. https://doi.org/10.1080/01621459.2012.682825
}
\seealso{
\code{\link{weights_lasso}}
}
\author{
Yangzhuoran Yang
}
