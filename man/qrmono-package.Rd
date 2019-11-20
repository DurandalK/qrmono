\name{qrmono-package}
\alias{qrmono-package}
\alias{qrmono}
\docType{package}
\title{
    Nonparametric Constrained Quantile Regression Curve Fitting
}
\description{
    This package provides a fast way to fit constrained quantile regression curve using P-spline technique and proximal algorithm. More details about this fitting method can be found in Dong and Feng (2019).
}
\details{
    Package \code{qrmono} contains the following functions for implementation.
    \describe{
    \item{qproxy}{Proximal algorithm for constraint quantile regression curve fitting using P-spline.}
    \item{qproxy_gacv}{Smoothing parameter determination using GACV criterion.}
    \item{proxyz}{A off-the-shelf version of \code{qproxy} without initialization needed.}
    }
}
\author{
Chen Dong, Xingdong Feng

Maintainer: Chen Dong <dongchen39@hotmail.com>
}
\references{
    Eilers P H C, Marx B D. (1996). Flexible Smoothing with B-splines and Penalties. \emph{Statistical Science}, \bold{11(2)}: 89-121.
    
    Yuan M. (2006). GACV for quantile smoothing splines. \emph{Computational Statistics & Data Analysis}, \bold{50(3)}: 813-829. 
    
    Parikh N, Boyd S. (2013). Proximal Algorithms. \emph{Foundations and TrendsR in Optimization}, \bold{1(3)}: 123-231.
    
    Dong C, Feng X D. (2019). Monotone nonparametric quantile regression and its bootstrap. Working paper.
}
\keyword{qrmono}
\examples{
  library(qrmono)
}
