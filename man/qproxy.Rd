\name{qproxy}
\alias{qproxy}
\title{
    Monotone P-spline quantile regression curve fitting using proximal algorithm
}
\description{
    \code{qproxy} uses QICD in Peng and Wang (2015) to estimate high dimensional quantile regression model and apply different criterions to select the penalty coefficient.
}
\usage{
qproxy(X, Y, tau = 0.5, lambda, weights = NULL, beta_init = NULL, 
       constraint = c("none", "increase", "decrease"), PeltyD = NULL, 
       pord = 2, method = c("Quad", "MM"), delta = NULL, eps = NULL,
       iter_out = 1e6, iter_in = 1e2, threshold = 1e-6)
}
\arguments{
    \item{X}{\emph{n} x \emph{p} numeric matrix, \emph{p} is spline basis matrix dimension.}
    \item{Y}{\emph{n} x 1 numeric vector, regression response.}
    \item{tau}{numeric value, current quantile level, should lie in (0, 1) and the default value is 0.5.}
    \item{lambda}{numeric vector, possible penalty level, should be (>0).}
    \item{weights}{\emph{n} x 1 numeric vector, pre-sepecfic data weights, default value can be void.}
    \item{beta_init}{\emph{p} x 1 numeric vector, pre-estimated quantile regression coefficients vector, default value can be void.}
    \item{constraint}{character string, could either be \code{"none"}, \code{"increase"} or \code{"decrease"}.}
    \item{PeltyD}{\emph{k} x \emph{p} matrix, Penalty matrix, under penalized-spline method, k = p - pord, default value is a differential matrix.}
    \item{pord}{integer, penalty term degree, the default value is 2 for p-spline.}
    \item{method}{character string, Approximation scenario, could either be \code{"Quad"}-quadratic approximation or \code{"MM"}-MM algorithm.}
    \item{delta}{numeric, Quad method approximation tolerance.}
    \item{eps}{numeric, MM algrithm approximation tolerance.}
    \item{iter_out}{integer, maximum iteration number in outer loop, default value is 1e6.}
    \item{iter_in}{integer, maximum iteration numer in inner loop, default value is 1e2.}
    \item{threshold}{numeric, threshold value for iteration algorithm convergence, default value is 1e-6.}

}
\details{
    During regression curve fitting, one can sepecify \code{constraint} to be one of \code{"none""}, \code{"increase""} and \code{"decrease"}. \code{"none""} stands for not applying any constraint on quantile regression curve fitting; \code{"increase""} stands for fitting regression curve assuming it is a monotone increasing function and \code{"decrease"} means the opposite.
    To perform proximal algorithm one can choose to use either quadratic function approximation at zero point (\code{method="Quad""}) or MM (Majorize-minimize) alogrithm (\code{method="MM""}), more details can be found at Dong and Feng (2019).
    A inital staring point \code{beta_init} can be specificied, if not, a least square regression estimate will be treated as staring point.
}
\value{
    \item{beta_est}{\emph{p} x 1 numeric vector, constrained quantile regression coefficients estimates.}
}
\references{
    Eilers P H C, Marx B D. (1996). Flexible Smoothing with B-splines and Penalties. \emph{Statistical Science}, \bold{11(2)}: 89-121.
    
    Parikh N, Boyd S. (2013). Proximal Algorithms. \emph{Foundations and TrendsR in Optimization}, \bold{1(3)}: 123-231.
    
    Dong C, Feng X D. (2019). Monotone nonparametric quantile regression and its bootstrap. Working paper.
}
\seealso{
    \code{\link{qproxy_gacv}}
}
\examples{
library(splines)

x <- runif(50, 0, 1)
y <- -3 * sin((x - 0.5) * pi) + rnorm(50)

kts <- c(rep(min(x) - 0.1, 2), seq(min(x), max(x), 20), rep(x + 0.1, 2))
B <- spline.des(kts, x, 2)$design

coef <- qproxy(B, y, tau = 0.5, lambda = 2, constraint = "decrease")
yfit <- drop(B \%*\% coef)
}
\keyword{Constrained quantile regression curve}
\keyword{Proximal algorithm}
