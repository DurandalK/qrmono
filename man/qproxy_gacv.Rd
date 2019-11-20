\name{qproxy_gacv}
\alias{qproxy_gacv}
\title{
    Smoothing parameter selection using GACV 
}
\description{
    \code{qproxy.gacv} helps determine the smoothing parameter by resorting to GACV proposed in (Yuan, 2004).
}
\usage{
qproxy_gacv(X, Y, tau = 0.5, weights = NULL, lambda = seq(0.01, 2, length.out = 50), 
            constraint = c("none", "increase", "decrease"), method = c("Quad", "MM"),
            PeltyD = NULL, pord = 2, delta_w = 1e-6, ...)
}
\arguments{
    \item{X}{\emph{n} x \emph{p} numeric matrix, \emph{p} is spline basis matrix dimension.}
    \item{Y}{\emph{n} x 1 numeric vector, regression response.}
    \item{tau}{numeric value, current quantile level, should lie in (0, 1) and the default value is 0.5.}
    \item{lambda}{\emph{l} x 1 vector, possible penalty coefficients vector to be examined.}
    \item{weights}{\emph{n} x 1 numeric vector, pre-sepecfic data weights, default value can be void.}
    \item{constraint}{character string, could either be \code{"none"}, \code{"increase"} or \code{"decrease"}.}
    \item{method}{character string, Approximation scenario, could either be \code{"Quad"}-quadratic approximation or \code{"MM"}-MM algorithm.}
    \item{PeltyD}{\emph{k} x \emph{p} matrix, Penalty matrix, under penalized-spline method, k = p - pord, default value is a differential matrix.}
    \item{pord}{integer, penalty term degree, the default value is 2 for p-spline.}
    \item{delta_w}{numeric value, hat-matrix approximation tolerance.}
    \item{...}{Further arguments to be passed down to function \code{qproxy.gacv}, more info can be found in help documentation of \code{qproxy} or \code{qproxyz}.}
}
\value{
    \item{lambda.min}{numeric value, ideal smoothing parameter.}
    \item{gacv}{\emph{l} x 1 vector numeric, gacv value for each input lambda.}
}
\references{
    Yuan M. (2006). GACV for quantile smoothing splines. \emph{Computational Statistics & Data Analysis}, \bold{50(3)}: 813-829.
}
\seealso{
    \code{\link{qproxy}}
}
\examples{
library(splines)

x <- runif(50, 0, 1)
y <- -3 * sin((x - 0.5) * pi) + rnorm(50)

kts <- c(rep(min(x) - 0.1, 2), seq(min(x), max(x), 20), rep(x + 0.1, 2))
B <- spline.des(kts, x, 2)$design

lam.min <- qproxy_gacv(B, y, tau = 0.5, lambda = seq(0.1, 5, length.out = 50), constraint = "increase")$lambda.min

coef <- qproxy(B, y, tau = 0.5, lambda = lam.min, constraint = "decrease")
yfit <- drop(B \%*\% coef)

}
\keyword{GACV}
