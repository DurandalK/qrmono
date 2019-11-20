# qrmono
This package provides a fast way to fit constrained quantile regression curve using P-spline technique and proximal algorithm. More details about this fitting method can be found in Dong and Feng (2019).
# Installation
```
#install.packages("devtools")
library(devtools)
install_github("DurandalK/qrmono")
```  
# Example
```
library(qrmono)
library(splines)

x <- runif(50, 0, 1)
y <- -3 * sin((x - 0.5) * pi) + rnorm(50)

kts <- c(rep(min(x) - 0.1, 2), seq(min(x), max(x), 20), rep(x + 0.1, 2))
B <- spline.des(kts, x, 2)$design

lam.min <- qproxy_gacv(B, y, tau = 0.5, lambda = seq(0.1, 5, length.out = 50), constraint = "decrease")$lambda.min

coef <- qproxy(B, y, tau = 0.5, lambda = lam.min, constraint = "decrease")
yfit <- drop(B %*% coef)
```
# References
Eilers P H C, Marx B D. (1996). Flexible Smoothing with B-splines and Penalties. Statistical Science, 11(2): 89-121.
    
Yuan M. (2006). GACV for quantile smoothing splines. Computational Statistics & Data Analysis, 50(3): 813-829. 
    
Parikh N, Boyd S. (2013). Proximal Algorithms. Foundations and TrendsR in Optimization, 1(3): 123-231.
    
Dong C, Feng X D. (2019). Monotone nonparametric quantile regression and its bootstrap. Working paper.

# Development
This R package is maintained by Chen Dong (dongchen39@hotmail.com).
