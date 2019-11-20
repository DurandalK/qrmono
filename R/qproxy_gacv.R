qproxy_gacv <- function(X, Y, tau = 0.5, weights = NULL, lambda = seq(0.01, 2, length.out = 50), 
                        constraint = c("none", "increase", "decrease"), 
                        method = c("Quad", "MM"), PeltyD = NULL, pord = 2, 
                        delta_w = 1e-6, ...){
    
    # Input: 
    #       Y, (n x 1) vector,
    #           response;
    #       X, (n x p) matrix,
    #           predictors;
    #       tau, numeric,
    #           current quantile level, should lie in (0, 1);
    #       weights, (n x 1) vector,
    #           weights, can be void;
    #       lambda, (l x 1) vector,
    #           possible penalty coefficients;
    #       constraint, char string,
    #           could be "none", "increase" or "decrease"; 
    #       method, char string,
    #           Approximation scenario, could either be "Quad"-quadratic function or "MM"-MM algorithm;
    #       delta_w, numeric,
    #           hat-matrix approximation tolerance;
    #       PeltyD, (k x p) matrix,
    #           Penalty matrix, under penalized-spline method, k = p - pord;
    #       pord, integer,
    #           penalty term degree; 
    #       ..., other related parameters in qproxy.
    #
    # Output:
    #       lambda.min, numeric,
    #           ideal smoothing parameter;
    #       gacv, (l x 1) vector,
    #           gacv value for all lambdas;
    #
    # Refer:
    #       Yuan M. (2006). GACV for quantile smoothing splines.
    #       Dong C, Feng X D. (2019). Monotone nonparametric quantile regression and its bootstrap.
    #
    
    l <- length(lambda)
    n <- length(Y)
    p <- ncol(X)
    
    if(is.null(weights)){
        weights <- rep(1, n)
    }
    
    if(is.null(PeltyD)){
        PeltyD <- diff(diag(p), differences = pord)
    }
    gacv.val <- rep(0, l)
    
    for(i in 1 : l){
        coef <- qproxy(X, Y, tau = tau, weights = weights, lambda = lambda[i], constraint = constraint, method = method, ...)
        resid <- Y - drop(X %*% coef)
        W <- ifelse(abs(resid) > delta_w, (tau - (resid < 0)) / (2 * resid), (tau - (resid < 0)) * abs(resid) / delta_w)
        W <- diag(W)
        H <- X %*% solve(t(X) %*% W %*% X + lambda[i] * (t(PeltyD) %*% PeltyD)) %*% t(X) %*% W
        gacv.val[i] <- crossprod(checkloss(resid, tau), rep(1 / (n - sum(diag(H))), n)) / n
    }
    idx.min <- which.min(gacv.val)
    return(list(lambda.min = lambda[idx.min], gacv = gacv.val))
}