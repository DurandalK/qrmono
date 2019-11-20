qproxy <- function(X, Y, tau = 0.5, lambda, weights = NULL, beta_init = NULL,
                   constraint = c("none", "increase", "decrease"),
                   PeltyD = NULL, pord = 2, method = c("Quad", "MM"), 
                   delta = NULL, eps = NULL, iter_out = 1e6, 
                   iter_in = 1e2, threshold = 1e-6){
    
    # Input: 
    #       Y, (n x 1) vector,
    #           response;
    #       X, (n x p) matrix,
    #           predictors;
    #       tau, numeric,
    #           current quantile level, should lie in (0, 1);
    #       lambda, numeric,
    #           penalty coefficient;
    #       weights, (n x 1) vector,
    #           weights, can be void;
    #       beta_init, (p x 1) vector,
    #           pre-estimated coefficients vector, can be void;
    #       constraint, char string,
    #           could be "none", "increase" or "decrease"; 
    #       PeltyD, (k x p) matrix,
    #           Penalty matrix, under penalized-spline method, k = p - pord;
    #       pord, integer,
    #           penalty term degree; 
    #       method, char string,
    #           Approximation scenario, could either be "Quad"-quadratic approximation or "MM"-MM algorithm;
    #       delta, numeric,
    #           Quad method approximation tolerance;
    #       eps, numeric,
    #           MM algrithm approximation tolerance;
    #       iter_out, integer,
    #           maximum iteration number in outer loop;
    #       iter_in, integer,
    #           maximum iteration numer in inner loop;
    #       threshold, numeric,
    #           threshold value for iteration algorithm convergence.
    #
    # Output:
    #       beta_est, (p x 1) vector,
    #           constrained quantile regression coefficients estimates.
    #
    # Refer:
    #       Dong C, Feng X D. (2019). Monotone nonparametric quantile regression 
    #           and its bootstrap. Working paper.
    #
    
    n <- nrow(X)
    p <- ncol(X)
    
    if(is.null(weights)){
        weights <- rep(1, n)
    }
    
    if(is.null(beta_init)){
        beta_init <- numeric(0)
    }
    
    if(constraint == "increase"){
        model <- 1
    }else if(constraint == "decrease"){
        model <- 2
    }else{
        model <- 0
    }
    
    if(length(method) > 1){
        method <- "Quad"    
    }
    
    if(is.null(eps)){
        eps <- threshold / n / 20
    }
    
    if(is.null(delta)){
        delta <- 1e-6
    }
    
    if(is.null(PeltyD)){
        PeltyD <- diff(diag(p), differences = pord)
    }
    
    if(method == "MM"){
        beta_est <- proxy_xi(Y, X, tau, beta_init, weights, model, eps, lambda, PeltyD, iter_out, iter_in, threshold)
    }else{
        beta_est <- proxy_delta(Y, X, tau, beta_init, weights, model, delta, lambda, PeltyD, iter_out, iter_in, threshold)
    }
    
    return(beta_est)
}