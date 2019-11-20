#include <Rcpp.h>
#include "nr3.h"
#include "ludcmp.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector proxy_delta(NumericVector Y, NumericMatrix X, double tau, NumericVector beta_in, NumericVector weights, int model, double delta, double lambda, NumericMatrix PeltyD, int maxint1, int maxint2, double threshold){
    int p = X.ncol(), n = X.nrow(), pd = PeltyD.nrow(), i, j, k, count1, count2, tl = (model != 0) ? (p + 1) : p;
    double step_size = 0.5, l, f_y, f_z, sum1, sum2, sum3;
    NumericVector beta_out(p), resid_y(n), resid_z(n), resid_new(n), rho_y(n), rho_z(n), rho_new(n), theta_init(tl);
        
    MatDoub C(p, tl), XC(n, tl), S(tl, tl), DC(pd, tl), SI(tl, tl);
    VecDoub theta_old(tl), theta_new(tl), theta_y(tl), theta_z(tl), mu(tl), f_prime(tl), vecs(n);
    
    NumericVector coef_initial(NumericVector Y, NumericMatrix X, double tau, NumericVector beta_in, int model, double lambda, NumericMatrix PeltyD);
    NumericVector checkloss(NumericVector res, double tau);
    
    C.assign(p, tl, 0.0);
    
    if(model == 1){
        // Set monotone constraint matrix C
        for(i = 0; i < p; i++){
            C[i][0] = 1.0;
            for(j = 1; j < i + 1; j++){
                C[i][j] = -1.0;
            }
        }
    }else if(model == 2){
        for(i = 0; i < p; i++){
            C[i][0] = -1.0;
            for(j = 1; j < i + 1; j++){
                C[i][j] = 1.0;
            }
        }
    }else{
        for(i = 0; i < p; i++){
            C[i][i] = 1.0;
        }
    }
    
    // Store XC to save computation time
    for(i = 0; i < n; i++){
        for(j = 0; j < tl; j++){
            XC[i][j] = 0.0;
            for(k = 0; k < p; k++){
                XC[i][j] += X(i, k) * C[k][j];
            }
        }
    }
    // Store S to save computation time
    for(i = 0; i < pd; i++){
        for(j = 0; j < tl; j++){
            DC[i][j] = 0.0;
            for(k = 0; k < p; k++){
                DC[i][j] += PeltyD(i, k) * C[k][j];
            }
        }
    }
    for(i = 0; i < tl; i++){
        for(j = i; j < tl; j++){
            S[i][j] = 0.0;
            for(k = 0; k < pd; k++){
                S[i][j] += DC[k][i] * DC[k][j];
            }
            S[j][i] = S[i][j];
        }
    }
    
    // coefs initialization
    theta_init = coef_initial(Y, X, tau, beta_in, model, lambda, PeltyD);
    for(k = 0; k < tl; k++){
        theta_old[k] = 0.0;
        theta_new[k] = theta_init[k];
    }
    
    for(l = 1.0, count1 = 0; count1 <= maxint1; count1++){
        // intermediate coefs
        for(k = 0; k < tl; k++){
            theta_y[k] = theta_new[k] + count1 / (count1 + 3.0) * (theta_new[k] - theta_old[k]);
        }
        // obtain current residuals
        for(i = 0; i < n; i++){
            resid_y[i] = Y[i];
            for(k = 0; k < tl; k++)
                resid_y[i] -= XC[i][k] * theta_y[k];  
        }
        // obtain current function value and derivative
        rho_y = checkloss(resid_y, tau);
        f_y = 0.0;
        for(i = 0; i < n; i++){
            f_y += rho_y[i] * weights[i];
            if(abs(resid_y[i]) > delta){
                vecs[i] = tau - ((resid_y[i] < 0.0) ? 1.0 : 0.0);
            }else{
                vecs[i] = 2.0 * abs(resid_y[i]) / delta * (tau - ((resid_y[i] < 0.0) ? 1.0 : 0.0));
            }
        }
        f_y /= n;
        for(k = 0; k < tl; k++){
            f_prime[k] = 0.0;
            for(i = 0; i < n; i++){
                f_prime[k] -= XC[i][k] * weights[i] * vecs[i];
            }
            f_prime[k] /= n;
        }
        
        // backtracking line search
        for(count2 = 0; count2 <= maxint2; count2++, l *= step_size){
            for(k = 0; k < tl; k++){
                mu[k] = theta_y[k] - l * f_prime[k];
            }
            if(model != 0){
                for(k = 0; k < tl; k++){
                    if(mu[k] > 0.0)
                        mu[k] = 0.0;
                }
            }
            for(i = 0; i < tl; i++){
                for(j = i; j < tl; j++){
                    SI[i][j] = S[i][j] * 2.0 * l * lambda;
                    SI[j][i] = SI[i][j];
                }
                SI[i][i] += 1.0;
            }
            
            LUdcmp M(SI);
            M.solve(mu, theta_z);
            
            for(i = 0; i < n; i++){
                resid_z[i] = Y[i];
                for(k = 0; k < tl; k++){
                    resid_z[i] -= XC[i][k] * theta_z[k];
                }
            }
            rho_z = checkloss(resid_z, tau);
            f_z = 0.0;
            for(i = 0; i < n; i++){
                f_z += rho_z[i] * weights[i];
            }
            f_z /= n;
            sum1 = 0.0;
            sum2 = 0.0;
            for(k = 0; k < tl; k++){
                sum1 += f_prime[k] * (theta_z[k] - theta_y[k]);
                sum2 += pow(theta_z[k] - theta_y[k], 2.0);
            }
            if(f_z <= f_y + sum1 + sum2 / (2.0 * l))
                break;
        }
        
        // condition to jump off
        for(i = 0; i < n; i++){
            resid_new[i] = Y[i];
            for(k = 0; k < tl; k++)
                resid_new[i] -= XC[i][k] * theta_new[k];  
        }
        rho_new = checkloss(resid_new, tau);
        sum3 = 0.0;
        for(i = 0; i < n; i++){
            sum3 += (rho_z[i] - rho_new[i]) * weights[i];
        }
        if(abs(sum3) / n < threshold)
            break;
        for(k = 0; k < tl; k++){
            theta_old[k] = theta_new[k];
            theta_new[k] = theta_z[k];
        }
    }
    
    for(j = 0; j < p; j++){
        beta_out[j] = 0.0;
        for(k = 0; k < tl; k++){
            beta_out[j] += C[j][k] * theta_z[k];
        }
    }
    return beta_out;
}

// [[Rcpp::export]]
NumericVector proxy_xi(NumericVector Y, NumericMatrix X, double tau, NumericVector beta_in, NumericVector weights, int model, double eps, double lambda, NumericMatrix PeltyD, int maxint1, int maxint2, double threshold){
    int p = X.ncol(), n = X.nrow(), pd = PeltyD.nrow(), i, j, k, count1, count2, tl = (model != 0) ? (p + 1) : p;
    double step_size = 0.5, l, f_y, f_z, sum1, sum2, sum3;
    NumericVector beta_out(p), resid_y(n), resid_z(n), resid_new(n), rho_y(n), rho_z(n), rho_new(n), vecs1(n), vecs2(n), theta_init(tl);
    
    VecDoub theta_old(tl), theta_new(tl), theta_y(tl), theta_z(tl), mu(tl), v_eps(n), Qxi_prime(tl);
    MatDoub C(p, tl), XC(n, tl), S(tl, tl), DC(pd, tl), SI(tl, tl);
    
    NumericVector coef_initial(NumericVector Y, NumericMatrix X, double tau, NumericVector beta_in, int model, double lambda, NumericMatrix PeltyD);
    NumericVector Q_xi(NumericVector r, NumericVector rk, double tau, double eps);
    
    C.assign(p, tl, 0.0);
    
    if(model == 1){
        // Set monotone constraint matrix C
        for(i = 0; i < p; i++){
            C[i][0] = 1.0;
            for(j = 1; j < i + 1; j++){
                C[i][j] = -1.0;
            }
        }
    }else if(model == 2){
        for(i = 0; i < p; i++){
            C[i][0] = -1.0;
            for(j = 1; j < i + 1; j++){
                C[i][j] = 1.0;
            }
        }
    }else{
        for(i = 0; i < p; i++){
            C[i][i] = 1.0;
        }
    }
    
    // Store XC to save computation time
    for(i = 0; i < n; i++){
        for(j = 0; j < tl; j++){
            XC[i][j] = 0.0;
            for(k = 0; k < p; k++){
                XC[i][j] += X(i, k) * C[k][j];
            }
        }
    }
    // Store S to save computation time
    for(i = 0; i < pd; i++){
        for(j = 0; j < tl; j++){
            DC[i][j] = 0.0;
            for(k = 0; k < p; k++){
                DC[i][j] += PeltyD(i, k) * C[k][j];
            }
        }
    }
    for(i = 0; i < tl; i++){
        for(j = i; j < tl; j++){
            S[i][j] = 0.0;
            for(k = 0; k < pd; k++){
                S[i][j] += DC[k][i] * DC[k][j];
            }
            S[j][i] = S[i][j];
        }
    }
    
    // coefs initialization
    theta_init = coef_initial(Y, X, tau, beta_in, model, lambda, PeltyD);
    for(k = 0; k < tl; k++){
        theta_old[k] = 0.0;
        theta_new[k] = theta_init[k];
    }
    
    for(l = 1.0, count1 = 0; count1 <= maxint1; count1++){
        // intermediate coefs
        for(k = 0; k < tl; k++){
            theta_y[k] = theta_new[k] + count1 / (count1 + 3.0) * (theta_new[k] - theta_old[k]);
        }
        // obtain current residuals
        for(i = 0; i < n; i++){
            resid_y[i] = Y[i];
            for(k = 0; k < tl; k++)
                resid_y[i] -= XC[i][k] * theta_y[k];  
        }
        // obtain current majorize function value
        rho_y = Q_xi(resid_y, resid_y, tau, eps);
        f_y = 0.0;
        for(i = 0; i < n; i++){
            f_y += rho_y[i] * weights[i];
            v_eps[i] = 1.0 - 2.0 * tau - resid_y[i] / (eps + abs(resid_y[i])); 
        }
        f_y /= n;
        for(k = 0; k < tl; k++){
            Qxi_prime[k] = 0.0;
            for(i = 0; i < n; i++){
                Qxi_prime[k] += XC[i][k] * v_eps[i] * weights[i];
            }
            Qxi_prime[k] /= (2.0 * n);
        }
        
        // backtracking line search
        for(count2 = 0; count2 <= maxint2; count2++, l *= step_size){
            for(k = 0; k < tl; k++){
                mu[k] = theta_y[k] - l * Qxi_prime[k];
            }
            if(model != 0){
                for(k = 0; k < tl; k++){
                    if(mu[k] > 0.0)
                        mu[k] = 0.0;
                }
            }
            for(i = 0; i < tl; i++){
                for(j = i; j < tl; j++){
                    SI[i][j] = S[i][j] * 2.0 * l * lambda;
                    SI[j][i] = SI[i][j];
                }
                SI[i][i] += 1.0;
            }
            
            LUdcmp M(SI);
            M.solve(mu, theta_z);
            
            for(i = 0; i < n; i++){
                resid_z[i] = Y[i];
                for(k = 0; k < tl; k++){
                    resid_z[i] -= XC[i][k] * theta_z[k];
                }
            }
            rho_z = Q_xi(resid_z, resid_y, tau, eps);
            f_z = 0.0;
            for(i = 0; i < n; i++){
                f_z += rho_z[i] * weights[i];
            }
            f_z /= n;
            sum1 = 0.0;
            sum2 = 0.0;
            for(k = 0; k < tl; k++){
                sum1 += Qxi_prime[k] * (theta_z[k] - theta_y[k]);
                sum2 += pow(theta_z[k] - theta_y[k], 2.0);
            }
            if(f_z <= (f_y + sum1 + sum2 / (2.0 * l)))
                break;
        }
        for(i = 0; i < n; i++){
            resid_new[i] = Y[i];
            for(k = 0; k < tl; k++){
                resid_new[i] -= XC[i][k] * theta_new[k];
            }
        }
        
        // condition to jump off
        vecs1 = Q_xi(resid_new, resid_new, tau, eps);
        vecs2 = Q_xi(resid_z, resid_new, tau, eps);
        sum3 = 0.0;
        for(i = 0; i < n; i++){
            sum3 += (vecs1[i] - vecs2[i]) * weights[i];
        }
        if(abs(sum3) / n < threshold)
            break;
        for(k = 0; k < tl; k++){
            theta_old[k] = theta_new[k];
            theta_new[k] = theta_z[k];
        }
    }
    for(j = 0; j < p; j++){
        beta_out[j] = 0.0;
        for(k = 0; k < tl; k++){
            beta_out[j] += C[j][k] * theta_z[k];
        }
    }
    return beta_out;
}

// [[Rcpp::export]]
NumericVector checkloss(NumericVector res, double tau){
    /*
    * Check-loss function
    * 
    * Input: 
    *      res, (n x 1) vector,
    *          residuals;
    *      tau, numeric,
    *          quantile level, should be in (0, 1).
    *  
    *  Output:
    *      out, (n x 1) vector,
    *          check-loss function value.
    */
    double pos, neg;
    int l = res.size(), i;
    NumericVector out(l);
    
    for(i = 0; i < l; i++){
        pos = (abs(res[i]) + res[i]) / 2;
        neg = (abs(res[i]) - res[i]) / 2;
        out[i] = (tau * pos + (1 - tau) * neg);
    }
    return out;
}

NumericVector coef_initial(NumericVector Y, NumericMatrix X, double tau, NumericVector beta_in, int model, double lambda, NumericMatrix PeltyD){
    int p = X.ncol(), n = X.nrow(), i, j, k1, k2, pd = PeltyD.nrow(), tl = (model != 0) ? (p + 1) : p;
    double num_val;
    
    NumericVector theta_init(tl);
    
    VecDoub sol(p);
    MatDoub trans_mat(p - 1, p);
    trans_mat.assign(p - 1, p, 0.0);
    
    if(!Rf_isNull(beta_in) && (Rf_length(beta_in) > 1)){
        for(j = 0; j < p; j++){
            sol[j] = beta_in[j];
        }
    }else{
        MatDoub XXD(p, p);
        VecDoub XY(p);
        
        for(i = 0; i < p; i++){
            for(j = i; j < p; j++){
                XXD[i][j] = 0.0;
                for(k1 = 0; k1 < n; k1++){
                    XXD[i][j] += X(k1, i) * X(k1, j);
                }
                num_val = 0.0;
                for(k2 = 0; k2 < pd; k2++){
                    num_val += PeltyD(k2, i) * PeltyD(k2, j);
                }
                XXD[i][j] += lambda * num_val;
                XXD[j][i] = XXD[i][j];
            }
        }
        
        for(j = 0; j < p; j++){
            XY[j] = 0.0;
            for(k1 = 0; k1 < n; k1++){
                XY[j] += X(k1, j) * Y[k1]; 
            }
        }
        
        LUdcmp M(XXD);
        M.solve(XY, sol);
    }
    
    if(model != 0){
        if(model == 1){
            for(i = 0; i < p - 1; i++){
                trans_mat[i][i] = 1.0;
                trans_mat[i][i + 1] = -1.0;
            }
            if(sol[0] <= 0.0){
                theta_init[0] = sol[0];
                theta_init[1] = 0.0;
            }else{
                theta_init[0] = 0.0;
                theta_init[1] = -sol[0];
            }
        }else{
            for(i = 0; i < p - 1; i++){
                trans_mat[i][i] = -1.0;
                trans_mat[i][i + 1] = 1.0;
            }
            if(sol[0] <= 0.0){
                theta_init[0] = 0.0;
                theta_init[1] = sol[0];
            }else{
                theta_init[0] = -sol[0];
                theta_init[1] = 0.0;
            }
        }
        for(k1 = 2; k1 < tl; k1++){
            theta_init[k1] = 0.0;
            for(k2 = 0; k2 < p; k2++){
                theta_init[k1] += trans_mat[k1 - 2][k2] * sol[k2];
            }
            theta_init[k1] = (theta_init[k1] <= 0.0) ? theta_init[k1] : 0.0;
        }
    }else{
        for(j = 0; j < p; j++){
            theta_init[j] = sol[j];
        }
    }
    return theta_init;
}

NumericVector Q_xi(NumericVector r, NumericVector rk, double tau, double eps){
    int n = r.size(), i;
    NumericVector out(n);
    
    for(i = 0; i < n; i++){
        out[i] = pow(r[i], 2.0) / (eps + abs(rk[i])) + (4.0 * tau - 2.0) * r[i];
        out[i] /= 4.0;
    }
    return out;
}


// [[Rcpp::export]]
NumericVector proxy_inital(NumericVector Y, NumericMatrix X, double tau, NumericVector weights, int model, double delta, double lambda, NumericMatrix PeltyD, int maxint1, int maxint2, double threshold){
    int p = X.ncol(), n = X.nrow(), pd = PeltyD.nrow(), i, j, k, count1, count2, tl = (model != 0) ? (p + 1) : p;
    double l, step_size = 0.5, num_val, f_y, sum1, sum2, sum3;
    NumericVector beta_out(p), resid_y(n), resid_upd(n), resid_new(n), rho_y(n), rho_upd(n), rho_new(n);
    
    MatDoub C(p, tl), XC(n, tl), S(tl, tl), DC(pd, tl);
    VecDoub theta_old(tl), theta_new(tl), theta_y(tl), theta_upd(tl), f_prime(tl), vecs1(n), vecs2(tl);
    
    NumericVector checkloss(NumericVector res, double tau);
    
    C.assign(p, tl, 0.0);
    
    if(model == 1){
        // Set monotone constraint matrix C
        for(i = 0; i < p; i++){
            C[i][0] = 1.0;
            for(j = 1; j < i + 1; j++){
                C[i][j] = -1.0;
            }
        }
    }else if(model == 2){
        for(i = 0; i < p; i++){
            C[i][0] = -1.0;
            for(j = 1; j < i + 1; j++){
                C[i][j] = 1.0;
            }
        }
    }else{
        for(i = 0; i < p; i++){
            C[i][i] = 1.0;
        }
    }
    
    // Store XC to save computation time
    for(i = 0; i < n; i++){
        for(j = 0; j < tl; j++){
            XC[i][j] = 0.0;
            for(k = 0; k < p; k++){
                XC[i][j] += X(i, k) * C[k][j];
            }
        }
    }
    // Store S to save computation time
    for(i = 0; i < pd; i++){
        for(j = 0; j < tl; j++){
            DC[i][j] = 0.0;
            for(k = 0; k < p; k++){
                DC[i][j] += PeltyD(i, k) * C[k][j];
            }
        }
    }
    for(i = 0; i < tl; i++){
        for(j = i; j < tl; j++){
            S[i][j] = 0.0;
            for(k = 0; k < pd; k++){
                S[i][j] += DC[k][i] * DC[k][j];
            }
            S[j][i] = S[i][j];
        }
    }
    
    // coefs initialization
    for(k = 0; k < tl; k++){
        theta_old[k] = 0.0;
        theta_new[k] = 0.0;
    }
    
    l = 1.0;
    for(count1 = 0; count1 <= maxint1; count1++){
        // intermediate coefs
        for(k = 0; k < tl; k++){
            theta_y[k] = theta_new[k] + count1 / (count1 + 3.0) * (theta_new[k] - theta_old[k]);
        }
        // obtain current residuals
        for(i = 0; i < n; i++){
            resid_y[i] = Y[i];
            for(k = 0; k < tl; k++)
                resid_y[i] -= XC[i][k] * theta_y[k];  
        }
        // obtain current function value and derivative
        rho_y = checkloss(resid_y, tau);
        num_val = 0.0;
        for(i = 0; i < tl; i++){
            vecs2[i] = 0.0;
            for(k = 0; k < tl; k++){
                vecs2[i] += theta_y[k] * S[k][i]; 
            }
            num_val += vecs2[i] * theta_y[i];
        }
        f_y = lambda * num_val;
        for(i = 0; i < n; i++){
            f_y += rho_y[i] * weights[i] / n;
            if(abs(resid_y[i]) > delta){
                vecs1[i] = tau - ((resid_y[i] < 0.0) ? 1.0 : 0.0);
            }else{
                vecs1[i] = 2.0 * abs(resid_y[i]) / delta * (tau - ((resid_y[i] < 0.0) ? 1.0 : 0.0));
            }
        }
        for(k = 0; k < tl; k++){
            f_prime[k] = 0.0;
            for(i = 0; i < n; i++){
                f_prime[k] += XC[i][k] * weights[i] * vecs1[i];
            }
            f_prime[k] = -f_prime[k] / n + 2.0 * lambda * vecs2[k];
        }
        
        for(count2 = 0; count2 <= maxint2; count2++){
            for(k = 0; k < tl; k++){
                theta_upd[k] = theta_y[k] - l * f_prime[k];
            }
            if(model != 0){
                for(k = 0; k < tl; k++){
                    if(theta_upd[k] > 0.0)
                        theta_upd[k] = 0.0;
                }
            }
            // jump out condition
            for(i = 0; i < n; i++){
                resid_upd[i] = Y[i];
                for(k = 0; k < tl; k++)
                    resid_upd[i] -= XC[i][k] * theta_upd[k];  
            }
            rho_upd = checkloss(resid_upd, tau);
            sum1 = 0.0;
            for(i = 0; i < n; i++){
                sum1 += rho_upd[i] * weights[i];      
            }
            sum2 = 0.0;
            sum3 = 0.0;
            for(k = 0; k < tl; k++){
                sum2 += f_prime[k] * (theta_upd[k] - theta_y[k]);
                sum3 += pow(theta_upd[k] - theta_y[k], 2.0);
            }
            if(sum1 / n <= f_y + sum2 + sum3 / (2.0 * l)) 
                break;
            l *= step_size;
        }
        
        for(i = 0; i < n; i++){
            resid_new[i] = Y[i];
            for(k = 0; k < tl; k++)
                resid_new[i] -= XC[i][k] * theta_new[k];  
        }
        rho_new = checkloss(resid_new, tau);
        sum1 = 0.0;
        for(i = 0; i < n; i++){
            sum1 += (rho_upd[i] - rho_new[i]) * weights[i];
        }
        if(abs(sum1) / n < threshold) 
            break;
        for(k = 0; k < tl; k++){
            theta_old[k] = theta_new[k];
            theta_new[k] = theta_upd[k];
        }
    }
    for(j = 0; j < p; j++){
        beta_out[j] = 0.0;
        for(k = 0; k < tl; k++){
            beta_out[j] += C[j][k] * theta_upd[k];
        }
    }
    return beta_out;
}