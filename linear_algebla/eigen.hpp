#ifndef EIGEN_LIB
#define EIGEN_LIB 1

#include "direct_method.hpp"

pair<double,vd> power_method(vvd A){
    
    int n = si(A);
    vd y(n),x(n,0.0);x[0] = 1.0;
    double lambda = 0.0;
    
    double d = 10.0;
    int loop = 0;
    while(d >= 0.00001){
        y.assign(n,0);
        rep(i,n) rep(j,n) y[i] += A[i][j]*x[j];
        double yy = 0.0;
        rep(i,n) yy += y[i]*y[i];
        yy = sqrt(yy);
        
        double lambda1 = lambda;
        lambda = 0.0;
        rep(i,n) lambda += x[i]*y[i];
        d = abs(lambda-lambda1);
        rep(i,n) x[i] = y[i]/yy;
        loop++;
    }
    //cout << "loop = " csp loop << endl;
    
    return make_pair(lambda,x);
}


pair<double,vd> power_method_shift(vvd A, double p){
    
    int n = si(A);
    vd y(n),x(n,0.0);x[0] = 1.0;
    double lambda = 0.0;
    rep(i,n) A[i][i] -= p;
    
    double d = 10.0;
    int loop = 0;
    while(d >= 0.00001){
        y.assign(n,0);
        rep(i,n) rep(j,n) y[i] += A[i][j]*x[j];
        double yy = 0.0;
        rep(i,n) yy += y[i]*y[i];
        yy = sqrt(yy);
        
        double lambda1 = lambda;
        lambda = 0.0;
        rep(i,n) lambda += x[i]*y[i];
        d = abs(lambda-lambda1);
        rep(i,n) x[i] = y[i]/yy;
        loop++;
    }
    
    //cout << "loop = " csp loop << endl;
    
    return make_pair(lambda+p,x);
}


pair<double,vd> inv_power_method(vvd A){
    
    vvd invA = inv_LU(A);
    auto [lambda,x] = power_method(invA);
    
    return make_pair(1.0/lambda,x);
}


#endif
