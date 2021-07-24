#ifndef ITERATIVE_METHOD_HPP
#define ITERATIVE_METHOD_HPP 1

#include "direct_method.hpp"


vd Jacobi(vvd A, vd b){
   
    int n = si(A);
    vd x(n,0.0);
    
    double d = 10.0;
    int loop = 0;
    while (d >= 0.000001) {
        loop++;
        vd x1 = x;
        rep(i,n) {
            x[i] = b[i];
            rep(j,n) if(i != j) x[i] -= A[i][j]*x1[j];
            x[i] /= A[i][i];
        }
        d = 0.0;
        rep(i,n) d += (x[i]-x1[i])*(x[i]-x1[i]);
        d = sqrt(d);
        //if(loop >= 20) break;
    }
    
    cout << "loop = " csp loop << endl;
    
    return x;
}


vd Gauss_Seidel(vvd A, vd b){
   
    int n = si(A);
    vd x(n,0.0);
    
    double d = 10.0;
    int loop = 0;
    while (d >= 0.000001) {
        loop++;
        vd x1 = x;
        rep(i,n) {
            x[i] = b[i];
            rep(j,n) if(i != j) x[i] -= A[i][j]*x[j];
            x[i] /= A[i][i];
        }
        d = 0.0;
        rep(i,n) d += (x[i]-x1[i])*(x[i]-x1[i]);
        d = sqrt(d);
        //if(loop >= 20) break;
    }
    
    cout << "loop = " csp loop << endl;
    
    return x;
}


vd SOR(vvd A, vd b, double omega = 1.2){
   
    int n = si(A);
    vd x(n,0.0);
    
    double d = 10.0;
    int loop = 0;
    while (d >= 0.000001) {
        loop++;
        vd x1 = x;
        vd xx(n);
        rep(i,n) {
            xx[i] = b[i];
            rep(j,n) if(i != j) xx[i] -= A[i][j]*x1[j];
            xx[i] /= A[i][i];
            x[i] += omega*(xx[i]-x[i]);
        }
        
        d = 0.0;
        rep(i,n) d += (x[i]-x1[i])*(x[i]-x1[i]);
        d = sqrt(d);
        //if(loop >= 20) break;
    }
    
    cout << "loop = " csp loop << endl;
    
    return x;
}



#endif
