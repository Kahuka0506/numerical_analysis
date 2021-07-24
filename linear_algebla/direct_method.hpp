#ifndef DIRECT_METHOD_HPP
#define DIRECT_METHOD_HPP 1

pair<vvd,vvd> LU_decomposition(vvd A){
    
    int n = si(A);
    vvd L(n,vd(n,0.0)), U(n,vd(n,0.0));
    rep(i,n){
        L[i][i] = A[i][i];
        U[i][i] = 1.0;
        reps(j,i+1,n) L[j][i] = A[j][i],U[i][j] = A[i][j]/L[i][i];
        reps(j,i+1,n) reps(k,i+1,n) A[j][k] -= L[j][i]*U[i][k];
    }
    
    return make_pair(L,U);
}

vd solve_LU(vvd A, vd b){
    
    int n = si(A);
    auto [L,U] = LU_decomposition(A);
    
    vd y(n,0.0);
    rep(i,n){
        rep(j,i) b[i] -= L[i][j]*y[j];
        y[i] = b[i]/L[i][i];
    }
    
    vd x(n,0.0);
    per(i,n){
        reps(j,i+1,n) y[i] -= U[i][j]*x[j];
        x[i] = y[i];
    }
    
    return x;
}

vvd inv_LU(vvd A){
    
    int n = si(A);
    auto [L,U] = LU_decomposition(A);
    
    vvd inv_A(n,vd(n,0.0));
    vd e(n,0.0);
    rep(ii,n){
        e.assign(n,0.0);
        e[ii] = 1.0;

        vd y(n,0.0);
        rep(i,n){
            rep(j,i) e[i] -= L[i][j]*y[j];
            y[i] = e[i]/L[i][i];
        }
        vd x(n,0.0);
        per(i,n){
            reps(j,i+1,n) y[i] -= U[i][j]*x[j];
            x[i] = y[i];
        }
        
        rep(i,n) inv_A[i][ii] = x[i];
    }
    
    return inv_A;
}



#endif
