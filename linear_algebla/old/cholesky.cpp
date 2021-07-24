


void cholesky_factrization(vector<vector<double>>& A, vector<vector<double>>& L){

    int N = int(A[0].size());
    
    L.assign(N,vector<double>(N,0.0));

    rep(i,N){

        rep(j,N){
            if(i == j){
                L[i][i] = A[i][j];
                rep(k,j) L[i][j] -=  L[i][k]*L[i][k];
                L[i][i] = sqrt(L[i][i]);
            }else if(i > j){
                L[i][j] = A[i][j];
                rep(k,j) L[i][j] -= L[i][k]*L[j][k];
                L[i][j] /= L[j][j];
            }
        }


    }

    return;

}


void solve_cholesky(vector<vector<double>>& L, vector<double>& b, vector<double>& x){

    int N = int(x.size());

    vector<double> y(N,0.0);
    rep(i,N){
        y[i] = b[i];
        rep(j,i) y[i] -= y[j]*L[i][j];
        y[i] /= L[i][i];
    }

    for(int i = N-1; i >= 0; i--){
        x[i] = y[i];
        reps(j,i+1,N) x[i] -= x[j]*L[j][i];
        x[i] /= L[i][i];
    }

    return;

}

void cholesky(vector<vector<double>> &A, vector<double> &b, vector<double> &x){
    
    int N = int(A[0].size());
    
    vector<vector<double>> L(N,vector<double>(N,0.0));
    
    x.resize(N);
    cholesky_factrization(A,L);
    
    vector<double> y(N,0.0);
    rep(i,N){
        y[i] = b[i];
        rep(j,i) y[i] -= y[j]*L[i][j];
        y[i] /= L[i][i];
    }

    for(int i = N-1; i >= 0; i--){
        x[i] = y[i];
        reps(j,i+1,N) x[i] -= x[j]*L[j][i];
        x[i] /= L[i][i];
    }

    return;
}



