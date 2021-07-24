
void LU_factorization(vector<vector<double>> &A, vector<vector<double>> &L, vector<vector<double>> &U){

    int N = int(A[0].size());
    
    vector<vector<double>> AA = A;
    L.assign(N,vector<double>(N,0.0));
    U.assign(N,vector<double>(N,0.0));

    rep(m,N){

        L[m][m] = AA[m][m];
        U[m][m] = 1;
        reps(i,m+1,N) L[i][m] = AA[i][m], U[m][i] = AA[m][i]/AA[m][m];
        reps(i,m+1,N) reps(j,m+1,N) {
            AA[i][j] -= L[i][m]*U[m][j];
        }
    }

    //outmat(L);
    //cout << endl;
    //outmat(U);

}
void LU_method(vector<vector<double>> &L, vector<vector<double>> &U, vector<double> &b, vector<double> &x){

    int N = int(A[0].size());
    
    x.resize(N);

    vector<double> y(N);
    rep(i,N){
        y[i] = b[i];
        rep(j,i) y[i] -= L[i][j]*y[j];
        y[i] /= L[i][i];
    }

    for (int i = N-1; i >= 0; i--) {
        x[i] = y[i];
        reps(j,i+1,N) x[i] -= U[i][j]*x[j];
    }

    //outve(x);

}
void invers(vector<vector<double>>& A, vector<vector<double>>& A_inv){

    int N = int(A[0].size());
    
    A_inv = A;
    vector<vector<double>> L,U;
    LU_factorization(A,L,U);
    vector<double> b(N,0.0);
    vector<double> x(N);
    rep(i,N){
        if(i != 0) b[i-1] = 0;
        b[i] = 1;
        LU_method(L,U,b,x);
        rep(j,N) A_inv[j][i] = x[j];
    }

}


void LU(vector<vector<double>> &A, vector<double> &b, vector<double> &x){
    vector<vector<double>> L,U;

    int N = int(A[0].size());
    
    vector<vector<double>> AA = A;
    L.assign(N,vector<double>(N,0.0));
    U.assign(N,vector<double>(N,0.0));

    rep(m,N){
        L[m][m] = AA[m][m];
        U[m][m] = 1;
        reps(i,m+1,N) L[i][m] = AA[i][m], U[m][i] = AA[m][i]/AA[m][m];
        reps(i,m+1,N) reps(j,m+1,N) AA[i][j] -= L[i][m]*U[m][j];
    }

    x.resize(N);
    vector<double> y(N);
    rep(i,N){
        y[i] = b[i];
        rep(j,i) y[i] -= L[i][j]*y[j];
        y[i] /= L[i][i];
    }

    for (int i = N-1; i >= 0; i--) {
        x[i] = y[i];
        reps(j,i+1,N) x[i] -= U[i][j]*x[j];
    }

    return;
}
