
int Gauss_Seidel(vector<vector<double>> &A, vector<double> &b, vector<double> &x){

    int N = int(A[0].size());
    
    double d = 10000000.0;
    int roop = 0;
    double b_nolm = 0;
    rep(i,N) b_nolm += b[i]*b[i];
    b_nolm = sqrt(b_nolm);


    while (d >= 0.00001) {

        vector<double> old_x = x;
        roop++;
        rep(i,N){
            x[i] = b[i];
            rep(j,N) {
                if(j == i) continue;
                x[i] -= A[i][j]*x[j];
            }
            x[i] /= A[i][i];
        }

        d = 0;
        vector<double> Ax(N,0);
        rep(i,N) rep(j,N) Ax[i] += A[i][j]*x[j];
        rep(i,N) d += (b[i]-Ax[i])*(b[i]-Ax[i]);
        d = sqrt(d);
        d /= b_nolm;

        if(roop % 10000 == 0) cout << roop csp d << endl;
    }

    //cout << "Gauss_Seidel" << endl;
    //cout << roop << endl;
    //outve(x);

    return roop;


}

