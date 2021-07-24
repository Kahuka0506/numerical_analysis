

int Jacobi(vector<vector<double>> &A, vector<double> &b, vector<double> &x){

    int N = int(A[0].size());
    
    vector<double> x_old;
    double d = 10000000.0;
    int roop = 0;
    double b_nolm = 0;
    rep(i,N) b_nolm += b[i]*b[i];
    b_nolm = sqrt(b_nolm);


    vector<double> Ax(N,0);

    while (d >= 0.00001) {
        roop++;
        x_old = x;

        rep(i,N){
            x[i] = b[i] - (Ax[i] - A[i][i]*x_old[i]);
            x[i] /= A[i][i];
            Ax[i] = 0;
        }

        d = 0;
        rep(i,N) rep(j,N) Ax[i] += A[i][j]*x[j];
        rep(i,N) d += (b[i]-Ax[i])*(b[i]-Ax[i]);
        d = sqrt(d);
        d /= b_nolm;

        //cout << roop csp d << endl;
    }


    //cout << "Jacobi" << endl;
    //cout << roop << endl;
    //outve(x);

    return roop;
}


