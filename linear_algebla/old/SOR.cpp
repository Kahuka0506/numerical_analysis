
int SOR(vector<vector<double>> &A, vector<double> &b, vector<double> &x, double omega){

    int N = int(A[0].size());
    
    double d = 10000000.0;
    int roop = 0;
    double b_nolm = 0;
    rep(i,N) b_nolm += b[i]*b[i];
    b_nolm = sqrt(b_nolm);

    double dd = 0.0;

    while (d >= 0.00001){
        roop++;
        dd = 0.0;
        rep(i,N){
            double a = b[i];
            rep(j,N){
                if(i == j) continue;
                a -= A[i][j]*x[j];
            }
            a /= A[i][i];
            dd += (omega*(a - x[i]))*(omega*(a - x[i]));
            x[i] = x[i] + omega*(a - x[i]);
        }

        d = 0;
        vector<double> Ax(N,0);
        rep(i,N) rep(j,N) Ax[i] += A[i][j]*x[j];
        rep(i,N) d += (b[i]-Ax[i])*(b[i]-Ax[i]);
        d = sqrt(d);
        d /= b_nolm;

        dd = sqrt(dd);
        //d = dd;

        //cout << roop csp d << endl;
    }

    //cout << "SOR" csp omega << endl;
    //cout << roop << endl;
    //outve(x);


    return roop;

}
