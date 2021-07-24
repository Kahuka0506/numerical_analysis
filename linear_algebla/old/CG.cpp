

int CG(vector<vector<double>> &A, vector<double> &b, vector<double> &x){

    int N = int(A[0].size());
    
    vector<double> r = b;
    rep(i,N) rep(j,N) r[i] -= A[i][j]*x[j];
    vector<double> p = r;

    double alpa = 0;
    double alpa_a = 0;
    double alpa_b = 0;

    double d = 1000000.0;

    double beta = 0;
    double beta_a = 0;
    double beta_b = 0;

    int roop = 0;

    vector<double> Ap(N,0);

    double b_nolm = 0;
    rep(i,N) b_nolm += b[i]*b[i];
    b_nolm = sqrt(b_nolm);


    while (d >= 0.00001) {
        roop++;
        d = 0;

        rep(i,N) Ap[i] = 0;
        rep(i,N) rep(j,N) Ap[i] += A[i][j]*p[j];

        alpa_a = 0,alpa_b = 0;
        beta_a = 0, beta_b = 0;

        rep(i,N) alpa_a += r[i]*r[i];
        rep(i,N) alpa_b += p[i]*Ap[i];
        alpa = alpa_a/alpa_b;

        rep(i,N) beta_b += r[i]*r[i];
        rep(i,N) x[i] = x[i] + alpa*p[i],r[i] = r[i] - alpa*Ap[i];
        rep(i,N) beta_a += r[i]*r[i];

        beta = beta_a/beta_b;

        rep(i,N) p[i] = r[i] + beta*p[i];


        rep(i,N) d += r[i]*r[i];
        d = sqrt(d);
        d /= b_nolm;
    }

    //cout << "CG" << endl;
    //cout << roop << endl;
    //outve(x);


    return roop;
}


