
int BCG(vector<vector<double>> &A, vector<double> &b, vector<double> &x){

    int N = int(A[0].size());
    
    vector<double> r = b;
    rep(i,N) rep(j,N) r[i] -= A[i][j]*x[j];
    vector<double> p = r;
    vector<double> r_chil = r;
    vector<double> p_chil = r;

    int roop = 0;
    double d = 100000.0;

    double b_nolm = 0;
    rep(i,N) b_nolm += b[i]*b[i];
    b_nolm = sqrt(b_nolm);

    double alpa,alpa_a,alpa_b;
    double beta,beta_a,beta_b;
    vector<double> Ap(N);
    vector<double> ATp(N);

    while (d >= 0.00001) {

        roop++;
        rep(i,N) Ap[i] = 0,ATp[i] = 0;
        rep(i,N) rep(j,N) Ap[i] += A[i][j]*p[j],ATp[i] += A[j][i]*p_chil[j];

        alpa_a = 0.0,alpa_b = 0.0;
        beta_a = 0.0,beta_b = 0.0;

        rep(i,N) alpa_a+=r_chil[i]*r[i], alpa_b+=p_chil[i]*Ap[i];
        alpa = alpa_a/alpa_b;


        rep(i,N) beta_b += r_chil[i]*r[i];
        rep(i,N) {
            x[i] = x[i] + alpa*p[i];
            r[i] = r[i] - alpa*Ap[i];
            r_chil[i] = r_chil[i] - alpa*ATp[i];
        }
        rep(i,N) beta_a += r_chil[i]*r[i];

        beta = beta_a/beta_b;


        rep(i,N) p[i] = r[i] + beta*p[i];
        rep(i,N) p_chil[i] = r_chil[i] + beta*p_chil[i];

        d = 0;
        rep(i,N) d += r[i]*r[i];
        d = sqrt(d);
        d /= b_nolm;

        //cout << roop csp d << endl;
    }



    //cout << "BCG" << endl;
    //cout << roop << endl;
    //outve(x);


    return roop;
}


