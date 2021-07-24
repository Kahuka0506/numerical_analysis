


int CGS(vector<vector<double>> &A, vector<double> &b, vector<double> &x){

    int N = int(A[0].size());
    
    vector<double> r = b;
    rep(i,N) rep(j,N) r[i] -= A[i][j]*x[j];
    vector<double> p = r;
    vector<double> r_0 = r;
    vector<double> r_chil = r;
    vector<double> r_chil_0 = r;
    vector<double> p_chil = r;

    int roop = 0;
    double d = 100000.0;

    double b_nolm = 0;
    rep(i,N) b_nolm += b[i]*b[i];
    b_nolm = sqrt(b_nolm);

    double alpa=1,alpa_a,alpa_b,alpa_old;
    double beta=0,beta_a,beta_b,beta_old;
    vector<double> Ap(N),At(N);

    vector<double> t(N,0),t_old,u(N,0),u_old,w(N,0),w_old,z(N,0),z_old,y(N,0);

    double eta,zeta;

    while (d >= 0.00001) {

        roop++;

        t_old=t, u_old=u, w_old=w, z_old=z;
        beta_old=beta, alpa_old=alpa;

        rep(i,N) Ap[i] = 0;
        rep(i,N) rep(j,N) Ap[i] += A[i][j]*p[j];

        alpa_a = 0,alpa_b = 0;
        beta_a = 0,beta_b = 0;

        rep(i,N) alpa_a+=r_chil_0[i]*r[i], alpa_b+=r_chil_0[i]*Ap[i];
        alpa = alpa_a/alpa_b;

        rep(i,N) {
            y[i] = t_old[i]-r[i]-alpa*(w_old[i]-Ap[i]);
            t[i] = r[i] - alpa*Ap[i];
            At[i] = 0;
        }

        eta = alpa*beta_old/alpa_old;
        zeta = alpa;

        rep(i,N) rep(j,N) At[i] += A[i][j]*t[j];

        rep(i,N) beta_b += r_chil_0[i]*r[i];
        beta_b *= zeta;

        rep(i,N){
            u[i] = zeta*Ap[i] + eta*(t_old[i]-r[i]+beta_old*u_old[i]);
            z[i] = zeta*r[i] + eta*z_old[i] - alpa*u[i];
            x[i] = x[i] + alpa*p[i] + z[i];
            r[i] = t[i] - eta*y[i] - zeta*At[i];
        }

        rep(i,N) beta_a += r_chil_0[i]*r[i];
        beta_a *= alpa;
        beta = beta_a/beta_b;

        rep(i,N){
            w[i] = At[i] + beta*Ap[i];
            p[i] = r[i] + beta*(p[i]-u[i]);
        }

        d = 0;
        rep(i,N) d += r[i]*r[i];
        d = sqrt(d);
        d /= b_nolm;

        //cout << roop csp d << endl;
    }



    //cout << "CGS" << endl;
    //cout << roop << endl;
    //outve(x);


    return roop;
}

