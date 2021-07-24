

int BiCGSTAB(ve<ve<double>> &A, ve<double> &b, ve<double> &x){

    int N = si(A[0]);
    //x.assign(N,0);
    ve<double> r = b;
    rep(i,N) rep(j,N) r[i] -= A[i][j]*x[j];
    ve<double> p = r;
    ve<double> r_0 = r;
    ve<double> r_chil = r;
    ve<double> r_chil_0 = r;
    ve<double> p_chil = r;

    int roop = 0;
    double d = 100000.0;

    double b_nolm = 0;
    rep(i,N) b_nolm += b[i]*b[i];
    b_nolm = sqrt(b_nolm);

    double alpa=1.0,alpa_a,alpa_b,alpa_old=1.0;
    double beta=0.0,beta_a,beta_b,beta_old;
    ve<double> Ap(N),At(N);

    ve<double> t(N,0),t_old,u(N,0),u_old,w(N,0),w_old,z(N,0),z_old,y(N,0);

    double eta=0,zeta,zeta_a,zeta_b;

    while (d >= 0.000001) {

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

        rep(i,N) rep(j,N) At[i] += A[i][j]*t[j];

        zeta_a = 0,zeta_b = 0;
        rep(i,N){
            zeta_a += t[i]*At[i];
            zeta_b += At[i]*At[i];
        }
        zeta = zeta_a/zeta_b;


        rep(i,N) beta_b += r_chil_0[i]*r[i];
        beta_b *= zeta;

        rep(i,N){
            u[i] = zeta*Ap[i] + eta*(t_old[i]-r[i]+beta_old*u_old[i]);
            z[i] = zeta*r[i] + eta*z_old[i] - alpa*u[i];
            x[i] = x[i] + alpa*p[i]+z[i];
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



    //cout << "BiCGSTAB" << endl;
    //cout << roop << endl;
    //outve(x);


    return roop;
}
