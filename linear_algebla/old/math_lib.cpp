#include<iostream>
#include<string>
#include<cstdio>
#include<cstring>
#include<vector>
#include<cmath>
#include<algorithm>
#include<functional>
#include<fstream>
#include<iomanip>
#include<queue>
#include<deque>
#include<ciso646>
#include<random>
#include<map>
#include<set>
#include<complex>
#include<bitset>
#include<stack>
#include<unordered_map>
#include<unordered_set>
#include<utility>
#include<cassert>

#include"math_lib.h"

#define rep(i,n) for(int i=0; i<(n); i++)
#define reps(i,s,n) for(int i=(s); i<(n); i++)
#define all(v) v.begin(),v.end()
#define outve(v) for(auto i : v) cout << i << " ";cout << endl
#define outmat(v) for(auto i : v){for(auto j : i) cout << j << " ";cout << endl;}
#define in(n,v) for(int i=0; i<(n); i++){cin >> v[i];}
#define fi first
#define se second
#define pb push_back
#define si(v) int(v.size())
#define csp << " " <<
#define outset(n) cout << fixed << setprecision(n);
using namespace std;
using ll = long long;
using vi = vector<int>;
using vvi = vector<vector<int>>;
using vd = vector<double>;
using vvd = vector<vector<double>>;
template<typename T> using ve = vector<T>;
template<typename T> bool chmax(T &a, T b) {if(a < b) {a = b;return 1;}return 0;}
template<typename T> bool chmin(T &a, T b) {if(a > b) {a = b;return 1;}return 0;}


void LU_factorization(vvd &A, vvd &L, vvd &U){

    int N = si(A[0]);
    vvd AA = A;
    L.assign(N,vd(N,0.0));
    U.assign(N,vd(N,0.0));

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
void LU_method(vvd &L, vvd &U, vd &b, vd &x){

    int N = si(L[0]);
    x.resize(N);

    vd y(N);
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
void invers(vvd& A, vvd& A_inv){

    int N = si(A[0]);
    A_inv = A;
    vvd L,U;
    LU_factorization(A,L,U);
    vd b(N,0.0);
    vd x(N);
    rep(i,N){
        if(i != 0) b[i-1] = 0;
        b[i] = 1;
        LU_method(L,U,b,x);
        rep(j,N) A_inv[j][i] = x[j];
    }

}


void LU(vvd &A, vd &b, vd &x){
    vvd L,U;

    int N = si(A[0]);
    vvd AA = A;
    L.assign(N,vd(N,0.0));
    U.assign(N,vd(N,0.0));

    rep(m,N){
        L[m][m] = AA[m][m];
        U[m][m] = 1;
        reps(i,m+1,N) L[i][m] = AA[i][m], U[m][i] = AA[m][i]/AA[m][m];
        reps(i,m+1,N) reps(j,m+1,N) AA[i][j] -= L[i][m]*U[m][j];
    }

    x.resize(N);
    vd y(N);
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



void cholesky_factrization(vvd& A, vvd& L){

    int N = si(A[0]);
    L.assign(N,vd(N,0.0));

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


void solve_cholesky(vvd& L, vd& b, vd& x){

    int N = si(x);

    vd y(N,0.0);
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

void cholesky(vvd &A, vd &b, vd &x){
    
    int N = si(A[0]);
    vvd L(N,vd(N,0.0));
    
    x.resize(N);
    cholesky_factrization(A,L);
    
    vd y(N,0.0);
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




int Jacobi(ve<ve<double>> &A, ve<double> &b, ve<double> &x){

    int N = si(A[0]);
    //x.assign(N,0);
    ve<double> x_old;
    double d = 10000000.0;
    int roop = 0;
    double b_nolm = 0;
    rep(i,N) b_nolm += b[i]*b[i];
    b_nolm = sqrt(b_nolm);


    ve<double> Ax(N,0);

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

int Gauss_Seidel(ve<ve<double>> &A, ve<double> &b, ve<double> &x){

    int N = si(A[0]);
    //x.assign(N,0);
    double d = 10000000.0;
    int roop = 0;
    double b_nolm = 0;
    rep(i,N) b_nolm += b[i]*b[i];
    b_nolm = sqrt(b_nolm);


    while (d >= 0.00001) {

        ve<double> old_x = x;
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
        ve<double> Ax(N,0);
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


int SOR(ve<ve<double>> &A, ve<double> &b, ve<double> &x, double omega){

    int N = si(A[0]);
    //x.assign(N,0);
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
        ve<double> Ax(N,0);
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


int CG(ve<ve<double>> &A, ve<double> &b, ve<double> &x){

    int N = si(A[0]);
    //x.assign(N,0);
    ve<double> r = b;
    rep(i,N) rep(j,N) r[i] -= A[i][j]*x[j];
    ve<double> p = r;

    double alpa = 0;
    double alpa_a = 0;
    double alpa_b = 0;

    double d = 1000000.0;

    double beta = 0;
    double beta_a = 0;
    double beta_b = 0;

    int roop = 0;

    ve<double> Ap(N,0);

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





int BCG(ve<ve<double>> &A, ve<double> &b, ve<double> &x){

    int N = si(A[0]);
    //x.assign(N,0);
    ve<double> r = b;
    rep(i,N) rep(j,N) r[i] -= A[i][j]*x[j];
    ve<double> p = r;
    ve<double> r_chil = r;
    ve<double> p_chil = r;

    int roop = 0;
    double d = 100000.0;

    double b_nolm = 0;
    rep(i,N) b_nolm += b[i]*b[i];
    b_nolm = sqrt(b_nolm);

    double alpa,alpa_a,alpa_b;
    double beta,beta_a,beta_b;
    ve<double> Ap(N);
    ve<double> ATp(N);

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




int CGS(ve<ve<double>> &A, ve<double> &b, ve<double> &x){

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

    double alpa=1,alpa_a,alpa_b,alpa_old;
    double beta=0,beta_a,beta_b,beta_old;
    ve<double> Ap(N),At(N);

    ve<double> t(N,0),t_old,u(N,0),u_old,w(N,0),w_old,z(N,0),z_old,y(N,0);

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
