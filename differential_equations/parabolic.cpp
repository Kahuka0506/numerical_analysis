#include <bits/stdc++.h>
#define rep(i,n) for(int i=0; i<n; i++)
#define reps(i,s,n) for(int i=s; i<n; i++)
#define per(i,n) for(int i=n-1; i>=0; i--)
#define pers(i,n,s) for(int i=n-1; i>=s; i--)
#define all(v) v.begin(),v.end()
#define fi first
#define se second
#define pb push_back
#define si(v) int(v.size())
#define lb(v,n) lower_bound(all(v),n)
#define lbi(v,n) int(lower_bound(all(v),n) - v.begin())
#define ub(v,n) upper_bound(all(v),n)
#define ubi(v,n) int(upper_bound(all(v),n) - v.begin())
 
#define mod 1000000007
#define infi 1010000000
#define infl 1100000000000000000
 
#define outve(v) for(auto i : v) cout << i << " ";cout << endl
#define outmat(v) for(auto i : v){for(auto j : i) cout << j << " ";cout << endl;}
#define in(n,v) for(int i=0; i<(n); i++){cin >> v[i];}
#define IN(n,m,v) rep(i,n) rep(j,m){cin >> v[i][j];}
#define cyes cout << "Yes" << endl
#define cno cout << "No" << endl
#define cYES cout << "YES" << endl
#define cNO cout << "NO" << endl
#define csp << " " <<
#define outset(n) cout << fixed << setprecision(n);
 
using namespace std;
using ll = long long;
using ull = unsigned long long;
using ld = long double;
using vi = vector<int>;
using vvi = vector<vector<int>>;
using vd = vector<double>;
using vvd = vector<vector<double>>;
using vl = vector<ll>;
using vvl = vector<vector<ll>>;
using vs = vector<string>;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
template<typename T> using ve = vector<T>;
template<typename T> using pq2 = priority_queue<T>;
template<typename T> using pq1 = priority_queue<T,vector<T>,greater<T>>;
template<typename T> bool chmax(T &a, T b) {if(a < b) {a = b;return 1;}return 0;}
template<typename T> bool chmin(T &a, T b) {if(a > b) {a = b;return 1;}return 0;}


#include "../linear_algebla/direct_method.hpp"

class parabolic {
    int n;
    double lambda;
    double dx,dt;
    double u00,u01;
    
public:
    vd x;
    vd u;
    
    parabolic(){
        n = 10, u00 = 0.0, u01 = 1.0;
        lambda = 0.3;
        dx = 1.0/(double)n;
        dt = 0.01;
        
        x.assign(n+1,0);
        x[0] = 0.0;
        reps(i,1,n) x[i+1] = x[i] + dx;
        u.assign(n+1,0);
        u[0] = u00;
        u.back() = u01;
    }
    
    
    void solve_explicit(){
        vd u1;
        u.assign(n+1,0);
        u[0] = u00;
        u.back() = u01;
        
        rep(_,10){
            u1 = u;
            reps(i,1,n){
                u[i] = dt*lambda/dx/dx*(u1[i+1] - 2.0*u1[i] + u1[i-1]) + u1[i];
            }
        }
        
        cout << "Explicit" << endl;
        cout << "stability " csp lambda*dt/dx/dx << endl;
        outve(u);
        cout << endl;
        
        
        
        FILE* fp;
        fp = fopen("./out/data_parabolic_explicit.txt","w");
        rep(i,n+1) fprintf(fp,"%.5lf %.5lf\n", x[i],u[i]);
        fclose(fp);
    }
    
    
    
    
    void solve_implicit(){
        
        u.assign(n+1,0);
        u[0] = u00;
        u.back() = u01;
        vvd A(n-1,vd(n-1,0.0));
        vd b(n-1,0.0);
        double r = dt*lambda/dx/dx;
        
        vd ut(n-1,0.0);
        rep(i,n-1){
            A[i][i] = 1.0 + 2.0*r;
            if(i != 0) A[i][i-1] = -r;
            if(i != n-2) A[i][i+1] = -r;
        }
        
        rep(_,10){
            b = ut;
            b[0] += r*u00;
            b[n-2] += r*u01;
            
            ut = solve_LU(A,b);
        }
        reps(i,1,n) u[i] = ut[i-1];
        
        cout << "Implicit" << endl;
        cout << "stable" << endl;
        outve(u);
        cout << endl;
        
        
        FILE* fp;
        fp = fopen("./out/data_parabolic_implicit.txt","w");
        rep(i,n+1) fprintf(fp,"%.5lf %.5lf\n", x[i],u[i]);
        fclose(fp);
    }
    
    
    
    void solve_crank_nocolson(){
        
        
        u.assign(n+1,0);
        u[0] = u00;
        u.back() = u01;
        vvd A(n-1,vd(n-1,0.0));
        vvd A1(n-1,vd(n-1,0.0));
        vd b(n-1,0.0);
        double r = dt*lambda/dx/dx;
        
        vd ut(n-1,0.0);
        rep(i,n-1){
            A[i][i] = 2.0 + 2.0*r;
            if(i != 0) A[i][i-1] = -r;
            if(i != n-2) A[i][i+1] = -r;
            
        
            A1[i][i] = 2.0 - 2.0*r;
            if(i != 0) A1[i][i-1] = r;
            if(i != n-2) A1[i][i+1] = r;
        }
        
        rep(_,10){
            b.assign(n-1,0.0);
            rep(i,n-1) rep(j,n-1) b[i] += A1[i][j]*ut[j];
            b[0] += r*u00 + r*u00;
            b[n-2] += r*u01 + r*u01;
            
            ut = solve_LU(A,b);
        }
        reps(i,1,n) u[i] = ut[i-1];
        
        cout << "Crank Nicolson" << endl;
        cout << "stable" << endl;
        outve(u);
        cout << endl;
        
        FILE* fp;
        fp = fopen("./out/data_parabolic_crank_nocolson.txt","w");
        rep(i,n+1) fprintf(fp,"%.5lf %.5lf\n", x[i],u[i]);
        fclose(fp);
    }
    
    
};


void solve(){
    
    outset(5);
    
    parabolic A;
    A.solve_explicit();
    A.solve_implicit();
    A.solve_crank_nocolson();
    
    
}



int main(){
    
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();

    return 0;
}
