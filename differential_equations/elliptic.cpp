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

#include "../linear_algebla/iterative_method.hpp"

class elliptic {
    int n;
    double dx;
    double u0 = 0.0, u1 = 0.0;
    
public:
    
    vd x;
    vd u;
    
    elliptic(){
        n = 50;
        dx = 2.0/(double)n;
        u.assign(n+1,0.0);
        x.assign(n+1,0.0);
        x[0] = -1.0;
        rep(i,n) x[i+1] = x[i] + dx;
        
        vvd A(n-1,vd(n-1,0.0));
        vd b(n-1);
        
        rep(i,n-1){
            A[i][i] = -2.0;
            if(i != 0) A[i][i-1] = 1.0;
            else b[i] -= u0;
            
            if(i!=n-2) A[i][i+1] = 1.0;
            else b[i] -= u1;
            
            if(abs(x[i+1]) <= dx) b[i] -= dx-abs(x[i+1]);
            //b[i] = sin(x[i+1]);
        }
        
        vd xx = Gauss_Seidel(A,b);
        reps(i,1,n) u[i] = xx[i-1];
        u[0] = u0;
        u.back() = u1;
        
        outve(u);
    }
    
    void plot_data(){
        FILE* fp;
        fp = fopen("./out/data_elliptic.txt","w");
        rep(i,n+1) fprintf(fp,"%.5lf %.5lf\n",x[i],u[i]);
        fclose(fp);
    }
};


void solve(){
    
    elliptic A;
    A.plot_data();
    
    
}



int main(){
    
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();

    return 0;
}
