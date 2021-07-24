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

class least_square {
    vd x,y;
    int n,m;
    
public:
    
    vd res;
    least_square(vd x_, vd y_, int m_):x(x_),y(y_),m(m_){
        n = si(x);
        res.assign(m,0.0);
        vvd A(n,vd(m));
        rep(i,n){
            A[i][0] = 1.0;
            reps(j,1,m) A[i][j] = A[i][j-1]*x[i];
        }
        vvd AA(m,vd(m,0.0));
        rep(i,m) rep(j,m) rep(k,n) AA[i][j] += A[k][i]*A[k][j];
        vd Ab(m,0.0);
        rep(i,m) rep(j,n) Ab[i] += A[j][i]*y[j];
        res = solve_LU(AA,Ab);
    }
    
    
    
    void plot_data(){
        double xs = *min_element(all(x));
        double xt = *max_element(all(x));
        
        double dx = 0.01;
        
        FILE* fp;
        fp = fopen("data_least_squares.txt","w");
        double yy,xi;
        while (xt > xs) {
            xi = 1.0;
            yy = 0.0;
            rep(i,n) yy += xi*res[i], xi *= xs;
            fprintf(fp,"%.5lf %.5lf\n",xs,yy);
            xs += dx;
        }
        fclose(fp);
    }
    
};

/*
vd least_square(vd x, vd y, int m){
    
    int n = si(x);
    vd res(m,0.0);
    vvd A(n,vd(m));
    rep(i,n){
        A[i][0] = 1.0;
        reps(j,1,m) A[i][j] = A[i][j-1]*x[i];
    }
    vvd AA(m,vd(m,0.0));
    rep(i,m) rep(j,m) rep(k,n) AA[i][j] += A[k][i]*A[k][j];
    vd Ab(m,0.0);
    rep(i,m) rep(j,n) Ab[i] += A[j][i]*y[j];
    
    res = solve_LU(AA,Ab);
    
    return res;
}

*/

void solve(){
    
    int N,M;
    cin >> N >> M;
    
    vd x(N),y(N);
    rep(i,N) cin >> x[i] >> y[i];
    
    //vd ans = least_square(x,y,M);
    //outve(ans);
    
    least_square LSM(x,y,M);
    outve(LSM.res);
    LSM.plot_data();
    
}



int main(){
    
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();

    return 0;
}
