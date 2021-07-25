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

//d^2x/dt^2 + adx/dt + bx = 0

class runge_kutta {
    double dt;
    vd x;
    ve<pair<double,double>> t_x;
    int n;
    
public:
    
    runge_kutta(){
        dt = 0.01;
        n = 2;
        x.assign(n,0.0);
        x[0] = 1.0;
        vd k1(n,2),k2(n,2),k3(n,2),k4(n,2);
        
        double a = 0.3, b = 4.0;
        vvd A = {{0.0,1.0},{-b,-a}};
        
        double t = 0.0;
        while (t < 20.0) {
            t_x.pb({t,x[0]});
            
            k1[0] = A[0][0]*x[0] + A[0][1]*x[1];
            k1[1] = A[1][0]*x[0] + A[1][1]*x[1];
        
            k2[0] = A[0][0]*(x[0]+k1[0]*dt/2.0) + A[0][1]*(x[1]+k1[1]*dt/2.0);
            k2[1] = A[1][0]*(x[0]+k1[0]*dt/2.0) + A[1][1]*(x[1]+k1[1]*dt/2.0);
        
            k3[0] = A[0][0]*(x[0]+k2[0]*dt/2.0) + A[0][1]*(x[1]+k2[1]*dt/2.0);
            k3[1] = A[1][0]*(x[0]+k2[0]*dt/2.0) + A[1][1]*(x[1]+k2[1]*dt/2.0);
        
            k4[0] = A[0][0]*(x[0]+k3[0]*dt) + A[0][1]*(x[1]+k3[1]*dt);
            k4[1] = A[1][0]*(x[0]+k3[0]*dt) + A[1][1]*(x[1]+k3[1]*dt);
            
            x[0] += dt/6.0*(k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]);
            x[1] += dt/6.0*(k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1]);
            
            t += dt;
        }
    }
    
    
    void plt_data(){
        FILE* fp;
        fp = fopen("./out/data_ode.txt","w");
        rep(i,si(t_x)) fprintf(fp,"%.4lf %.4lf\n", t_x[i].fi, t_x[i].se);
        fclose(fp);
        return;
    }
};



void solve(){
    
    runge_kutta A;
    A.plt_data();
    
    
}



int main(){
    
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();

    return 0;
}
