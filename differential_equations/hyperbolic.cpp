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


class hyperbollic {
    int n;
    double dx,dt;
    double c;
    
public:
    vd x,u,u1,u2;
    
    hyperbollic(){
        n = 30;
        c = 1.0;
        dx = 2.0/(double)n;
        dt = 0.0002;
        x.assign(n+1,0.0);
        rep(i,n) x[i+1] = x[i] + dx;
        u.assign(n+1,0);
        u1.assign(n+1,0);
        
        rep(i,n+1) u[i] = 0.05*x[i]*(2.0-x[i]);
        rep(i,n+1) u1[i] = u[i] + dt*dt/2.0*c*c*(0.05*(-2.0));
        
        rep(_,3001){
            u2 = u1;
            u1 = u;
            reps(i,1,n){
                u[i] = (dt/dx)*(dt/dx)*c*c*(u1[i+1]-2.0*u1[i]+u1[i-1]) + 2.0*u1[i] - u2[i];
            }
            
            
            if(_%1000 == 0){
                FILE* fp;
                string s = "./out/data_hyperbollic_" + to_string(_) + ".txt";
                fp = fopen(s.c_str(),"w");
                rep(i,n+1) fprintf(fp,"%.5lf %.5lf\n", x[i], u[i]);
                fclose(fp);
            }
        }
        outve(u);
    }
};


void solve(){
    
    hyperbollic A;
    
    
}



int main(){
    
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();

    return 0;
}
