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




vector<double> LU(vector<vector<double>> A, vector<double> b){
    
    int n = int(b.size());
    vector<double> x(n,0.0);

    
    vector<vector<double>> L(n,vector<double>(n,0));
    vector<vector<double>> U(n,vector<double>(n,0));
    rep(i,n){
        reps(j,i,n) L[j][i] = A[j][i];
        reps(j,i,n) U[i][j] = A[i][j]/L[i][i];
        
        reps(j,i+1,n) reps(k,i+1,n) A[j][k] -= L[j][i]*U[i][k];
    }
    
    vector<double> y = x;
    rep(i,n){
        rep(j,i) b[i] -= L[i][j]*y[j];
        y[i] = b[i]/L[i][i];
    }
    per(i,n){
        reps(j,i+1,n) y[i] -= U[i][j]*x[j];
        x[i] = y[i];
    }
    
    return x;
}






class spline_interpolation {
    
    int N,sz;
    
    
public:
    vector<double> a,b,c,d,u,v;
    vector<vector<double>> A;
    vector<double> x,y;
    
    
    spline_interpolation(vector<double> x_, vector<double> y_):x(x_),y(y_){
        sz = int(x.size());
        N = sz-1;
        
        A.resize(N,vector<double>(N));
        u.resize(N);
        v.resize(N);
        
        for (int i = 0; i < N; i++) {
            double h0 = x[i+1]-x[i];
            double h1 = x[i+2]-x[i+1];
            A[i][i] = 2.0*(h0 + h1);
            if(i != 0) A[i][i-1] = h0;
            if(i != N-1) A[i][i+1] = h1;
        
            v[i] = 6.0*((y[i+2]-y[i+1])/h1 - (y[i+1]-y[i])/h0);
        }
        
        
        u = LU(A,v);
        u.insert(u.begin(),0.0);
        u.pb(0.0);
        
        for (int i = 0; i < N; i++) {
            a.push_back((u[i+1]-u[i])/(x[i+1]-x[i])/6.0);
            b.push_back(u[i]/2.0);
            c.push_back((y[i+1]-y[i])/(x[i+1]-x[i])-(x[i+1]-x[i])*(2.0*u[i]+u[i+1])/6.0);
            d.push_back(y[i]);
        }
    }
    
    
    void show(){
        for (int i = 0; i < N; i++) {
            printf("%.5lf %.5lf %.5lf %.5lf\n",a[i],b[i],c[i],d[i]);
        }
    }
    
    
    void plot_data(){
        double xs = *min_element(x.begin(),x.end());
        double xt = *max_element(x.begin(),x.end());
        
        double dx = 0.01;
        
        FILE* fp;
        fp = fopen("data_spline_interpolation.txt","w");
        
        double yy,xx;
        
        while (xt > xs) {
            int idx = int(upper_bound(x.begin(),x.end(),xs)-x.begin())-1;
            xx = xs-x[idx];
            yy = a[idx]*xx*xx*xx + b[idx]*xx*xx + c[idx]*xx + d[idx];
            
            fprintf(fp,"%.5lf %.5lf\n",xs,yy);
            xs += dx;
        }
        fclose(fp);
    }
    
};





void solve(){
    
    int N;
    cin >> N;
    vector<double> x(N),y(N);
    for (int i = 0; i < N; i++) {
        cin >> x[i] >> y[i];
    }
    
    spline_interpolation spline(x,y);
    spline.show();
    spline.plot_data();
    
}



int main(){
    
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();

    return 0;
}
