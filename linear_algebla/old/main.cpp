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






int main(){

    ios::sync_with_stdio(false);
    cin.tie(nullptr);


    int N;
    cin >> N;
    vvd A(N,vd(N));
    rep(i,N) rep(j,N) cin >> A[i][j];
    
    vd b(N);
    rep(i,N) cin >> b[i];
    
    vd x(N);
    
    LU(A,b,x);
    outve(x);
    
    cholesky(A,b,x);
    outve(x);
    
    Jacobi(A,b,x);
    outve(x);
    
    return 0;
}
