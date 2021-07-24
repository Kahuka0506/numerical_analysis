#ifndef MATH_LIB_H
#define MATH_LIB_H

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

using namespace std;
using ll = long long;
using vi = vector<int>;
using vvi = vector<vector<int>>;
using vd = vector<double>;
using vvd = vector<vector<double>>;



void LU_factorization(vvd &A, vvd &L, vvd &U);
void LU_method(vvd &L, vvd &U, vd &b, vd &x);
void LU_invers(vvd& A, vvd& A_inv);
void LU(vvd &A, vd &b, vd &x);
void cholesky_factrization(vvd& A, vvd& L);
void solve_cholesky(vvd& L, vd& b, vd& x);
void cholesky(vvd& A, vd& b, vd& x);

int Jacobi(vector<vector<double>> &A, vector<double> &b, vector<double> &x);
int Gauss_Seidel(vector<vector<double>> &A, vector<double> &b, vector<double> &x);
int SOR(vector<vector<double>> &A, vector<double> &b, vector<double> &x, double omega);


int CG(vector<vector<double>> &A, vector<double> &b, vector<double> &x);
int BCG(vector<vector<double>> &A, vector<double> &b, vector<double> &x);
int CGS(vector<vector<double>> &A, vector<double> &b, vector<double> &x);
int BiCGSTAB(vector<vector<double>> &A, vector<double> &b, vector<double> &x);



#endif
