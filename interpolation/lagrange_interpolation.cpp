#include <bits/stdc++.h>
using namespace std;





class lagrange_interpolation {
    int N;
    vector<double> x,y;
    
public:
    lagrange_interpolation(vector<double> x_, vector<double> y_):x(x_),y(y_){
        N = int(x.size());
    }
    
    void plot_data(){
        double xs = *min_element(x.begin(),x.end());
        double xt = *max_element(x.begin(),x.end());
        
        double dx = 0.01;
        
        FILE* fp;
        fp = fopen("data_lagrange_interpolation.txt","w");
        
        double Lk,L;
        
        while (xt > xs) {
            L = 0.0;
            for(int i = 0; i < N; i++) {
                Lk = 1.0;
                for(int j = 0; j < N; j++){
                    if(i == j) continue;
                    Lk *= (xs-x[j])/(x[i]-x[j]);
                }
                L += Lk*y[i];
            }
            
            fprintf(fp,"%.5lf %.5lf\n",xs,L);
   
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
    
    lagrange_interpolation L(x,y);
    L.plot_data();
    
    
}



int main(){
    
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();

    return 0;
}
