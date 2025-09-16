#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <string>

using namespace std;

inline double y_of_x(double x) {
    return (x <= 0.0 || x >= 1.0) ? 0.0 : 0.1*sqrt(x*(1.0-x));
}

bool gauss(vector<vector<double>>& A, int rhs_col) {
    int n = A.size();
    for (int col=0, row=0; col<n && row<n; ++col, ++row) {
        int piv = row;
        for (int r=row+1; r<n; ++r)
            if (fabs(A[r][col]) > fabs(A[piv][col])) piv=r;
        if (fabs(A[piv][col])<1e-14) return false;
        swap(A[piv], A[row]);
        double div=A[row][col];
        for (int c=col;c<=rhs_col;++c) A[row][c]/=div;
        for (int r=0;r<n;++r) if(r!=row){
            double f=A[r][col];
            for (int c=col;c<=rhs_col;++c) A[r][c]-=f*A[row][c];
        }
    }
    return true;
}

struct DubletSolver {
    int N;
    double XS, XF;
    vector<double> T;
    vector<double> M;

    void build_endpoints() {
        const double PI=3.1415926585;
        T.resize(N+1);
        for(int i=0;i<=N;++i){
            double fract=0.5*(1-cos(PI*i/N));
            T[i]=XS+(XF-XS)*fract;
        }
    }

    void build_and_solve() {
        const double PI=3.1415926585;
        vector<vector<double>> A(N, vector<double>(N+1,0.0));
        for(int i=0;i<N;++i){
            double x1=0.5*(T[i]+T[i+1]);
            double y1=y_of_x(x1);
            double fac1=atan2(T[0]-x1,y1);
            for(int j=0;j<N;++j){
                double fac2=atan2(T[j+1]-x1,y1);
                A[i][j]=(fac2-fac1)/PI;
                fac1=fac2;
            }
            A[i][N]=1.0;
        }
        if(!gauss(A,N)) throw runtime_error("Singular system");
        M.resize(N);
        for(int i=0;i<N;++i) M[i]=A[i][N];
    }

    void press(double x, double& U, double& Cp) const {
        double yb=y_of_x(x);
        U=1.0; double V=0.0;
        double vf1=1.0/((T[0]-x)*(T[0]-x)+yb*yb);
        double uf1=(T[0]-x)*vf1;
        for(int j=0;j<N;++j){
            double dx=T[j+1]-x;
            double vf2=1.0/(dx*dx+yb*yb);
            double uf2=dx*vf2;
            U+=M[j]*(uf2-uf1);
            V-=M[j]*yb*(vf2-vf1);
            vf1=vf2; uf1=uf2;
        }
        Cp=1.0-U*U-V*V;
    }
};

int main() {
    cout << "\n=== PROGRAM DUBLET (Airfoil by Doublet Distribution) ===\n";
    DubletSolver S;
    cout << "N = "; cin >> S.N;
    while(true){
        cout << "XS, XF = "; cin >> S.XS >> S.XF;
        try{
            S.build_endpoints();
            S.build_and_solve();
        }catch(exception& e){
            cout << "Error: "<<e.what()<<"\n"; continue;
        }
        double U0,Cp0,U1,Cp1;
        S.press(0.0,U0,Cp0);
        S.press(1.0,U1,Cp1);
        cout<<"U at x=0: "<<U0<<"\n";
        cout<<"U at x=1: "<<U1<<"\n";
        cout<<"Accept results (Y/N)? "; string ans; cin>>ans;
        if(ans[0]=='Y'||ans[0]=='y') break;
    }
    const double PI=3.1415926585;
    cout<<"\nDOUBLET STRENGTH DISTRIBUTION\n";
    for(int i=0;i<=S.N;++i){
        double Mi=(i<S.N)?S.M[i]:0.0;
        cout<<setw(10)<<S.T[i]<<setw(10)<<Mi+PI<<"\n";
    }
    cout<<"\nBODY SHAPE\n";
    for(int i=0;i<S.N;++i){
        double xx=0.5*(S.T[i]+S.T[i+1]);
        cout<<setw(10)<<xx<<setw(10)<<y_of_x(xx)<<"\n";
    }
    cout<<"\nNPRINT = "; int NPRINT; cin>>NPRINT;
    if (NPRINT < 2) NPRINT = 2;
    for(int i=0;i<NPRINT;++i){
        double x=(double)i/(NPRINT-1);
        double U,Cp; S.press(x,U,Cp);
        cout<<setw(10)<<x<<setw(10)<<Cp<<"\n";
    }
    return 0;
}
