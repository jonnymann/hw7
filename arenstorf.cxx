#include <iostream>
#include <math>


using namespace std;

void f(double* x, double* y, double* kx, double* ky, double mu){
    double r = sqrt(pow(x[0]+mu,2)+pow(y[0],2));
    double s = sqrt(pow(x[0]+1-mu,2)+pow(y[0],2));
    kx[0]=x[1];
    kx[1]=x[0]-2*y[1]-(1-mu)*(x[0]+mu)/pow(r,3)-mu*(x[0]-1+mu)/pow(s,3);
    ky[0]=y[1];
    ky[1]=y[0]-2*x[1]-(1-mu)*y[0]/pow(y[0],3)-mu*y[0]/pow(s,3);
    
}

int main(){
    const int N = 1000;
    double step = 1e-3;
    double mu = 0.0012277471;
    const int dim = 2;
    
    double kx1[dim],kx2[dim],kx3[dim],kx4[dim],kx5[dim],kx6[dim],kx7[dim];
    double ky1[dim],ky2[dim],ky3[dim],ky4[dim],ky5[dim],ky6[dim],ky7[dim];
    double x0[dim],y0[dim],xn4[dim],yn4[dim],xn5[dim],yn5[dim]; 

    /* Erde-Mond+leichte Masse:
        x[0]=0.994;
        x[1]=0;
        y[0]=0;
        y[1]=-2.00158510637908;
    // (T=17.065216560157)
        */
    
    
        f(x0,y0,kx1,ky1,mu);
        xn4[0]=x0[0]+0.2*step*kx1[0];
        xn4[1]=x0[1]+0.2*step*kx1[1];
        yn4[0]=y0[0]+0.2*step*ky1[0];
        yn4[1]=y0[1]+0.2*step*ky1[1];
    
        f(xn4,yn4,kx2,ky2,mu);
        xn4[0]=x0[0]+3/40*step*kx1[0]+9/40*step*kx2[0];
        xn4[1]=x0[1]+3/40*step*kx1[1]+9/40*step*kx2[1];
        yn4[0]=y0[0]+3/40*step*ky1[0]+9/40*step*ky2[0];
        yn4[1]=y0[1]+3/40*step*ky1[1]+9/40*step*ky2[1];

        f(xn4,yn4,kx3,ky3,mu);
        xn4[0]=x0[0]+44/45*step*kx1[0]-56/15*step*kx2[0]+32/9*step*kx3[0];
        xn4[1]=x0[1]+44/45*step*kx1[1]-56/15*step*kx2[1]+32/9*step*kx3[1];
        yn4[0]=y0[0]+44/45*step*ky1[0]-56/15*step*ky2[0]+32/9*step*ky3[0];
        yn4[1]=y0[1]+44/45*step*ky1[1]-56/15*step*ky2[1]+32/9*step*ky3[1];

        f(xn4,yn4,kx4,ky4,mu);
        xn4[0]=x0[0]+19372/6561*step*kx1[0]-25360/2187*step*kx2[0]+64448/6561*step*kx3[0]-212/729*step*kx4[0];
        xn4[1]=x0[1]+19372/6561*step*kx1[1]-25360/2187*step*kx2[1]+64448/6561*step*kx3[1]-212/729*step*kx4[1];
        yn4[0]=y0[0]+19372/6561*step*ky1[0]-25360/2187*step*ky2[0]+64448/6561*step*ky3[0]-212/729*step*ky4[0];
        yn4[1]=y0[1]+19372/6561*step*ky1[1]-25360/2187*step*ky2[1]+64448/6561*step*ky3[1]-212/729*step*ky4[1];

        f(xn4,yn4,kx5,ky5,mu);
        xn4[0]=x0[0]+9017/3168*step*kx1[0]-355/33*step*kx2[0]+46732/5247*step*kx3[0]-49/176*step*kx4[0]-5103/18656*step*kx5[0];
        xn4[1]=x0[1]+9017/3168*step*kx1[1]-355/33*step*kx2[1]+46732/5247*step*kx3[1]-49/176*step*kx4[1]-5103/18656*step*kx5[1];
        yn4[0]=y0[0]+9017/3168*step*ky1[0]-355/33*step*ky2[0]+46732/5247*step*ky3[0]-49/176*step*ky4[0]-5103/18656*step*ky5[0];
        yn4[1]=y0[1]+9017/3168*step*ky1[1]-355/33*step*ky2[1]+46732/5247*step*ky3[1]-49/176*step*ky4[1]-5103/18656*step*ky5[1];

        f(xn4,yn4,kx6,ky6,mu);
        xn4[0]=x0[0]+35/384*step*kx1[0]+500/1113*step*kx3[0]-125/192*step*kx4[0]-2187/6784*step*kx5[0]+11/84*step*kx6[0];
        xn4[1]=x0[1]+35/384*step*kx1[1]+500/1113*step*kx3[1]-125/192*step*kx4[1]-2187/6784*step*kx5[1]+11/84*step*kx6[1];
        yn4[0]=y0[0]+35/384*step*ky1[0]+500/1113*step*ky3[0]-125/192*step*ky4[0]-2187/6784*step*ky5[0]+11/84*step*ky6[0];
        yn4[1]=y0[1]+35/384*step*ky1[1]+500/1113*step*ky3[1]-125/192*step*ky4[1]-2187/6784*step*ky5[1]+11/84*step*ky6[1];        
    
        f(xn4,yn4,kx7,ky7,mu);
        xn5[0]=x0[0]+5179/57600*step*kx1[0]+7571/16695*step*kx3[0]-393/640*step*kx4[0]-92097/339200*step*kx5[0]+187/2100*step*kx6[0]+1/40*step*kx7[0];
        xn5[1]=x0[1]+5179/57600*step*kx1[1]+7571/16695*step*kx3[1]-393/640*step*kx4[1]-92097/339200*step*kx5[1]+187/2100*step*kx6[1]+1/40*step*kx7[1];
        yn5[0]=y0[0]+5179/57600*step*ky1[0]+7571/16695*step*ky3[0]-393/640*step*ky4[0]-92097/339200*step*ky5[0]+187/2100*step*ky6[0]+1/40*step*ky7[0];
        yn5[1]=y0[1]+5179/57600*step*ky1[1]+7571/16695*step*ky3[1]-393/640*step*ky4[1]-92097/339200*step*ky5[1]+187/2100*step*ky6[1]+1/40*step*ky7[1];  
        
        return 0;
}