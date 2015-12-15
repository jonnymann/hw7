#include <iostream>
#include <cmath>


using namespace std;

void f(double* y,double* k, double mu){
    double r = sqrt(pow(y[0]+mu,2)+pow(y[2],2));
    double s = sqrt(pow(y[0]-1+mu,2)+pow(y[2],2));
    k[0]=y[1];
    k[1]=y[0]+2.*y[3]-(1.-mu)*(y[0]+mu)/pow(r,3)-mu*(y[0]-1.+mu)/pow(s,3);
    k[2]=y[3];
    k[3]=y[2]-2.*y[1]-(1.-mu)*y[2]/pow(r,3)-mu*y[2]/pow(s,3);
    
}

void RK4Step(double* y0, double* yn4, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double mu, double step, int dim){
	
	f(y0,k1,mu);
	for(int i=0;i<dim;i++)
		yn4[i]=y0[i]+0.2*step*k1[i];
    
    f(yn4,k2,mu);
	for(int i=0;i<dim;i++)
		yn4[i]=y0[i]+3./40*step*k1[i]+9./40*step*k2[i];
       
    f(yn4,k3,mu);
	for(int i=0;i<dim;i++)
		yn4[i]=y0[i]+44./45*step*k1[i]-56./15*step*k2[i]+32./9*step*k3[i];

    f(yn4,k4,mu);
	for(int i=0;i<dim;i++)
		yn4[i]=y0[i]+19372./6561*step*k1[i]-25360./2187*step*k2[i]+64448./6561*step*k3[i]-212./729*step*k4[i];

    f(yn4,k5,mu);
	for(int i=0;i<dim;i++)
		yn4[i]=y0[i]+9017./3168*step*k1[i]-355./33*step*k2[i]+46732./5247*step*k3[i]+49./176*step*k4[i]-5103./18656*step*k5[i];

    f(yn4,k6,mu);
	for(int i=0;i<dim;i++)
		yn4[i]=y0[i]+35./384*step*k1[i]+500./1113*step*k3[i]+125./192*step*k4[i]-2187./6784*step*k5[i]+11./84*step*k6[i];
}



int main(){
    const int N = 1000;
    double step = 1e-5;
    const int dim = 4;
	double t=0, dis;

    double k1[dim],k2[dim],k3[dim],k4[dim],k5[dim],k6[dim],k7[dim];
    double y0[dim],yn4[dim],yn5[dim]; 

    double mu = 0.012277471;
	double T0=17.065216560157;
	double eps = 1e-5;
	y0[0]=0.994;
    y0[1]=0;
    y0[2]=0;
    y0[3]=-2.00158510637908;
	while (t <T0)
	{
		RK4Step(y0,yn4,k1,k2,k3,k4,k5,k6,mu,step,dim);		// Runge Kutta Stufe 4

		f(yn4,k7,mu);
		for(int i=0;i<dim;i++)								// Runge Kutta Stufe 5
			yn5[i]=y0[i]+5179./57600*step*k1[i]+7571./16695*step*k3[i]+393./640*step*k4[i]-92097./339200*step*k5[i]+187./2100*step*k6[i]+1./40*step*k7[i];  
		
		dis=abs(yn5[0]-yn4[0]);
		for(int i=1; i<dim;i++){
			if(abs(yn5[i]-yn4[i])>dis)
				dis=abs(yn5[i]-yn4[i]);							// Berechne Distanz
		}

		step=step*pow(eps/dis,0.2);							// Neue Schrittweite

		RK4Step(y0,yn4,k1,k2,k3,k4,k5,k6,mu,step,dim);		// Runge Kutta 4 mit neuer Schrittweite
		cout << t << '\t' << step << '\t' << yn4[0] << '\t' << yn4[2] << endl;					// Ausgabe

		for(int i=0; i<dim; i++)							// Vorbereitung auf nächsten Schleifendurchgang
			y0[i]=yn4[i];
		t=t+step;
	}
	return 0;
}