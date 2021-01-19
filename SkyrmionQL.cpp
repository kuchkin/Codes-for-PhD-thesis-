#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string> 
#include <sstream>
#include <cmath>
#include <cstring>

using namespace std;

const double pi = 3.1415926535897932384626433832795;

const int N = 1e6;
double L = 8.0; //Box size in units of Ld
double dr = L/(N+1.);
double h = 0.65; //External magnetic field in units of Bd
double u = 0.0; //Anisotropy in units of Ld and Bd
double To = pi, Tn = 0.; // skyrmion boundary conditions
double tolerance = 1e-10, dif = 1;


double *T = new double [N], *T0 = new double [N], *a = new double [N], *b = new double [N];


void CalculateF1(double *a){
	for(int i = 0; i < N; i++){
		double r = (i+1)*dr;
		a[i] = -(cos(2*a[i])/r - 4*pi*sin(2*a[i]) + 4*pi*pi*r*(h*cos(a[i]) + 2*u*cos(2*a[i])));
	}
}

void CalculateF2(double *b){
	for(int i = 0; i < N; i++){
		double r = (i+1)*dr;
		b[i] = 0.5*sin(2*b[i])/r - 4*pi*sin(b[i])*sin(b[i]) + \
		4*pi*pi*r*(h*sin(b[i]) + u*sin(2*b[i])) + a[i]*b[i];
	}
}

double Energy(double *t){
	double E = 0;
	E += 0.5*dr*((0.5*(t[1]-To)/dr)*(0.5*(t[1]-To)/dr)+sin(t[0])*sin(t[0])/(dr*dr)) + \
		2*pi*(dr*0.5*(t[1]-To)/dr) + \
		4*pi*pi*dr*(2*h*sin(t[0]/2)*sin(t[0]/2) + u*sin(t[0])*sin(t[0]));

	for(int i = 1; i < N-1; i++){
		double r = (i+1)*dr;
		E += 0.5*r*((0.5*(t[i+1]-t[i-1])/dr)*(0.5*(t[i+1]-t[i-1])/dr)+sin(t[i])*sin(t[i])/(r*r)) + \
		2*pi*(r*0.5*(t[i+1]-t[i-1])/dr + 0.5*sin(2*t[i])) + \
		4*pi*pi*r*(2*h*sin(t[i]/2)*sin(t[i]/2) + u*sin(t[i])*sin(t[i]));
	}

	double r = N*dr;
	E += 0.5*r*((0.5*(Tn-t[N-2])/dr)*(0.5*(Tn-t[N-2])/dr)+sin(t[N-1])*sin(t[N-1])/(r*r)) + \
		2*pi*(r*0.5*(Tn-t[N-2])/dr + 0.5*sin(2*t[N-1])) + \
		4*pi*pi*r*(2*h*sin(t[N-1]/2)*sin(t[N-1]/2) + u*sin(t[N-1])*sin(t[N-1]));
	E = 2*pi*E*dr;
	return(E);
}

void SaveTheta(double *t){
	printf("Saving data...\n");
	fstream o;
	o.open("T.csv", ios::out);
	o<<0<<","<<To<<std::endl;
	for(int i = 0; i < N; i++){
		double r = (i+1)*dr;
		o<<r<<","<<t[i]<<std::endl;
	}
	o<<L<<","<<Tn<<std::endl;
	o.close();
	printf("Data are saved to T.csv\n");
}

double Difference(double *t0, double *t){
	double res = 0.0;
	for(int i = 0; i < N; i++)
		res += (t0[i]-t[i])*(t0[i]-t[i])*dr;
	return res;
}

void Solve(double *a, double *b, double *t){
	for (int i = 0; i < N; i++){
		double r = (i+1)*dr;
		a[i] -= 2*r/(dr*dr);
	}
	b[0] -= To*0.5/dr;
	b[N-1] -= Tn*1.5/dr;

	for(int i = 1; i < N; i++){
		double r = (i+1)*dr;
		a[i] = a[i] - (r/dr-0.5)*(r/dr-0.5)/(dr*dr*a[i-1]);
		b[i] = b[i] - b[i-1]*(r/dr-0.5)/(dr*a[i-1]);
	}
	t[N-1] = b[N-1]/a[N-1];
	for(int i = N-2; i > -1; i--){
		double r = (i+1)*dr;
		t[i] = (b[i]-(r/dr+0.5)*t[i+1]/dr)/a[i];
	}
	for(int i = 0; i < N; i++){
		a[i] = t[i];
		b[i] = t[i];
	}

}

int main ()
{
	for (int i = 0; i<N; i++){
		double r = (i+1)*dr;
		T0[i] = 2.*atan(exp(-2*pi*sqrt(h+2*u)*r)/r);
		a[i] = T0[i]; b[i] = T0[i];
	}
	double Eprev = Energy(T0), Enow;
	printf("Energy of the inital state = %.15f\n",Eprev);

	while(dif>tolerance){
		CalculateF1(a);
		CalculateF2(b);
		Solve(a,b,T);
		Enow = Energy(T);		
		if(Enow>Eprev){
			printf("Can not reach the tolernace %0.1e, try to increase discretization or improve function T0!\n",tolerance);
			dif = tolerance/10;
		}
		else{
			Eprev = Enow;
			dif = Difference(T0,T);
			for(int i = 0; i<N; i++)
				T0[i] = T[i];
			printf("Energy = %0.15f, Error = %0.15f\n",Enow,dif);
		}	
	}
	printf("Done!\n");
	SaveTheta(T);

	delete T0;
	delete T;
	delete a;
	delete b;

	return 0;
}
