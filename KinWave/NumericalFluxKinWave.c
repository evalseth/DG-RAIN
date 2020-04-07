#include <stdio.h>
#include <math.h>
#include "globals.h"

extern int sign(double x);
double upwindKinWave(double H_L, double H_R, double S0_L, double S0_R, double n)
{
	// upwind flux 
	
	// /*** Pure  advection ****/
	// double F_L = 1.0/n*sqrt(S0_L)*H_L;
	// double F_R = 1.0/n*sqrt(S0_L)*H_R;
	// /*******************************/
	

	double F_L = 1.0/n*sqrt(S0_L)*pow(H_L,5.0/3);
	double F_R = 1.0/n*sqrt(S0_R)*pow(H_R,5.0/3);

	double s = (F_L - F_R)/(H_L-H_R);
	double Fhat = 0.5*(F_L + F_R) + 0.5*sign(s)*(F_L - F_R);

	//printf("H_L = %3.14f \t H_R = %3.14f \n", H_L, H_R);
	//printf("Fhat = %3.14f\n", Fhat);

	return Fhat;

}
