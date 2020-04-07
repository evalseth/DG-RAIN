#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <nubspline.h>
#include "SimulationSteps.h"

void eval_NUBspline_1d_d_vg(NUBspline_1d_d * restrict spline, double x, 
	double* restrict val, double* restrict grad);

void eval_NUBspline_1d_d_vgl(NUBspline_1d_d * restrict spline, double x, 
	double* restrict val, double* restrict grad, double* restrict lapl);

void GetBathymetryData(double* x, double* z, int NumEdges)
{
	for(int i = 0; i < NumEdges; i++)
	{
		double currx = x[i];
		//z[i] = 0.01*currx;
		//z[i] = 0.01*currx*currx;
		z[i] = 0.25*log(currx+1);
		//z[i] = sin(currx);
	}

}

void InterpolateWithSplines(int NumPointsToInterpolate, int NumPointsToEvaluate, double* x, double* NodalX, double *z, double* S0, double* S0Deriv)
{
#ifndef NUBSPLINE_EVAL_STD_D_H
	printf("header file not included\n");
#endif

	double newX[NumPointsToInterpolate];
	double x_0 = 0;
	double h_el = 2*PI/(NumPointsToInterpolate-1);
	newX[0] = x_0;
	for (int i = 1; i < NumPointsToInterpolate; i++)
	{
		newX[i] = newX[i-1] + h_el;
	}

	NUgrid* myGrid = create_general_grid(newX, NumPointsToInterpolate);
	//NUgrid* myGrid = create_general_grid(x, NumPointsToInterpolate);
	BCtype_d myBC;
	myBC.lCode = NATURAL;
	myBC.rCode = NATURAL;
	//myBC.lCode = DERIV2;
	//myBC.rCode = DERIV2;
	//myBC.lVal = 0.25;
	//myBC.rVal = -0.25/pow(1+(2*PI),2);
	NUBspline_1d_d *mySplineZ = create_NUBspline_1d_d(myGrid, myBC, z);
	for(int i =0; i < NumPointsToEvaluate; i++)
	{
		double dummyz;
		eval_NUBspline_1d_d_vgl(mySplineZ, NodalX[i], &dummyz, &S0[i], &S0Deriv[i]);
		//eval_NUBspline_1d_d_vg(mySplineZ, NodalX[i], &dummyz, &S0[i]);
		//printf("z = %3.14f \t Real = %3.14f\n", dummyz, 0.25*log(1+NodalX[i]));
		//printf("S0 = %3.14f \t Real = %3.14f\n", S0[i], 0.25/(1+NodalX[i]));
		//printf("DerivS0 = %3.14f \t Real = %3.14f\n", S0Deriv[i], -0.25/(pow(1+NodalX[i],2)));
	}
//	exit(1);
	
	/*
	for (int i = 0; i < 99; i++)
	{
		double xval1 = 0 + i*(2*PI)/99;
		double dummyz;
		double S0Val;
		double DerivS0Val;
		eval_NUBspline_1d_d_vgl(mySplineZ, xval1, &dummyz, &S0Val, &DerivS0Val);  
		printf("%3.14f\n", DerivS0Val);
	}
	*/
//	exit(1);

}

/*int main()
{
	double x_0 = 0.0; 
	double x_1 = 2*PI;

	int NumEl = 50;
	int NumPoints = NumEl+1;
	double h_el = (x_1-x_0)/NumEl;

	double* x = calloc(NumPoints,sizeof(double));
	double* z = calloc(NumPoints,sizeof(double));
	double* S0 = calloc(NumPoints,sizeof(double));

	x[0] = x_0;
	x[NumEl] = x_1;
	for (int i = 1; i < NumEl; i++)
	{
		x[i] = x_0 + h_el*i;
	}

	GetBathymetryData(x, z, NumPoints);
	InterpolateWithSplines(NumPoints, x, z, S0); 
	
	for (int i = 0; i < NumPoints; i++)
	{
		printf("%3.14f \t %3.14f \t %3.14f \n", x[i], z[i], S0[i]);
	}
}
*/

