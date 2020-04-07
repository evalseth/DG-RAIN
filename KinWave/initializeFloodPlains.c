#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "globals.h"

double dt;
double *KinRHS;
double *H;

void initializeFloodPlains()
{
	H = malloc(TotalNumNodes*sizeof(double));
	KinRHS = malloc(TotalNumNodes*sizeof(double));
	for (int i =0; i < NumEl; i++)
	{
		//H[i] =0.5;
		for (int j =0; j < Np; j++)
		{
			double currx = NodalX[i*Np+j];
			H[i*Np+j+1] = sin(currx)+2;
			KinRHS[i*Np+j+1] = 0.0;
		}
	}
	H[0] = H[1];
	H[TotalNumNodes-1] = H[TotalNumNodes-2];
	KinRHS[0] = 0.0;
	KinRHS[TotalNumNodes-1] = KinRHS[TotalNumNodes-2];
}

