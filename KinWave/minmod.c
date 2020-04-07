#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globals.h"

/************************ Slope limiter reduces convergence rate to 1.5 ******************/


extern void* xcalloc(int items, int size);
extern void GetLGLWeights(int N, double *w);

// a function that returns -1 if a <0, 0 if a = 0 and 1 if a >1
int sign(double a)
{
	if (a < 0)
		return (-1);
	if (a == 0)
		return (0);
	else
		return (1);
}

void minmod()
{
	double* avgHt = xcalloc(NumEl, sizeof(double));
	
	double LGLWeight[Np];

	GetLGLWeights(P, LGLWeight);

	for (int k =0; k < NumEl; k++)
	{
		int nodeNum = EdgToNode[k]+1; 
		double avgVal = 0;
		for (int i =0; i < Np; i++)
		{
			avgVal += LGLWeight[i]*H[nodeNum+i];
		}
		avgHt[k] = 0.5*avgVal;
		//printf("avgHt = %3.13f\n", avgHt[k]);
	}

	double eps0 = 1e-8;
	
	for (int k =0; k < NumEl; k++)
	{
		int nodeNum = EdgToNode[k] + 1;
		double a, b, sigma, Htilda[2];
		if (k > 0)
			a = 2*(avgHt[k] - avgHt[k-1])/(dh[k]+dh[k-1]); 
		else
			a = 0;

		if (k < NumEl - 1)
			b = 2*(avgHt[k+1] - avgHt[k])/(dh[k] + dh[k+1]);
		
		else 
			b = 0;

		sigma = sign(a)*0.5*(1+sign(a*b))*fmin(fabs(a),fabs(b));
		
		Htilda[0] = avgHt[k] - 0.5*sigma*dh[k];
		Htilda[1] = avgHt[k] + 0.5*sigma*dh[k];

		if (fabs(Htilda[0] - H[nodeNum]) > eps0 || fabs(Htilda[1] - H[nodeNum+Np-1]) > eps0)
		{
			//printf("H0 = %3.13f \t H1 = %3.13f\n", H[nodeNum], H[nodeNum+Np-1]);
			H[nodeNum] = Htilda[0];
			for (int j =0; j < Np-1; j++)
			{
				H[nodeNum+j] = Htilda[0] + j*(Htilda[1] - Htilda[0])/(Np-1);
			}
			H[nodeNum+Np-1] = Htilda[1];
			//printf("Htilda0 = %3.13f\t Htilda1 = %3.13f\n", Htilda[0], Htilda[1]);
		}

	}

	free(avgHt);
//	exit(1);

}
