#include <stdio.h>
#include <math.h>
#include "Math_Functions.h"
#include "Constitutive_Equations.h"
#include "Globals.h"
#include "Channels_KinRoutes.h"

extern void GetLGLWeights(int N, double* w);
extern void* xcalloc(int numItems, int size);

void minmodChan(struct channel* Chan)
{
	int NumEl = Chan->NumEl;
	double *dh = Chan->dh;
	int Np = Chan->Np;
	int P = Chan->P;
	int NumNodes = Chan->NumNodes;

	// average value of Zeta and Q
	double* avgZeta = xcalloc(NumEl, sizeof(double));
	double* avgQ = xcalloc(NumEl, sizeof(double));
	double** Zeta = xcalloc(NumEl, sizeof(double*));
	
	for (int k = 0; k < NumEl; k++)
		Zeta[k] = xcalloc(Np, sizeof(double));

	double LGLWeight[Np];
	GetLGLWeights(P, LGLWeight);

	for (int i=0; i < NumEl; ++i)
	{
		double avgZetaVal = 0; 
		double avgQVal = 0;
		for (int j = 0; j < Np; j++)
		{
			int nodeNum = i*Np+j;
			double Aval = Chan->A[nodeNum+1];
			double bval = Chan->NodalB[nodeNum];
			double zval = Chan->NodalZ[nodeNum];
			double m1val = Chan->Nodalm1[nodeNum];
			double m2val = Chan->Nodalm2[nodeNum];
			double Hval = getH(Aval, bval, m1val, m2val);
			double Qval = Chan->Q[nodeNum+1];
			double zetaVal = Hval - zval;
			Zeta[i][j] = zetaVal;

			avgZetaVal += LGLWeight[j]*zetaVal;
			avgQVal += LGLWeight[j]*Qval;
		}
		avgZeta[i] = 0.5*avgZetaVal;
		avgQ[i] = 0.5*avgQVal;

	}

	double eps0 = 1e-8;
	for (int i =0; i < NumEl; ++i)
	{	
		double a1, a2, b1, b2, sigma1, sigma2;
		double Zetatilda[2], Qtilda[2];
		if (i > 0)
		{
			a1 = 2*(avgZeta[i] - avgZeta[i-1])/(dh[i] + dh[i-1]);
			a2 = 2*(avgQ[i] - avgQ[i-1])/(dh[i]+dh[i-1]);
		}
		else
		{	
			a1 = 0;
			a2 = 0;
		}
		if (i < NumEl-1)
		{	
			b1 = 2*(avgZeta[i+1] - avgZeta[i])/(dh[i]+dh[i+1]);
			b2 = 2*(avgQ[i+1] - avgQ[i])/(dh[i]+dh[i+1]);
		}
		else
		{
			b1 = 0; 
			b2 = 0;
		}

		sigma1 = sign(a1)*0.5*(1+sign(a1*b1))*fmin(fabs(a1),fabs(b1));
		sigma2 = sign(a2)*0.5*(1+sign(a2*b2))*fmin(fabs(a2),fabs(b2));

		Zetatilda[0] = avgZeta[i] - 0.5*sigma1*dh[i];
		Zetatilda[1] = avgZeta[i] + 0.5*sigma1*dh[i];
		Qtilda[0] = avgQ[i] - 0.5*sigma2*dh[i];
		Qtilda[1] = avgQ[i] + 0.5*sigma2*dh[i];

		// limit only if limiting is required
		if (fabs(Zetatilda[0]-Zeta[i][0]) > eps0 || fabs(Zetatilda[1] - Zeta[i][Np-1]) > eps0)
		{
			for (int j = 0; j < Np; j++)
			{
				Zeta[i][j] = Zetatilda[0] + j*(Zetatilda[1] - Zetatilda[0])/(Np-1);
				double hval = Zeta[i][j] + Chan->NodalZ[i*Np+j];
				double bval = Chan->NodalB[i*Np+j];

				double m1val = Chan->Nodalm1[i*Np+j];
				double m2val = Chan->Nodalm2[i*Np+j];
				Chan->A[i*Np+j+1] = hval*bval + 0.5*m1val*hval*hval+ 0.5*m2val*hval*hval;
			}
		}
	
		// limit only if limiting is required
		if (fabs(Qtilda[0]-Chan->Q[i*Np+1]) > eps0 || fabs(Qtilda[1] - Chan->Q[i*Np+Np]) > eps0)
		{
			for (int j = 0; j < Np; j++)
			{
				Chan->Q[i*Np+j+1] = Qtilda[0] + j*(Qtilda[1] - Qtilda[0])/(Np-1);
			}
			
		}

	}


/*	for (int i = 0; i < NumNodes; ++i)
	{
		int modify;
		#ifdef WDON
		modify = (height[i][0] > H0)*(height[i][1]>H0)*((Zeta[i][0]+zEl[i][0]) > 0)*((Zeta[i][1]+zEl[i][1])>0);
		#else
		modify = 1;
		#endif

		for (int k =0; k < 2; ++k)
		{
			if (modify)
			{
				double zeta = Zeta[i][k];
				height[i][k] = zeta + zEl[i][k];
				A[i*2+k+1] = height[i][k] * bEl[i][k];
				Q[i*2+k+1] = QEl[i][k];
			}
		}
	}
*/
	Chan->A[0] = Chan->A[1];
	Chan->Q[0] = Chan->Q[1];
	Chan->A[NumNodes-1] = Chan->A[NumNodes-2];
	Chan->Q[NumNodes-1] = Chan->Q[NumNodes-2];

	for (int k = 0; k < NumEl; k++)
		free(Zeta[k]);

	free(Zeta);
	free(avgZeta);
	free(avgQ);
}

void minmodHeightChan(struct channel* Chan)
{
	int NumEl = Chan->NumEl;
	double *dh = Chan->dh;
	int Np = Chan->Np;
	int P = Chan->P;
	int NumNodes = Chan->NumNodes;

	// average value of Zeta and Q
	double* avgHeight = xcalloc(NumEl, sizeof(double));
	double* avgQ = xcalloc(NumEl, sizeof(double));
	double** Height = xcalloc(NumEl, sizeof(double*));
	
	for (int k = 0; k < NumEl; k++)
		Height[k] = xcalloc(Np, sizeof(double));

	double LGLWeight[Np];
	GetLGLWeights(P, LGLWeight);

	for (int i=0; i < NumEl; ++i)
	{
		double avgHeightVal = 0; 
		double avgQVal = 0;
		for (int j = 0; j < Np; j++)
		{
			int nodeNum = i*Np+j;
			double Aval = Chan->A[nodeNum+1];
			double bval = Chan->NodalB[nodeNum];
			double zval = Chan->NodalZ[nodeNum];
			double m1val = Chan->Nodalm1[nodeNum];
			double m2val = Chan->Nodalm2[nodeNum];
			double Hval = getH(Aval, bval, m1val, m2val);
			double Qval = Chan->Q[nodeNum+1];
			Height[i][j] = Hval;

			avgHeightVal += LGLWeight[j]*Hval;
			avgQVal += LGLWeight[j]*Qval;
		}
		avgHeight[i] = 0.5*avgHeightVal;
		avgQ[i] = 0.5*avgQVal;

	}

	double eps0 = 1e-8;
	for (int i =0; i < NumEl; ++i)
	{	
		double a1, a2, b1, b2, sigma1, sigma2;
		double Htilda[2], Qtilda[2];
		if (i > 0)
		{
			a1 = 2*(avgHeight[i] - avgHeight[i-1])/(dh[i] + dh[i-1]);
			a2 = 2*(avgQ[i] - avgQ[i-1])/(dh[i]+dh[i-1]);
		}
		else
		{	
			a1 = 0;
			a2 = 0;
		}
		if (i < NumEl-1)
		{	
			b1 = 2*(avgHeight[i+1] - avgHeight[i])/(dh[i]+dh[i+1]);
			b2 = 2*(avgQ[i+1] - avgQ[i])/(dh[i]+dh[i+1]);
		}
		else
		{
			b1 = 0; 
			b2 = 0;
		}

		sigma1 = sign(a1)*0.5*(1+sign(a1*b1))*fmin(fabs(a1),fabs(b1));
		sigma2 = sign(a2)*0.5*(1+sign(a2*b2))*fmin(fabs(a2),fabs(b2));

		Htilda[0] = avgHeight[i] - 0.5*sigma1*dh[i];
		Htilda[1] = avgHeight[i] + 0.5*sigma1*dh[i];
		Qtilda[0] = avgQ[i] - 0.5*sigma2*dh[i];
		Qtilda[1] = avgQ[i] + 0.5*sigma2*dh[i];

		// limit only if limiting is required
		if (fabs(Htilda[0]-Height[i][0]) > eps0 || fabs(Htilda[1] - Height[i][Np-1]) > eps0)
		{
			for (int j = 0; j < Np; j++)
			{
				double hval = Htilda[0] + j*(Htilda[1] - Htilda[0])/(Np-1);
				Height[i][j] = hval;
				double bval = Chan->NodalB[i*Np+j];

				double m1val = Chan->Nodalm1[i*Np+j];
				double m2val = Chan->Nodalm2[i*Np+j];
				Chan->A[i*Np+j+1] = hval*bval + 0.5*m1val*hval*hval+ 0.5*m2val*hval*hval;
			}
		}
	
		// limit only if limiting is required
		if (fabs(Qtilda[0]-Chan->Q[i*Np+1]) > eps0 || fabs(Qtilda[1] - Chan->Q[i*Np+Np]) > eps0)
		{
			for (int j = 0; j < Np; j++)
			{
				Chan->Q[i*Np+j+1] = Qtilda[0] + j*(Qtilda[1] - Qtilda[0])/(Np-1);
			}
			
		}

	}


/*	for (int i = 0; i < NumNodes; ++i)
	{
		int modify;
		#ifdef WDON
		modify = (height[i][0] > H0)*(height[i][1]>H0)*((Zeta[i][0]+zEl[i][0]) > 0)*((Zeta[i][1]+zEl[i][1])>0);
		#else
		modify = 1;
		#endif

		for (int k =0; k < 2; ++k)
		{
			if (modify)
			{
				double zeta = Zeta[i][k];
				height[i][k] = zeta + zEl[i][k];
				A[i*2+k+1] = height[i][k] * bEl[i][k];
				Q[i*2+k+1] = QEl[i][k];
			}
		}
	}
*/
	Chan->A[0] = Chan->A[1];
	Chan->Q[0] = Chan->Q[1];
	Chan->A[NumNodes-1] = Chan->A[NumNodes-2];
	Chan->Q[NumNodes-1] = Chan->Q[NumNodes-2];

	for (int k = 0; k < NumEl; k++)
		free(Height[k]);

	free(Height);
	free(avgHeight);
	free(avgQ);
}



/************************ Slope limiter reduces convergence rate to 1.5 ******************/


void minmodKinRoute(struct kinematicRoute* kinRoute)
{
	int NumEl = kinRoute->NumEl;
	int NumNodes = kinRoute->NumNodes;
	int P = kinRoute->P;
	int Np = kinRoute->Np;

	double* avgHt = xcalloc(NumEl, sizeof(double));

	double *dh = kinRoute->dh;

	double LGLWeight[Np];

	GetLGLWeights(P, LGLWeight);

	for (int k =0; k < NumEl; k++)
	{
		int nodeNum = k*Np+1; 
		double avgVal = 0;
		for (int i =0; i < Np; i++)
		{
			avgVal += LGLWeight[i]*kinRoute->H[nodeNum+i];
		}
		avgHt[k] = 0.5*avgVal;
		//printf("avgHt = %3.13f\n", avgHt[k]);
	}

	double eps0 = 1e-8;
	
	for (int k =0; k < NumEl; k++)
	{
		int nodeNum = k*Np + 1;
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

		if (fabs(Htilda[0] - kinRoute->H[nodeNum]) > eps0 || fabs(Htilda[1] - kinRoute->H[nodeNum+Np-1]) > eps0)
		{
			//printf("H0 = %3.13f \t H1 = %3.13f\n", H[nodeNum], H[nodeNum+Np-1]);
			kinRoute->H[nodeNum] = Htilda[0];
			for (int j =0; j < Np-1; j++)
			{
				kinRoute->H[nodeNum+j] = Htilda[0] + j*(Htilda[1] - Htilda[0])/(Np-1);
			}
			kinRoute->H[nodeNum+Np-1] = Htilda[1];
			//printf("Htilda0 = %3.13f\t Htilda1 = %3.13f\n", Htilda[0], Htilda[1]);
		}

	}

	free(avgHt);

	kinRoute->H[0] = kinRoute->H[1];
	kinRoute->H[NumNodes-1] = kinRoute->H[NumNodes-2];

}
