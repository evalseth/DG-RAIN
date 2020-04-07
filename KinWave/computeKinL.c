#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "globals.h"

extern double upwindKinWave(double H_L, double H_R, double S0_L, double S0_R, double n);
extern void *xcalloc(int items, int size);



void computeKinL(double t)
{
	double *Fhat = xcalloc(NumEdges,sizeof(double));

	/************ Compute the numerical flux at the edges *********************/
	for (int i =0; i < NumEdges; i++)
	{
			double H_L, H_R, S0_L, S0_R, nVal;
			int nodeNum = EdgToNode[i];
			H_L = H[nodeNum];
			H_R = H[nodeNum+1];
	
			//printf("nodeNum = %d\n", nodeNum);
			if (i < NumEdges-1)
				S0_L = S0[i*Np];
			else
				S0_L = S0[i*Np-1];
			S0_R = S0_L;
		
			//printf("S0 = %3.14f\n",S0_L);

		/*	if (i < NumEdges-1)
				S0_R = (z[i+1] - z[i])/dh[i];
			if (i > 0)
				S0_L = (z[i] - z[i-1])/dh[i-1];
			else
				S0_L = S0_R;
			*/

			nVal = n[i];
			
			Fhat[i] = upwindKinWave(H_L, H_R, S0_L, S0_R, nVal);
		
	}

	for (int k = 0; k < NumEl; k++)
	{
		double h_el = dh[k];
		gsl_vector *F = gsl_vector_alloc(Np); 
		//double S0Val = (z[k+1] - z[k])/h_el;
		double nval = n[k];					// for now assume n is constant. Will fix this later
		int begNode = EdgToNode[k]+1;
		for (int i =0; i < Np; i++)
		{
			//gsl_vector_set(F,i,1./nval*sqrt(S0)*H[begNode+i]);
			//gsl_vector_set(F,i,1./nval*sqrt(S0Val)*pow(H[begNode+i],5.0/3));
			gsl_vector_set(F,i,1./nval*sqrt(S0[k*Np+i])*pow(H[begNode+i],5.0/3));

		}

		double val1 = Fhat[k];
		double val2 = Fhat[k+1];
		gsl_vector *localFhat = gsl_vector_alloc(2);
		gsl_vector_set(localFhat, 0, val1);
		gsl_vector_set(localFhat, 1, val2);

		gsl_vector *localRHS = gsl_vector_calloc(Np);
		
		// calculate the volume integral
		gsl_blas_dgemv(CblasNoTrans, 2.0/h_el, VolMat, F, 1.0, localRHS);

/*		printf("**** Volume Integral **********\n");
		for(int i =0; i < Np; i++)
		{
				printf("%3.13f\n", gsl_vector_get(localRHS,i));
		}
*/
/*		printf("***** local Fhat *********\n");
		for (int i =0; i < Np; i++)
			printf("%3.13f\n", gsl_vector_get(localFhat,i));
*/

		// calculate the surface integral
		gsl_vector *SurfPart = gsl_vector_calloc(Np);
		
		gsl_blas_dgemv(CblasNoTrans, 2.0/h_el, LIFT, localFhat, 1.0, SurfPart); 

/*		printf("*************** LIFT **********\n");
			gsl_vector_set(F,i,1./nval*sqrt(S0[k*Np+i])*pow(H[begNode+i],5.0/3));
		for (int i = 0; i < Np; i++)
		{
			for(int j = 0; j < 2; j++)
				printf("%3.14f\t", gsl_matrix_get(LIFT,i,j));
			printf("\n");
		}
*/
/*		printf("***** Surface Integral *********\n");
		for(int i =0; i < Np; i++)
		{
				printf("%3.13f\n", gsl_vector_get(SurfPart, i));
		}

		printf("**************\n");
*/
		// subtract the surface integral from the volume integral
		gsl_vector_add(localRHS, SurfPart);
	
		// for manufactured solution
/*		gsl_vector *s = gsl_vector_alloc(Np);
		for (int i = 0; i < Np; i++)
		{
			double currx = NodalX[k*Np+i];
			double intermVal = cbrt(sin(currx)*cos(t));
			intermVal = pow(intermVal, 2.0);
			double derivS = S0Deriv[k*Np+i];
			double Sval = S0[k*Np+i];
			double sval = 0.5/(nval*sqrt(Sval))*derivS*pow(sin(currx)*cos(t)+2, 5.0/3) + 5*sqrt(Sval)/(3*nval)*pow(sin(currx)*cos(t)+2, 2.0/3)*cos(currx)*cos(t) - sin(currx)*sin(t);
			//double sval = -sin(currx)*sin(t) - 0.5/nval*pow(1+currx, -2)*pow(sin(currx)*cos(t)+2,5.0/3) + 0.5/nval*pow(1+currx,-1)*5.0/3*pow(sin(currx)*cos(t)+2,2.0/3)*cos(currx)*cos(t);
			//double sval = 5.0/(3*nval)*sqrt(S0)*pow(sin(currx)*cos(t)+2, 2.0/3)*cos(currx)*cos(t) - sin(currx)*sin(t);
			//double sval = 1.0/nval*sqrt(S0)*cos(currx)*cos(t) - sin(currx)*sin(t);
			gsl_vector_set(s,i,sval);
		}

	gsl_vector_add(localRHS, s);
	gsl_vector_free(s);
*/

/*	for (int i =0; i < Np; i++)
		printf("localRHS = %3.14f \t STerm = %3.14f\n",gsl_vector_get(localRHS,i), gsl_vector_get(s,i));
*/	

		gsl_vector_free(F);
		gsl_vector_free(localFhat);
		gsl_vector_free(SurfPart);

		for(int i =0; i < Np; i++)
		{
			KinRHS[begNode+i] = gsl_vector_get(localRHS, i);
			//printf("KinRHS = %3.14f\n", KinRHS[begNode+i]);
		}
		
		gsl_vector_free(localRHS);

	}

	free(Fhat);
//	exit(1);
}
