#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "Channels_KinRoutes.h"
#include "Math_Functions.h"
#include "Globals.h"
#include "Constitutive_Equations.h"
#include "Numerical_Flux.h"
#include "Infiltration.h"

/*****************************************************************************/
/***************************************************************************/
/************ Compute the spatial discretization of the equation ***********/
/***************************************************************************/


extern const double g;
extern void* xcalloc(int items, int size);
extern void GetLGLWeights(int N, double *w);


void computeChanL(struct channel* Chan, double time, int channelNumber, double* RHSA, double* RHSQ)
{
	int NumEdges = Chan->NumEdges;
	int NumEl = Chan->NumEl;
	int NumNodes = Chan->NumNodes;
	int Np = Chan->Np;

	/*
	if (time > 49067.2)
	{
		for (int i = 0; i < NumNodes+2; i++)
			printf("A = %lf \t Q = %lf \n", Chan->A[i], Chan->Q[i]);
	}
	*/
	#ifdef WDON
/**** Compute the mass over an element to use later for ensuring ********/
/**** positive mass with the flux in a wetting and drying treatment *************/
	double* mass = xcalloc(NumEl, sizeof(double));
	double LGLWeight[Np];
	GetLGLWeights(P, LGLWeight);

	for (int i = 0; i < NumEl; ++i)
	{
		double avgArea = 0;
		int begNode = i*Np + 1;
		for (int j =0; j < Np; j++)
		{
			avgArea += LGLWeight[j]*A[begNode+1];
		}
		mass[i] = avgArea;
	}
#endif
		
/************** Compute the numerical flux at the faces *****************/
	double* Fhat1L = xcalloc(NumEdges, sizeof(double));
	double* Fhat2L = xcalloc(NumEdges, sizeof(double));
	double* Fhat1R = xcalloc(NumEdges, sizeof(double));
	double* Fhat2R = xcalloc(NumEdges, sizeof(double));

	for (int i=0; i < NumEdges; ++i)
	{
		double A_L, Q_L, A_R, Q_R, m1val, m2val;

		int leftNode = i*Np;
		int rightNode = leftNode+1;
	
		m1val = Chan->m1[i];
		m2val = Chan->m2[i];
		double bval = Chan->b[i];

		double tmpF[2];
		
		#ifdef WDON
		// Check to see if the elements separated by this boundary are both dry
		if (i > 0 && i < NumEdges-1)
		{
			if (WD[i-1] == 0 && WD[i] == 0)
			{
				// Reflection flux for the left element
				A_L = Chan->A[leftNode];
				A_R = A_L;
				Q_L = Chan->Q[leftNode];
				Q_R = -Q_L;
				//LF(tmpF, A_L, A_R, Q_L, Q_R,bval, m1val, m2val, 0);
				RoeFluxChan(tmpF, A_L, A_R, Q_L, Q_R,bval, m1val, m2val,0);
				Fhat1L[i] = tmpF[0];
				Fhat2L[i] = tmpF[1];
				
				// Reflection flux for the right element
				A_R = Chan->A[rightNode];
				A_L = A_R;
				Q_R = Chan->Q[rightNode];
				Q_L = -Q_R;
				//LF(tmpF, A_L, A_R, Q_L, Q_R, bval, m1val, m2val, 0);
				RoeFluxChan(tmpF, A_L, A_R, Q_L, Q_R, bval, m1val, m2val, 0);
				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];

				if (isnan(tmpF[0]) || isnan(tmpF[1]))
				{
					printf("both elements dry flux not a number, edge %d time %lf \n",i, time);
					printf("A_L = %lf, A_R = %lf, Q_L = %lf, Q_R = %lf\n", A_L, A_R, Q_L, Q_R);
					exit(EXIT_FAILURE);
				}
				
			}

		}

		// if the elements are not both dry
		if ((i == 0) || (i == NumEl) || WD[i-1] == 1 || WD[i] == 1)
		{
			A_L = Chan->A[leftNode];
			A_R = Chan->A[rightNode];
	
			Q_L = Chan->Q[leftNode];
			Q_R = Chan->Q[rightNode];

			RoeFluxChan(tmpF, A_L, A_R, Q_L,Q_R,bval, m1val, m2val, g);
			//LF(tmpF, A_L, A_R, Q_L,Q_R,bval, m1val, m2val, g);
	
			if (isnan(tmpF[0]) || isnan(tmpF[1]))
			{
				printf(" both elements wet flux not a number, edge %d \n",i);
				exit(EXIT_FAILURE);
			}
			
			if (i==0 || WD[i-1] == 1)
			{
				Fhat1L[i] = tmpF[0];
				Fhat2L[i] = tmpF[1];	

			}
			else
			{
				double newtmpF[2];
				//LF(newtmpF, A_L, A_R, Q_L, Q_R, bval, m1val, m2val, 0);
				RoeFluxChan(newtmpF, A_L, A_R, Q_L, Q_R, bval, m1val, m2val, 0);
				Fhat1L[i] = newtmpF[0];
				Fhat2L[i] = newtmpF[1];
				if (isnan(newtmpF[0]) || isnan(newtmpF[1]))
				{
					printf("left element dry flux not a number, edge %d \n",i);
					exit(EXIT_FAILURE);
				}
	
			}
			if (i == NumEl || WD[i]==1)
			{
				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];
			}
			else
			{
				//LF(tmpF, A_L,A_R,Q_L,Q_R,b[i], m1val, m2val, 0);
				RoeFluxChan(tmpF, A_L,A_R,Q_L,Q_R,bval, m1val, m2val, 0);
				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];		
				if (isnan(tmpF[0]) || isnan(tmpF[1]))
				{
					printf("right element dry flux not a number, edge %d \n",i);
					exit(EXIT_FAILURE);
				}
	
			}
		}

		#else
		A_L = Chan->A[leftNode];
		A_R = Chan->A[rightNode];
	
		Q_L = Chan->Q[leftNode];
		Q_R = Chan->Q[rightNode];

		RoeFluxChan(tmpF, A_L, A_R, Q_L,Q_R,bval,m1val,m2val,g);
		//LF(tmpF, A_L, A_R, Q_L,Q_R,bval,m1val,m2val,g);
	
		if (isnan(tmpF[0]) || isnan(tmpF[1]))
		{
			printf("both elements dry flux not a number, edge %d time %lf \n",i, time);
			printf("A_L = %lf, A_R = %lf, Q_L = %lf, Q_R = %lf\n", A_L, A_R, Q_L, Q_R);

			
			//printf("A_L = %lf A_R = %lf Q_L = %lf Q_R = %lf, b = %lf, m1 = %lf, m2 = %lf \n", A_L, A_R, Q_L, Q_R, Chan->b[i], m1val, m2val);
			exit(EXIT_FAILURE);
		}
			
		Fhat1L[i] = tmpF[0];
		Fhat2L[i] = tmpF[1];	

		Fhat1R[i] = tmpF[0];
		Fhat2R[i] = tmpF[1];
		
		#endif

		/**************************************************************************************************/
		double u_L = Q_L/A_L;
		double u_R = Q_R/A_R;
		//approximate c_L (exact for rectangular channels, but not for trapezoidal
		double c_L = sqrt(g*A_L/Chan->b[i]);
		double c_R = sqrt(g*A_R/Chan->b[i]);

		if (i==0)
			Chan->max_lambda = fmax((fabs(u_L)+c_L), (fabs(u_R)+c_R));			
		// Compute maximum eigenvalue for the next time step
		double current_max =fmax((fabs(u_L) + c_L), (fabs(u_R) + c_R));
		Chan->max_lambda = max(Chan->max_lambda,current_max);
	}

	#ifdef WDON
	for (int i =1; i < NumEdges-1; ++i)
	{
		int leftNode = i*Np;
		int rightNode = leftNode+1;
		double maxBetaOverAlpha = 2;
		// Check to see if this flux might possibly result in negative mass
		// If the mass will be negative on the left side
	  if (Fhat1L[i]*maxBetaOverAlpha*dt > mass[i-1])
		{
	    double A_L = Chan->A[leftNode];
			double A_R = A_L;
			double Q_L = Chan->Q[leftNode];
	    double Q_R = -Q_L;
	    double tmpF[2];
			double localG = g;
				
			if (WD[i] == 0)
				localG = 0;
			
			//LF(tmpF, A_L, A_R, Q_L, Q_R, bval,m1val,m2val,localG);
			RoeFluxChan(tmpF, A_L, A_R, Q_L, Q_R, bval,m1val,m2val,localG);
	      		Fhat1L[i] = tmpF[0];
	      		Fhat2L[i] = tmpF[1];
			printf("element %d will be dry\n",i);
	  }
		
		// If the mass will be negative on the right side
		if (-Fhat1R[i]*maxBetaOverAlpha*dt > mass[i])
		{
			double A_R = Chan->A[rightNode];
			double A_L = A_R;
			double Q_R = Chan->Q[rightNode];
			double Q_L = -Q_R;
			double tmpF[2];
			double localG = g;
			if (WD[i] == 0)
				localG = 0;
			//LF(tmpF, A_L, A_R, Q_L, Q_R, bval,m1val,m2val,localG);
			RoeFluxChan(tmpF, A_L, A_R, Q_L, Q_R, bval, m1val, m2val,localG);
			Fhat1R[i] = tmpF[0];
			Fhat2R[i] = tmpF[1];
			printf("element %d will be dry \n",i);
				
		}
	}
	free(mass);
	#endif


	for (int k=0; k < NumEl; ++k) 
	{
		
		double h = Chan->dh[k];	
		gsl_vector *F1 = gsl_vector_alloc(Np);
		gsl_vector *F2 = gsl_vector_alloc(Np);
		gsl_vector *ST21 = gsl_vector_alloc(Np);
		gsl_vector *ST22 = gsl_vector_alloc(Np);
		gsl_vector *ST23 = gsl_vector_alloc(Np);
		
		//ST11 contains the qL term that comes from overland flow
		gsl_vector *ST11 = gsl_vector_calloc(Np);

		int begNode = k*Np;
		for (int i = 0; i < Np; i++)
		{
			double Aval = Chan->A[begNode+i+1];
			double Qval = Chan->Q[begNode+i+1];
			double qL = Chan->qL[begNode+i];
			double bval = Chan->NodalB[begNode+i];
			double S0 = Chan->dz[begNode+i];
			double m1val = Chan->Nodalm1[begNode+i];
			double m2val = Chan->Nodalm2[begNode+i];
			double dm1val = Chan->dm1[begNode+i];
			double dm2val = Chan->dm2[begNode+i];
			double dbval = Chan->db[begNode+i];
			double nval = Chan->NodalnFriction[begNode+i];
			double I1val = getI1(Aval, bval, m1val, m2val);
			double I2val = getI2(Aval, bval, dbval, m1val, dm1val, m2val, dm2val);
			double Sfval = getS_f(Aval, Qval, bval, m1val, m2val, nval);
			double beta = 1.0;
			double localG;
			#ifdef WDON
			if(WD[k]==1)
				localG = g;
			else
				localG = 0;
			#else
				localG = g;
			#endif
			
			gsl_vector_set(F1, i, Qval);
			gsl_vector_set(F2, i, beta*Qval*Qval/Aval + localG*I1val); 
			gsl_vector_set(ST21, i, localG*I2val);
			gsl_vector_set(ST22, i, localG*Aval*S0);
			gsl_vector_set(ST23, i, localG*Aval*Sfval);
			gsl_vector_set(ST11, i, qL);

		}

		//gsl_blas_dgemv(CblasNoTrans, 2.0/h, InvM, ST, 1.0, ST11);
		
		gsl_vector *localFhat1 = gsl_vector_alloc(2);
		gsl_vector_set(localFhat1, 0, Fhat1R[k]);
		gsl_vector_set(localFhat1, 1, Fhat1L[k+1]);

		gsl_vector *localFhat2 = gsl_vector_alloc(2);
		gsl_vector_set(localFhat2, 0, Fhat2R[k]);
		gsl_vector_set(localFhat2, 1, Fhat2L[k+1]);

		// cacluate the volume integral of the flux
		gsl_vector *localRHS1 = gsl_vector_calloc(Np);
		gsl_vector *localRHS2 = gsl_vector_calloc(Np);
		gsl_blas_dgemv(CblasNoTrans, 2.0/h, VolMat, F1, 1.0, localRHS1);
		gsl_blas_dgemv(CblasNoTrans, 2.0/h, VolMat, F2, 1.0, localRHS2);

		// calculate the surface integral of the flux
		gsl_vector *SurfPart1 = gsl_vector_calloc(Np);
		gsl_vector *SurfPart2 = gsl_vector_calloc(Np);
		gsl_blas_dgemv(CblasNoTrans, 2.0/h, LIFT, localFhat1, 1.0, SurfPart1);
		gsl_blas_dgemv(CblasNoTrans, 2.0/h, LIFT, localFhat2, 1.0, SurfPart2);

		// calculate the RHS
		gsl_vector_add(localRHS1, SurfPart1);
		gsl_vector_add(localRHS1, ST11);
		gsl_vector_add(localRHS2, SurfPart2);
		gsl_vector_add(localRHS2, ST21);
		gsl_vector_add(localRHS2, ST22);
		gsl_vector_sub(localRHS2, ST23);

		for(int i = 0; i < Np; i++)
		{
			RHSA[begNode+i+1] = gsl_vector_get(localRHS1, i);
/*			if (isnan(RHSA[begNode+i+1]))
			{
				printf("SurfPart\n");
				for(int k = 0; k < Np; k++)
				{
					printf("%lf \n", gsl_vector_get(SurfPart1, k));
	
				}

				printf("ST11\n");
				for(int k = 0; k < Np; k++)
				{
					printf("%lf \n", gsl_vector_get(ST11, k));
	
				}

			}
*/
			RHSQ[begNode+i+1] = gsl_vector_get(localRHS2, i);
		}

		gsl_vector_free(F1);
		gsl_vector_free(F2);
		gsl_vector_free(localFhat1);
		gsl_vector_free(localFhat2);
		gsl_vector_free(SurfPart1);
		gsl_vector_free(SurfPart2);
		gsl_vector_free(ST11);
		gsl_vector_free(ST21);
		gsl_vector_free(ST22);
		gsl_vector_free(ST23);
		gsl_vector_free(localRHS1);
		gsl_vector_free(localRHS2);


	}

	// The value at the ghost nodes will be the same as the values at the other side of the face

	RHSA[0] = RHSA[1];
	RHSA[NumNodes + 1] = RHSA[NumNodes];

	RHSQ[0] = RHSQ[1];
	RHSQ[NumNodes + 1] = RHSQ[NumNodes];

	free(Fhat1L);
	free(Fhat1R);
	free(Fhat2L);
	free(Fhat2R);

}


void computeKinL(struct kinematicRoute* kinRoute, double time, int kinRouteNumber, double dt, double* RHS)
{

	int Np = kinRoute->Np;
	int NumEdges =  kinRoute->NumEdges;
	int NumEl = kinRoute->NumEl;
	int NumNodes = kinRoute->NumNodes;

	double *Fhat = xcalloc(NumEdges,sizeof(double));
	
	/************ Compute the numerical flux at the edges *********************/
	for (int i =0; i < NumEdges; i++)
	{
			double H_L, H_R, S0_L, S0_R, nVal;
			int leftNode = i*Np;
			int rightNode = leftNode+1;
			H_L = kinRoute->H[leftNode];
			H_R = kinRoute->H[rightNode];
	
			//printf("nodeNum = %d\n", nodeNum);
			if (i < NumEdges-1)
				S0_R = kinRoute->dz[i*Np];
			else
				S0_R = kinRoute->dz[i*Np-1];
			
			S0_L = S0_R;
		
			nVal = kinRoute->nFriction[i];
			
			Fhat[i] = upwindKinRoute(H_L, H_R, S0_L, S0_R, nVal);
			//printf("Fhat = %lf S0 = %lf nf = %lf\n", Fhat[i], S0_L, nVal);
		
	}

	for (int k = 0; k < NumEl; k++)
	{
		double h_el = kinRoute->dh[k];
		gsl_vector *F = gsl_vector_alloc(Np); 
		gsl_vector *ST = gsl_vector_calloc(Np);
		gsl_vector *ST1 = gsl_vector_calloc(Np);
		
		// infiltration parameters for sandy clay loam
//		double K = 0.15*0.0328/3600;
//		double theta_e = 0.330;
//		double psi = 21.85*0.0328;
//		double s_e = 0.0;
//		double psi_delta_theta = psi*(1-s_e)*theta_e;
//
//		double cumulF = calculateInfiltration(psi, s_e, theta_e, K, time);
//		double f = K*(psi_delta_theta/cumulF + 1);
//
		for (int i =0; i < Np; i++)
		{
			double H = kinRoute->H[k*Np+i+1];
			double nval = kinRoute->NodalnFriction[k*Np+i];
			double S0 = kinRoute->dz[k*Np+i];
			gsl_vector_set(F,i,1./nval*sqrt(S0)*pow(H,5.0/3));
			
			// parking lot case
			if ( time < 1800)
				gsl_vector_set(ST, i, 2.0/12/3600);
				//gsl_vector_set(ST, i, 0.000014); 	
			//if (time <= 60.0)
			//	gsl_vector_set(ST, i, 0.001);
			
			//if (H > f || gsl_vector_get(ST, i) > f)
			//	gsl_vector_set(ST1, i, -f);
			//else 
			//	gsl_vector_set(ST1, i, -gsl_vector_get(ST, i));

		}
		

		double val1 = Fhat[k];
		double val2 = Fhat[k+1];
		gsl_vector *localFhat = gsl_vector_alloc(2);
		gsl_vector_set(localFhat, 0, val1);
		gsl_vector_set(localFhat, 1, val2);

		gsl_vector *localRHS = gsl_vector_calloc(Np);
		
		// calculate the volume integral
		gsl_blas_dgemv(CblasNoTrans, 2.0/h_el, VolMat, F, 1.0, localRHS);

		// calculate the surface integral
		gsl_vector *SurfPart = gsl_vector_calloc(Np);
		
		gsl_blas_dgemv(CblasNoTrans, 2.0/h_el, LIFT, localFhat, 1.0, SurfPart); 

		// subtract the surface integral from the volume integral
		gsl_vector_add(localRHS, SurfPart);

		// add source term to the right hand side
		gsl_vector_add(localRHS, ST);
		//gsl_vector_add(localRHS, ST1);
	
		for(int i =0; i < Np; i++)
		{
			RHS[k*Np+i+1] = gsl_vector_get(localRHS, i);
		}
	
		RHS[0] = RHS[1];
		RHS[NumNodes+1] = RHS[NumNodes];

		gsl_vector_free(F);
		gsl_vector_free(localFhat);
		gsl_vector_free(SurfPart);
		gsl_vector_free(localRHS);
		gsl_vector_free(ST);
		gsl_vector_free(ST1);

	}

	free(Fhat);
}


