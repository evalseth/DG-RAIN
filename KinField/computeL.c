#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "Watershed.h"
#include "Math_Functions.h"
#include "Globals.h"
#include "Constitutive_Equations.h"
#include "Numerical_Flux.h"
#include "ManufacturedSolution.h"
#include "Infiltration.h"
/*****************************************************************************************//**
 * @file computeL.c
 *
 * This file contains code to evaluate the right hand side of the following discrete ODE obtained
 * from the DG discretization of the 2-D shallow water equations:
 * 
 * \f$\frac{\partial \hat{\mathbf{w}}}{\partial t} = M^{-1}L(\hat{\mathbf{w}},t)\f$
 *
 * *******************************************************************************************/


/***************************************************************************/
/************ Compute the spatial discretization of the equation ***********/
/***************************************************************************/
extern void GetLGLWeights(int N, double *w);

void computeChanL(struct channel* Chan, double time, double *dt, int channelNumber, int stage, double* RHSA, double* RHSQ)
{
	int fluxType = 1;
	int NumEdges = Chan->NumEdges;
	int NumEl = Chan->NumEl;
	int NumNodes = Chan->NumNodes;
	int Np = Chan->Np;
	int P = Chan->P;

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
			avgArea += LGLWeight[j]*Chan->A[begNode+j];
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
		//if (i > 0 && i < NumEdges-1)
		{
			if (Chan->WD[i-1] == 0 && Chan->WD[i] == 0)
			{
				// Reflection flux for the left element
				A_L = Chan->A[leftNode];
				A_R = A_L;
				Q_L = Chan->Q[leftNode];
				Q_R = -Q_L;
				
				if (fluxType == 1)
					RoeFluxChan(tmpF, A_L, A_R, Q_L, Q_R,bval, m1val, m2val,0);
				else if (fluxType == 2)
					LFChan(tmpF, A_L, A_R, Q_L, Q_R,bval, m1val, m2val, 0);
		
				else
				{
					printf("Unknown flux type. Exiting now \n");
					exit(1);
				}

				Fhat1L[i] = tmpF[0];
				Fhat2L[i] = tmpF[1];

				// Reflection flux for the right element
				A_R = Chan->A[rightNode];
				A_L = A_R;
				Q_R = Chan->Q[rightNode];
				Q_L = -Q_R;
				
				if (fluxType == 1)
					RoeFluxChan(tmpF, A_L, A_R, Q_L, Q_R, bval, m1val, m2val, 0);
				else if (fluxType == 2)
					LFChan(tmpF, A_L, A_R, Q_L, Q_R, bval, m1val, m2val, 0);
				else
				{
					printf("Unknown flux type. Exiting now \n");
					exit(1);
				}
				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];

				if (isnan(tmpF[0]) || isnan(tmpF[1]))
				{
					printf("both elements dry flux not a number, edge %d \n",i);
					exit(EXIT_FAILURE);
				}

			}

		}

		// if the elements are not both dry
		if ((i == 0) || (i == NumEl) || Chan->WD[i-1] == 1 || Chan->WD[i] == 1)
		{
			A_L = Chan->A[leftNode];
			A_R = Chan->A[rightNode];

			Q_L = Chan->Q[leftNode];
			Q_R = Chan->Q[rightNode];

			if (fluxType == 1)
				RoeFluxChan(tmpF, A_L, A_R, Q_L,Q_R,bval, m1val, m2val, g);
			else if (fluxType == 2)
				LFChan(tmpF, A_L, A_R, Q_L,Q_R,bval, m1val, m2val, g);

			else
			{
				printf("Unknown flux type. Exiting now \n");
				exit(1);
			}
			
			if (isnan(tmpF[0]) || isnan(tmpF[1]))
			{
				printf("A_L = %lf A_R = %lf Q_L = %lf Q_R = %lf\n", A_L, A_R, Q_L, Q_R);
				printf(" both elements wet flux not a number, edge %d \n",i);
				exit(EXIT_FAILURE);
			}

			if (i==0 || Chan->WD[i-1] == 1)
			{
				Fhat1L[i] = tmpF[0];
				Fhat2L[i] = tmpF[1];	

			}
			else
			{
				//printf("left element dry, edge %d\n", i);
				double newtmpF[2];
				if (fluxType == 1)
					RoeFluxChan(newtmpF, A_L, A_R, Q_L, Q_R, bval, m1val, m2val, 0);
				else if (fluxType == 2)
					LFChan(newtmpF, A_L, A_R, Q_L, Q_R, bval, m1val, m2val, 0);
				else
				{
					printf("Unknown flux type. Exiting now \n");
					exit(1);
				}
				Fhat1L[i] = newtmpF[0];
				Fhat2L[i] = newtmpF[1];
				if (isnan(newtmpF[0]) || isnan(newtmpF[1]))
				{
					printf("left element dry flux not a number, edge %d \n",i);
					exit(EXIT_FAILURE);
				}

			}
			if (i == NumEl || Chan->WD[i]==1)
			{
				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];
			}
			else
			{
				//printf("right element dry, edge %d\n",i);
				if (fluxType == 1)
					RoeFluxChan(tmpF, A_L,A_R,Q_L,Q_R,bval, m1val, m2val, 0);
				else if (fluxType == 2)
					LFChan(tmpF, A_L,A_R,Q_L,Q_R,bval, m1val, m2val, 0);
				else
				{
					printf("Unknown flux type. Exiting now \n");
					exit(1);
				}
				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];		
				if (isnan(tmpF[0]) || isnan(tmpF[1]))
				{
					printf("right element dry flux not a number, edge %d \n",i);
					printf("A_L = %lf A_R = %lf Q_L = %lf Q_R = %lf bval = %lf m1val = %lf m2val = %lf\n", A_L, A_R, Q_L, Q_R, bval, m1val, m2val);
					printf("F = %lf\n", tmpF[0]);
					exit(EXIT_FAILURE);
				}

			}
		}

#else
		A_L = Chan->A[leftNode];
		A_R = Chan->A[rightNode];

		Q_L = Chan->Q[leftNode];
		Q_R = Chan->Q[rightNode];


		if (fluxType == 1)
			RoeFluxChan(tmpF, A_L, A_R, Q_L,Q_R,bval,m1val,m2val,g);
		else if (fluxType == 2)
			LFChan(tmpF, A_L, A_R, Q_L,Q_R,bval,m1val,m2val,g);
		else
		{
			printf("Unknown flux type. Exiting now \n");
			exit(1);
		}
		
		if (isnan(tmpF[0]) || isnan(tmpF[1]))
		{
			printf("flux not a number, edge %d \n",i);

			printf("A_L = %3.16f A_R = %3.16f Q_L = %3.16f Q_R = %31.6f, b = %lf, m1 = %lf, m2 = %lf \n", A_L, A_R, Q_L, Q_R, Chan->b[i], m1val, m2val);

			printf("Fhat = %lf Ghat = %lf\n", tmpF[0], tmpF[1]);
			exit(EXIT_FAILURE);
		}

		Fhat1L[i] = tmpF[0];
		Fhat2L[i] = tmpF[1];	

		Fhat1R[i] = tmpF[0];
		Fhat2R[i] = tmpF[1];

			//printf("edge %d A_L = %3.16f A_R = %3.16f Q_L = %3.16f Q_R = %3.16f, b = %lf, m1 = %lf, m2 = %lf \n",i, A_L, A_R, Q_L, Q_R, Chan->b[i], m1val, m2val);
/*		if (i == (Chan->NumEdges-1))
		{
			printf("A_L = %lf A_R = %lf Q_L = %lf Q_R = %lf, b = %lf, m1 = %lf, m2 = %lf \n", A_L, A_R, Q_L, Q_R, Chan->b[i], m1val, m2val);
			printf("Fhat = %lf Ghat = %lf \n", Fhat1L[i], Fhat2L[i]);
		}
*/
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
		int cont = 1;
		while(cont)
		{
			if (Fhat1L[i]*maxBetaOverAlpha*(*dt) > mass[i-1])
			{
				double A_L = Chan->A[leftNode];
				double A_R = A_L;
				double Q_L = Chan->Q[leftNode];
				double Q_R = -Q_L;
				double bval = Chan->NodalB[leftNode];
				double m1val = Chan->Nodalm1[leftNode];
				double m2val = Chan->Nodalm2[leftNode];

				double tmpF[2];
				double localG = g;

				if (Chan->WD[i] == 0)
					localG = 0;

				if (fluxType == 1)
					RoeFluxChan(tmpF, A_L, A_R, Q_L, Q_R, bval,m1val,m2val,localG);
				else if (fluxType == 2)
					LFChan(tmpF, A_L, A_R, Q_L, Q_R, bval,m1val,m2val,localG);
				else
				{
					printf("Unknown flux type. Exiting now \n");
					exit(1);
				}
				Fhat1L[i] = tmpF[0];
				Fhat2L[i] = tmpF[1];
				printf("element %d will be dry\n",i);
				(*dt) = 0.5*(*dt);
			}
			else
				cont = 0;
		}

		cont = 1;
		while(cont)
		{
			// If the mass will be negative on the right side
			if (-Fhat1R[i]*maxBetaOverAlpha*(*dt) > mass[i])
			{
				double A_R = Chan->A[rightNode];
				double A_L = A_R;
				double Q_R = Chan->Q[rightNode];
				double Q_L = -Q_R;
				double bval = Chan->NodalB[rightNode];
				double m1val = Chan->Nodalm1[rightNode];
				double m2val = Chan->Nodalm2[rightNode];

				double tmpF[2];
				double localG = g;
				if (Chan->WD[i] == 0)
					localG = 0;

				if (fluxType == 1)
					RoeFluxChan(tmpF, A_L, A_R, Q_L, Q_R, bval, m1val, m2val,localG);
				else if (fluxType == 2)
					LFChan(tmpF, A_L, A_R, Q_L, Q_R, bval,m1val,m2val,localG);
				else
				{
					printf("Unknown flux type. Exiting now \n");
					exit(1);
				}

				Fhat1R[i] = tmpF[0];
				Fhat2R[i] = tmpF[1];
				printf("element %d will be dry \n",i);

				*dt = 0.5*(*dt);

			}
			else
				cont = 0;
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
			//double qL = 0;
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
			double beta = Chan->beta[begNode+i];
			double localG;
#ifdef WDON
			if(Chan->WD[k]==1)
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

	free(Fhat1L);
	free(Fhat1R);
	free(Fhat2L);
	free(Fhat2R);

}

void computeLKinematicEls(double time, double* RHS, int fp)
{
	int NumEl = FloodplainList[fp]->NumEl;
	// calculate and store the flux (Q) at the downstream edge of each kinematic element
	// these fluxes will be used to calculate the upwinded flux values
	double *Fhat = xcalloc(NumEl, sizeof(double));
	for (int i = 0; i < NumEl; ++i)
	{
		struct kinematicEl *kinEl = KinematicElList[i];
		if (kinEl->isActive)
		{
			int Np = kinEl->Np;
			double myA = kinEl->A[Np-1];
			double myS0 = kinEl->dz[Np-1];
			double myNf = kinEl->NodalnFriction[Np-1];
			double myWeq = kinEl->weq;
			double myH = myA/myWeq;

			// using manning's N.
			Fhat[i] = myWeq*sqrt(myS0)*pow(myH, 5.0/3)/myNf;

			// using chezy's relationship. nf stands for C, chezy's coefficient
			//Fhat[i] = myWeq*sqrt(myS0*myH)*myH*myNf;
		}
	}

	// now calculate the right hand side
	for (int i = 0; i < NumEl; ++i)
	{
		struct kinematicEl *kinEl = KinematicElList[i];
		if (kinEl->isActive)
		{
			//printf("active element = %d el1 = %d el2 = %d\n", i, kinEl->el1, kinEl->el2);
			double h_el = kinEl->dh;
			int Np = kinEl->Np;

			gsl_vector *F = gsl_vector_calloc(Np);

			gsl_vector *ST = gsl_vector_calloc(Np);			 // will eventually store R-I (Rainfall - Infiltration)
			gsl_vector *ST1 = gsl_vector_calloc(Np); 		 // Infiltration rate

			double weq = kinEl->weq;

			for (int j = 0; j < Np; ++j)
			{
				double A = kinEl->A[j];
				double nf = kinEl->NodalnFriction[j];
				double S0 = kinEl->dz[j];
				//printf(" i = %d A = %1.13f S0 = %lf\n", i, A, S0);
				double H = A/weq;

				// using manning's N
				double myF = weq*sqrt(S0)*pow(H,5.0/3)/nf;

				// using Chezy's relationship (nf stands for C)
				//double myF = weq*sqrt(S0*H)*H*nf;

				if (S0 < 0 )
					printf("z1 = %lf z2 = %lf\n", kinEl->z1, kinEl->z2);
				if (isnan(myF))
				{
					printf("myF is nan; A = %lf, nf = %lf, S0 = %lf, H = %lf weq = %lf, kinElNum = %d\n", A, nf, S0, H, weq, i);
					exit(1);
				}
				gsl_vector_set(F, j, myF);

				// parking lot case
				if (time  < 1800)
				{
					double rainfall = weq*2.0/12/3600;
					gsl_vector_set(ST, j, rainfall);
				}

				// infiltration
				// solve for Green-Ampt Infiltration
				double K = kinEl->K[j];
				// assume s_e = 0; 
				double s_e =  0;
				double theta_e = kinEl->EffPor[j];
				double psi = kinEl->suctionHead[j];
				double F = calculateInfiltration(psi, s_e, theta_e, K, time);
				double psi_delta_theta = psi*(1-s_e)*theta_e;
				double f = K*(psi_delta_theta/F + 1);
				gsl_vector_set(ST1, j, -f);

			
				// second parking lot case
				//if (time  <= 180)
				//	gsl_vector_set(ST, j, weq*2.0/12/3600);
				//else if (time <= 360)
				//	gsl_vector_set(ST, j, weq*4.0/12/3600);

			}
			// Now calculate total upstream flux. 
			// If the element is at the boundary of the watershed, flux coming in will be 0.
			int numUpstreamEls = kinEl->numUpstreamEls;
			double Fhat1 = 0.0;
			//if (numUpstreamEls > 1)
			//	printf(" i = %d UpEl1 = %d UpEl2 = %d\n", i, kinEl->upstreamEls[0], kinEl->upstreamEls[1]);
			for (int j = 0; j < numUpstreamEls; ++j)
			{
				int upstreamEl = kinEl->upstreamEls[j];
				//printf("upstream el = %d\t", upstreamEl);
				Fhat1 += Fhat[upstreamEl];
			}
			//printf("downstream el = %d\n", kinEl->el2);
			double Fhat2 = Fhat[i];
			gsl_vector *localFhat = gsl_vector_alloc(2);
			gsl_vector_set(localFhat, 0, Fhat1);
			gsl_vector_set(localFhat, 1, Fhat2);

			//printf("el %d fhat1 = %lf fhat2 = %lf\n", i, Fhat1, Fhat2);

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
			gsl_vector_add(localRHS, ST1);

			for (int j = 0; j < Np; j++)
			{
				RHS[i*Np+j] = gsl_vector_get(localRHS, j);
				if (isnan(RHS[i*Np+j]))
				{
					printf("SurfPart = %lf \n", gsl_vector_get(SurfPart,j));
				}
			}
			gsl_vector_free(F);
			gsl_vector_free(localFhat);
			gsl_vector_free(SurfPart);
			gsl_vector_free(localRHS);
			gsl_vector_free(ST);
			gsl_vector_free(ST1);

		}

	}
	free(Fhat);

}


/***********************************************************************************************//**
 * Function for evaluating the right hand side of the discrete ODE obtained from the DG discretization
 * of the 2-D shallow water equations
 * @param [in] junc a pointer to the junction structure corresponding to the junction that is 
 * currently being worked on 
 * @param [in] time  a double representing the current time of the simulation
 * @param [out] RHSZeta a pointer to an array of size NumEl x 3 in which the right hand side for
 * the water surface elevation will be stored
 * @param [out] RHSQx a pointer to an array of size NumEl x 3 in which the right hand side for the
 * momentum in the x-direction is stored
 * @param [out] RHSQy a pointer to an array of size NumEl x 3 in which the right hand side for the
 * momentum in the y-direction is stored
 *
 ***************************************************************************************************/

void compute2DL(struct TwoDRegion *junc, double time, double *RHSZeta, double *RHSQx, double *RHSQy)
{

	/********************************* Compute Roe's flux at the face **********************************************/
	int NumEdges = junc->TotalNumEdges;
	int NumEl = junc->NumEl;
	int Nfp = junc->P + 1;
	int Np = junc->Np;

	double **Fhat1dotn, **Fhat2dotn, **Fhat3dotn;
	Fhat1dotn = malloc(NumEl*sizeof(double*));
	Fhat2dotn = malloc(NumEl*sizeof(double*));
	Fhat3dotn = malloc(NumEl*sizeof(double*));
	for (int i = 0; i < NumEl; i++)
	{
		Fhat1dotn[i] = xcalloc(3*Nfp, sizeof(double));
		Fhat2dotn[i] = xcalloc(3*Nfp, sizeof(double));
		Fhat3dotn[i] = xcalloc(3*Nfp, sizeof(double));
	}

	for (int i=0; i<NumEdges; ++i)
	{
		int el1 = junc->EdgtoEls[i*2];
		int el2 = junc->EdgtoEls[i*2+1];

		int ledg1 = junc->GlobaltoLocalEdg[i*2];
		int ledg2 = junc->GlobaltoLocalEdg[i*2+1];

		// take the normal from element 1 side
		double nx = junc->nx[el1*3+ledg1];
		double ny = junc->ny[el1*3+ledg1];

		double tx = -ny;
		double ty = nx;

		for (int j = 0; j < Nfp; j++)
		{
			double zeta_in, Qx_in, Qy_in;
			double zeta_ex, Qx_ex, Qy_ex;
			double z_in, z_ex;			// should be the same as z_ex
			
			int lv1 = junc->GlobalEdgPosNegNodes[i][2*j];
			int lv2 = junc->GlobalEdgPosNegNodes[i][2*j+1];

			int pos1 = junc->PosInFVec[i][2*j];
			int pos2 = junc->PosInFVec[i][2*j+1];

			double Fn_in[3], Fn_ex[3];

			zeta_in = junc->zeta[el1][lv1];
			Qx_in = junc->Qx[el1][lv1];
			Qy_in = junc->Qy[el1][lv1];
			z_in = junc->NodalZ[el1][lv1];
		
			double Q_T_in = Qx_in*tx + Qy_in*ty;

			// need to change this later so that the boundary condition is given per node on the edge
			// if the edge is an exterior edge connected to a channel 
			int bdrypres = junc->BdryPrescribed[i];
			// manufactured solution
			if (bdrypres == 555)
			//if (bdrypres == 555 || (el1 == el2))
			{
				double xval = junc->NodalX[el1][lv1];
				double yval = junc->NodalY[el1][lv1];
				zeta_ex = getmanH(xval, yval, time);
				Qx_ex = getQx(xval, yval, time);
				Qy_ex = getQy(xval, yval, time);
			}
			// if inflow or outflow boundary, i.e. if bdrypres = 1 or 2
			//for floodplains
			else if (bdrypres == 111)	
			{
				zeta_ex = junc->bzeta[i];
				Qx_ex = Qx_in;
				Qy_ex = Qy_in;
			}

			else if (bdrypres == 222)
			{
				zeta_ex = zeta_in;
				double Q_N_ex = junc->bQn[i];
				double Q_T_ex = Q_T_in;
				double denom = 1./(nx*ty-ny*tx);
				Qx_ex = (ty*Q_N_ex - ny*Q_T_ex)*denom;
				Qy_ex = (-tx*Q_N_ex + nx*Q_T_ex)*denom;
			}

			else if (bdrypres == 333)
			{
				zeta_ex = junc->bQn[i];
				double Q_N_ex = junc->bQn[i];
				double Q_T_ex = Q_T_in;
				double denom = 1./(nx*ty-ny*tx);
				Qx_ex = (ty*Q_N_ex - ny*Q_T_ex)*denom;
				Qy_ex = (-tx*Q_N_ex + nx*Q_T_ex)*denom;
		
			}
			// for junctions
			else if (bdrypres == 1 || bdrypres == 2)
			{
				zeta_ex = junc->bzeta[i];
				double Q_N_ex = junc->bQn[i];
				// for the first iteration
				if (zeta_ex == -10000.0 && Q_N_ex == -10000.0)
				{
					zeta_ex = zeta_in;
					Qx_ex = Qx_in;
					Qy_ex = Qy_in;
					z_ex = z_in;
				}
				else
				{
					double Q_T_ex = Q_T_in;
					double denom = 1./(nx*ty-ny*tx);
					Qx_ex = (ty*Q_N_ex - ny*Q_T_ex)*denom;
					Qy_ex = (-tx*Q_N_ex + nx*Q_T_ex)*denom;
					z_ex = z_in;
				}
			}
			// if edge i is a boundary edge but not connected to a channel, implement no flux boundary condition
			else if (el1 == el2)
			{
				zeta_ex  = junc->zeta[el2][lv2];
				// Compute the velocity in the normal direction
				double Q_N_in = Qx_in*nx + Qy_in*ny;
				double Q_T_in = Qx_in*tx + Qy_in*ty;

				// Reflect the velocity in the normal direction
				double Q_N_ex = -Q_N_in;
				double Q_T_ex = Q_T_in;

				// Compute the x and y components of the external state flow
				double denom = 1./(nx*ty - ny*tx);
				Qx_ex = (ty*Q_N_ex - ny*Q_T_ex)*denom;
				Qy_ex = (-tx*Q_N_ex + nx*Q_T_ex)*denom;

				z_ex = z_in;

			}

			// if the edge is not a boundary edge
			else
			{
				zeta_ex = junc->zeta[el2][lv2];
				Qx_ex = junc->Qx[el2][lv2];
				Qy_ex = junc->Qy[el2][lv2];
				z_ex = junc->NodalZ[el2][lv2];
			}

			#ifdef WDON
			// Check to see if both of the elements separated by this edge are dry
			if (junc->WD[el1] == 0 && junc->WD[el2] == 0)
			{
				// Reflection flux for the interior element
				double zeta_ex_ref = zeta_in;
				double Q_N_in = Qx_in*nx + Qy_in*ny;
				double Q_T_in = Qx_in*tx + Qy_in*ty;
			
				double Q_N_ex = -Q_N_in;
				double Q_T_ex = Q_T_in;
				
				double denom = 1./(nx*ty - ny*tx);
				double Qx_ex_ref = (ty*Q_N_ex - ny*Q_T_ex)*denom;
				double Qy_ex_ref = (-tx*Q_N_ex + nx*Q_T_ex)*denom;
			
				double max_lam_in = RoeFluxJunc(zeta_in, zeta_ex_ref, Qx_in,Qx_ex_ref,
					Qy_in, Qy_ex_ref, z_in, nx, ny, 0, Fn_in);
			
				// Reflection flux for the exterior element
				double zeta_in_ref = zeta_ex;
				Q_N_ex = Qx_ex*nx + Qy_ex*ny;
				Q_T_ex = Qx_ex*tx + Qy_ex*ty;
			
				Q_N_in = -Q_N_ex;
				Q_T_in = Q_T_ex;
				
				double Qx_in_ref = (ty*Q_N_in - ny*Q_T_in)*denom;
				double Qy_in_ref = (-tx*Q_N_in + nx*Q_T_ex)*denom;
				double max_lam_ex = RoeFluxJunc(zeta_in_ref, zeta_ex, Qx_in_ref, Qx_ex, Qy_in_ref, Qy_ex, z_in, nx, ny, 0, Fn_ex);
			
				if (i ==0)
					junc->max_lambda = max(max_lam_in, max_lam_ex);
				else
				{
					junc->max_lambda = fmax(junc->max_lambda, max_lam_in);
					junc->max_lambda = fmax(junc->max_lambda, max_lam_ex);
				}
			
			}
			else // if the elements aren't both dry
			{
				double Fhatdotn[3];
				double current_max_lam= RoeFluxJunc(zeta_in, zeta_ex, Qx_in, Qx_ex, Qy_in, Qy_ex, z_in, nx, ny, g,Fhatdotn);
				if (isnan(Fhatdotn[0]) || isnan(Fhatdotn[1]) || isnan(Fhatdotn[2]))
				{
					printf("2D both elements wet flux not a number, edge %d, time %e \n",i,time);
					exit(EXIT_FAILURE);
				}
				
				if (junc->WD[el1] == 1)
				{
					for (int j =0; j < 3; ++j)
					{
						Fn_in[j] = Fhatdotn[j];
					}
				}
				else
				{
					double max_lam_in = RoeFluxJunc(zeta_in, zeta_ex, Qx_in, Qx_ex, Qy_in, Qy_ex,
						z_in, nx, ny, 0, Fn_in);
					if (isnan(Fn_in[0]) || isnan(Fn_in[1]) || isnan(Fn_in[2]))
					{
						printf("2D interior element dry flux not a number, edge %d \n",i);
						exit(EXIT_FAILURE);
					}
					current_max_lam = max(current_max_lam, max_lam_in);
				}
		
				if (junc->WD[el2] == 1)
				{
					for (int j = 0; j < 3; ++j)
					{
						Fn_ex[j] = Fhatdotn[j];
					}
				}
				else
				{
					double max_lam_ex = RoeFluxJunc(zeta_in, zeta_ex, Qx_in, Qx_ex,
						Qy_in, Qy_ex, z_in, nx, ny, 0, Fn_ex);
					if (isnan(Fn_ex[0]) || isnan(Fn_ex[1]) || isnan(Fn_ex[2]))
					{
						printf("2D exterior element dry flux not a number, edge %d \n",i);
						exit(EXIT_FAILURE);
					}
					current_max_lam = max(current_max_lam, max_lam_ex);
				}
					
				if (i == 0)
					junc->max_lambda = current_max_lam;
				else
					junc->max_lambda = fmax(junc->max_lambda, current_max_lam);
			}
			
		#else

			double Fhatdotn[3];
			double current_max_lam= RoeFluxJunc(zeta_in, zeta_ex, Qx_in, Qx_ex, Qy_in, Qy_ex, z_in, nx, ny, g, Fhatdotn);
			//printf("zeta_in = %e, zeta_ex = %e, z_edge = %e, BdryPrescribed = %d \n", zeta_in, zeta_ex, z_in, bdrypres);
			//printf("Fhat1dotn = %e\n", Fhatdotn[0]);
			//printf("Qx_in = %e, Qy_in = %e, Qx_ex = %e, Qy_ex = %e\n", Qx_in, Qy_in, Qx_ex, Qy_ex);
			if (isnan(Fhatdotn[0]) || isnan(Fhatdotn[1]) || isnan(Fhatdotn[2]))
			{
				printf("2D numerical flux not a number, edge %d \n",i);
				printf("zeta_in = %e, zeta_ex = %e, z_edge = %e, BdryPrescribed = %d \n", zeta_in, zeta_ex, z_in, bdrypres);
				exit(EXIT_FAILURE);
			}

			for (int j = 0; j < 3; ++j)
			{
				Fn_in[j] = Fhatdotn[j];
				Fn_ex[j] = Fhatdotn[j];
			}

			if (i == 0)
				junc->max_lambda = current_max_lam;
			else
				junc->max_lambda = max(junc->max_lambda, current_max_lam);

		#endif

			// store Fhatdotn for the two elements connected by this edge
			Fhat1dotn[el1][pos1] = Fn_in[0];
			Fhat2dotn[el1][pos1] = Fn_in[1];
			Fhat3dotn[el1][pos1] = Fn_in[2];

			// store the value for the exterior element, only if the edge is not a boundary edge
			if (el2 != el1)
			{
				Fhat1dotn[el2][pos2] = -Fn_ex[0];
				Fhat2dotn[el2][pos2] = -Fn_ex[1];
				Fhat3dotn[el2][pos2] = -Fn_ex[2];
			}

		} // end node loop


	}	// end edge loop


	for (int k=0; k<NumEl; ++k)
	{
		gsl_vector *F1x = gsl_vector_alloc(Np);
		gsl_vector *F2x = gsl_vector_alloc(Np);
		gsl_vector *F3x = gsl_vector_alloc(Np);
		gsl_vector *F1y = gsl_vector_alloc(Np);
		gsl_vector *F2y = gsl_vector_alloc(Np);
		gsl_vector *F3y = gsl_vector_alloc(Np);
		gsl_vector *ST21 = gsl_vector_calloc(Np);
		gsl_vector *ST22 = gsl_vector_calloc(Np);
		gsl_vector *ST31 = gsl_vector_calloc(Np);
		gsl_vector *ST32 = gsl_vector_calloc(Np);

		double localG = g;
		#ifdef WDON
			if(junc->WD[k] == 0)
				localG = 0;
		#endif

		for (int i = 0; i < Np; i++)
		{
			double zeta = junc->zeta[k][i];
			double Qx = junc->Qx[k][i];
			double Qy = junc->Qy[k][i];
			double z = junc->NodalZ[k][i];
			double H = zeta + z;
			double f2x = Qx*Qx/H + 0.5*localG*(H*H - z*z);
			double f2y = Qx*Qy/H;
			double f3x = f2y;
			double f3y = Qy*Qy/H + 0.5*localG*(H*H - z*z);
			gsl_vector_set(F1x, i, Qx);
			gsl_vector_set(F1y, i, Qy);
			gsl_vector_set(F2x, i, f2x);
			gsl_vector_set(F2y, i, f2y);
			gsl_vector_set(F3x, i, f3x);
			gsl_vector_set(F3y, i, f3y);


			double u = Qx/H;
			double v = Qy/H;
			double nfric = junc->NodalnFriction[k][i];
			double tau = g*nfric*nfric*sqrt(u*u + v*v)/pow(H,4.0/3);
			double dzx = junc->Nodaldzx[k][i];
			double dzy = junc->Nodaldzy[k][i];
			//printf("dzx = %lf , dzy = %lf\n", dzx, dzy);
			gsl_vector_set(ST21, i, localG*zeta*dzx);
			gsl_vector_set(ST22, i, -tau*Qx);
			gsl_vector_set(ST31, i, localG*zeta*dzy);
			gsl_vector_set(ST32, i, -tau*Qy);
			
			// Source terms for manufactured solution
			//double xval = junc->NodalX[k][i];
			//double yval = junc->NodalY[k][i];

			//double st2 = getS2(xval, yval, time);
			//double st3 = getS3(xval, yval, time);
			//gsl_vector_set(ST21, i, st2);
			//gsl_vector_set(ST31, i, st3);
			
			
		}

		// Calculate the volume integral of the flux
		double rx = junc->rx[k];
		double ry = junc->ry[k];
		double sx = junc->sx[k];
		double sy = junc->sy[k];

		gsl_vector *localRHS1 = gsl_vector_calloc(Np);
		gsl_vector *localRHS2 = gsl_vector_calloc(Np);
		gsl_vector *localRHS3 = gsl_vector_calloc(Np);


		gsl_vector *lR12 = gsl_vector_calloc(Np);
		gsl_vector *lR13 = gsl_vector_calloc(Np);
		gsl_vector *lR14 = gsl_vector_calloc(Np);
		gsl_vector *lR22 = gsl_vector_calloc(Np);
		gsl_vector *lR23 = gsl_vector_calloc(Np);
		gsl_vector *lR24 = gsl_vector_calloc(Np);
		gsl_vector *lR32 = gsl_vector_calloc(Np);
		gsl_vector *lR33 = gsl_vector_calloc(Np);
		gsl_vector *lR34 = gsl_vector_calloc(Np);
	
		gsl_blas_dgemv(CblasNoTrans, rx, Drw, F1x, 0.0, localRHS1);
		gsl_blas_dgemv(CblasNoTrans, sx, Dsw, F1x, 0.0, lR12);
		gsl_blas_dgemv(CblasNoTrans, ry, Drw, F1y, 0.0, lR13);
		gsl_blas_dgemv(CblasNoTrans, sy, Dsw, F1y, 0.0, lR14);

	
		gsl_vector_add(localRHS1, lR12);
		gsl_vector_add(localRHS1, lR13);
		gsl_vector_add(localRHS1, lR14);

		gsl_blas_dgemv(CblasNoTrans, rx, Drw, F2x, 0.0, localRHS2);
		gsl_blas_dgemv(CblasNoTrans, sx, Dsw, F2x, 0.0, lR22);
		gsl_blas_dgemv(CblasNoTrans, ry, Drw, F2y, 0.0, lR23);
		gsl_blas_dgemv(CblasNoTrans, sy, Dsw, F2y, 0.0, lR24);

		gsl_vector_add(localRHS2, lR22);
		gsl_vector_add(localRHS2, lR23);
		gsl_vector_add(localRHS2, lR24);


		gsl_blas_dgemv(CblasNoTrans, rx, Drw, F3x, 0.0, localRHS3);
		gsl_blas_dgemv(CblasNoTrans, sx, Dsw, F3x, 0.0, lR32);
		gsl_blas_dgemv(CblasNoTrans, ry, Drw, F3y, 0.0, lR33);
		gsl_blas_dgemv(CblasNoTrans, sy, Dsw, F3y, 0.0, lR34);

	
		gsl_vector_add(localRHS3, lR32);
		gsl_vector_add(localRHS3, lR33);
		gsl_vector_add(localRHS3, lR34);

		
		// calculate the surface integral of the flux
		gsl_vector *localFhat1dotn = gsl_vector_alloc(3*Nfp);
		gsl_vector *localFhat2dotn = gsl_vector_alloc(3*Nfp);
		gsl_vector *localFhat3dotn = gsl_vector_alloc(3*Nfp);

		for (int i = 0; i < 3; i++)
		{
			double edgJac = junc->edgJac[k*3+i];
			for (int j = 0; j < Nfp; j++)
			{
					int index = i*Nfp+j;
					gsl_vector_set(localFhat1dotn, index, edgJac*Fhat1dotn[k][index]);
					gsl_vector_set(localFhat2dotn, index, edgJac*Fhat2dotn[k][index]);
					//printf("edgJac = %lf Fhat = %lf\n", edgJac, Fhat2dotn[k][index]);
					gsl_vector_set(localFhat3dotn, index, edgJac*Fhat3dotn[k][index]);
			}
		}

		double jac = junc->jac[k];
		
		gsl_vector *SurfPart1 = gsl_vector_calloc(Np);
		gsl_vector *SurfPart2 = gsl_vector_calloc(Np);
		gsl_vector *SurfPart3 = gsl_vector_calloc(Np);
	
		double fac = 1.0/jac;
		gsl_blas_dgemv(CblasNoTrans, fac, LIFT2D, localFhat1dotn, 0.0, SurfPart1);
		gsl_blas_dgemv(CblasNoTrans, fac, LIFT2D, localFhat2dotn, 0.0, SurfPart2);
		gsl_blas_dgemv(CblasNoTrans, fac, LIFT2D, localFhat3dotn, 0.0, SurfPart3);
		
		// subtract this from the RHS
		gsl_vector_sub(localRHS1, SurfPart1);
		gsl_vector_sub(localRHS2, SurfPart2);
		gsl_vector_sub(localRHS3, SurfPart3);

		// add source and sink terms
		//gsl_vector_add(localRHS1, ST11);
		gsl_vector_add(localRHS2, ST21);
		gsl_vector_add(localRHS2, ST22);
		gsl_vector_add(localRHS3, ST31);
		gsl_vector_add(localRHS3, ST32);


		int begNode = k*Np;
		for(int i = 0; i < Np; i++)
		{
			RHSZeta[begNode+i] = gsl_vector_get(localRHS1,i);
			RHSQx[begNode+i] = gsl_vector_get(localRHS2,i);
			RHSQy[begNode+i] = gsl_vector_get(localRHS3, i);
		}

	
		// free all gsl vectors
		gsl_vector_free(localRHS1);
		gsl_vector_free(localRHS2);
		gsl_vector_free(localRHS3);
		gsl_vector_free(localFhat1dotn);
		gsl_vector_free(localFhat2dotn);
		gsl_vector_free(localFhat3dotn);
		gsl_vector_free(SurfPart1);
		gsl_vector_free(SurfPart2);
		gsl_vector_free(SurfPart3);
		gsl_vector_free(lR12);
		gsl_vector_free(lR13);	
		gsl_vector_free(lR14);
		gsl_vector_free(lR22);
		gsl_vector_free(lR23);
		gsl_vector_free(lR24);
		gsl_vector_free(lR32);
		gsl_vector_free(lR33);
		gsl_vector_free(lR34);
		gsl_vector_free(F1x); 
		gsl_vector_free(F2x); 
		gsl_vector_free(F3x); 
		gsl_vector_free(F1y); 	
		gsl_vector_free(F2y); 
		gsl_vector_free(F3y); 
		gsl_vector_free(ST21);
		gsl_vector_free(ST22);
		gsl_vector_free(ST31);
		gsl_vector_free(ST32);
	
	}

	// free allocated space
	for (int i = 0; i < NumEl; i++)
	{
		free(Fhat1dotn[i]);
		free(Fhat2dotn[i]);
		free(Fhat3dotn[i]);
	}
	free(Fhat1dotn);
	free(Fhat2dotn);
	free(Fhat3dotn);

}
