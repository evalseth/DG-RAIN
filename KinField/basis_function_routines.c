#include <stdio.h>
#include<stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "Watershed.h"
#include "Globals.h"

gsl_matrix *LIFT;
gsl_matrix *VolMat;
gsl_matrix *MassMatrix;

extern void InterpolateWithGSLSplines(int NumPointsToInterpolate, int NumPointsToEvaluate, double*x,
		double*y,	double* NodalX, double* NodalY, double* func, double* NodalFunc, double* NodalDeriv);

// Evaluate normalized Jacobi Polynomial of type (alpha,beta) > -1 
// at points x for order N and returns P of length Nc (= length(x))
void JacobiP(double *x, double alpha, double beta, int N, int Nc, double *P)
{
	double aold=0.0, anew=0.0, bnew=0.0, h1=0.0;
	double gamma0=0.0, gamma1=0.0;
	double ab=alpha+beta, ab1=alpha+beta+1.0, a1=alpha+1.0, b1=beta+1.0;
	double PL[N+1][Nc];
	double prow[Nc], x_bnew[Nc];

	// Initial values P_0(x) and P_1(x)
	gamma0 = pow(2.0,ab1)/(ab1)*tgamma(a1)*tgamma(b1)/tgamma(ab1);

	if (0==N) 
	{
		for (int i = 0; i < Nc; i++)
			P[i]  = 1.0/sqrt(gamma0);  
		return;
	}
	else
	{
		gamma1 = (a1)*(b1)/(ab+3.0)*gamma0;
		for (int i = 0; i < Nc; i++)
		{
			PL[0][i] = 1.0/sqrt(gamma0);
			prow[i] = ((ab+2.0)*x[i]/2.0 + (alpha-beta)/2.0) / sqrt(gamma1);
		}
	}
	if (1 ==N)
	{
		for (int i = 0; i < Nc; i++)
			P[i] = prow[i];
		return;
	}
	else
	{
		for (int i = 0; i < Nc; i++)
			PL[1][i] = prow[i];
	}

	// Repeat value in recurrence
	aold = 2.0/(2.0+ab)*sqrt((a1)*(b1)/(ab+3.0));

	// Forward recurrence using the symmetry of the recurrence.
	for (int i=1; i<=(N-1); ++i) 
	{
		h1 = 2.0*i+ab;
		anew = 2.0/(h1+2.0)*sqrt((i+1)*(i+ab1)*(i+a1)*(i+b1)/(h1+1.0)/(h1+3.0));
		bnew = - (alpha*alpha-beta*beta)/h1/(h1+2.0);
		for (int j = 0; j < Nc; j++)
		{
			x_bnew[j] = x[j]-bnew;
			PL[i+1][j]= 1.0/anew*( -aold*PL[i-1][j] + x_bnew[j]*PL[i][j]);
		}
		aold =anew;
	}

	for (int i=0; i < Nc; i++)
		P[i] = PL[N][i];

}

void GetLGLPoints(int N, double *x)
{
	if (N == 0)
	{
		printf("Support for constants doesn't exist\n");
		exit(1);
	}
	else if (N == 1)
	{
		x[0] = -1.0;
		x[1] = 1.0;
	}
	else if (N==2)
	{
		x[0] = -1.0;
		x[1] = 0.0;
		x[2] = 1.0;
	}
	else if (N==3)
	{
		x[0] = -1.0;
		x[1] = -0.4472;
		x[2] = 0.4472;
		x[3] = 1.0;
	}
	else if (N==4)
	{
		x[0] = -1.0;
		x[1] = -0.6547;
		x[2] = 0;
		x[3] = 0.6547;
		x[4] = 1.0;
	}
	else
	{
		printf("Currently only have up to fourth order polynomial encoded for LGL points\n");
		exit(1);
	}
}

void GetLGLWeights(int N, double *w)
{
	if (N == 0)
	{
		printf("Support for constants doesn't exist\n");
		exit(1);
	}
	else if (N == 1)
	{
		w[0] = 1.0;
		w[1] = 1.0;
	}
	else if (N==2)
	{
		w[0] = 0.3333333333333;
		w[1] = 1.3333333333333;
		w[2] = 0.3333333333333;
	}
	else if (N==3)
	{
		w[0] = 0.16666666666667;
		w[1] = 0.83333333333333;
		w[2] = 0.83333333333333;
		w[3] = 0.16666666666667;
	}
	else if (N==4)
	{
		w[0] = 0.1;
		w[1] = 0.5444444;
		w[2] = 0.7111111;
		w[3] = 0.5444444;
		w[4] = 0.1;
	}
	else
	{
		printf("Currently only have up to fourth order polynomial encoded for LGL points\n");
		exit(1);
	}

}

void calculateNodalCoordinates(int N, int NumEl, double* X, double* Xnodes)
{
	double r[N+1];
	GetLGLPoints(N, r);
	for(int i = 0; i < NumEl; i++)
	{
		for (int j = 0; j < N+1; j++)
		{
			Xnodes[i*(N+1)+j] = X[i] + 0.5*(r[j]+1)*(X[i+1]-X[i]);
		}
	}

}

void LegendrePoly(double *x, int Np, int N, double* P)
{
	// Purpose: Evaluates the normalized Legendre polynomials 
	// 					at points x for order N and returns the values in the
	// 					array P that is Np long

	double aold = 0.0, anew=0.0;

	double **PL = xcalloc(N+1, sizeof(double*));
	for (int i =0; i < N+1; i++)
	{
		PL[i] =  xcalloc(Np, sizeof(double));
	}

	// Initial values P_0(x) and P_1(x)
	if (N==0)
	{
		double constVal = 1.0/sqrt(2);
		for (int j = 0; j < Np; j++)
		{
			P[j] = constVal;
		}

		free(PL[0]);
		free(PL); 
		return;
	}
	else
	{
		for (int j = 0; j < Np; j++)
		{
			PL[0][j] = 1.0/sqrt(2);
		}
	}

	double *prow = xcalloc(Np, sizeof(double));
	double coeff = sqrt(1.5);
	for (int i = 0; i < Np; i++)
	{
		prow[i] = coeff*x[i];
	}

	if(N==1)
	{
		for (int i = 0; i < Np; i++)
		{
			P[i] = prow[i];
		}
		free(prow);
		free(PL[0]);
		free(PL[1]);
		free(PL);
		return;
	}
	else
	{
		for (int j =0; j < Np; j++)
		{
			PL[1][j] = coeff*x[j];
		}

	}

	// Repeat value in recurrence
	aold = sqrt(1.0/3);

	// Forward recurrence
	for (int i = 2; i <= N; i++)
	{
		anew = sqrt(i*i/((2.0*i+1)*(2*i-1)));
		for (int j = 0; j < Np; j++)
		{
			PL[i][j] = (x[j]*PL[i-1][j]-aold*PL[i-2][j])/anew;
		}
		aold = anew;
	}

	for (int j = 0; j < Np; j++)
	{
		P[j] = PL[N][j];
	}

	for (int i =0; i < N; i++)
	{
		free(PL[i]);
	}
	free(PL);

	return;
}

void Vandermonde1D(double *x, int Np, int N, gsl_matrix* V1D)
{
	// Purpose: Evaluates the transpose of the Vandermonde matrix V_{ij} = phi_i(r_j);

	double V1DT[Np];
	for(int i = 0; i < N+1; i++)
	{
		LegendrePoly(x, Np, i, V1DT);
		for (int j = 0; j < Np; j++)
		{
			gsl_matrix_set(V1D, j, i, V1DT[j]);
			//gsl_matrix_set(*V1D, j, i, V1DT[j]);
		}
	}
}

void GradJacobiP(double *x, double alpha, double beta, int N, int Np, double *dP)
{
	if (0 == N)
	{
		for (int i = 0; i < Np; i++)
			dP[i] = 0;
	}
	else
	{
		double JP[Np];
		JacobiP(x, alpha+1,beta+1, N-1, Np, JP);
		for (int i = 0; i < Np; i++)
			dP[i] = sqrt(N*(N+alpha+beta+1))*JP[i];
	}
}

void GradLegendrePoly(double *x, int Np, int N, double* DP)
{
	// Purpose: Evaluates the derivatives of normalized Legendre polynomial  
	// 					at points x for order N and returns the values in the
	// 					array P that is Np long
	//

	double aold, anew;

	double **DPL = xcalloc(N+1, sizeof(double*));
	for (int i =0; i < N+1; i++)
	{
		DPL[i] =  xcalloc(Np, sizeof(double));
	}

	// Initial values P_0(x) and P_1(x)
	if (N==0)
	{
		for (int j = 0; j < Np; j++)
		{
			DP[0] = 0;
		}
		free(DPL[0]);
		free(DPL); 
		return;
	}
	else
	{
		for (int j = 0; j < Np; j++)
		{
			DPL[0][j] = 0;
		}
	}

	double *dprow = xcalloc(Np, sizeof(double));
	double deriv = sqrt(1.5);
	for (int i = 0; i < Np; i++)
	{
		dprow[i] = deriv;
	}

	if(N==1)
	{
		for (int i = 0; i < Np; i++)
		{
			DP[i] = deriv;
		}
		free(dprow);
		free(DPL[0]);
		free(DPL[1]);
		free(DPL);
		return;
	}
	else
	{
		for (int j =0; j < Np; j++)
		{
			DPL[1][j] = deriv;
		}

	}

	// Repeat value in recurrence
	aold = sqrt(1.0/3);

	// Forward recurrence
	for (int i = 2; i <= N; i++)
	{
		anew = sqrt(i*i/((2.0*i+1)*(2*i-1)));
		double P[Np];
		LegendrePoly(x, Np, i-1, P);
		for (int j = 0; j < Np; j++)
		{
			DPL[i][j] = 1.0/(anew)*(P[j] + x[j]*DPL[i-1][j] - aold*DPL[i-2][j]);
		}
		aold = anew;
	}

	for (int j = 0; j < Np; j++)
	{
		DP[j] = DPL[N][j];
	}

	for (int i =0; i < N; i++)
	{
		free(DPL[i]);
	}
	free(DPL);

	return;

}

void GradVandermonde1D(double *r, int Np, int N, gsl_matrix* GV)
{
	// Purpose: Evaluates the gradient of modal basis (i) at (r)
	// 					at order N and return them in GV whose size is
	// 					ixNp
	//
	double *GVT = xcalloc(Np,sizeof(double));
	for (int i =0; i < N+1; i++)
	{
		GradLegendrePoly(r, Np, i, GVT);
		for (int j = 0; j < Np; j++)
			gsl_matrix_set(GV,j,i, GVT[j]);
		//gsl_matrix_set(*GV,j,i,GVT[j]);
	}
	free(GVT);

}

void calculateMassMatrix(int Np, gsl_matrix* V, gsl_matrix* MassMatrix)
{
	gsl_matrix *InvM = gsl_matrix_alloc(Np,Np);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V, V, 0.0, InvM);

	int s = 0;	
	gsl_permutation *p = gsl_permutation_alloc(Np);
	gsl_linalg_LU_decomp(InvM, p, &s);
	gsl_linalg_LU_invert(InvM, p, MassMatrix);

	gsl_permutation_free(p);
	gsl_matrix_free(InvM);

}

void VolIntMat1D(double *x, int N, gsl_matrix* V, gsl_matrix *Vr, gsl_matrix** VolMat)
{

	// S = M*Dr
	// inv(M)*S = Dr. But we need VolMat = (inv)M*S'
	// inv(M)*S' = VV'*Dr'*M' = V*Vr'*M;

	// calculate M
	gsl_matrix *M = gsl_matrix_calloc(N+1, N+1);
	calculateMassMatrix(N+1, V, M);

	// calculate VolMat
	gsl_matrix *tmp = gsl_matrix_calloc(N+1, N+1);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Vr, M, 0.0, tmp);

	*VolMat = gsl_matrix_calloc(N+1, N+1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, tmp, 0.0, *VolMat);

	gsl_matrix_free(M);
	gsl_matrix_free(tmp);

}

void Lift1D(int Np, gsl_matrix *V, gsl_matrix **LIFT)
{
	*LIFT = gsl_matrix_alloc(Np,2);

	gsl_matrix *EMAT = gsl_matrix_calloc(Np, 2);
	gsl_matrix_set(EMAT, 0, 0, 1.0);
	gsl_matrix_set(EMAT, Np-1, 1, -1.0);

	gsl_matrix *tmp = gsl_matrix_alloc(Np,2);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, V, EMAT, 0.0, tmp);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, tmp, 0.0, *LIFT);

	gsl_matrix_free(EMAT);
	gsl_matrix_free(tmp);

}




void setup_1D_domains()
{
	int P = 1;
	int Np = P+1;

	for (int i = 0; i < NumChannels; i++)
	{
		ChannelList[i]->P = P;
		ChannelList[i]->Np = Np;

		int NumEl = ChannelList[i]->NumEl;
		int NumNodes = Np*ChannelList[i]->NumEl;
		ChannelList[i]->NumNodes = NumNodes;
		ChannelList[i]->NodalX = xcalloc(NumNodes, sizeof(double));
		ChannelList[i]->NodalY = xcalloc(NumNodes, sizeof(double));
		ChannelList[i]->NodalZ = xcalloc(NumNodes, sizeof(double));
		ChannelList[i]->NodalB = xcalloc(NumNodes, sizeof(double));
		ChannelList[i]->NodalDepth = xcalloc(NumNodes, sizeof(double));
		ChannelList[i]->Nodalm1 = xcalloc(NumNodes, sizeof(double));
		ChannelList[i]->Nodalm2 = xcalloc(NumNodes, sizeof(double));
		ChannelList[i]->NodalnFriction = xcalloc(NumNodes, sizeof(double));
		ChannelList[i]->dz = xcalloc(NumNodes, sizeof(double));
		ChannelList[i]->db = xcalloc(NumNodes, sizeof(double));
		ChannelList[i]->dm1 = xcalloc(NumNodes, sizeof(double));
		ChannelList[i]->dm2 = xcalloc(NumNodes, sizeof(double));

		calculateNodalCoordinates(P, NumEl, ChannelList[i]->x, ChannelList[i]->NodalX);
		calculateNodalCoordinates(P, NumEl, ChannelList[i]->y, ChannelList[i]->NodalY);
		calculateNodalCoordinates(P, NumEl, ChannelList[i]->z, ChannelList[i]->NodalZ);
		calculateNodalCoordinates(P, NumEl, ChannelList[i]->nFriction, ChannelList[i]->NodalnFriction);
		calculateNodalCoordinates(P, NumEl, ChannelList[i]->b, ChannelList[i]->NodalB);
		//calculateNodalCoordinates(P, NumEl, ChannelList[i]->depth, ChannelList[i]->NodalDepth);
		calculateNodalCoordinates(P, NumEl, ChannelList[i]->m1, ChannelList[i]->Nodalm1);
		calculateNodalCoordinates(P, NumEl, ChannelList[i]->m2, ChannelList[i]->Nodalm2);

		for (int j = 0; j < NumEl; j++)
		{
			double dzval =  (ChannelList[i]->z[j+1] - ChannelList[i]->z[j])/ChannelList[i]->dh[j];
			double dbval =  (ChannelList[i]->b[j+1] - ChannelList[i]->b[j])/ChannelList[i]->dh[j];
			double dm1val =  (ChannelList[i]->m1[j+1] - ChannelList[i]->m1[j])/ChannelList[i]->dh[j];
			double dm2val =  (ChannelList[i]->m2[j+1] - ChannelList[i]->m2[j])/ChannelList[i]->dh[j];

			for (int k = 0; k <Np; k++)
			{
				ChannelList[i]->dz[j*Np+k] = dzval;
				ChannelList[i]->db[j*Np+k] = dbval;
				ChannelList[i]->dm1[j*Np+k] = dm1val;
				ChannelList[i]->dm2[j*Np+k] = dm2val;
			}
		}
		//InterpolateWithGSLSplines(ChannelList[i]->NumEdges, ChannelList[i]->NumNodes, ChannelList[i]->x, ChannelList[i]->y, ChannelList[i]->NodalX, ChannelList[i]->NodalY, ChannelList[i]->z, ChannelList[i]->NodalZ, ChannelList[i]->dz);

		//InterpolateWithGSLSplines(ChannelList[i]->NumEdges, NumNodes, ChannelList[i]->x, ChannelList[i]->y, ChannelList[i]->NodalX, ChannelList[i]->NodalY, ChannelList[i]->b, ChannelList[i]->NodalB, ChannelList[i]->db);
		//InterpolateWithGSLSplines(ChannelList[i]->NumEdges, NumNodes, ChannelList[i]->x, ChannelList[i]->y, ChannelList[i]->NodalX, ChannelList[i]->NodalY, ChannelList[i]->m1, ChannelList[i]->Nodalm1, ChannelList[i]->dm1);
		//InterpolateWithGSLSplines(ChannelList[i]->NumEdges, NumNodes, ChannelList[i]->x, ChannelList[i]->y, ChannelList[i]->NodalX, ChannelList[i]->NodalY, ChannelList[i]->m2, ChannelList[i]->Nodalm2, ChannelList[i]->dm2);
		//InterpolateWithGSLSplines(ChannelList[i]->NumEdges, NumNodes, ChannelList[i]->x, ChannelList[i]->y, ChannelList[i]->NodalX, ChannelList[i]->NodalY, ChannelList[i]->nFriction, ChannelList[i]->NodalnFriction, NULL);

		//printf("****************After**********\n");
		//for (int j = 0; j < NumEl; j++)
		//{
		//	for (int k = 0; k <Np; k++)
		//		{
		//			printf("z = %3.16f dz = %3.16f\n", ChannelList[i]->NodalZ[j*Np+k], ChannelList[i]->dz[j*Np+k]);
		//		}
		//}

	}

	for (int fp = 0; fp < NumFloodplains; fp++)
	{
		int NumEl = FloodplainList[fp]->NumEl;
		for (int i = 0; i < NumEl; i++)
		{
			struct kinematicEl *myKinEl = KinematicElList[i];
			if (myKinEl->isActive)
			{
				myKinEl->P = P;
				myKinEl->Np = Np;
				myKinEl->NodalX = xcalloc(Np, sizeof(double));
				myKinEl->NodalY = xcalloc(Np, sizeof(double));
				myKinEl->NodalZ = xcalloc(Np, sizeof(double));
				myKinEl->NodalnFriction = xcalloc(Np,sizeof(double));
				myKinEl->dz = xcalloc(Np, sizeof(double));

				double x[2] = {myKinEl->x1, myKinEl->x2};
				double y[2] = {myKinEl->y1, myKinEl->y2};
				double z[2] = {myKinEl->z1, myKinEl->z2};
				double nf[2] = {myKinEl->nf1, myKinEl->nf2};

				calculateNodalCoordinates(P, 1, x, myKinEl->NodalX);
				calculateNodalCoordinates(P, 1, y, myKinEl->NodalY);
				calculateNodalCoordinates(P, 1, z, myKinEl->NodalZ);
				calculateNodalCoordinates(P, 1, nf, myKinEl->NodalnFriction);

				for (int j = 0; j < Np; j++)
				{
					myKinEl->dz[j] = (myKinEl->z2 - myKinEl->z1)/myKinEl->dh;
				}
				//		
				//		//// THESE FUNCTION CALLS DO NOT WORK BECAUSE WE CAN'T DO CUBIC INTERPOLATION WITH TWO POINTS.
				//		//// SO JUST DO LINEAR INTERPOLATION FOR NOW
				//		//InterpolateWithGSLSplines(2, Np, x, y, KinematicElList[i]->NodalX, KinematicElList[i]->NodalY, z, KinematicElList[i]->NodalZ, KinematicElList[i]->dz);
				//		//InterpolateWithGSLSplines(2, Np, x, y, KinematicElList[i]->NodalX, KinematicElList[i]->NodalY, nf, KinematicElList[i]->NodalnFriction, NULL);
			}

		}
	}

	// calculate LIFT, volume matrix and mass matrix
	double r[Np];
	GetLGLPoints(P,r);

	gsl_matrix *V = gsl_matrix_alloc(Np, Np);
	gsl_matrix *GV = gsl_matrix_alloc(Np, Np);

	Vandermonde1D(r,Np,P,V);
	GradVandermonde1D(r, Np, P, GV);
	VolIntMat1D(r, P, V, GV, &VolMat);
	Lift1D(Np, V, &LIFT);

	MassMatrix = gsl_matrix_alloc(Np, Np);
	calculateMassMatrix(Np, V, MassMatrix);

	gsl_matrix_free(V);
	gsl_matrix_free(GV);


}




