#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "Watershed.h"
#include "Globals.h"
#include "Math_Functions.h"

int **Fmask;
gsl_matrix *LIFT2D;
gsl_matrix *Drw;
gsl_matrix *Dsw;
gsl_matrix *MassMatrix2D;

extern void GetLGLPoints(int N, double *x);
extern void Vandermonde1D(double *x, int Np, int N, gsl_matrix* V1D);
extern void LegendrePoly(double *x, int Np, int N, double* P);
extern void JacobiP(double *x, double alpha, double beta, int N, int Nc, double *P);
extern void GradJacobiP(double *x, double alpha, double beta, int N, int Np, double *dP);
void calculateMassMatrix(int Np, gsl_matrix* V, gsl_matrix* MassMatrix);

// Compute scaled warp function at order N based on rout 
// interpolation node
void warpFactor(int N, double *rout, int Nr, double *warp)
{
	// Compute LGL and equidistant node distribution
	double LGLr[N+1];
	GetLGLPoints(N, &LGLr[0]);

	double req[N+1];
	double h = 2.0/N;
	req[0] = -1.0;
	for (int i = 1; i < N+1; i++)
	{
		req[i] = req[i-1] + h;
	}

	// Compute V based on req
	gsl_matrix *Veq = gsl_matrix_alloc(N+1,N+1);
	Vandermonde1D(req, N+1, N, Veq);

	// invert V
	gsl_permutation *p = gsl_permutation_alloc(N+1);
	//gsl_matrix *Veqcopy = gsl_matrix_alloc(N+1, N+1);
	//gsl_matrix_memcpy(Veqcopy, Veq);

	int s;
	gsl_matrix *VeqInv = gsl_matrix_alloc(N+1, N+1);
	gsl_linalg_LU_decomp(Veq, p, &s);
	gsl_linalg_LU_invert(Veq, p, VeqInv);

	// First evaluate Lengendre Polynomial at rout
	gsl_matrix *LPmat = gsl_matrix_alloc(N+1,Nr);
	for (int i = 0; i < N+1; i++)
	{
		double P[Nr];
		LegendrePoly(rout, Nr, i, &P[0]);
		for (int j = 0; j < Nr; j++)
		{
			gsl_matrix_set(LPmat, i, j, P[j]); 
		}
	}

	// Evaluate Lagrange polynomial at rout
	gsl_matrix *Lmat = gsl_matrix_calloc(N+1, Nr);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, VeqInv, LPmat, 0.0, Lmat);

	// Compute warp factor
	for (int i = 0; i < Nr; i++)
	{
		warp[i] = 0;
		for (int j = 0; j < N+1; j++)
		{
			double Lji = gsl_matrix_get(Lmat, j,i);
			warp[i] += Lji*(LGLr[j]-req[j]);
		}

		// scale factor
		double zerof = (fabs(rout[i])< 1.0 - 1.0e-10);
		double sf = 1.0 - (zerof*rout[i])*(zerof*rout[i]);
		warp[i] = warp[i]/sf + warp[i]*(zerof-1);
	}

	gsl_matrix_free(Veq);
	gsl_permutation_free(p);
	gsl_matrix_free(VeqInv);
	gsl_matrix_free(LPmat);
	gsl_matrix_free(Lmat);

}

// Compute nodes in equilateral triangle for polynomial or order N
void getEquiNodes(int N, double *x, double *y)
{
	if (N == 0)
	{
		printf("Support for constants doesn't exist\n");
		exit(1);
	}

	double alpopt[4] = {0.0000, 0.0000, 1.4152, 0.1001};

	// set optimized parameter alpha depending on order N
	double alpha = alpopt[N-1];

	//total number of nodes
	int Np = (N+1)*(N+2)/2;

	double L1[Np], L2[Np], L3[Np];
	double blend1[Np], blend2[Np], blend3[Np];
	double r1[Np], r2[Np], r3[Np];
	int sk = 0;
	for (int n = 1; n < N+2; n++)
	{
		for (int m = 1; m < N+3-n; m++)
		{
			// Create equidistributed nodes on equilateral triangle
			L1[sk] = (n-1.0)/N;
			L3[sk] = (m-1.0)/N;
			L2[sk] = 1.0 - L1[sk]-L3[sk];

			r1[sk] = L3[sk]-L2[sk];
			r2[sk] = L1[sk]-L3[sk];
			r3[sk] = L2[sk]-L1[sk];

			x[sk] = -L2[sk]+L3[sk];
			y[sk] = (-L2[sk] - L3[sk] + 2*L1[sk])/sqrt(3.0);

			// compute blending function at each node for each edge
			blend1[sk] = 4*L2[sk]*L3[sk];
			blend2[sk] = 4*L1[sk]*L3[sk];
			blend3[sk] = 4*L1[sk]*L2[sk];

			sk ++;
		}
	}

	// amount of warp for each node, for each edge
	double warpf1[Np], warpf2[Np], warpf3[Np];
	warpFactor(N, &r1[0], Np, &warpf1[0]);
	warpFactor(N, &r2[0], Np, &warpf2[0]);
	warpFactor(N, &r3[0], Np, &warpf3[0]);

	sk = 0;
	for (int n = 1; n < N+2; n++)
	{
		for(int m = 1; m < N+3-n; m++)
		{
			// combine blend and warp
			double warp1 = blend1[sk]*warpf1[sk]*(1+(alpha*L1[sk])*(alpha*L1[sk]));
			double warp2 = blend2[sk]*warpf2[sk]*(1+(alpha*L2[sk])*(alpha*L1[sk]));
			double warp3 = blend3[sk]*warpf3[sk]*(1+(alpha*L3[sk])*(alpha*L1[sk]));

			// Accumulate deformations associated with each edge
			x[sk] += 1*warp1 + cos(2*PI/3)*warp2 + cos(4*PI/3)*warp3;
			y[sk] += 0*warp1 + sin(2*PI/3)*warp2 + sin(4*PI/3)*warp3;
			sk++;
		}
	}
}

// Map points from equilateral triangle to standard triangle
void xytors(int Np, double *x, double *y, double *r, double *s)
{
	for (int i = 0; i < Np; i++)
	{
		double L1 = (sqrt(3.0)*y[i]+1)/3;
		double L2 = (-3*x[i] -sqrt(3.0)*y[i]+2)/6;
		double L3 = (3.0*x[i] -sqrt(3.0)*y[i]+2)/6;

		r[i] = -L2 + L3 - L1;
		s[i] = -L2 - L3 + L1;

	}
}

// Map points from the standard triangle to the coordinates used
// in basis functions
void rstoab(int Np, double *r, double *s, double *a, double *b)
{
	for (int i = 0; i < Np; i++)
	{
		if (s[i] != 1)
			a[i] = 2*(1+r[i])/(1-s[i])-1;
		else
			a[i] = -1;
		b[i] = s[i];
	}
}

// Obtain nodal coordinates in the physical domain
void interpolateVertsToNodalCoordinates(int N, int NumEl, double* VX, int *EltoVert, double **X)
{
	int Np = (N+1)*(N+2)/2;

	// Compute nodal set
	double xeq[Np], yeq[Np] ,r[Np], s[Np];
	getEquiNodes(N, xeq, yeq);
	xytors(Np, xeq, yeq, r, s);

	// Interpolate values from vertices to nodes
	for (int i = 0; i < NumEl; i++)
	{
		int v1 = EltoVert[i*3];
		int v2 = EltoVert[i*3+1];
		int v3 = EltoVert[i*3+2];
		for (int j = 0; j < Np; j++)
		{
			double lambda1 = 0.5*(s[j]+1);
			double lambda2 = -0.5*(r[j]+s[j]);
			double lambda3 = 0.5*(r[j]+1);
			X[i][j] = lambda2*VX[v1] + lambda3*VX[v2] + lambda1*VX[v3];
		}

	}
}

// Evaluates the transpose of the 2D Vandermonde Matrix V_{ij} = phi_i(r_j, s_j); Np = length(r)
void Vandermonde2D(double *r, double *s, int Np, int N, gsl_matrix* V2D)
{
	// Transfer to (a,b) coordinates
	double a[Np], b[Np];
	rstoab(Np, r, s, a, b); 

	// build the Vandermonde matrix
	int sk = 0;
	for (int i = 0; i < N+1; i++)
	{
		for (int j = 0; j < N+1-i; j++)
		{
			double h1[Np];
			double h2[Np];
			JacobiP(a, 0, 0, i, Np, &h1[0]);
			JacobiP(b, 2*i+1, 0, j, Np, &h2[0]);

			for (int k = 0; k < Np; k++)
			{
				double simplex2DP = sqrt(2.0)*h1[k]*h2[k]*pow(1-b[k],i);
				gsl_matrix_set(V2D, k, sk, simplex2DP);
			}
			sk++;
		}
	}
}

// Evaluate the derivatives of the modal basis (id, jd) on the 2D simplex at (a,b)
void GradSimplex2DP(double *a, double *b, int id, int jd, int Np, double *dmodedr, double *dmodeds)
{
	double fa[Np], dfa[Np], gb[Np], dgb[Np];

	JacobiP(a, 0, 0, id, Np, fa);
	GradJacobiP(a, 0, 0, id, Np, dfa);
	JacobiP(b, 2*id+1, 0, jd, Np, gb);
	GradJacobiP(b, 2*id+1, 0, jd, Np, dgb);

	// s-derivative
	for (int i = 0; i < Np; i++)
	{
		// r derivative
		dmodedr[i] = dfa[i]*gb[i];

		// s derivative
		dmodeds[i] = dfa[i]*(gb[i]*(0.5*(1+a[i])));

		double tmp = dgb[i]*pow(0.5*(1-b[i]), id);
		if (id > 0)
		{
			dmodedr[i] *= pow(0.5*(1-b[i]), id-1);
			dmodeds[i] *= pow(0.5*(1-b[i]), id-1);
			tmp -= 0.5*id*gb[i]*(pow(0.5*(1-b[i]), id-1));
		}
		dmodeds[i] += fa[i]*tmp;

		// normalize
		dmodedr[i] *= pow(2, id+0.5);
		dmodeds[i] *= pow(2, id+0.5);
	}

}

// Evaluates the gradient of the modal basis (i,j) at (r,s) at order N
void GradVandermonde2D(double *r, double *s, int Np, int N, gsl_matrix* Vr, gsl_matrix* Vs)
{

	// find tensor-product coordinates
	double a[Np], b[Np];
	rstoab(Np, r,s, a, b);

	double dmodedr[Np];
	double dmodeds[Np];

	int sk = 0;
	for (int i = 0; i < N+1; i++)
	{
		for (int j = 0; j < N+1-i; j++)
		{
			GradSimplex2DP(a, b, i, j, Np, dmodedr, dmodeds);
			for (int k = 0; k < Np; k++)
			{
				gsl_matrix_set(Vr, k, sk, dmodedr[k]);
				gsl_matrix_set(Vs, k, sk, dmodeds[k]);
			}
			sk++;
		}
	}

}

void Dmatrices2D(double *r, double *s, int N, int Np, gsl_matrix* V, gsl_matrix *Dr, gsl_matrix *Ds)
{
	int num_cols = (N+1)*(N+2)/2; // Np should equal num_cols
	gsl_matrix *Vr = gsl_matrix_alloc(Np, Np);
	gsl_matrix *Vs = gsl_matrix_alloc(Np, Np);
	GradVandermonde2D(r, s, Np, N, Vr, Vs);

	gsl_permutation *p = gsl_permutation_alloc(Np);
	gsl_matrix *Vcopy = gsl_matrix_alloc(Np, num_cols);
	gsl_matrix_memcpy(Vcopy, V);

	gsl_matrix *Vinv = gsl_matrix_alloc(Np, Np);

	int signum;
	gsl_linalg_LU_decomp(Vcopy, p, &signum);
	gsl_linalg_LU_invert(Vcopy, p, Vinv);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Vr, Vinv, 0.0, Dr);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Vs, Vinv, 0.0, Ds);

	gsl_permutation_free(p);
	gsl_matrix_free(Vcopy);
	gsl_matrix_free(Vinv);
	gsl_matrix_free(Vr);
	gsl_matrix_free(Vs);
}

// compute surface to volume lift terms for DG formulation and also populate Fmast matrix that contains indices of the edge points
void Lift2D(int N, int Np, gsl_matrix *V, double *r, double *s, int **Fmask, gsl_matrix *LIFT2D)
{
	int Nfp = N+1;

	// construct Fmask
	int f1 = 0, f2 = 0, f3 = 0;
	for (int i = 0; i < Np; i++)
	{
		if (fabs(s[i]+1) < 1e-8)
		{
			Fmask[f1][0] = i;
			f1++;
		}
		if (fabs(r[i]+s[i]) < 1e-8)
		{
			Fmask[f2][1] = i;
			f2++;
		}
		if (fabs(r[i]+1) < 1e-8)
		{
			Fmask[f3][2] = i;
			f3++;
		}
	}
	if (f1 > (Nfp) || f2 > (Nfp) || f3 > (Nfp))
		printf("Something is wrong. More than N+1 nodes found on an edge\n");

	gsl_matrix *Emat = gsl_matrix_calloc(Np, 3*Nfp);

	// calculate the coordinates on the standard triangle at edges
	double face1R[N+1], face2R[N+1], face3S[N+1];
	for (int i = 0; i < N+1; i++)
	{
		face1R[i] = r[Fmask[i][0]];
		face2R[i] = r[Fmask[i][1]];
		face3S[i] = s[Fmask[i][2]];
	}

	// edge 1
	gsl_matrix *V1D =gsl_matrix_alloc(N+1,N+1);
	Vandermonde1D(face1R, N+1, N, V1D);
	gsl_matrix *massEdge1 = gsl_matrix_alloc(N+1, N+1);
	calculateMassMatrix(N+1, V1D, massEdge1);	


	for (int i = 0; i < N+1; i++)
	{
		for (int j = 0; j < Nfp; j++)
		{
			gsl_matrix_set(Emat, Fmask[i][0], j, gsl_matrix_get(massEdge1, i, j));
		}
	}

	// edge 2
	Vandermonde1D(face2R, N+1, N, V1D);
	gsl_matrix *massEdge2 = gsl_matrix_alloc(N+1, N+1);
	calculateMassMatrix(N+1, V1D, massEdge2);	

	for (int i = 0; i < N+1; i++)
	{
		for (int j = 0; j < Nfp; j++)
		{
			gsl_matrix_set(Emat, Fmask[i][1], Nfp+j, gsl_matrix_get(massEdge2, i, j));
		}
	}

	// edge 3
	Vandermonde1D(face3S, N+1, N, V1D);
	gsl_matrix *massEdge3 = gsl_matrix_alloc(N+1, N+1);
	calculateMassMatrix(N+1, V1D, massEdge3);	

	for (int i = 0; i < N+1; i++)
	{
		for (int j = 0; j < Nfp; j++)
		{
			gsl_matrix_set(Emat, Fmask[i][2], 2*Nfp+j, gsl_matrix_get(massEdge3, i, j));
		}
	}

	// inv(mass matrix) *\I_n (L_i, L_j)_{edge_n}
	gsl_matrix *tmp = gsl_matrix_alloc(Np, 3*Nfp);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, V, Emat, 0.0, tmp);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, tmp, 0.0, LIFT2D);


	gsl_matrix_free(tmp);
	gsl_matrix_free(massEdge1);
	gsl_matrix_free(massEdge2);
	gsl_matrix_free(massEdge3);
	gsl_matrix_free(Emat);
	gsl_matrix_free(V1D);

}

void calculate_2D_volume_integral_operators(int Np, int N, double *r, double *s, gsl_matrix *V, gsl_matrix *Drw, gsl_matrix *Dsw)
{
	gsl_matrix *Vr = gsl_matrix_alloc(Np, Np);
	gsl_matrix *Vs = gsl_matrix_alloc(Np, Np);
	GradVandermonde2D(r, s, Np, N, Vr, Vs);

	gsl_matrix *M = gsl_matrix_alloc(Np,Np);
	calculateMassMatrix(Np, V, M);

	gsl_matrix *tmp = gsl_matrix_alloc(Np,Np);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V, Vr, 0.0, tmp);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp, M, 0.0, Drw);

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V, Vs, 0.0, tmp);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp, M, 0.0, Dsw);

	gsl_matrix_free(tmp);
	gsl_matrix_free(M);
	gsl_matrix_free(Vr);
	gsl_matrix_free(Vs);

}

void prepare_2D_structure(struct TwoDRegion* currDomain, int P, int Np, int Nfp, int** Fmask, gsl_matrix* Dr, gsl_matrix* Ds)
{
	int myNumEl = currDomain->NumEl;
	int NumEdges = currDomain->TotalNumEdges;

	currDomain->P = P;
	currDomain->Np = Np;

	currDomain->NodalX = malloc(myNumEl*sizeof(double*));
	currDomain->NodalY = malloc(myNumEl*sizeof(double*));
	currDomain->NodalZ = malloc(myNumEl*sizeof(double*));
	currDomain->NodalnFriction = malloc(myNumEl*sizeof(double*));
	currDomain->Nodaldzx = malloc(myNumEl*sizeof(double*));
	currDomain->Nodaldzy = malloc(myNumEl*sizeof(double*));

	for (int el = 0; el < myNumEl; el++)
	{
		currDomain->NodalX[el] = xcalloc(Np,sizeof(double));
		currDomain->NodalY[el] = xcalloc(Np,sizeof(double));
		currDomain->NodalZ[el] = xcalloc(Np,sizeof(double));
		currDomain->NodalnFriction[el] = xcalloc(Np,sizeof(double));
		currDomain->Nodaldzx[el] = xcalloc(Np, sizeof(double));
		currDomain->Nodaldzy[el] = xcalloc(Np, sizeof(double));
	}

	interpolateVertsToNodalCoordinates(P, myNumEl, currDomain->Vx, currDomain->EltoVert, currDomain->NodalX);
	interpolateVertsToNodalCoordinates(P, myNumEl, currDomain->Vy, currDomain->EltoVert, currDomain->NodalY);
	interpolateVertsToNodalCoordinates(P, myNumEl, currDomain->Vz, currDomain->EltoVert, currDomain->NodalZ);
	interpolateVertsToNodalCoordinates(P, myNumEl, currDomain->VnFriction, currDomain->EltoVert, currDomain->NodalnFriction);


	// calculate and store the jacobian, normal vectors and surface jacobian
	// Also calculate and store rx, ry, sx and sy
	currDomain->jac = xcalloc(myNumEl, sizeof(double));
	currDomain->edgJac = xcalloc(myNumEl*3, sizeof(double));
	currDomain->nx = xcalloc(myNumEl*3,sizeof(double));
	currDomain->ny = xcalloc(myNumEl*3,sizeof(double));
	currDomain->rx = xcalloc(myNumEl, sizeof(double));
	currDomain->ry = xcalloc(myNumEl, sizeof(double));
	currDomain->sx = xcalloc(myNumEl, sizeof(double));
	currDomain->sy = xcalloc(myNumEl, sizeof(double));

	for (int el = 0; el < myNumEl; el++)
	{
		int v1 = currDomain->EltoVert[el*3];
		int v2 = currDomain->EltoVert[el*3+1];
		int v3 = currDomain->EltoVert[el*3+2];
		//printf("el = %d, v1 = %d, v2 = %d, v3 = %d\n", el, v1, v2, v3);
		double x1 = currDomain->Vx[v1];
		double x2 = currDomain->Vx[v2];
		double x3 = currDomain->Vx[v3];
		double y1 = currDomain->Vy[v1];
		double y2 = currDomain->Vy[v2];
		double y3 = currDomain->Vy[v3];

		double xr = 0.5*(x2-x1);
		double yr = 0.5*(y2-y1);
		double xs = 0.5*(x3-x1);
		double ys = 0.5*(y3-y1);

		double jac = xr*ys - xs*yr;
		currDomain->jac[el] = jac;
		if ( jac < 0)
		{
			printf("The jacobian is negative - %lf. Check the vertex numbering for element %d.\n", jac, el);
			exit(1);
		}

		currDomain->rx[el] = ys/jac;
		currDomain->ry[el] = -xs/jac;
		currDomain->sx[el] = -yr/jac;
		currDomain->sy[el] = xr/jac;


		double nx[3], ny[3];
		nx[0] = yr; 
		ny[0] = -xr;
		nx[1] = ys - yr;
		ny[1] = -xs + xr;
		nx[2] = -ys;
		ny[2] = xs;

		for (int edg = 0; edg < 3; edg++)
		{
			double sj = sqrt(nx[edg]*nx[edg] + ny[edg]*ny[edg]);
			currDomain->nx[el*3+edg] = nx[edg]/sj;
			currDomain->ny[el*3+edg] = ny[edg]/sj;
			currDomain->edgJac[el*3+edg] = sj;
		}

		// calculate dz
		gsl_vector *myz = gsl_vector_alloc(Np);

		for (int k = 0; k < Np; k++)
			gsl_vector_set(myz, k, currDomain->NodalZ[el][k]);

		gsl_vector *dzrx = gsl_vector_alloc(Np);
		gsl_vector *dzsx = gsl_vector_alloc(Np);
		gsl_vector *dzry = gsl_vector_alloc(Np);
		gsl_vector *dzsy = gsl_vector_alloc(Np);

		gsl_blas_dgemv(CblasNoTrans, currDomain->rx[el], Dr, myz, 0.0, dzrx);
		gsl_blas_dgemv(CblasNoTrans, currDomain->sx[el], Ds, myz, 0.0, dzsx);
		gsl_blas_dgemv(CblasNoTrans, currDomain->ry[el], Dr, myz, 0.0, dzry);
		gsl_blas_dgemv(CblasNoTrans, currDomain->sy[el], Ds, myz, 0.0, dzsy);


		for (int k = 0; k < Np; k++)
		{
			currDomain->Nodaldzx[el][k] = gsl_vector_get(dzrx,k) + gsl_vector_get(dzsx,k);
			currDomain->Nodaldzy[el][k] = gsl_vector_get(dzry,k) + gsl_vector_get(dzsy,k);
		}

		gsl_vector_free(dzrx);
		gsl_vector_free(dzsx);
		gsl_vector_free(dzry);
		gsl_vector_free(dzsy);

	}

	// store GlobalEdgPosNegnodes
	currDomain->GlobalEdgPosNegNodes = malloc(NumEdges*sizeof(int*));
	currDomain->PosInFVec = malloc(NumEdges*sizeof(int*));
	for (int edg = 0; edg < NumEdges; edg++)
	{
		currDomain->GlobalEdgPosNegNodes[edg] = xcalloc(Nfp*2, sizeof(int));
		currDomain->PosInFVec[edg] = xcalloc(Nfp*2, sizeof(int));
		int el1 = currDomain->EdgtoEls[edg*2];
		int el2 = currDomain->EdgtoEls[edg*2+1];
		int local_edg1 = currDomain->GlobaltoLocalEdg[edg*2];
		int local_edg2 = currDomain->GlobaltoLocalEdg[edg*2+1];
		double x1[Nfp], x2[Nfp], y1[Nfp], y2[Nfp];
		for (int s = 0; s < Nfp; s++)
		{
			int lv1 = Fmask[s][local_edg1];
			int lv2 = Fmask[s][local_edg2];
			x1[s] = currDomain->NodalX[el1][lv1];
			y1[s] = currDomain->NodalY[el1][lv1];
			x2[s] = currDomain->NodalX[el2][lv2];
			y2[s] = currDomain->NodalY[el2][lv2];
		}
		double dis[Nfp][Nfp];
		for (int s = 0; s	< Nfp; s++)
		{
			double tol = 1e-8;
			int found = 0;
			for(int k = 0; k < Nfp; k++)
			{
				dis[s][k] = sqrt((x1[s] - x2[k])*(x1[s]-x2[k]) + (y1[s]-y2[k])*(y1[s]-y2[k]));
				if(dis[s][k] < tol)
				{
					currDomain->PosInFVec[edg][2*s] = Nfp*local_edg1+s;
					currDomain->PosInFVec[edg][2*s+1] = Nfp*local_edg2+k;
					currDomain->GlobalEdgPosNegNodes[edg][2*s] = Fmask[s][local_edg1];
					currDomain->GlobalEdgPosNegNodes[edg][2*s+1] = Fmask[k][local_edg2];
					found = 1;
					break;
				}
			}
			if (!found)
				printf("Unable to find corresponding positive node\n");
		}

	}

	// create VtoNode
	currDomain->VtoNode = malloc(currDomain->NumEl*sizeof(int*));
	for (int el = 0; el < currDomain->NumEl; el++)
	{
		currDomain->VtoNode[el] = xcalloc(3, sizeof(int));
		for (int i = 0; i < 3; i++)
		{
			int v = currDomain->EltoVert[el*3+i];
			double v_x = currDomain->Vx[v];
			double v_y = currDomain->Vy[v];
			for (int j = 0; j < Np; j++)
			{
				double n_x = currDomain->NodalX[el][j];
				double n_y = currDomain->NodalY[el][j];
				double dist = sqrt((n_x - v_x)*(n_x -v_x) + (n_y - v_y)*(n_y - v_y));
				if (dist < 1e-8)
				{
					currDomain->VtoNode[el][i] = j;
					break;
				}
			}
		}
	}



}

void setup_2D_domains()
{
	int P = 1;
	int Np = (P+1)*(P+2)/2;
	int Nfp = P+1;

	// Compute nodal set
	double xeq[Np], yeq[Np] ,r[Np], s[Np];
	getEquiNodes(P, xeq, yeq);
	xytors(Np, xeq, yeq, r, s);

	// Compute Vandermonde matrix
	gsl_matrix *V2D = gsl_matrix_alloc(Np, Np);
	Vandermonde2D(r, s, Np, P, V2D);

	// calculate 2D mass matrix
	MassMatrix2D = gsl_matrix_alloc(Np, Np);
	calculateMassMatrix(Np, V2D, MassMatrix2D);

	
	// calculate LIFT and Volume integral matrices
	Fmask = malloc((Nfp)*sizeof(int*));
	for (int i = 0; i < P+1; i++)
		Fmask[i] = malloc(3*sizeof(int));
	LIFT2D = gsl_matrix_alloc(Np, 3*Nfp);
	Lift2D(P, Np, V2D, r, s, Fmask, LIFT2D);


	// calculate Dr and Ds to calculate dz
	gsl_matrix *Dr = gsl_matrix_alloc(Np,Np);
	gsl_matrix *Ds = gsl_matrix_alloc(Np,Np);
	Dmatrices2D(r, s, P, Np, V2D, Dr, Ds);

	Drw = gsl_matrix_alloc(Np, Np);
	Dsw = gsl_matrix_alloc(Np,Np);
	calculate_2D_volume_integral_operators(Np, P, r, s, V2D, Drw, Dsw);

	gsl_matrix_free(V2D);

	for (int i = 0; i < NumJunctions; i++)
	{
		struct TwoDRegion *junc = JunctionList[i];
		prepare_2D_structure(junc, P, Np, Nfp, Fmask, Dr, Ds);
	}

	for (int i = 0; i < NumFloodplains; i++)
	{
		struct TwoDRegion *floodplain = FloodplainList[i];
		prepare_2D_structure(floodplain, P, Np, Nfp, Fmask, Dr, Ds);

	}

}

