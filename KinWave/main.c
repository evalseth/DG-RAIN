#include <stdio.h>
#include<gsl/gsl_matrix.h>
#include "SimulationSteps.h"
#include "globals.h"


extern void GetBathymetryData(double* x, double* z, int NumEdges);

double* X;
double* NodalX;
double* dh;
double *z;
double *S0;
double *S0Deriv;
double *n;
int *EdgToNode;
int NumEl;
int NumEdges;
int TotalNumNodes;
int P;
int Np;
gsl_matrix *LIFT;
gsl_matrix *VolMat;
gsl_matrix *MassMatrix;

void initialize_grid()
{
	double x_0 = 0.0; 
	double x_1 = 2*PI;

	NumEdges = NumEl+1;
	EdgToNode = malloc(NumEdges*sizeof(int));
	
	P = 2;
	Np = P+1;
	TotalNumNodes = Np*NumEl+2;

	double h_el = (x_1-x_0)/NumEl;

	X = malloc((NumEdges)*sizeof(double));
	dh = malloc(NumEl*sizeof(double));
	z = malloc(NumEdges*sizeof(double));
	S0 = malloc((TotalNumNodes-2)*sizeof(double));
	S0Deriv = malloc((TotalNumNodes-2)*sizeof(double));
	n = malloc(NumEdges*sizeof(double));
	X[0] = x_0;
	X[NumEl] = x_1;
	dh[0] = h_el;
	n[0] = 0.3;

	EdgToNode[0] = 0;
	for (int i = 1; i < NumEl; i++)
	{
		X[i] = x_0 + h_el*i;
		dh[i] = h_el;
		EdgToNode[i] = EdgToNode[i-1]+Np;
		n[i] = 0.3;

	}
	EdgToNode[NumEl] = EdgToNode[NumEl-1]+Np;
	n[NumEl] = 0.3;

	NodalX = malloc(Np*NumEl*sizeof(double));
	CalculateNodalCoordinates(P, NumEl, X, NodalX);
	GetBathymetryData(X, z, NumEl+1);

	int NumPointsToInterpolate = 9;
	InterpolateWithSplines(NumPointsToInterpolate, TotalNumNodes-2, X, NodalX, z, S0, S0Deriv);
	//InterpolateWithSplines(NumEl+1, TotalNumNodes-2, X, NodalX, z, S0, S0Deriv);

	
/*	for (int i = 0; i < 100; i++)
	{
		double xval1 = 0 + i*(2*PI)/99;
		printf("%3.14f \t %3.14f\n", xval1, S0[i]);
	}
*/
	//for (int i = 0; i < TotalNumNodes-2; i++)
		//printf("%3.14f \t %3.14f \t %3.14f\n", X[i], -z[i], S0[i]);

}


int main(int argc, char **argv)
{
	printf("Enter the number of elements:\n");
	scanf("%d", &NumEl);

	initialize_grid();

	double FinalTime =2*PI;

	calculateLIFTVolMat(P, Np, &LIFT, &VolMat, &MassMatrix);
	
	initializeFloodPlains();

	time_evolution(FinalTime);

	//calculateL2Error();

	for(int i =0; i <TotalNumNodes; i++)
		printf("%3.14f\n", H[i]);
	
}
