#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>

#include "globals.h"

double log2(double x)
{
	return (log(x)/log(2));
}

void calculateL2Error()
{
	double error_sq = 0.0;
	gsl_vector *diff = gsl_vector_alloc(Np);
	for (int i = 0; i < NumEl; i++)
	{
		double el_error;
		for (int j = 0; j < Np; j++)
		{
			double currX = NodalX[i*Np+j];
			//double TrueSol = 2;
			double TrueSol = sin(currX)+2;
			double numSol = H[i*Np+j+1];
			printf("currX = %3.14f \t TrueSol = %3.14f \t NumSol = %3.14f \n", currX, TrueSol, numSol);
			gsl_vector_set(diff, j, TrueSol-numSol);
		}
	gsl_vector *tmp = gsl_vector_alloc(Np);
	gsl_blas_dgemv(CblasNoTrans, 0.5*dh[i], MassMatrix, diff, 0.0, tmp);
	gsl_blas_ddot(diff, tmp, &el_error);
	
	error_sq += el_error;
	}

	double error = sqrt(error_sq);
	printf("%d \t %3.14f\n", NumEl, error);
	printf("%3.14f \t %3.14f\n", log2(NumEl), log2(error));
	gsl_vector_free(diff);

}


