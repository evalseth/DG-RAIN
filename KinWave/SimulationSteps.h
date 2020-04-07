#ifndef SIMULATION_STEPS

#define SIMULATION_STEPS

#include <gsl/gsl_matrix.h>

#define PI 3.1415926535

extern void CalculateNodalCoordinates(int N, int NumEl, double* X, double* Xnodes);
extern void calculateLIFTVolMat(int N, int Np, gsl_matrix **LIFT, gsl_matrix **VolMat, gsl_matrix** MassMatrix);
extern void initializeFloodPlains();
extern void time_evolution(double FinalTime);
extern void calculateL2Error();
extern void InterpolateWithSplines(int NumPointsToInterpolate, int NumPointsToEvaluate, double* x, double* NodalX, double* z, double* S0, double* S0Deriv);

#endif
