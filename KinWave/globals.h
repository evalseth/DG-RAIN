#ifndef GLOBALSi

#define GLOBALS

#include <gsl/gsl_matrix.h>

extern double *H;
extern double *KinRHS;
extern double *X;
extern double *NodalX;
extern double *n;
extern double *z;
extern double *S0;
extern double *S0Deriv;
extern double* dh;
extern int *EdgToNode; 	// element i contains the left node number corresponding to edge i
extern int NumEl;
extern int P;
extern int Np;
extern int NumEdges;
extern int TotalNumNodes;
extern gsl_matrix *LIFT;
extern gsl_matrix *VolMat;
extern gsl_matrix *MassMatrix;

extern double dt;

#endif
