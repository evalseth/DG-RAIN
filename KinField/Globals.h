#ifndef GLOBALS

#define GLOBALS

#include <gsl/gsl_matrix.h>
#include "Watershed.h"

extern const double g;

extern int *Valence;				 // valence of each kinematic element

extern int NumChannels; 
extern int NumKinematicEls;
extern int NumJunctions;
extern int NumFloodplains;

extern struct kinematicEl **KinematicElList;
extern struct channel **ChannelList;
extern struct TwoDRegion **JunctionList;
extern struct TwoDRegion **FloodplainList;

extern gsl_matrix *LIFT;
extern gsl_matrix *VolMat;
extern gsl_matrix *MassMatrix;

extern int **Fmask;
extern gsl_matrix *LIFT2D;
extern gsl_matrix *Drw;
extern gsl_matrix *Dsw;
extern gsl_matrix *MassMatrix2D;

#ifdef WDON
extern const double H0;
extern const double VELZERO;
#endif

extern void *xcalloc(int items, int size);

#endif
