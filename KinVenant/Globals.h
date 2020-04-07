#ifndef GLOBALS
#define GLOBALS

#include <gsl/gsl_matrix.h>
#include "Channels_KinRoutes.h"

extern int NumChannels;
extern int NumKinRoutes;
extern struct channel** ChannelList;
extern struct kinematicRoute** KinRouteList;

extern gsl_matrix *LIFT;
extern gsl_matrix *VolMat;
extern gsl_matrix *MassMatrix;
extern gsl_matrix *InvM;

#endif
