/************************************************************************************//**
*
* @file oneTimeStep.h
* 
* This file contains function prototypes for the functions that are exectued in order
* to step forward in time.
* ***************************************************************************************/ 


#ifndef ONE_TIME_STEP

#define ONE_TIME_STEP

#include "Channels_KinRoutes.h"

/* @cond FUNCTION_PROTOTPYES */
extern void computeChanL(struct channel* Chan, double time, int channelNumber, double* RHSA, double* RHSQ);
extern void minmodChan(struct channel* Chan);
extern void minmodHeightChan(struct channel* Chan);
extern void wetDryStatus1D(struct channel *Chan);
extern void PDop1D(struct channel *Chan);

extern void computeKinL(struct kinematicRoute* KinRoute, double time, int kinRoutenumber, double dt, double* RHS);
extern void minmodKinRoute(struct kinematicRoute* KinRoute);

extern void coupleKinRoutesWithChannels(double dt, double time);
extern void boundary_conditions(double time);
extern void outputDataToFile(double time);
extern double calculateTotalWater(struct channel *Chan);

/* @uncond */
#endif
