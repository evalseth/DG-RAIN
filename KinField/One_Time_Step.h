/************************************************************************************//**
*
* @file oneTimeStep.h
* 
* This file contains function prototypes for the functions that are exectued in order
* to step forward in time.
* ***************************************************************************************/ 


#ifndef ONE_TIME_STEP

#define ONE_TIME_STEP

#include "Watershed.h"

/* @cond FUNCTION_PROTOTPYES */
extern void computeLKinematicEls(double time, double*RHS, int fp);
extern void minmodKinematicEls(int fp);


extern void computeChanL(struct channel* Chan, double time, double* dt, int channelNumber, int stage, double* RHSA, double* RHSQ);
extern void minmodChan(struct channel* Chan);
extern void minmodHeightChan(struct channel* Chan);
#ifdef WDON
extern void wetDryStatusChan(struct channel *Chan);
extern void PDopChan(struct channel *Chan);
#endif

extern void compute2DL(struct TwoDRegion* junc, double time, double* RHSZeta, double* RHSQx, double* RHSQy);
extern void SlopeLimiterJunc(struct TwoDRegion* junc);
#ifdef WDON
extern void wetDryStatusJunc(struct TwoDRegion* junc);
extern void PDopJunc(struct TwoDRegion* junc);
#endif

extern void couple_kinEls_with_channels(double time, double dt);
extern void couple_channels_with_junctions();

extern void boundary_conditions(double time);


extern void outputFloodplainDataToFile(double time, double finalTime, double recordTimeIntervals);
extern void outputFloodplainNodalDataToFile(double time, double finalTime, double recordTimeIntervals);
extern void outputDataToFile(double time, double finalTime, double recordTimeIntervals);
extern double calculateTotalWater(struct channel *Chan);
extern void append_to_file(char* fileName, double time, double data);
extern void output2DNodalError(double time);
extern void calculate_l2_error_2d(double time);

/* @uncond */
#endif
