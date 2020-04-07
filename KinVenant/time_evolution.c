/************************************************************************//**
 * @file time_evolution.c
 *
 * This file contains code to advance the solutions one step forward in time
 *
 * *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "Channels_KinRoutes.h"
#include "Globals.h"
#include "Math_Functions.h"
#include "One_Time_Step.h"
#include "Constitutive_Equations.h"

/* @cond GLOBAL_VARIABLES */
extern const double g;
/* @endcond */

extern void* xcalloc(int NumEl, int size);

/***************************************************************************//**
 *
 * Function to advance a channel structure one step forward in time using 2-stage
 * RKDG scheme
 * @param[in] Chan channel structure that we are currently working on
 * @param[in] time current time of the simulation
 * @param[in] dt size of the time step
 * @param[in] channelNumber identity (number) of the channel structure that we
 * are currently working on
 *
 * *********************************************************************************/
void oneTimeStepChannel(struct channel *Chan, double time, double dt, int channelNumber)
{
	int TotalNumNodes = Chan->NumNodes+2;
	double* RHSA = xcalloc(TotalNumNodes, sizeof(double));
	double* RHSQ = xcalloc(TotalNumNodes, sizeof(double));

	// Calculate wet/dry status of elements
#ifdef WDON
	wetDryStatus1D(Chan);
#endif

	// Intermediate RK step
	computeChanL(Chan, time, channelNumber, RHSA, RHSQ);


	double* oldA = xcalloc(TotalNumNodes, sizeof(double));
	double* oldQ = xcalloc(TotalNumNodes, sizeof(double));
	memcpy(oldA, Chan->A, TotalNumNodes*sizeof(double));
	memcpy(oldQ, Chan->Q, TotalNumNodes*sizeof(double));

	// w_1 = w + dt*L 
	for (int i = 0; i < TotalNumNodes; i++)
	{
		Chan->A[i] = oldA[i]+dt*RHSA[i];
		Chan->Q[i] = oldQ[i]+dt*RHSQ[i];

	}

/*	int NumNodes = Chan->NumNodes;
	if (time > 49067.2)
	{
		for (int i = 0; i < NumNodes+2; i++)
			printf("oldA = %lf \t oldQ = %lf \n", oldA[i], oldQ[i]);
	}
	if (time > 49067.2)
	{
		for (int i = 0; i < NumNodes+2; i++)
			printf("RHSA = %lf \t RHSQ = %lf \n", RHSA[i], RHSQ[i]);
	}
*/

	// Apply minmod slope limiter on w_1
	minmodHeightChan(Chan);

/*	if (time > 49067.2)
	{
		for (int i = 0; i < NumNodes+2; i++)
			printf("A = %lf \t Q = %lf \n", Chan->A[i], Chan->Q[i]);
	}
*/

	boundary_conditions(time);

#ifdef WDON
	// apply the positive-depth operaton on w_1
	PDop1D(Chan);

	// Caculate wet/dry status again
	wetDryStatus1D(Chan);

#endif 

	//	internal_BC();

	// Compute w
	computeChanL(Chan, time, channelNumber, RHSA, RHSQ);

	// w_1 = (w + w_1 + dt*L(w_1))/2
	for (int i = 0; i < TotalNumNodes; i++)
	{
		Chan->A[i] = 0.5*(oldA[i]+Chan->A[i]+dt*RHSA[i]);
		Chan->Q[i] = 0.5*(oldQ[i]+Chan->Q[i]+dt*RHSQ[i]);

	}

	// Apply minmod slope limiter to w
	minmodHeightChan(Chan);

#ifdef WDON
	PDop1D(Chan);
#endif

	free(RHSA);
	free(RHSQ);
	free(oldA);
	free(oldQ);

	boundary_conditions(time);

}

/***************************************************************************//**
 *
 * Function to advance a kinematicRoute structure one step forward in time using 2-stage
 * RKDG scheme
 * @param[in] kinRoute the kinematicRoute structure that we are currently working on
 * @param[in] time current time of the simulation
 * @param[in] dt size of the time step
 *
 * *********************************************************************************/

void oneTimeStepKinRoute(struct kinematicRoute *kinRoute, double time, double dt, int kinRouteNumber)
{
	int TotalNumNodes = kinRoute->NumNodes+2;
	double* oldH = xcalloc(TotalNumNodes, sizeof(double));
	memcpy(oldH, kinRoute->H, TotalNumNodes*sizeof(double));

	double* RHS = xcalloc(TotalNumNodes, sizeof(double));

	// Intermediate RK step
	computeKinL(kinRoute, time, kinRouteNumber, dt, RHS);

	//printf("before minmod\n");
	for(int i = 0; i < TotalNumNodes; i++)
	{
		kinRoute->H[i] = kinRoute->H[i] + dt*RHS[i];
		//printf("H = %3.13f\n", kinRoute->H[i]);
	}

	//printf("after minmod\n");
	// slope limiter
	minmodKinRoute(kinRoute);

	//for (int i = 0; i < TotalNumNodes; ++i)
	//{
	//	printf("H = %3.13f\n", kinRoute->H[i]);
	//}

	// Second RK Stage
	computeKinL(kinRoute, time, kinRouteNumber, dt, RHS);

	//printf("after second L\n");
	for(int i = 0; i < TotalNumNodes; i++)
	{
		kinRoute->H[i] = 0.5*(oldH[i] + kinRoute->H[i] + dt*RHS[i]);

	//	printf("H = %3.13f\n", kinRoute->H[i]);
	}

	// slope limiter
	minmodKinRoute(kinRoute);

	//printf("after second minmod\n");
	//for (int i = 0; i < TotalNumNodes; ++i)
	//{
	//	printf("H = %3.13f\n", kinRoute->H[i]);
	//}


}

/*******************************************************************//**
 *
 * Function to evolve channels as well as kinRoutes through time
 * until the final simulation time. It also outputs the data at a certain
 * time
 *
 * *********************************************************************/
void time_evolution(double FinalTime)
{

	double time = 0;
	//double dt = 1e-5;
	double dt = 0.05;
	int Nstep = 0;

	outputDataToFile(time);

	while (time < FinalTime)
	{
		time = time + dt;
		Nstep++;
		
		double maxLambda = 0;

		// step forward in time
		for(int i=0; i<NumKinRoutes; ++i)
		{
			oneTimeStepKinRoute(KinRouteList[i],time,dt,i);
			maxLambda = fmax(maxLambda, KinRouteList[i]->max_lambda);
		}

		coupleKinRoutesWithChannels(dt,time);
		
		for(int i=0; i<NumChannels; ++i)
		{
			oneTimeStepChannel(ChannelList[i],time,dt,i);
			maxLambda = fmax(maxLambda, ChannelList[i]->max_lambda);
		}

		if (Nstep % 5000 == 0 || time == FinalTime)
		{
			outputDataToFile(time);
		
			double totalWater = calculateTotalWater(ChannelList[0]);

			printf("Total Water at time %lf  = %lf\n", time, totalWater);
		}

		double minLength = 99999999;	
		for (int i = 0; i < NumChannels; ++i)
		{
			minLength = fmin(minLength, ChannelList[i]->mindh);
		}

		for (int i = 0; i < NumKinRoutes; ++i)
		{
			minLength = fmin(minLength, KinRouteList[i]->mindh);
		}



		boundary_conditions(time);

		/**************** Set next time step ****************/
		// set next time step
		//printf("max_lambda = %e\n", maxlambda);
		//printf("mindx = %e\n", mindx);
	
		//printf("Step size = %lf\n", dt);
		//dt = 0.1*minLength/maxLambda;
		//if (time < 60 && (time + dt) > 60)
		//	dt = 60 - time;
		//if ((time + dt) > FinalTime)
		//	dt = FinalTime - time;
		//printf("Completed time = %lf\n", time);

		
	} // end while loop

	printf("collectedWater = %lf\n", ChannelList[0]->collectedQL);

}
