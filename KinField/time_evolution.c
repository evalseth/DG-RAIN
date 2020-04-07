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
#include "Watershed.h"
#include "Globals.h"
#include "One_Time_Step.h"
#include "Constitutive_Equations.h"
#include "Math_Functions.h"


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
	wetDryStatusChan(Chan);
#endif

	// Intermediate RK step
	computeChanL(Chan, time, &dt, channelNumber, 0 ,RHSA, RHSQ);

	double* oldA = xcalloc(TotalNumNodes, sizeof(double));
	double* oldQ = xcalloc(TotalNumNodes, sizeof(double));
	memcpy(oldA, Chan->A, TotalNumNodes*sizeof(double));
	memcpy(oldQ, Chan->Q, TotalNumNodes*sizeof(double));

	// w_1 = w + dt*L 
	for (int i = 0; i < TotalNumNodes; i++)
	{
		Chan->A[i] = oldA[i]+dt*RHSA[i];
		Chan->Q[i] = oldQ[i]+dt*RHSQ[i];
		//if (channelNumber == 5)
		//{
		//	printf("after first computeL , RHSA = %lf A = %lf \n", RHSA[i], Chan->A[i]);
		//	printf("after first computeL , RHSQ = %lf Q = %lf \n", RHSQ[i], Chan->Q[i]);
		//}

	}

	// Apply minmod slope limiter on w_1
	//minmodHeightChan(Chan);
	minmodChan(Chan);
	
	//if (channelNumber == 5)
	//{
	//	for (int i = 0; i < TotalNumNodes; i++)
	//	{
	//		printf("after first minmod, A = %lf \n", Chan->A[i]);
	//		printf("after first minmod, Q = %lf \n", Chan->Q[i]);
	//	}
	//}


	//boundary_conditions(time);

#ifdef WDON
	// apply the positive-depth operaton on w_1
	PDopChan(Chan);

	// Caculate wet/dry status again
	wetDryStatusChan(Chan);

#endif 

	// Compute w
	computeChanL(Chan, time, &dt, channelNumber, 1,RHSA, RHSQ);

	// w_1 = (w + w_1 + dt*L(w_1))/2
	for (int i = 0; i < TotalNumNodes; i++)
	{
		Chan->A[i] = 0.5*(oldA[i]+Chan->A[i]+dt*RHSA[i]);
		Chan->Q[i] = 0.5*(oldQ[i]+Chan->Q[i]+dt*RHSQ[i]);
		//if (channelNumber == 5)
		//{
		//	printf("after second computeL , RHSA = %lf A = %lf \n", RHSA[i], Chan->A[i]);
		//	printf("after second computeL , RHSQ = %lf Q = %lf \n", RHSQ[i], Chan->Q[i]);
		//}

	}

	// Apply minmod slope limiter to w
	//minmodHeightChan(Chan);
	minmodChan(Chan);
	
	//if (channelNumber == 5)
	//{
	//	for (int i = 0; i < TotalNumNodes; i++)
	//	{
	//		printf("after second minmod, A = %lf \n", Chan->A[i]);
	//		printf("after second minmod, Q = %lf \n", Chan->Q[i]);
	//	}
	//}

	//for (int i = 0; i < TotalNumNodes; i++)
	//	printf("A = %lf \n", Chan->A[i]);
#ifdef WDON
	PDopChan(Chan);
#endif

	free(RHSA);
	free(RHSQ);
	free(oldA);
	free(oldQ);

	//if (channelNumber == 5)
	//	exit(1);
	//boundary_conditions(time);

}

/***************************************************************************//**
*
* Routine to advance a junction structure one step forward in time using 2-stage
* RKDG scheme
* @param[in] junc the junction tructure that we are currently working on
* @param[in] time current time of the simulation
* @param[in] dt size of the time step
*
* *********************************************************************************/

void oneTimeStep2D(struct TwoDRegion *junc, double time, double dt)
{

	int NumEl = junc->NumEl;
	int Np = junc->Np;
	int NumNodes= Np*NumEl;
	double *RHSZeta = xcalloc(NumNodes, sizeof(double));
	double *RHSQx = xcalloc(NumNodes, sizeof(double));
	double *RHSQy = xcalloc(NumNodes, sizeof(double));

	#ifdef WDON
	wetDryStatusJunc(junc);
	#endif

	// Intermediate RK step
	//printf("Intermediate RK step\n");
	compute2DL(junc, time, RHSZeta, RHSQx, RHSQy);

	double *oldZeta = xcalloc(NumNodes, sizeof(double));
	double *oldQx = xcalloc(NumNodes, sizeof(double));
	double *oldQy = xcalloc(NumNodes, sizeof(double));

	for (int i = 0; i < NumEl; i++)
	{
		int begNode = i*Np;
		for (int j = 0; j < Np; j++)
		{
			oldZeta[begNode+j] = junc->zeta[i][j];
			oldQx[begNode+j] = junc->Qx[i][j];
			oldQy[begNode+j] = junc->Qy[i][j];
		}
		
	}

	// w_1 = w + dt*L
	for (int i = 0; i < NumEl; i++)
	{
		int begNode = i*Np;
		for (int j = 0; j < Np; j++)
		{
			junc->zeta[i][j] += dt*RHSZeta[begNode+j];
			junc->Qx[i][j] += dt*RHSQx[begNode+j];
			junc->Qy[i][j] += dt*RHSQy[begNode+j];
			//printf("juncZeta = %lf\n", junc->zeta[i][j]);
		}
	}

	// Apply slope limiter on w_1
	//SlopeLimiterJunc(junc);
	//SlopeLimiterHeightJunc(junc);

	#ifdef WDON
	PDopJunc(junc);
	wetDryStatusJunc(junc);
	#endif

	// Compute w
	//printf("Second RK Step\n");
	compute2DL(junc, time, RHSZeta, RHSQx, RHSQy);

	// w_1 = (w + w_1 + dt*L(w_1))/2;
	for (int i = 0; i < NumEl; i++)
	{
		int begNode = i*Np;
		for (int j = 0; j < Np; j++)
		{
			junc->zeta[i][j] = 0.5*(oldZeta[begNode+j] + junc->zeta[i][j] + dt*RHSZeta[begNode+j]);
			junc->Qx[i][j] = 0.5*(oldQx[begNode+j] + junc->Qx[i][j] + dt*RHSQx[begNode+j]);
			junc->Qy[i][j] = 0.5*(oldQy[begNode+j] + junc->Qy[i][j] + dt*RHSQy[begNode+j]);
			//printf("juncZeta = %lf \n", junc->zeta[i][j]);
		}
	}

	// Apply slope limiter to w
	//printf("After slope limiting \n");
	//SlopeLimiterJunc(junc);
	//SlopeLimiterHeightJunc(junc);

/*	for (int i = 0; i < NumEl; i++)
	{
		int begNode = i*Np;
		for (int j = 0; j < Np; j++)
		{

			printf("juncZeta = %lf \n", junc->zeta[i][j]);
		}
	}
*/
	#ifdef WDON
	PDopJunc(junc);
	#endif

	free(RHSZeta);
	free(RHSQx);
	free(RHSQy);
	free(oldZeta);
	free(oldQx);
	free(oldQy);
	
}



/***************************************************************************//**
 *
 * Function to advance all of the kinematic elements one step forward in time 
 * using 2-stage RKDG scheme 
 * @param[in] time current time of the simulation
 * @param[in] dt size of the time step
 *
 * *********************************************************************************/
double oneTimeStepKinematicEls(double time, double dt, int fp)
{
	int NumEl = FloodplainList[fp]->NumEl;
	int Np = KinematicElList[0]->Np;
	int TotalNumNodes = Np*NumEl ;
	double* RHS = xcalloc(TotalNumNodes, sizeof(double));

	// Intermediate RK step
	computeLKinematicEls(time, RHS, fp);
	
	double* oldA = xcalloc(TotalNumNodes, sizeof(double));

	// w_1 = w+dt*L;
	//printf("before minmod\n");
	for (int i = 0; i < NumEl; ++i)
	{
		struct kinematicEl* kinEl = KinematicElList[i];
		if (kinEl->isActive)
		{
			for (int j = 0; j < Np; ++j) 
			{
				oldA[i*Np+j] = kinEl->A[j];
				kinEl->A[j] = kinEl->A[j] + dt*RHS[i*Np+j];
				//printf("H = %3.13f\n", kinEl->A[j]/kinEl->weq);

			}
		}
	}

	// Apply minmod slope limiter on w_1
	minmodKinematicEls(fp);

	//printf("After minmod\n");
	//for (int i = 0; i < NumEl; ++i)
	//{
	//	struct kinematicEl* kinEl = KinematicElList[i];
	//	if (kinEl->isActive)
	//	{
	//		int Np = kinEl->Np;
	//		for (int j = 0; j < Np; ++j)
	//		{
	//			printf("H = %3.13f\n", kinEl->A[j]/kinEl->weq);
	//		}
	//		
	//	}
	//}


	// Second RK Stage
	computeLKinematicEls(time, RHS, fp);

	//printf("After second L\n");
	// w_1 = (w + w_1 + dt*L(w_1))/2
	for (int i = 0; i < NumEl; ++i)
	{
		struct kinematicEl* kinEl = KinematicElList[i];
		if (kinEl->isActive)
		{
			for (int j = 0; j < Np; ++j)
			{
				kinEl->A[j] = 0.5*(oldA[i*Np+j] + kinEl->A[j] + dt*RHS[i*Np+j]);

				//printf("H = %3.13f\n", kinEl->A[j]/kinEl->weq);
			}
		}
	}

	//printf("After second minmod\n");
	// Apply minmod slope limiter to w
	minmodKinematicEls(fp);
	//for (int i = 0; i < NumEl; ++i)
	//{
	//	struct kinematicEl* kinEl = KinematicElList[i];
	//	if (kinEl->isActive)
	//	{
	//		int Np = kinEl->Np;
	//		for (int j = 0; j < Np; ++j)
	//		{
	//			printf("H = %3.13f\n", kinEl->A[j]/kinEl->weq);
	//		}
	//		
	//	}
	//}


	double max_lambda = 0;
	for (int i = 0; i < NumEl; ++i)
	{
		struct kinematicEl* kinEl = KinematicElList[i];
		int Np = kinEl->Np;
		for (int j = 0; j < Np; ++j)
		{
			double A = kinEl->A[j];
			double S0 = kinEl->dz[j];
			double nf = kinEl->NodalnFriction[j];
			double weq = kinEl->weq;
			double u = sqrt(S0)/nf *pow(A/weq, 2.0/3);
			max_lambda = fmax(max_lambda, fabs(u));
		}
	}

	free(RHS);
	free(oldA);

	//boundary_conditions(time);
	
	return max_lambda;

}


/*******************************************************************//**
 *
 * Function to evolve channels, junctions and kinematic elements through time
 * until the final simulation time. It also outputs the data at certain
 * times.
 *
 * *********************************************************************/
void time_evolution(double FinalTime)
{

	double time = 0;
	double dt = 1e-5;
	int Nstep = 0;
	double nextRecordTime = 0.0;
	const double recordTimeInterval = 20.0;

	while (time <= FinalTime)
	{
		if ( time >= nextRecordTime || time == FinalTime)
		{
			double collectedWater = calculateTotalWater(ChannelList[0]);
			char* fileName = "CollectedWaterInChannel.dat";
			//append_to_file(fileName, time/60, collectedWater); 
			//outputFloodplainNodalDataToFile(time, FinalTime, recordTimeInterval);
			//outputFloodplainDataToFile(time, FinalTime, recordTimeInterval);
			//calculate_l2_error_2d(time);
			//output2DNodalError(time);
			outputDataToFile(time, FinalTime, recordTimeInterval);
			nextRecordTime += recordTimeInterval;
		}

		double maxLambda = 0;
	
		// advance the kinematic elements by one time step
		for (int i = 0; i < NumFloodplains; i++)
		{
			maxLambda = fmax(maxLambda, oneTimeStepKinematicEls(time, dt, i));
		}

		//couple kinematic elements with channels
		couple_kinEls_with_channels(time, dt);
	
		// advance channels by one time step
		for(int i=0; i<NumChannels; ++i)
		{
			//printf("Channel %d\n", i);
			oneTimeStepChannel(ChannelList[i],time,dt,i);
			maxLambda = fmax(maxLambda, ChannelList[i]->max_lambda);
		}
		
		// advance junctions by one time step
		for (int  i = 0; i < NumJunctions; ++i)
		{
			//printf("Evolving junction %d\n", i);
			oneTimeStep2D(JunctionList[i], time, dt);
			maxLambda = fmax(maxLambda, JunctionList[i]->max_lambda);
		}
		
		// couple channels with junctions
		couple_channels_with_junctions();


		double minLength = 99999999;	
		for (int i = 0; i < NumChannels; ++i)
		{
			minLength = fmin(minLength, ChannelList[i]->mindh);
		}

		for (int k = 0; k < NumFloodplains; k++)
		{
			int NumEl = FloodplainList[k]->NumEl;
			for (int i = 0; i < NumEl; ++i)
			{
				if (KinematicElList[i]->isActive)
					minLength = fmin(minLength, KinematicElList[i]->dh);
			}
			//oneTimeStep2D(FloodplainList[k], time, dt);
			//maxLambda = fmax(maxLambda, FloodplainList[k]->max_lambda);
			//minLength = fmin(minLength, FloodplainList[k]->minEdgLength);
		}

		for (int i = 0; i < NumJunctions; ++i)
		{
			minLength = fmin(minLength, JunctionList[i]->minEdgLength);
		}


		boundary_conditions(time);

		/**************** Set next time step ****************/
		//printf("****************************************\n");
		//printf("Reached time %lf \n", time);
		//printf("For the next time step:\n");
		//printf("minLength = %lf\n", minLength);
		
		//printf("maxLambda = %lf\n", maxLambda);
		//printf("minLength = %lf\n", minLength);
	
		dt = 0.3*minLength/maxLambda;
		if (time == FinalTime)
			dt = 1;
		else if ((time + dt) > FinalTime)
			dt = FinalTime - time;
		//printf("****************************************\n");
	
		printf("dt = %.13f\n", dt);
		printf("time = %lf\n", time);

		time = time + dt;
		Nstep++;
	} // end while loop

}
