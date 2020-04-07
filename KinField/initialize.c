/************************************************************************************************//**
* @file initialize_channels.c
*
* This file contains code to initialize the channel and junction structures created in create_channel_network.c
* and their member fields.
* It also initializes the kinematic elements created in process_watershed.c
*
* *************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Watershed.h"
#include "Globals.h"
#include "Math_Functions.h"
#include "ManufacturedSolution.h"

/***********************************************************************************************//**
* This function allocates necessary space for all the member fields of the channel structure. Then
* it assigns initial values for the quantities A and Q. It also assigns an initial wet/dry state
* for each element of the %channel.
* 
* *************************************************************************************************/ 

void initialize_channels()
{	
	for (int i=0; i<NumChannels; ++i)
	{
		int NumNodes = ChannelList[i]->NumNodes;
	
		ChannelList[i]->A = xcalloc(NumNodes+2, sizeof(double));
		ChannelList[i]->Q = xcalloc(NumNodes+2, sizeof(double));
		ChannelList[i]->beta = xcalloc(NumNodes,sizeof(double));
		ChannelList[i]->qL = xcalloc(NumNodes, sizeof(double));
	
		#ifdef WDON
		ChannelList[i]->WD = malloc(ChannelList[i]->NumEl*sizeof(int));
		#endif

		for (int k=0; k<NumNodes; ++k)
		{
			
			double z_node = ChannelList[i]->NodalZ[k];
			double b_node = ChannelList[i]->NodalB[k];
			double m1val = ChannelList[i]->Nodalm1[k];
			double m2val = ChannelList[i]->Nodalm2[k];
			//printf("znode = %lf ,bnode = %lf, m1val = %lf, m2val = %lf\n", z_node, b_node, m1val, m2val);
			
			double zeta = 1.0;
			double height = zeta + z_node;
			//double height = 1.0;
			if (b_node == 0)
			{
				printf("width zero at node %d of channel %d \n", k, i);
				exit(EXIT_FAILURE);
			}

			ChannelList[i]->A[k+1] = height*b_node + 0.5*m1val*height*height + 0.5*m2val*height*height ;
			ChannelList[i]->Q[k+1] = 0;
			ChannelList[i]->beta[k] = 1.0;
		}

		ChannelList[i]->A[0] = ChannelList[i]->A[1];
		ChannelList[i]->Q[0] = ChannelList[i]->Q[1];
		ChannelList[i]->A[NumNodes+1] = ChannelList[i]->A[NumNodes];
		ChannelList[i]->Q[NumNodes+1] = ChannelList[i]->Q[NumNodes];

		ChannelList[i]->collectedQL = 0;

		#ifdef WDON
		if (k < NumNodes-1)
			ChannelList[i]->WD[k] = -1;
		#endif
	}

}

/*********************************************************************************************//**
* This function allocates necessary space for all the member fields of the junction structure.
* Then it assigns initial values for the quantities zeta, Qx and Qy. It also assigns an intiial
* wet/dry state for each element of the %junction.
*
* ************************************************************************************************/ 
void initialize_junctions()
{

	for (int i =0; i<NumJunctions; ++i)
	{
		int NumEl = JunctionList[i]->NumEl;
		int Np = JunctionList[i]->Np;

		JunctionList[i]->zeta = malloc(NumEl*sizeof(double*));
		JunctionList[i]->Qx = malloc(NumEl*sizeof(double*));
		JunctionList[i]->Qy = malloc(NumEl*sizeof(double*));


		#ifdef WDON
		JunctionList[i]->WD = malloc(NumEl*sizeof(int));
		#endif
	
		double zeta = 1.0;
		//double height= 1.0;

		for (int k=0; k<NumEl; ++k)
		{
			JunctionList[i]->zeta[k] = xcalloc(Np, sizeof(double));
			JunctionList[i]->Qx[k] = xcalloc(Np, sizeof(double));
			JunctionList[i]->Qy[k] = xcalloc(Np, sizeof(double));
			for (int j=0; j <Np; ++j)
			{
				double z_val = JunctionList[i]->NodalZ[k][j];
				JunctionList[i]->zeta[k][j] = zeta;
				//JunctionList[i]->zeta[k][j] = height - z_val;
				JunctionList[i]->Qx[k][j] = 0;
				JunctionList[i]->Qy[k][j] =0;
			}
			
			#ifdef WDON
			JunctionList[i]->WD[k] = -1;
			#endif
		}

	}

}

/*********************************************************************************************//**
* This function assigns initial values for H in kinematic elements
* ************************************************************************************************/ 
void initialize_kinematicEls()
{
	int NumEl = FloodplainList[0]->NumEl;
	for (int i = 0; i < NumEl; ++i)
	{		
		struct kinematicEl *myKinEl = KinematicElList[i]; 
		// if the kinematic element is active
		if (myKinEl->isActive)
		{
			int numNodes = myKinEl->Np;
			myKinEl->A = xcalloc(numNodes, sizeof(double));
			myKinEl->K = xcalloc(numNodes, sizeof(double));
			myKinEl->EffPor = xcalloc(numNodes, sizeof(double));
			myKinEl->suctionHead = xcalloc(numNodes, sizeof(double));
			for (int j = 0; j < numNodes; ++j)
			{
				myKinEl->A[j] = 1e-7*myKinEl->weq;
				 // infiltration parameters for sand
				myKinEl->K[j] = 11.78*0.0328/3600; // 11.78 cm/hr
				myKinEl->EffPor[j] = 0.417;
				myKinEl->suctionHead[j] = 4.95*0.0328; 	// 4.95 cm

			}
		}
	}

}

/*************************************************************************//**
* This function assign initial values for zeta, Qx and Qy in Floodplains
***************************************************************************/
void initialize_floodplains()
{
	for (int i = 0; i < NumFloodplains; i++)
	{
		int NumEl = FloodplainList[i]->NumEl;
		int Np = FloodplainList[i]->Np;

		FloodplainList[i]->zeta = malloc(NumEl*sizeof(double*));
		FloodplainList[i]->Qx = malloc(NumEl*sizeof(double*));
		FloodplainList[i]->Qy = malloc(NumEl*sizeof(double*));

		double zeta = 0.1;
		//double height= 1.0;

		for (int k=0; k<NumEl; ++k)
		{
			FloodplainList[i]->zeta[k] = xcalloc(Np, sizeof(double));
			FloodplainList[i]->Qx[k] = xcalloc(Np, sizeof(double));
			FloodplainList[i]->Qy[k] = xcalloc(Np, sizeof(double));
			for (int j=0; j <Np; ++j)
			{
				//double xval = FloodplainList[i]->NodalX[k][j];
				//double yval = FloodplainList[i]->NodalY[k][j];
				//double h = getmanH(xval, yval, 0);
				//double Qx = getQx(xval, yval, 0);
				//double Qy = getQy(xval, yval, 0);
				//FloodplainList[i]->zeta[k][j] = h;
				//FloodplainList[i]->Qx[k][j] = Qx;
				//FloodplainList[i]->Qy[k][j] =Qy;

				double z_val = FloodplainList[i]->NodalZ[k][j];
				FloodplainList[i]->zeta[k][j] = zeta;
				FloodplainList[i]->Qx[k][j] = 0;
				FloodplainList[i]->Qy[k][j] =0;
				
			}
			
		}
		
	}

}
