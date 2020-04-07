/************************************************************************************************//**
* @file initialize.c
*
* This file contains code to initialize the channel structures and kinematicRoute structures 
* created in create_flow_paths.c and their member fields
*
* *************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Math_Functions.h"
#include "Channels_KinRoutes.h"
#include "Globals.h"

/***********************************************************************************************//**
* This function allocates necessary space for all the member fields of the channel structure. Then
* it assigns initial values for the quantities A and Q.
* **************************************************************************************************/ 
extern void* xcalloc(int NumEl, int size);

void initialize_channels()
{	
	for (int i=0; i<NumChannels; i++)
	{
		int NumNodes = ChannelList[i]->NumNodes;
		
		ChannelList[i]->A = xcalloc(NumNodes+2,sizeof(double));
		ChannelList[i]->Q = xcalloc(NumNodes+2,sizeof(double));
		ChannelList[i]->beta = xcalloc(NumNodes,sizeof(double));
		ChannelList[i]->qL = xcalloc(NumNodes, sizeof(double));

		for (int k=0; k<NumNodes+1; k++)
		{
			double z_node = ChannelList[i]->NodalZ[k];
			double b_node = ChannelList[i]->NodalB[k];
			double m1val = ChannelList[i]->Nodalm1[k];
			double m2val = ChannelList[i]->Nodalm2[k];
			//double zeta = 2;
			//double height = zeta + z_node;
			double height = 1e-8;
			ChannelList[i]->beta[k] = 1.0;
			ChannelList[i]->A[k+1] = height*b_node + 0.5*m1val*height*height + 0.5*m2val*height*height ;
			ChannelList[i]->Q[k+1] = 0;
		}
		ChannelList[i]->A[0] = ChannelList[i]->A[1];
		ChannelList[i]->Q[0] = ChannelList[i]->Q[1];
		ChannelList[i]->A[NumNodes+1] = ChannelList[i]->A[NumNodes];
		ChannelList[i]->Q[NumNodes+1] = ChannelList[i]->Q[NumNodes];
		
		ChannelList[i]->collectedQL = 0;
	}

}

/***********************************************************************************************//**
* This function allocates necessary space for all the member fields of the KinRoute structure. Then
* it assigns initial values for the quantities H,
* **************************************************************************************************/
void initialize_kinRoutes()
{
	for(int i = 0; i < NumKinRoutes; i++)
	{
		int NumNodes = KinRouteList[i]->NumNodes;

		KinRouteList[i]->H = xcalloc(NumNodes+2, sizeof(double));
	
		for(int k =0; k < NumNodes+2; k++)
		{
			KinRouteList[i]->H[k] = 1e-8;
		}

	}

}


