
/****************** Boundary Conditions ************************************/	

#include <math.h>
#include "Globals.h"

//extern const double g;

void boundary_conditions(double time)
{

	// keep the height constant at the upstream end of kinematic routes
	for (int i = 0; i < NumKinRoutes; i++)
	{
		int NumNodes = KinRouteList[i]->NumNodes;
		KinRouteList[i]->H[NumNodes+1] = KinRouteList[i]->H[NumNodes];
		KinRouteList[i]->H[0] = 0.0000001;
		
	}

	// implement a no-flow boundary condition on the channel
	int chanNumNodes = ChannelList[0]->NumNodes;
	ChannelList[0]->A[0] = ChannelList[0]->A[1];
	ChannelList[0]->A[chanNumNodes+1] = ChannelList[0]->A[chanNumNodes];
	

	//ChannelList[0]->Q[0] = ChannelList[0]->Q[1];
	//ChannelList[0]->Q[chanNumNodes+1] = ChannelList[0]->Q[chanNumNodes];
	
	ChannelList[0]->Q[0] = -ChannelList[0]->Q[1];
	ChannelList[0]->Q[chanNumNodes+1] = -ChannelList[0]->Q[chanNumNodes];

	//ChannelList[0]->A[0] = 0.01;
	//ChannelList[0]->A[chanNumNodes+1] = 0.01;
	

}




