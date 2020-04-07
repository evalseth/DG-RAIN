
/****************** Boundary Conditions ************************************/	

#include <math.h>
#include "Watershed.h"
#include "Globals.h"


void boundary_conditions(double time)
{

	// keep the height constant at the upstream end of kinematic routes
	for (int fp = 0; fp < NumFloodplains; fp++)
	{
		int NumEl = FloodplainList[fp]->NumEl;

		for (int i = 0; i < NumEl; i++)
		{
			struct kinematicEl* kinEl = KinematicElList[i];
			if (kinEl->isActive && kinEl->numUpstreamEls == 0)
			{
				kinEl->A[0] = 1e-5*kinEl->weq;
			}
			
		}
	}

	for (int i = 0; i < NumChannels; i++)
	{
		int upstreamBType = ChannelList[i]->BCtype[0];
		int downstreamBType = ChannelList[i]->BCtype[1];

		if (upstreamBType == 0) 	// no flow
		{
			ChannelList[i]->A[0] = ChannelList[i]->A[1];
			ChannelList[i]->Q[0] = -ChannelList[i]->Q[1];
		}

		else if (upstreamBType == 1) // Water height imposed
		{
			double height = ChannelList[i]->BCvals[0];
			double z_node = ChannelList[i]->NodalZ[0];
			double b_node = ChannelList[i]->NodalB[0];
			double m1val = ChannelList[i]->Nodalm1[0];
			double m2val = ChannelList[i]->Nodalm2[0];

			ChannelList[i]->A[0] = height*b_node + 0.5*m1val*height*height + 0.5*m2val*height*height;
			ChannelList[i]->Q[0] = ChannelList[i]->Q[1];
		}

		else if (upstreamBType == 2) // discharge imposed
		{
			ChannelList[i]->A[0] = ChannelList[i]->A[1];
			ChannelList[i]->Q[0] = ChannelList[i]->BCvals[1];
		}

		else if (upstreamBType == 3) //both imposed
		{
			double height = ChannelList[i]->BCvals[0];
			double z_node = ChannelList[i]->NodalZ[0];
			double b_node = ChannelList[i]->NodalB[0];
			double m1val = ChannelList[i]->Nodalm1[0];
			double m2val = ChannelList[i]->Nodalm2[0];

			ChannelList[i]->A[0] = height*b_node + 0.5*m1val*height*height + 0.5*m2val*height*height;

			ChannelList[i]->Q[0] = ChannelList[i]->BCvals[1];
		}

		else if (upstreamBType == 4) // just inflow
		{
			ChannelList[i]->A[0] = ChannelList[i]->A[1];
			ChannelList[i]->Q[0] = ChannelList[i]->Q[1];
		}	

		int NodeNum = ChannelList[i]->NumNodes+1;
		if (downstreamBType == 0) 	// no flow
		{
			ChannelList[i]->A[NodeNum] = ChannelList[i]->A[NodeNum-1];
			ChannelList[i]->Q[NodeNum] = -ChannelList[i]->Q[NodeNum-1];
		}

		else if (downstreamBType == 1) // Water height imposed
		{
			double height = ChannelList[i]->BCvals[2];
			double z_node = ChannelList[i]->NodalZ[NodeNum-1];
			double b_node = ChannelList[i]->NodalB[NodeNum-1];
			double m1val = ChannelList[i]->Nodalm1[NodeNum-1];
			double m2val = ChannelList[i]->Nodalm2[NodeNum-1];

			ChannelList[i]->A[NodeNum] = height*b_node + 0.5*m1val*height*height + 0.5*m2val*height*height;
			ChannelList[i]->Q[NodeNum] = ChannelList[i]->Q[NodeNum-1];
		}

		else if (downstreamBType == 2) // discharge imposed
		{
			ChannelList[i]->A[NodeNum] = ChannelList[i]->A[NodeNum-1];
			ChannelList[i]->Q[NodeNum] = ChannelList[i]->BCvals[3];
		}

		else if ( downstreamBType == 3)
		{
			double height = ChannelList[i]->BCvals[2];
			double z_node = ChannelList[i]->NodalZ[NodeNum-1];
			double b_node = ChannelList[i]->NodalB[NodeNum-1];
			double m1val = ChannelList[i]->Nodalm1[NodeNum-1];
			double m2val = ChannelList[i]->Nodalm2[NodeNum-1];

			ChannelList[i]->A[NodeNum] = height*b_node + 0.5*m1val*height*height + 0.5*m2val*height*height;

			ChannelList[i]->Q[NodeNum] = ChannelList[i]->BCvals[3];
		}

		else if (downstreamBType == 4) // outflow
		{
			ChannelList[i]->A[NodeNum] = ChannelList[i]->A[NodeNum-1];
			ChannelList[i]->Q[NodeNum] = ChannelList[i]->Q[NodeNum-1];

		}

	}

}




