#include <stdio.h>
#include <math.h>
#include "Watershed.h"
#include "Globals.h"
#include "Math_Functions.h"
#include "Constitutive_Equations.h"

/**********************************************************************************************//**
* This routine couples the channels with the respective kinematic elements
* *************************************************************************************************/

// This is a one way coupling. Kinematic elements do not get input from channels.
// for now, just assume that there is only one floodplain
void couple_kinEls_with_channels(double time, double dt)
{
	static int Nstep = 0;
	double totalqL = 0;
	double totalArea = 0;
	//reset the lateral flow to zero in the beginning
	for (int i = 0; i < NumChannels; ++i)
	{
		int NumNodes = ChannelList[i]->NumNodes;
		for (int j = 0; j < NumNodes; ++j)
			ChannelList[i]->qL[j] = 0;
	}

	double *chanqL = xcalloc(NumChannels, sizeof(double));
	// begin coupling
	int NumEl = FloodplainList[0]->NumEl;
	for (int i = 0; i < NumEl; ++i)
	{
		struct kinematicEl* kinEl = KinematicElList[i];
		if (kinEl->isActive && kinEl->numDownstreamEls == 0)
		{
			int connectedChannel = kinEl->connChanNum;
			
			//printf("channel connected with KinEl %d = %d\n", i, connectedChannel);
			int Np = kinEl->Np;
			double A = kinEl->A[Np-1];
			double S0 = kinEl->dz[Np-1];
			double nf = kinEl->NodalnFriction[Np-1];
			double weq = kinEl->weq;
			//double H = A;
			double H = A/weq;
			double Htil = H - 1e-7;					// current wetting/drying treatment
			if (Htil < 0)
				Htil = 0;

			// Manning's relationship
			double u = sqrt(S0)*pow(Htil, 2.0/3)/nf;
			
			// Chezy's relationship
			//double u = sqrt(S0*Htil)*nf;

			double qL = Htil*u;
		
			chanqL[connectedChannel] += qL*weq;

						
		}
		
		
	}
	
	for (int i = 0; i < NumChannels; i++)
	{
		int ChanNumNodes = ChannelList[i]->NumNodes;

		for (int j = 0; j < ChanNumNodes; ++j)
		{
			ChannelList[i]->qL[j] = chanqL[i]/ChannelList[i]->channelLength;
		}
		//ChannelList[0]->collectedQL = ChannelList[0]->collectedQL + dt*totalqL;
		printf("Time = %lf totalqL = %lf\n", time, totalqL);

	}

	Nstep++;

}




/**********************************************************************************************//**
* This function couples the channels and junctions with each other by applying the coupling 
* conditions on all of the channels and the junctions
*
* *************************************************************************************************/

// assumes we are using P1 approximations in junctions
void couple_channels_with_junctions()
{
	// loop over the channels
	for (int i = 0; i < NumChannels; ++i)
	{

		int numInflowBC, numOutflowBC; 		// inflow and outflow pertain to junctions and not channels
		int NumInflowEdges = ChannelList[i]->NumInflowJunctionEdges;
		double totalInflowQn = 0;
		double totalInflowArea = 0;

		// get channel values at inflow (junction's outflow)
		double A_in = ChannelList[i]->A[1];
		double b_in = ChannelList[i]->NodalB[0];
		double m1_in = ChannelList[i]->Nodalm1[0];
		double m2_in = ChannelList[i]->Nodalm2[0];
		double H_in = getH(A_in, b_in, m1_in, m2_in);
		double Qn_in = ChannelList[i]->Q[1];
		double c_in = sqrt(g*A_in/b_in);
		double u_in = Qn_in/A_in;

		//printf("Width of channel %d at inflow = %lf\n", i, b_in);

		if (u_in + c_in < 0)
			numInflowBC = 0;
		else if (u_in - c_in < 0)
			numInflowBC = 1;
		else
			numInflowBC = 2;

		// get channels values at outflow (junction's inflow)
		int  NumNodes = ChannelList[i]->NumNodes;
		double A_out = ChannelList[i]->A[NumNodes];
		double b_out = ChannelList[i]->NodalB[NumNodes-1];
		double m1_out = ChannelList[i]->Nodalm1[NumNodes-1];
		double m2_out = ChannelList[i]->Nodalm2[NumNodes-1];
		double H_out = getH(A_out, b_out, m1_out, m2_out);
		double Qn_out = ChannelList[i]->Q[NumNodes];
		double c_out = sqrt(g*A_out/b_out);
		double u_out = Qn_out/A_out;
		double z_out = ChannelList[i]->NodalZ[NumNodes-1];

		//printf("Width of channel %d at outflow = %lf\n", i, b_out);
		if (u_out+c_out < 0)
			numOutflowBC = 2;
		else if (u_out - c_out < 0)
			numOutflowBC = 1;
		else
			numOutflowBC = 0;

		// loop over all the junction edges that are connected to the 
		// beginning of this channel
		for (int j = 0; j < NumInflowEdges; ++j)
		{
			int globalEdgNum = ChannelList[i]->InflowJunctionEdges[2*j];
			int juncNum = ChannelList[i]->InflowJunctionEdges[2*j+1];

			struct TwoDRegion *junc = JunctionList[juncNum];

			int el = junc->EdgtoEls[2*globalEdgNum];	
			int edg = junc->GlobaltoLocalEdg[globalEdgNum*2];
		
			int edgv1 = junc->EdgtoVert[globalEdgNum*2];
			int edgv2 = junc->EdgtoVert[globalEdgNum*2+1];
			double xv1 = junc->Vx[edgv1];
			double xv2 = junc->Vx[edgv2];
			double yv1 = junc->Vy[edgv1];
			double yv2 = junc->Vy[edgv2];
			double edgLength = sqrt((xv1-xv2)*(xv1-xv2)+(yv1-yv2)*(yv1-yv2));
			//printf("x1 = %lf, x2 = %lf, y1 = %lf, y2 = %lf \n", xv1, xv2, yv1, yv2);
			//printf("edgLength of edg %d of junction %d = %lf\n", globalEdgNum, juncNum, edgLength);

			double nx = junc->nx[el*3+edg];
			double ny = junc->ny[el*3+edg];

			// calculate average value on the edge of the junction
			int n1 = Fmask[0][edg];
			int n2 = Fmask[1][edg];
			double juncZeta = 0.5*(junc->zeta[el][n1] + junc->zeta[el][n2]);
			double juncZ = 0.5*(junc->NodalZ[el][n1] + junc->NodalZ[el][n2]);
			double juncH = juncZeta + juncZ;
			double Qx = 0.5*(junc->Qx[el][n1] + junc->Qx[el][n2]);
			double Qy = 0.5*(junc->Qy[el][n1] + junc->Qy[el][n2]);
		
			double juncQn = Qx*nx+Qy*ny;
	
			totalInflowArea += juncH*edgLength;
			totalInflowQn += juncQn*edgLength;	
		
			JunctionList[juncNum]->bzeta[globalEdgNum] = H_in - juncZ;
			JunctionList[juncNum]->bQn[globalEdgNum] = Qn_in/b_in;
	
		} // end inflow junction edges loop

		if (NumInflowEdges != 0)
		{
			ChannelList[i]->Q[0] = totalInflowQn;
			ChannelList[i]->A[0] = totalInflowArea;
		}
		
		
		int NumOutflowEdges = ChannelList[i]->NumOutflowJunctionEdges;
		double totalOutflowQn = 0;
		double totalOutflowArea = 0;


		// loop over all the edges of the junction that are connected to the end
		// of this channel
		for (int j = 0; j < NumOutflowEdges; ++j)
		{
			int globalEdgNum = ChannelList[i]->OutflowJunctionEdges[2*j];
			int juncNum = ChannelList[i]->OutflowJunctionEdges[2*j+1];
			
			struct TwoDRegion *junc = JunctionList[juncNum];

			int el = junc->EdgtoEls[2*globalEdgNum];	
			int edg = junc->GlobaltoLocalEdg[globalEdgNum*2];
		
			int edgv1 = junc->EdgtoVert[globalEdgNum*2];
			int edgv2 = junc->EdgtoVert[globalEdgNum*2+1];
			double xv1 = junc->Vx[edgv1];
			double xv2 = junc->Vx[edgv2];
			double yv1 = junc->Vy[edgv1];
			double yv2 = junc->Vy[edgv2];
			double edgLength = sqrt((xv1-xv2)*(xv1-xv2)+(yv1-yv2)*(yv1-yv2));
			//printf("x1 = %lf, x2 = %lf, y1 = %lf, y2 = %lf \n", xv1, xv2, yv1, yv2);
			//printf("edgLength for edg %d of junction %d= %lf\n", globalEdgNum, juncNum, edgLength);

			double nx = junc->nx[el*3+edg];
			double ny = junc->ny[el*3+edg];

			// calculate average value on the edge of the junction
			int n1 = Fmask[0][edg];
			int n2 = Fmask[1][edg];
			double juncZeta = 0.5*(junc->zeta[el][n1] + junc->zeta[el][n2]);
			double juncZ = 0.5*(junc->NodalZ[el][n1] + junc->NodalZ[el][n2]);
			double juncH = juncZeta + juncZ;
			double Qx = 0.5*(junc->Qx[el][n1] + junc->Qx[el][n2]);
			double Qy = 0.5*(junc->Qy[el][n1] + junc->Qy[el][n2]);
		
			double juncQn = Qx*nx+Qy*ny;
		
			totalOutflowArea += juncH*edgLength;
			totalOutflowQn += juncQn*edgLength;	

			
			JunctionList[juncNum]->bzeta[globalEdgNum] = H_out - juncZ;
			JunctionList[juncNum]->bQn[globalEdgNum] = -Qn_out/b_out;
				
		} // end outflow junction edges loop

		if (NumOutflowEdges != 0)
		{
		
			int NumNodes = ChannelList[i]->NumNodes;
			//printf("ChanZ = %lf\n", ChannelList[i]->NodalZ[NumNodes-1]);
	
		
				ChannelList[i]->Q[NumNodes+1] = -totalOutflowQn;
				ChannelList[i]->A[NumNodes+1] = totalOutflowArea;
		}

	} // end channel loop

}
