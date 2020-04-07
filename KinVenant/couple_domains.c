#include <stdio.h>
#include <math.h>
#include "Channels_KinRoutes.h"
#include "Globals.h"

void coupleKinRoutesWithChannels(double dt, double time)
{
	static int Nstep = 0;
	for (int i = 0; i < NumKinRoutes; i++)
	{
		int connectedChannel = KinRouteList[i]->connectedChan;
		int connectedChanEl = KinRouteList[i]->connectedChanEl;

		int NumNodes = KinRouteList[i]->NumNodes;
		double H = KinRouteList[i]->H[NumNodes];
		double S0 = KinRouteList[i]->dz[NumNodes-1];
		double n = KinRouteList[i]->NodalnFriction[NumNodes-1];
		double Htil = H - 0.0000001;
		if (Htil < 0)
			Htil = 0;
		double u;
		u = 1.0/n*sqrt(S0)*pow(Htil,2.0/3);
		double qL = u*Htil;				 // assumes unit width of kinematic routes
		int Np = ChannelList[connectedChannel]->Np;
		int begNode = connectedChanEl*Np;
		double elLen = ChannelList[connectedChannel]->dh[connectedChanEl];
		for (int j = 0; j < Np; j++)
		{
			ChannelList[connectedChannel]->qL[begNode+j] = qL/elLen;
		}
	
		ChannelList[connectedChannel]->collectedQL = ChannelList[connectedChannel]->collectedQL + dt*qL;

		// write out qL every half minute
		if (Nstep % 600 == 0)
		{
			// convert m^2/s to ft^2/s
			double qLFt = qL*10.764;
			// convert ft^2/s to ft^2/min
			qLFt = qLFt*60;

			FILE* Qfile = fopen("KinFlowQ.dat", "a");
			fprintf(Qfile, "%lf \t %lf \n", time/60, qL*60);
			//fprintf(Qfile, "%lf \t %lf \n", time/60, qLFt);
			fclose(Qfile);
	
		}

	}

	Nstep++;

}
