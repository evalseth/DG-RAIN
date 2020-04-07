#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Channels_KinRoutes.h"
#include "Globals.h"

void create_flow_paths(char* ChannelGrid, char* KinGrid)
{
	FILE* ChanGrid = fopen(ChannelGrid, "r");
	char line[100];
	fgets(line, sizeof(line), ChanGrid);
	sscanf(line, "%d", &NumChannels);

	// allocate space for an array to store channels
	ChannelList = malloc(NumChannels*sizeof(struct channel*));
	
	for (int i = 0; i < NumChannels; i++)
	{

		ChannelList[i] = malloc(sizeof(struct channel));
		
		char newline[100];
		int NumEl;
		fgets(newline, sizeof(newline),ChanGrid);
		sscanf(newline, "%d", &NumEl);
		int NumEdges = NumEl+1;
		
		ChannelList[i]->NumEdges = NumEdges;
		ChannelList[i]->NumEl = NumEl;
		ChannelList[i]->x = malloc(NumEdges*sizeof(double));
		ChannelList[i]->y = malloc(NumEdges*sizeof(double));
		ChannelList[i]->z = malloc(NumEdges*sizeof(double));
		ChannelList[i]->b = malloc(NumEdges*sizeof(double));
		ChannelList[i]->m1 = malloc(NumEdges*sizeof(double));
		ChannelList[i]->m2 = malloc(NumEdges*sizeof(double));
		ChannelList[i]->nFriction = malloc(NumEdges*sizeof(double));
		ChannelList[i]->dh = malloc(NumEl*sizeof(double));

		for (int j = 0; j < NumEdges; j++)
		{
			double x, y, z, b, m1, m2, n;
			int dummy;
			char grid[500];
			fgets(grid, sizeof(grid), ChanGrid);
			sscanf(grid, "%d %lf %lf %lf %lf %lf %lf %lf", &dummy, &x, &y, &z, &b, &m1, &m2, &n);
			ChannelList[i]->x[j] = x;
			ChannelList[i]->y[j] = y;
			ChannelList[i]->z[j] = z;
			ChannelList[i]->b[j] = b;
			ChannelList[i]->m1[j] = m1;
			ChannelList[i]->m2[j] = m2;
			ChannelList[i]->nFriction[j] = n;
		}

		double mindh = 9999999999999999;
		for (int j = 0; j < NumEl; j++)
		{
			double dh = sqrt((ChannelList[i]->x[j+1] - ChannelList[i]->x[j])*(ChannelList[i]->x[j+1]-ChannelList[i]->x[j]) + (ChannelList[i]->y[j+1] - ChannelList[i]->y[j])*  (ChannelList[i]->y[j+1] - ChannelList[i]->y[j]));
			ChannelList[i]->dh[j] = dh;
		mindh = fmin(mindh, dh);
		}
	
	ChannelList[i]->mindh = mindh;
	
	} // end loop for Channels

	fclose(ChanGrid);

	FILE* KinematicGrid = fopen(KinGrid, "r");
	char kinLine[100];
	fgets(kinLine, sizeof(kinLine), KinematicGrid);
	sscanf(kinLine, "%d", &NumKinRoutes);

	// allocate space for an array to store kinematic routes
	KinRouteList = malloc(NumKinRoutes*sizeof(struct kinematicRoute*));
	for (int i = 0; i < NumKinRoutes; i++)
	{
		char newline[500];
		int NumEl, NumEdges;
		fgets(newline, sizeof(newline), KinematicGrid);
		sscanf(newline, "%d", &NumEl);
		NumEdges = NumEl + 1;

		KinRouteList[i] = malloc(sizeof(struct kinematicRoute));
		KinRouteList[i]->NumEdges = NumEdges;
		KinRouteList[i]->NumEl = NumEl;
		KinRouteList[i]->x = malloc(NumEdges*sizeof(double));
		KinRouteList[i]->y = malloc(NumEdges*sizeof(double));
		KinRouteList[i]->z = malloc(NumEdges*sizeof(double));
		KinRouteList[i]->nFriction = malloc(NumEdges*sizeof(double));
		KinRouteList[i]->dh = malloc(NumEl*sizeof(double));

		for (int j = 0; j < NumEdges; j++)
		{
			double x, y, z, n;
			char grid[1000];
			int nodeNum;
			fgets(grid, sizeof(grid), KinematicGrid);
			sscanf(grid, "%d %lf %lf %lf %lf", &nodeNum, &x, &y, &z, &n);
			KinRouteList[i]->x[j] = x;
			KinRouteList[i]->y[j] = y;
			KinRouteList[i]->z[j] = z;
			KinRouteList[i]->nFriction[j] = n;

		}
		
		double mindh = 9999999999999999;
		for (int j = 0; j < NumEl; j++)
		{
			double dh = sqrt((KinRouteList[i]->x[j+1] - KinRouteList[i]->x[j])*(KinRouteList[i]->x[j+1]-KinRouteList[i]->x[j]) + (KinRouteList[i]->y[j+1] - KinRouteList[i]->y[j])*(KinRouteList[i]->y[j+1] - KinRouteList[i]->y[j]));
			KinRouteList[i]->dh[j] = dh;
		mindh = fmin(mindh, dh);
		}
	
		KinRouteList[i]->mindh = mindh;

	}

	// throw away a newline
	char grid[100];
	fgets(grid, sizeof(grid), KinematicGrid);
	for (int j = 0; j < NumKinRoutes; j++)
	{
		int kinRouteNum, channelNum, chanElNum;
		char newgrid[1000];
		fgets(newgrid, sizeof(newgrid), KinematicGrid);
		sscanf(newgrid, "%d %d %d", &kinRouteNum, &channelNum, &chanElNum);
		KinRouteList[kinRouteNum]->connectedChan = channelNum;
		KinRouteList[kinRouteNum]->connectedChanEl = chanElNum;

	}

	fclose(KinematicGrid);

}

