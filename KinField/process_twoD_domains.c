#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Watershed.h"

int NumJunctions;
int NumFloodplains;
struct TwoDRegion **JunctionList;
struct TwoDRegion **FloodplainList;

extern void* xcalloc(int items, int size);

int createEdgtoVertAndElCount(int* EltoVert, int currNumEl, int NumVerts, int* EdgtoVert, int* ElCount)
{
	int **NtoNodes = malloc(NumVerts*sizeof(int*));
	int *NodesCount = malloc(NumVerts*sizeof(int));

	for (int n = 0; n < NumVerts; ++n)
	{
		NtoNodes[n] = malloc(10*sizeof(int));
		ElCount[n] = 0;
		for (int j = 0; j < 10; ++j)
			NtoNodes[n][j] = -1;
	}

	for (int el = 0; el < currNumEl; ++el)
	{
		for (int n = 0; n < 3; ++n)
		{
			int node = EltoVert[el*3+n];
			ElCount[node] = ElCount[node]+1;
			int numConnNodes = 0;
			while(NtoNodes[node][numConnNodes] >= 0)
			{
				++numConnNodes;
			}

			int n1 = (n+1)%3;
			int n2 = (n+2)%3;
			int neighb0 = EltoVert[el*3+n1];
			int neighb1 = EltoVert[el*3+n2];
			int neighbfound = 0;
			int foundneighb[2] = {-1,-1};
			for (int j = 0; j < 10; ++j)
			{
				if (NtoNodes[node][j] < 0 || neighbfound == 2)
					break;
				if (NtoNodes[node][j] == neighb0)
				{
					neighbfound+=1;
					foundneighb[0] = 1;
				}
				if (NtoNodes[node][j] == neighb1)
				{
					neighbfound+=1;
					foundneighb[1] = 1;
				}

			}
			if (foundneighb[0] == -1)
			{
				NtoNodes[node][numConnNodes] = neighb0;
				++numConnNodes;
			}
			if (foundneighb[1] == -1)
			{
				NtoNodes[node][numConnNodes] = neighb1;
				++numConnNodes;
			}
			NodesCount[node] = numConnNodes;

		}
	}


	int EdgCount = 0;
	for (int n = 0; n < NumVerts; ++n)
	{
		int NumConnEdges = NodesCount[n];
		for (int k = 0; k < NumConnEdges; ++k)
		{
			int connNode = NtoNodes[n][k];
			if (connNode > n)
			{
				EdgtoVert[EdgCount*2] = n;
				EdgtoVert[EdgCount*2+1] = connNode;
				EdgCount += 1;

			}

		}

	}

	free(NodesCount);
	for (int i = 0; i < NumVerts; i++)
		free(NtoNodes[i]);
	free(NtoNodes);

	return EdgCount;

}


void read2DGridFile(char *gridFile, char* regionType)
{
	FILE *Coordinates = fopen(gridFile, "r");
	char juncLine[100];
	int NumRegions;
	fgets(juncLine, sizeof(juncLine), Coordinates);
	sscanf(juncLine, "%d", &NumRegions);

	struct TwoDRegion **currRegion;
	// allocate space for an array to store junctions
	if (strcmp(regionType,"Junction") == 0)
	{
		JunctionList = malloc(NumRegions*sizeof(struct TwoDRegion*));
		NumJunctions = NumRegions;
		currRegion = JunctionList;
	}
	else
	{
		FloodplainList = malloc(NumRegions*sizeof(struct TwoDRegion*));
		NumFloodplains = NumRegions;
		currRegion = FloodplainList;
	}
	for (int i = 0; i < NumRegions; ++i)
	{
		char newline[500];
		int currNumEl, NumVerts;
		fgets(newline, sizeof(newline), Coordinates);
		sscanf(newline, "%d\t%d", &currNumEl, &NumVerts);

		currRegion[i] = malloc(sizeof(struct TwoDRegion));
		currRegion[i]->NumEl = currNumEl;
		currRegion[i]->NumVerts = NumVerts;

		currRegion[i]->Vx = xcalloc(NumVerts, sizeof(double));
		currRegion[i]->Vy = xcalloc(NumVerts, sizeof(double));
		currRegion[i]->Vz = xcalloc(NumVerts, sizeof(double));
		currRegion[i]->VnFriction = xcalloc(NumVerts, sizeof(double));

		// read and store the coordinates
		for (int j = 0; j < NumVerts; ++j)
		{
			char grid1[500];
			int nodeNum;
			fgets(grid1, sizeof(grid1), Coordinates);
			//printf("grid1 = %s\n", grid1);
			sscanf(grid1, "%d %lf %lf %lf %lf", &nodeNum, &currRegion[i]->Vx[j], &currRegion[i]->Vy[j], &currRegion[i]->Vz[j], &currRegion[i]->VnFriction[j]);
		}

		// read and store the adjacency matrix
		currRegion[i]->EltoVert = malloc(3*currNumEl*sizeof(int));
		for (int j = 0; j < currNumEl; ++j)
		{
			char table[100];
			int elNum, NumVertsPerEl;
			fgets(table, sizeof(table), Coordinates);
			sscanf(table, "%d %d %d %d %d", &elNum, &NumVertsPerEl, &currRegion[i]->EltoVert[3*j], &currRegion[i]->EltoVert[j*3+1], &currRegion[i]->EltoVert[j*3+2]);
			// convert to 0 based indexing
			for (int s = 0; s < 3; s++)
			{
				currRegion[i]->EltoVert[3*j+s] -= 1;
				//printf("%lf \t", currRegion[i]->EltoVert[3*j+s]);
			}
			//printf("\n");
		}

		// Euler's relation: NumNodes - NumEdges + Numfaces = 1 -> NumEdges = NumNodes + NumEl - 1;
		int NumEdges = NumVerts + currNumEl - 1;
		currRegion[i]->TotalNumEdges = NumEdges;

		currRegion[i]->BdryPrescribed = malloc(NumEdges*sizeof(int));
		currRegion[i]->ChannelNumber = malloc(NumEdges*sizeof(int));
		if (strcmp(regionType, "Floodplain") == 0)
		{
			currRegion[i]->ChannelElNumber = malloc(NumEdges*sizeof(int));
		}

		// initialize arrays to store boundary conditions
		currRegion[i]->bzeta = malloc(NumEdges*sizeof(double));
		currRegion[i]->bQn = malloc(NumEdges*sizeof(double));
		
		// create EdgtoVert and ElCount
		currRegion[i]->EdgtoVert = malloc(2*NumEdges*sizeof(int));
		currRegion[i]->ElCount = xcalloc(NumVerts, sizeof(int));

		int EdgCount = createEdgtoVertAndElCount(currRegion[i]->EltoVert, currRegion[i]->NumEl, currRegion[i]->NumVerts, currRegion[i]->EdgtoVert, currRegion[i]->ElCount);

		if (EdgCount != currRegion[i]->TotalNumEdges)	
		{
			printf("EdgCount = %d\n", EdgCount);
			printf("Total number of edges is not correct\n");
			exit(EXIT_FAILURE);
		}
		
		// Read and mark edges and information about channels connected to them
		int NumInflowOutflow;
		char inOut[100];
		fgets(inOut, sizeof(inOut), Coordinates);
		sscanf(inOut, "%d", &NumInflowOutflow);
		int **InOutNodes = malloc(NumInflowOutflow*sizeof(int*));

		for (int j = 0; j < NumInflowOutflow; ++j)
		{
			InOutNodes[j] = malloc(5*sizeof(int));
			char bdry[1000];
			fgets(bdry, sizeof(bdry), Coordinates);
			if (strcmp(regionType, "Junction") == 0) 
				sscanf(bdry, "%d %d %d %d", &InOutNodes[j][0], &InOutNodes[j][1], &InOutNodes[j][2], &InOutNodes[j][3]); 
			else
				sscanf(bdry, "%d %d %d %d %d", &InOutNodes[j][0], &InOutNodes[j][1], &InOutNodes[j][2], &InOutNodes[j][3], &InOutNodes[j][4]); 
			InOutNodes[j][0] -=1; 
			InOutNodes[j][1] -=1; 
		}

		double minEdgLength = 99999999999;
		for (int e = 0; e < NumEdges; ++e)
		{
			currRegion[i]->bzeta[e] = -10000.0;
			currRegion[i]->bQn[e] = -10000.0;

			int n0 = currRegion[i]->EdgtoVert[e*2+0];
			int n1 = currRegion[i]->EdgtoVert[e*2+1];

			double edgLength = sqrt((currRegion[i]->Vx[n1]-currRegion[i]->Vx[n0])*(currRegion[i]->Vx[n1]-currRegion[i]->Vx[n0])+(currRegion[i]->Vy[n1]-currRegion[i]->Vy[n0])*(currRegion[i]->Vy[n1]-currRegion[i]->Vy[n0]));
			minEdgLength = fmin(minEdgLength,edgLength);

			int inflowOutflow = 0;
			for (int k = 0 ; k < NumInflowOutflow; ++k)
			{
				int b0 = InOutNodes[k][0];
				int b1 = InOutNodes[k][1];
				if (((b0==n0) && (b1 ==n1)) || ((b0==n1) && (b1==n0)))
				{
					currRegion[i]->BdryPrescribed[e] = InOutNodes[k][2];
					currRegion[i]->ChannelNumber[e] = InOutNodes[k][3];
					if (strcmp(regionType, "Floodplain") == 0)
						currRegion[i]->ChannelElNumber[e] = InOutNodes[k][4];
					inflowOutflow = 1;
					break;
				}

			}

			if (!inflowOutflow)
			{
				currRegion[i]->BdryPrescribed[e] = 0;	
				currRegion[i]->ChannelNumber[e] = -1;
			}

		}
		
		currRegion[i]->minEdgLength = minEdgLength;

		for (int j = 0; j < NumInflowOutflow; j++)
			free(InOutNodes[j]);

		// if floodplain, read additional information about external boundaries 
		if(strcmp(regionType, "Floodplain") == 0)
		{
			int NumPresBdryEdges;
			char BEarr[100];
			fgets(BEarr, sizeof(BEarr), Coordinates);
			sscanf(BEarr, "%d", &NumPresBdryEdges);
			for (int j = 0; j < NumPresBdryEdges; j++)
			{
				int v1, v2, btype;
				char vrts[200];
				fgets(vrts, sizeof(vrts), Coordinates);
				sscanf(vrts, "%d %d %d", &v1, &v2, &btype);
				v1 -= 1;
				v2 -= 1;

				for (int e = 0; e < NumEdges; ++e)
				{
					int n1 = currRegion[i]->EdgtoVert[e*2+0];
					int n2 = currRegion[i]->EdgtoVert[e*2+1];
					
					if( ((n1 == v1) && (n2 == v2)) || ((n2 == v1) && (n1 == v2)))
					{
						currRegion[i]->BdryPrescribed[e] = btype;
						char bc[500];
						/// 111 - only zeta prescribed
						/// 222 - only normal flow prescribed
						/// 333 - both zeta and normal flow prescribed
						/// 444 - just outflow
						/// 555 - manufactured solution
						if (btype == 111)
						{
							fgets(bc, sizeof(bc), Coordinates);
							sscanf(bc, "%lf", &currRegion[i]->bzeta[e]);

						}
						else if (btype == 222)
						{
							fgets(bc, sizeof(bc), Coordinates);
							sscanf(bc, "%lf", &currRegion[i]->bQn[e]);

						}

						else if (btype == 333)
						{
							fgets(bc, sizeof(bc), Coordinates);
							sscanf(bc, "%lf %lf", &currRegion[i]->bzeta[e], &currRegion[i]->bQn[e]);

						}

					}

				}


			}

		}

	}

	fclose(Coordinates);


}

void createEdgtoEls(int *EltoVert, int *EdgtoVert, int currNumEl, int NumEdges, int *EdgtoEls)
{
	for (int e = 0; e < NumEdges; ++e)
	{
		EdgtoEls[e*2] = -1;
		EdgtoEls[e*2+1] = -1;
	}

	for (int e = 0; e < NumEdges; ++e)
	{
		int v0 = EdgtoVert[e*2];
		int v1 = EdgtoVert[e*2+1];
		int elfound = 0;
		for (int el = 0; el < currNumEl; ++el)
		{
			for (int k = 0; k < 3; ++k)
			{
				int n0 = EltoVert[el*3+k];
				int ind1 = (k+1)%3;
				int n1 = EltoVert[el*3+ind1];
				if (((n0 == v0) && (n1 == v1))|| ((n0 == v1) &&(n1==v0)))
				{
					if (EdgtoEls[e*2] < 0)
						EdgtoEls[e*2] = el;
					else if (EdgtoEls[e*2+1] < 0)
						EdgtoEls[e*2+1] = el;
					elfound+=1;
					break;
				}
			}
			if (elfound==2)
				break;
			if ((el == (currNumEl-1)) && elfound == 1)
				EdgtoEls[e*2+1] = EdgtoEls[e*2];

		}
		if (elfound == 0)
		{
			printf("edge %d not an element edge \n",e);
			exit(EXIT_FAILURE);
		}

	}
}

void createGlobaltoLocalEdg(int* EltoVert, int* EdgtoVert, int *EdgtoEls, int NumEdges, int* GlobaltoLocalEdg)
{
	for (int e = 0; e < NumEdges; ++e)
	{
		int v0 = EdgtoVert[e*2];
		int v1 = EdgtoVert[e*2+1];
		int el1 = EdgtoEls[e*2];
		int el2 = EdgtoEls[e*2+1];

		int el1_v0 = EltoVert[el1*3];
		int el1_v1 = EltoVert[el1*3+1];
		int el1_v2 = EltoVert[el1*3+2];

		int el2_v0 = EltoVert[el2*3];
		int el2_v1 = EltoVert[el2*3+1];
		int el2_v2 = EltoVert[el2*3+2];

		int el1_local_edg;
		int el2_local_edg;

		if ((v0 == el1_v0 && v1 == el1_v1) || (v0 == el1_v1 && v1 == el1_v0))
			el1_local_edg = 0;
		else if ((v0 == el1_v1 && v1 == el1_v2) || (v0 == el1_v2 && v1 == el1_v1))
			el1_local_edg = 1;
		else if ((v0 == el1_v2 && v1 == el1_v0) || (v0 == el1_v0 && v1 == el1_v2))
			el1_local_edg = 2;
		else
		{
			printf("error in calculating interior element local edge number \n");
			exit(EXIT_FAILURE);
		}

		GlobaltoLocalEdg[e*2] = el1_local_edg;

		if ((v0 == el2_v0 && v1 == el2_v1) || (v0 == el2_v1 && v1 == el2_v0))
			el2_local_edg = 0;
		else if ((v0 == el2_v1 && v1 == el2_v2) || (v0 == el2_v2 && v1 == el2_v1))
			el2_local_edg = 1;
		else if ((v0 == el2_v2 && v1 == el2_v0) || (v0 == el2_v0 && v1 == el2_v2))
			el2_local_edg = 2;
		else
		{
			printf("error in calculating exterior element local edge number \n");
			exit(EXIT_FAILURE);
		}

		GlobaltoLocalEdg[e*2+1] = el2_local_edg;

	}

}

void createVtoEl(int* EltoVert, int* ElCount, int NumEl, int NumVerts, int** VtoEl)
{
	for (int n = 0; n < NumVerts; ++n)
	{	
		int NumConnEl = ElCount[n];
		VtoEl[n] = malloc(NumConnEl*sizeof(int));
		for (int i = 0; i < NumConnEl; ++i)
			VtoEl[n][i] = -1;
	}

	for (int el =0; el < NumEl; ++el)
	{
		for (int n =0; n < 3; ++n)
		{	
			int ind = el*3+n;
			int node = EltoVert[ind];
			for (int k = 0; k < ElCount[node]; ++k)
			{
				if (VtoEl[node][k] < 0)
				{
					VtoEl[node][k] = el;
					break;
				} 
			}
		}
	}
}
