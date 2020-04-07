//  main.c
//  FlowPaths
//
//  Created by Prapti Neupane on 6/9/15.
//  Copyright (c) 2015 Prapti Neupane. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <linux/limits.h>
#include <sys/stat.h>
#include <Python.h>
#include "Globals.h"
#include "Data_Structures.h"
#include "process_twoD_domains.h"

/*************** Globals variables prototyped in Globals.h **********************/
int NumKinematicEls;

struct kinematicEl **KinematicElList;
struct channel **ChannelList;

int *Valence;

int *FlowEdg; 	   // stores the global number of the flow edge
int *LocalFlowEdg; // stores the local number of the flow edge

/*************************************************************/

//int numChannelEdges;
//int **channel_array;

extern void create_channel_network(char* ChannelNodes, char* JunctionNodes);
void process_floodplain_attributes();
/********************************************************************************************/
/********************************************************************************************/

/**************** Routines to read in and process the initial mesh *************************/
void create_floodplain(char *FloodplainMesh)
{
	read2DGridFile(FloodplainMesh, "Floodplain");

	for (int i = 0; i < NumFloodplains; i++)
	{
		int NumEdges = FloodplainList[i]->TotalNumEdges;
		//FloodplainList[i]->EdgtoVert = malloc(2*NumEdges*sizeof(int));
		//FloodplainList[i]->ElCount = xcalloc(FloodplainList[i]->NumVerts, sizeof(int));
		//int EdgCount = createEdgtoVertAndElCount(FloodplainList[i]->EltoVert, FloodplainList[i]->NumEl, FloodplainList[i]->NumVerts, FloodplainList[i]->EdgtoVert, FloodplainList[i]->ElCount);
		//
		//if (EdgCount != FloodplainList[i]->TotalNumEdges)	
		//{
		//	printf("EdgCount = %d\n", EdgCount);
		//	printf("Total number of edges is not correct\n");
		//	exit(EXIT_FAILURE);
		//}

		// Create EdgtoEls
		FloodplainList[i]->EdgtoEls = malloc(2*NumEdges*sizeof(int));
		createEdgtoEls(FloodplainList[i]->EltoVert, FloodplainList[i]->EdgtoVert, FloodplainList[i]->NumEl, FloodplainList[i]->TotalNumEdges, FloodplainList[i]->EdgtoEls);

		// Create GlobaltoLocalEdg
		FloodplainList[i]->GlobaltoLocalEdg = xcalloc(2*NumEdges, sizeof(int));
		createGlobaltoLocalEdg(FloodplainList[i]->EltoVert, FloodplainList[i]->EdgtoVert, FloodplainList[i]->EdgtoEls, FloodplainList[i]->TotalNumEdges, FloodplainList[i]->GlobaltoLocalEdg);
		
		// create VToEl
		FloodplainList[i]->VtoEl = malloc(FloodplainList[i]->NumVerts*sizeof(int*));
		createVtoEl(FloodplainList[i]->EltoVert, FloodplainList[i]->ElCount, FloodplainList[i]->NumEl, FloodplainList[i]->NumVerts, FloodplainList[i]->VtoEl);		

	} // end floodplain loop

	process_floodplain_attributes();
}


void process_floodplain_attributes()
{
	for (int fp = 0; fp < NumFloodplains; fp++)
	{
		int NumEl = FloodplainList[fp]->NumEl;
		int NumEdges = FloodplainList[fp]->TotalNumEdges;

		// Create ElToEdg
		FloodplainList[fp]->EltoEdg = malloc(3*NumEl*sizeof(int));
		for (int el = 0; el < NumEl; el++)
		{
			FloodplainList[fp]->EltoEdg[el*3]=-1;
			FloodplainList[fp]->EltoEdg[el*3+1]=-1;
			FloodplainList[fp]->EltoEdg[el*3+2]=-1;
		}

		for (int e = 0; e < NumEdges; e++)
		{
			int el1 = FloodplainList[fp]->EdgtoEls[e*2];
			int el2 = FloodplainList[fp]->EdgtoEls[e*2+1];
			int edgNum = 0;
			while(FloodplainList[fp]->EltoEdg[el1*3+edgNum] > -1  && edgNum < 2)
				edgNum++;

			FloodplainList[fp]->EltoEdg[el1*3+edgNum] = e;

			if (el2 != el1)
			{
				int edgNum = 0;
				while(FloodplainList[fp]->EltoEdg[el2*3+edgNum] > -1 && edgNum < 2)
					edgNum++;
				FloodplainList[fp]->EltoEdg[el2*3+edgNum] = e;
			}

		}

		//arrange EltoEdg in the same order as vertices
		for (int el = 0; el < NumEl; el++)
		{
			int edg[3] = {FloodplainList[fp]->EltoEdg[el*3], FloodplainList[fp]->EltoEdg[el*3+1], FloodplainList[fp]->EltoEdg[el*3+2]};

			int node[3]= {FloodplainList[fp]->EltoVert[el*3], FloodplainList[fp]->EltoVert[el*3+1], FloodplainList[fp]->EltoVert[el*3+2]};

			for (int i= 0; i < 3; i++)
			{
				int mynode1 = FloodplainList[fp]->EdgtoVert[2*edg[i]];
				int mynode2 = FloodplainList[fp]->EdgtoVert[2*edg[i]+1];
				for (int j = 0; j < 3; j++)
				{
					if ((mynode1 == node[j%3] && mynode2 == node[(j+1)%3]) || (mynode2 == node[j%3] && mynode1 == node[(j+1)%3]))
					{

						FloodplainList[fp]->EltoEdg[el*3+j] = edg[i];
						break;
					}

				}

			}

		}
	}
}

void calculate_flow_edge(int fpNum, int* FlowEdg, int* LocalFlowEdg)
{
	int NumEl = FloodplainList[fpNum]->NumEl;
	int NumNodes = FloodplainList[fpNum]->NumVerts;

	// Initialize the Python Interpreter
	Py_Initialize();

	PyObject *moduleName, *module, *moduleFlowDirFunc, *moduleContents, *args;
	PyObject *input_x, *input_y, *input_z, *triangles;
	PyObject *localFlowDir;

		// Build the name object
	moduleName = PyString_FromString((char*) "calculate_flow_direction");
	if (!moduleName)
	{
		PyErr_Print();
	}

	// Load the module object
	module = PyImport_Import(moduleName);
	if (!module)
	{
		PyErr_Print();
	}

	moduleContents = PyModule_GetDict(module);
	if (!moduleContents)
		PyErr_Print();

	moduleFlowDirFunc = PyDict_GetItemString(moduleContents, (char*) "calculate_local_flow_dir");
	if (!moduleFlowDirFunc)
		PyErr_Print();
	
	if (PyCallable_Check(moduleFlowDirFunc))
	{
		input_x = PyList_New(NumNodes);
		input_y = PyList_New(NumNodes);
		input_z = PyList_New(NumNodes);
		triangles = PyList_New(NumEl);

		for (int k = 0; k < NumNodes; k++)
		{
			PyList_SetItem(input_x, k, PyFloat_FromDouble(FloodplainList[fpNum]->Vx[k]));
			PyList_SetItem(input_y, k, PyFloat_FromDouble(FloodplainList[fpNum]->Vy[k]));
			PyList_SetItem(input_z, k, PyFloat_FromDouble(FloodplainList[fpNum]->Vz[k]));
		}
	
		for (int k = 0; k < NumEl; k++)
		{
			PyObject *currTri = PyList_New(3);
			for (int j = 0; j < 3; j++)
			{
				PyList_SetItem(currTri, j, PyInt_FromLong(FloodplainList[fpNum]->EltoVert[k*3+j]));
			}
			PyList_SetItem(triangles, k, currTri);
		}

		args = PyTuple_New(4);
		if (!args)
			PyErr_Print();
		PyTuple_SetItem(args, 0, input_x);
		PyTuple_SetItem(args, 1, input_y);
		PyTuple_SetItem(args, 2, input_z);
		PyTuple_SetItem(args, 3, triangles);

		localFlowDir = PyObject_CallObject(moduleFlowDirFunc, args);
		PyErr_Print();

				for (int el = 0; el < NumEl; el++)
		{
			int flowEdg = (int) PyInt_AsLong(PyList_GetItem(localFlowDir, el));
			LocalFlowEdg[el] = flowEdg;
			FlowEdg[el] = FloodplainList[fpNum]->EltoEdg[el*3+flowEdg];
		}
	}
	else
		PyErr_Print();

	// clean up
	Py_Finalize();

}

/**************************************************************************************************/
/**************************************************************************************************/


void write_flowpaths_to_file(int fp, int* FlowEdg)
{
	
	FILE *fpfile;
	char fileName[PATH_MAX + 1];
	if(getcwd(fileName, sizeof(fileName))==NULL)
		printf("Error while getting current working directory\n");

	char* filePath = "/Output/FlowPaths.out";
	strcat(fileName, filePath);

	if (access("Output", F_OK) == -1)
		if (mkdir("Output", 0777) == -1)
		{
			perror("The following error occured while making the directory output");
		}

	fpfile = fopen(fileName, "w");
	fprintf(fpfile, "#FlowEdge written out by process_watershed.c\n");
	for (int i = 0; i < FloodplainList[fp]->NumEl; i++)
	{
		int edg = FlowEdg[i];
		int n1 = FloodplainList[fp]->EdgtoVert[edg*2];
		int n2 = FloodplainList[fp]->EdgtoVert[edg*2+1];
		double xval = 0.5*(FloodplainList[fp]->Vx[n1]+FloodplainList[fp]->Vx[n2]);
		double yval = 0.5*(FloodplainList[fp]->Vy[n1]+FloodplainList[fp]->Vy[n2]);
		fprintf(fpfile,"%lf %lf\n", xval, yval);
	}

	fclose(fpfile);

}

int is_channel_edge(int edg, int floodplainNum)
{
	if(FloodplainList[floodplainNum]->ChannelNumber[edg] > -1)
		return 1;
	else
		return 0;
	
	//int v1 = EdgtoVert[2*edg];
	//int v2 = EdgtoVert[2*edg+1];
	//for (int i = 0; i < numChannelEdges; i++)
	//{
	//	int cv1 = channel_array[i][0];
	//	int cv2 = channel_array[i][1];
	//	if ((cv1 == v1 && cv2 == v2) || (cv1 == v2 && cv2 == v1))
	//			return 1;
	//}
	//return 0;
}

void redirect_flow_to_channels(int floodplainNum, int* FlowEdg, int* LocalFlowEdg)
{
	for (int i = 0; i < FloodplainList[floodplainNum]->TotalNumEdges; i++)
	{
		if (is_channel_edge(i, floodplainNum))
		{
			int ind1 = i*2;
			int ind2 = i*2+1;
			int el1 = FloodplainList[floodplainNum]->EdgtoEls[ind1];
			int el2 = FloodplainList[floodplainNum]->EdgtoEls[ind2];
			FlowEdg[el1] = i;
			FlowEdg[el2] = i;
			int v1 = FloodplainList[floodplainNum]->EdgtoVert[ind1];
			int v2 = FloodplainList[floodplainNum]->EdgtoVert[ind2];
		
			for (int k = 0; k < 3; k++)
			{
				int el1v1 = FloodplainList[floodplainNum]->EltoVert[el1*3+k];
				int el1v2 = FloodplainList[floodplainNum]->EltoVert[el1*3+k+1];
				if (((el1v1 == v1 ) && (el1v2 == v2)) || ((el1v1 == v2 && el1v2 == v1)))
				{
					LocalFlowEdg[el1] = k;
					break;
				}
			}
			for (int k = 0; k < 3; k++)
			{
				int el2v1 = FloodplainList[floodplainNum]->EltoVert[el1*3+k];
				int el2v2 = FloodplainList[floodplainNum]->EltoVert[el1*3+k+1];
				if (((el2v1 == v1 ) && (el2v2 == v2)) || ((el2v1 == v2 && el2v2 == v1)))
				{
					LocalFlowEdg[el2] = k;
					break;
				}
			}
		}
	}
}


/***********************************************************************************/
/************************************************************************************/

/**********************************************************************
 * Routines that help in the creation kinematic (flow field) elements
 * and their storage in KinematicElList 
 * ********************************************************************/


void get_channel_num_and_edg_num(int edg, int floodplainNum, int* returnArray)
{
	returnArray[0] = FloodplainList[floodplainNum]->ChannelNumber[edg];
	returnArray[1] = FloodplainList[floodplainNum]->ChannelElNumber[edg];

	//int v1 = EdgtoVert[2*edg];
	//int v2 = EdgtoVert[2*edg+1];
	//for(int i = 0; i < numChannelEdges; i++)
	//{
	//	int cv1 = channel_array[i][0];
	//	int cv2 = channel_array[i][1];
	//	if ((cv1 == v1 && cv2 == v2) || (cv1 == v2 && cv2 == v1))
	//	{
	//		returnArray[0] = channel_array[i][2];
	//		returnArray[1] = channel_array[i][3];
	//		return;
	//	}    
	//}

}

void create_kinematic_elements(int fpNum, int* FlowEdg, int* LocalFlowEdg)
{
	int NumEl = FloodplainList[0]->NumEl;
	KinematicElList = (struct kinematicEl**) xcalloc(NumEl,sizeof(struct kinematicEl*));
	Valence = (int *) xcalloc(NumEl, sizeof(int));
	double *areas = (double*) xcalloc(NumEl, sizeof(double));

	double totalArea = 0;
	for (int i = 0; i < NumEl; i++)
	{
		KinematicElList[i] = (struct kinematicEl*) xcalloc(1,sizeof(struct kinematicEl));
		KinematicElList[i]->numDownstreamEls=0;
		KinematicElList[i]->numUpstreamEls=0;
		KinematicElList[i]->isActive = 1;

		//calculate elemental areas
		double x1 = FloodplainList[fpNum]->Vx[FloodplainList[fpNum]->EltoVert[3*i]];
		double x2 = FloodplainList[fpNum]->Vx[FloodplainList[fpNum]->EltoVert[3*i+1]];
		double x3 = FloodplainList[fpNum]->Vx[FloodplainList[fpNum]->EltoVert[3*i+2]];
		double y1 = FloodplainList[fpNum]->Vy[FloodplainList[fpNum]->EltoVert[3*i]];
		double y2 = FloodplainList[fpNum]->Vy[FloodplainList[fpNum]->EltoVert[3*i+1]];
		double y3 = FloodplainList[fpNum]->Vy[FloodplainList[fpNum]->EltoVert[3*i+2]];

		double area = 0.5*fabs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
		areas[i] = area;
		//printf("area = %lf\n", area);
		if (area == 0)
			printf("x1 = %lf x2 = %lf x3 = %lf y1 = %lf y2 = %lf y3 = %lf\n", x1, x2, x3, y1, y2, y3);
		totalArea+= area;
	}

	int numInactiveEls = 0;
	for (int i =0; i < NumEl; i++)
	{
		int flowEdg = FlowEdg[i];
		int el1 = FloodplainList[fpNum]->EdgtoEls[2*flowEdg];
		int el2 = FloodplainList[fpNum]->EdgtoEls[2*flowEdg+1];
		int nextEl = (i == el1) ? el2 : el1;

		// if flow is routed towards boundary, find another edge with the lowest z value
		if (nextEl == i)
		{

		  int local_flow_edg = LocalFlowEdg[i];
			int local_edg1 = (local_flow_edg+1)%3;
			int local_edg2 = (local_flow_edg+2)%3;
			int edg1 = FloodplainList[fpNum]->EltoEdg[3*i+local_edg1];
			int edg2 = FloodplainList[fpNum]->EltoEdg[3*i+local_edg2];
			int n11 = FloodplainList[fpNum]->EdgtoVert[edg1*2];
			int n12 = FloodplainList[fpNum]->EdgtoVert[edg1*2+1];
			int n21 = FloodplainList[fpNum]->EdgtoVert[edg2*2];
			int n22 = FloodplainList[fpNum]->EdgtoVert[edg2*2+1];
			double myz1 = 0.5*(FloodplainList[fpNum]->Vz[n11] + FloodplainList[fpNum]->Vz[n12]);
			double myz2 = 0.5*(FloodplainList[fpNum]->Vz[n21] + FloodplainList[fpNum]->Vz[n22]);

			if (myz1 < myz2)
			{
				LocalFlowEdg[i] = local_edg1;
				FlowEdg[i] = edg1;
			}
			else
			{
				LocalFlowEdg[i] = local_edg2;
				FlowEdg[i] = edg2;
			}

			flowEdg = FlowEdg[i];
			el1 = FloodplainList[fpNum]->EdgtoEls[2*flowEdg];
			el2 = FloodplainList[fpNum]->EdgtoEls[2*flowEdg+1];
			nextEl = (i == el1) ? el2 : el1;

		}

		int nextFlowEdg = FlowEdg[nextEl];

		// if the next element sends the flow back into this element
		if (!is_channel_edge(flowEdg, fpNum) && (nextFlowEdg == flowEdg))
		{
			int local_flow_edg = LocalFlowEdg[i];
			int local_edg1 = (local_flow_edg+1)%3;
			int local_edg2 = (local_flow_edg+2)%3;
			int edg1 = FloodplainList[fpNum]->EltoEdg[3*i+local_edg1];
			int edg2 = FloodplainList[fpNum]->EltoEdg[3*i+local_edg2];
			int n11 = FloodplainList[fpNum]->EdgtoVert[edg1*2];
			int n12 = FloodplainList[fpNum]->EdgtoVert[edg1*2+1];
			int n21 = FloodplainList[fpNum]->EdgtoVert[edg2*2];
			int n22 = FloodplainList[fpNum]->EdgtoVert[edg2*2+1];
			double myz1 = 0.5*(FloodplainList[fpNum]->Vz[n11] + FloodplainList[fpNum]->Vz[n12]);
			double myz2 = 0.5*(FloodplainList[fpNum]->Vz[n21] + FloodplainList[fpNum]->Vz[n22]);
			int rem_edg;
			int rem_local_edg;

			if (myz1 < myz2)
			{
				LocalFlowEdg[i] = local_edg1;
				FlowEdg[i] = edg1;
				rem_edg = edg2;
				rem_local_edg = local_edg2;
			}
			else
			{
				LocalFlowEdg[i] = local_edg2;
				FlowEdg[i] = edg2;
				rem_edg = edg1;
				rem_local_edg = local_edg1;
			}

			flowEdg = FlowEdg[i];
			el1 = FloodplainList[fpNum]->EdgtoEls[2*flowEdg];
			el2 = FloodplainList[fpNum]->EdgtoEls[2*flowEdg+1];
			nextEl = (i == el1) ? el2 : el1;

			// if this next edge is a boundary edge, then pick the last remaining edge
			// if the next element picked then will also send the flow back to this element, 
			// then for now, we just turn this element off.
			if (nextEl == i)
			{
				flowEdg = rem_edg;
				FlowEdg[i] = rem_edg;
				LocalFlowEdg[i] = rem_local_edg;
				el1 = FloodplainList[fpNum]->EdgtoEls[2*flowEdg];
				el2 = FloodplainList[fpNum]->EdgtoEls[2*flowEdg+1];
				nextEl = (i == el1) ? el2 : el1;

			}

		}
		
		nextFlowEdg = FlowEdg[nextEl];

		if (!is_channel_edge(flowEdg,fpNum) && KinematicElList[i]->isActive)
		{
			
			KinematicElList[i]->node1Area = areas[i];
			KinematicElList[i]->node2Area = areas[nextEl];

			int n11 = FloodplainList[fpNum]->EdgtoVert[flowEdg*2];
			int n12 = FloodplainList[fpNum]->EdgtoVert[flowEdg*2+1];
			int n21 = FloodplainList[fpNum]->EdgtoVert[nextFlowEdg*2];
			int n22 = FloodplainList[fpNum]->EdgtoVert[nextFlowEdg*2+1];

			double myx1 = 0.5*(FloodplainList[fpNum]->Vx[n11] + FloodplainList[0]->Vx[n12]);
			double myy1 = 0.5*(FloodplainList[fpNum]->Vy[n11] + FloodplainList[0]->Vy[n12]);
			double myx2 = 0.5*(FloodplainList[fpNum]->Vx[n21] + FloodplainList[0]->Vx[n22]);
			double myy2 = 0.5*(FloodplainList[fpNum]->Vy[n21] + FloodplainList[0]->Vy[n22]);
			double myz1 = 0.5*(FloodplainList[fpNum]->Vz[n11] + FloodplainList[0]->Vz[n12]);
			double myz2 = 0.5*(FloodplainList[fpNum]->Vz[n21] + FloodplainList[0]->Vz[n22]);
			double mynf1 = 0.5*(FloodplainList[fpNum]->VnFriction[n11]+FloodplainList[0]->VnFriction[n12]);
			double mynf2 = 0.5*(FloodplainList[fpNum]->VnFriction[n21]+FloodplainList[0]->VnFriction[n22]);

			//int n11 = EltoVert[i*3];
			//int n12 = EltoVert[i*3+1];
			//int n13 = EltoVert[i*3+2];

			//int n21 = EltoVert[nextEl*3];
			//int n22 = EltoVert[nextEl*3+1];
			//int n23 = EltoVert[nextEl*3+2];
			//
			//double denom = 1.0/3;
			//double myx1 = (x_coord[n11]+x_coord[n12]+x_coord[n13])*denom;
			//double myy1 = (y_coord[n11]+y_coord[n12]+y_coord[n13])*denom;
			//double myx2 = (x_coord[n21]+x_coord[n22]+x_coord[n23])*denom;
			//double myy2 = (y_coord[n21]+y_coord[n22]+y_coord[n23])*denom;
			//double myz1 = (z[n11]+z[n12]+z[n13])*denom;
			//double myz2 = (z[n21]+z[n22]+z[n23])*denom;
			//double mynf1 = (nFriction[n11]+nFriction[n12]+nFriction[n13])*denom;
			//double mynf2 = (nFriction[n21]+nFriction[n22]+nFriction[n23])*denom;


			//printf("nextEl4 = %d\n", nextEl);
			KinematicElList[i]->x1 = myx1;
			KinematicElList[i]->y1 = myy1;

			KinematicElList[i]->z1 = myz1; 
			KinematicElList[i]->nf1 = mynf1;

			KinematicElList[i]->x2 = myx2;
			KinematicElList[i]->y2 = myy2;
			KinematicElList[i]->z2 = myz2;
			KinematicElList[i]->nf2 = mynf2;

			KinematicElList[i]->dh = sqrt((myx2-myx1)*(myx2-myx1) + (myy2-myy1)*(myy2-myy1));

			if (KinematicElList[i]->dh == 0)
			{
				printf("Inactivating element %d because it is a sink \n", i);
				printf("myEl = %d nextEl = %d\n", i, nextEl);
				printf("flowEdg = %d nextFlowEdg = %d\n", flowEdg, nextFlowEdg);
				printf("n11 = %d n12 = %d\n", n11 ,n12);
				printf("midx = %lf, midy = %lf \n", myx1, myy1);
				KinematicElList[i]->isActive = 0;
				numInactiveEls++;
				continue;
			}

			KinematicElList[i]->el1 = i;
			KinematicElList[i]->el2 = nextEl;
			//printf("i = %d el1 = %d el2 = %d nextEl = %d\n", i, KinematicElList[i]->el1, KinematicElList[i]->el2, nextEl);
			Valence[i] += 1;
			Valence[nextEl] += 1;

			int ind = KinematicElList[nextEl]->numUpstreamEls;
			KinematicElList[nextEl]->upstreamEls[ind] = i;
			KinematicElList[nextEl]->numUpstreamEls +=1;

			if(is_channel_edge(nextFlowEdg, fpNum))
			{
				int ChannelInfo[2];
				get_channel_num_and_edg_num(nextFlowEdg, fpNum, &ChannelInfo[0]);
				KinematicElList[i]->connChanNum = ChannelInfo[0];
				KinematicElList[i]->connChanEl = ChannelInfo[1];

				KinematicElList[i]->numDownstreamEls = 0;
				if (KinematicElList[nextEl]->isActive)
				{
					//printf("Inactivating kinematic element because element %d is connected to a channel edge\n", nextEl);
					KinematicElList[nextEl]->isActive=0;
					numInactiveEls++;
				}

			}
			else
			{
				KinematicElList[i]->numDownstreamEls = 1;
				KinematicElList[i]->connChanNum = -1;
				KinematicElList[i]->connChanEl = -1;
			}

		}
		else if (is_channel_edge(flowEdg, fpNum) && KinematicElList[i]->isActive)
		{
			KinematicElList[i]->isActive = 0;
			numInactiveEls++;
		}

	}
	NumKinematicEls = NumEl - numInactiveEls;
	//printf("Number of inactive elements = %d\n", numInactiveEls);
	//printf("Number of flow field elements = %d\n" , NumKinematicEls);
	//printf("Num Els = %d\n", NumEl);
	//printf("Num Channel Edges = %d\n", numChannelEdges);

	//free(areas);

	double length = 0;
	for (int i = 0; i < NumEl; ++i)
	{
		if (KinematicElList[i]->isActive)
			length += KinematicElList[i]->dh;
	}

	printf("Total Area = %lf\n", totalArea);
	printf("Total Length = %lf\n", length);
	printf("weq = %lf\n", totalArea/length);


	for (int i = 0; i < NumEl; ++i)
	{
		struct kinematicEl *myKinEl = KinematicElList[i];
		if (myKinEl->isActive)
		{
			int el1 = myKinEl->el1;
			int el2 = myKinEl->el2;
			myKinEl->weq = totalArea/length;
			//if (myKinEl->numUpstreamEls == 2)
			//	printf("i = %d upstream el1 = %d upstream el2 = %d\n", i, myKinEl->upstreamEls[0], myKinEl->upstreamEls[1]);
			//if (myKinEl->numDownstreamEls == 1)
			//	printf("dowsntream el = %d\t", el2);
			//Printf("\n");

		}

	}
}

/********************************************************************************************/

// free arrays that haven't been freed up yet
void free_storage()
{
	//for (int i = 0; i < numChannelEdges; ++i)
	//	free(channel_array[i]);

	//free(channel_array);

	free(LocalFlowEdg);
	free(FlowEdg);
}

// so far only works for 1 floodplain domain
void process_watershed(char* FloodplainMesh, char* ChannelCoordinates, char* JunctionCoordinates)
{
	create_channel_network(ChannelCoordinates, JunctionCoordinates);
	create_floodplain(FloodplainMesh);
	for (int fp = 0; fp < NumFloodplains; fp++)
	{
		int NumEl = FloodplainList[fp]->NumEl;
		int* FlowEdg = xcalloc(NumEl, sizeof(int));
		int* LocalFlowEdg = xcalloc(NumEl, sizeof(int));

		calculate_flow_edge(fp, FlowEdg, LocalFlowEdg);
		redirect_flow_to_channels(fp, FlowEdg, LocalFlowEdg);
	
		create_kinematic_elements(fp, FlowEdg, LocalFlowEdg);
		write_flowpaths_to_file(fp, FlowEdg);

		free(LocalFlowEdg);
		free(FlowEdg);
	}
}
