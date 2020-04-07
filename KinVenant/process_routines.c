#include <stdio.h>
#include "Globals.h"
#include "Constitutive_Equations.h"
#include "Channels_KinRoutes.h"

extern void GetLGLWeights(int N, double* w);
extern void* xcalloc(int numItems, int size);

static int fileNumber = 0;

void outputDataToFile(double time)
{
	//Output data 
	FILE* file1;
	FILE* file2;

	for (int i=0; i<NumChannels; i++)
	{
		char filename1[100];
		char filename2[100];
		sprintf(filename1, "%3.3d_ZetaChannel%d.dat", fileNumber, i);
		sprintf(filename2, "%3.3d_Q%d.dat", fileNumber, i);
		file1 = fopen(filename1, "w");
		file2 = fopen(filename2, "w");

		fprintf(file1, "# H and Zeta at time %3.3f for Channel %d\n", time, i); 
		fprintf(file2, "#Q at time %3.3f for Channel %d\n", time, i);

		int NumNodes = ChannelList[i]->NumNodes;	

		for (int j = 0; j<NumNodes; j++)
		{	
			double bval = ChannelList[i]->NodalB[j];
			double zval = ChannelList[i]->NodalZ[j];
			double x_val = ChannelList[i]->NodalX[j];
			double y_val = ChannelList[i]->NodalY[j];
			double A = ChannelList[i]->A[j+1]; 
			double Q = ChannelList[i]->Q[j+1];
			double m1val = ChannelList[i]->Nodalm1[j];
			double m2val = ChannelList[i]->Nodalm2[j];
			double H = getH(A, bval, m1val, m2val); 
			double zeta = H - zval;
			fprintf(file1, "%3.3f\t%3.3f\t%3.3f\t%3.13f \t %3.13f\n", x_val, y_val, -zval, H, zeta);
			fprintf(file2, "%3.3f\t%3.3f\t%3.3f\t%3.13f\n", x_val, y_val, -zval, Q);

		}

		fclose(file1);
		fclose(file2);
	}

	FILE *file3;
	for (int i = 0; i < NumKinRoutes; i++)
	{
		char filename[100];
		sprintf(filename, "%3.3d_ZetaKinRoute%d.dat", fileNumber, i);
		file3 = fopen(filename, "w");
		fprintf(file3, "H and Zeta at time %3.3f for Kinematic Route %d \n", time, i);

		int NumNodes = KinRouteList[i]->NumNodes;
		for (int j = 0; j < NumNodes; j++)
		{
			double xval = KinRouteList[i]->NodalX[j];
			double yval = KinRouteList[i]->NodalY[j];
			double H = KinRouteList[i]->H[j+1];
			double zval = KinRouteList[i]->NodalZ[j];
			double zeta = H - zval;
			fprintf(file3, "%3.3f\t%3.3f\t%3.3f\t%3.13f \t %3.13f\n", xval, yval, -zval, H, zeta);
		}

		fclose(file3);
	}

	fileNumber++;
}

double calculateTotalWater(struct channel *Chan)
{
	int NumNodes = Chan->NumNodes;
	int NumEl = Chan->NumEl;
	int Np = Chan->Np;
	int P = Chan->P;

	double totalWater = 0;
	double LGLWeight[Np];
	GetLGLWeights(P, LGLWeight);
	for (int k = 0; k < NumEl; k++)
	{
		double water = 0;
		double height[Np];
		double dh = Chan->dh[k];
		for (int i = 0; i < Np; i++)
		{
			double m1val = Chan->Nodalm1[k*Np+i];
			double m2val = Chan->Nodalm2[k*Np+i];
			double bval = Chan->NodalB[k*Np+i];
			double A = Chan->A[k*Np+i+1];
			double h = getH(A, bval, m1val, m2val);
			water += LGLWeight[i]*(h-0.0000001);
		}
		totalWater += 0.5*dh*water;
	}

	return totalWater;

}

double ftTom(double ft)
{
	double meter = ft*0.3048;
	return meter;
}

double mToft (double m)
{
	double ft = m*3.28084;
	return ft;

}

