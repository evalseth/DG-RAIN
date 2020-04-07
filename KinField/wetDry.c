/*****************************************************************************************//**
* @file wetDry.c
*
* This file contains all the functions associated with the wetting and drying treatment.
*
* ****************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Watershed.h"

#ifdef WDON

/***********************************************************************************//**
* 
* A function that checks to see if a value is positive.
* @param a a double value to be tested
* @return 1 if positive, 0 if negative
* **************************************************************************************/
int BigTheta(double a)
{
	if (a > 0)
		return(1);
	else
		return(0);
}

/*********************************************************************************//**
*
* A function that assigns a wet or dry status to each element of a channel and stores
* it in channel::WD. 0 represents a dry element and 1 represents a wet element.
* @param[in] Chan pointer to the channel structure that we are currently working on
* **********************************************************************************/
void wetDryStatusChan(struct channel *Chan)
{
	for (int i =0; i < Chan->NumEl; ++i)
	{
		double A0 = Chan->A[2*i+1];
		double A1 = Chan->A[2*i+2];
		double b0 = Chan->b[i];
		double b1 = Chan->b[i+1];
		double h0 = A0/b0;
		double h1 = A1/b1;
		double avgH = (h0+h1)/2;

		int maxHeightIndex;
		double maxHeight;
		if (h0 > h1)
		{
			maxHeightIndex = i;
			maxHeight = h0;
		}
		else
		{
			maxHeightIndex = i+1;
			maxHeight = h1;
		}
		
		double ZetaMax = maxHeight - Chan->z[maxHeightIndex];
		double z0 = Chan->z[i];
		double z1 = Chan->z[i+1];
		double zmin = z0*(z0 < z1) + z1*(z1 < z0); 

		if (Chan->WD[i] == -1) 	// if uninitialized
			Chan->WD[i] = BigTheta(avgH-H0)*BigTheta(ZetaMax - (H0-zmin));
		else
			Chan->WD[i] = Chan->WD[i]*BigTheta(avgH - H0)+(1-Chan->WD[i])*BigTheta(avgH-H0)*BigTheta(ZetaMax - (H0 - zmin));

	}


}

/**************************************************************************************//**
*
* This function does some post processing to ensure that the water depth remains positive
* in the 1-D channels
* @param[in] Chan pointer to the channel structure that we are currently working on
*
* *****************************************************************************************/
void PDopChan(struct channel *Chan)
{

	// redistribute mass as necessary
	for (int i =0; i < Chan->NumEl; ++i)
	{
		double A0 = Chan->A[2*i+1];
		double A1 = Chan->A[2*i+2];
		double h0 = A0/Chan->b[i];
		double h1 = A1/Chan->b[i+1];
		double avgH = (h0+h1)/2;
		
//		printf("avgH = %e \n", avgH);

		if (h0 >= H0 && h1 >= H0)
		{	
			h0 = h0;
			h1 = h1;
		}		
	
		else if (avgH <= H0)
		{
			h0 = avgH;
			h1 = avgH;
//			printf("avgH is less than H0\n");
		}

		else
		{
//			printf("h0 = %e \t h1 = %e\n", h0, h1);
//			printf("avgH is greater than H0\n");
			if (h0 <= h1)
			{	
				double oldh0 = h0;
				h0 = H0;
				h1 = h1 - (h0 - oldh0);
			}

			else
			{
				double oldh1 = h1;
				h1 = H0;
				h0 = h0 - (h1-oldh1);
			}

		}

		A0 = h0*Chan->b[i];
		A1 = h1*Chan->b[i+1];
		Chan->A[2*i+1] = A0;
		Chan->A[2*i+2] = A1;
		if (isnan(A0) || isnan(A1))
		{
			printf("water height negative during wetting and drying treatment for el %d with avgH = %e\n", i, avgH);
			exit(EXIT_FAILURE);
		}
	
	}			

	// redistribute momentum as necessary
	for (int el = 0; el < Chan->NumEl; ++el)
	{
		int wetDryFlagNode[2];
		double u[2];
		for (int i =0; i < 2; ++i)
		{
			double Anode = Chan->A[2*el+i+1];
			double bnode = Chan->b[el+i];
			double hnode = Anode/bnode;
			wetDryFlagNode[i] = BigTheta(hnode - H0);
			
			double Qnode = Chan->Q[2*el+i+1];
			u[i] = Qnode/Anode;
	
		}
		if (wetDryFlagNode[0] == 0 && wetDryFlagNode[1] == 0)
		{
			u[0] = 0;
			u[1] = 0;


		}

		else
		{
			int npos = wetDryFlagNode[0] + wetDryFlagNode[1];
			//double delu = u[0]*(1-wetDryFlagNode[0]) + u[1]*(1-wetDryFlagNode[1]);
			double delu = u[0]*(1-wetDryFlagNode[0])*(fabs(Chan->Q[2*el+1]>VELZERO)) + u[1]*(1-wetDryFlagNode[1])*(fabs(Chan->Q[2*el+2])>VELZERO);
			for (int i =0; i < 2; ++i)
			{
				u[i] = wetDryFlagNode[i]*(u[i] + delu/npos);
			}
			

		}
		Chan->Q[2*el+1] = Chan->A[2*el+1]*u[0];
		Chan->Q[2*el+2] = Chan->A[2*el+2]*u[1];


		if (isnan(u[0]) || isnan(u[1]))
			printf("WD-modified velocity NAN for element %d \n",el);


	}

}

/*********************************************************************************//**
*
* A function that assigns a wet or dry status to each element of a junction and stores
* it in junction::WD. 0 represents a dry element and 1 represents a wet element.
* @param[in] junc pointer to the the junction structure that we are currently working on
* **********************************************************************************/
// only works for P=1 elements
void wetDryStatusJunc(struct TwoDRegion *junc)
{
	for (int el = 0; el < junc->NumEl; ++el)
	{
		double zetaEl[3] = {junc->zeta[el][0], junc->zeta[el][1], junc->zeta[el][2]};
		double zEl[3] = {junc->NodalZ[el][0], junc->NodalZ[el][1], junc->NodalZ[el][2]};

		double height[3] = {zetaEl[0] + zEl[0], zetaEl[1] + zEl[1], zetaEl[2] + zEl[2]};
		double avgH = (height[0]+height[1]+height[2])/3;

		int maxHeightIndex = 0;
		double maxHeight = height[0];

		for (int i =1; i < 3; ++i)
		{
			if (height[i] > maxHeight)
			{
				maxHeightIndex = i;
				maxHeight = height[i];
			} 
		}

		double ZetaMax = zetaEl[maxHeightIndex];
		
		double zmin = zEl[0];
		for (int i = 1; i < 3; ++i)
		{
			if (zEl[i] < zmin)
			{
				zmin = zEl[i];
			} 

		}

		if (junc->WD[el] == -1)		// if uninitialized
			junc->WD[el] = BigTheta(avgH - H0)*BigTheta(ZetaMax - (H0 - zmin));
		else
			junc->WD[el] = junc->WD[el]*BigTheta(avgH - H0) + (1-junc->WD[el])*BigTheta(avgH - H0)*BigTheta(ZetaMax - (H0 - zmin));


	}

}

/**************************************************************************************//**
*
* This function does some post processing to ensure that the water depth remains positive
* in the 2-D junctions
* @param[in] junc pointer to the junction structure that we are currently working on
*
* *****************************************************************************************/
void PDopJunc(struct TwoDRegion *junc)
{
	// redistribute mass and momentum as necessary
	for (int el = 0; el < junc->NumEl; ++el)
	{
		
		double zeta[3] = {junc->zeta[el][0], junc->zeta[el][1], junc->zeta[el][2]};
		double z[3] = {junc->NodalZ[el][0], junc->NodalZ[el][1], junc->NodalZ[el][2]};
		double avgH = 0;
		double Ht[3];
		for (int s = 0; s < 3; s++)
		{
			Ht[s] = zeta[s] + z[s];
			avgH += Ht[s];
		}
		avgH /= 3;

		// redistribute mass
		if (Ht[0] <= H0 || Ht[1] <= H0 || Ht[2] <= H0)
		{
			if (avgH <= H0)
			{
				for (int i =0; i < 3; ++i)
					Ht[i] = avgH;	
			}

			else
			{
				int minHtIndex = 0;
				double minHt = Ht[0];
				for (int i = 1; i < 3; ++i)
				{
					if (Ht[i] < minHt)
					{
						minHtIndex = i;
						minHt = Ht[i];
					} 

				}
				
				int nextMinHtIndex = 1;
				double nextMinHt = Ht[1];
				for (int i = 0; i < 3; ++i)
				{
					if (Ht[i] < nextMinHt && i != minHtIndex)
					{
						nextMinHtIndex = i;
						nextMinHt = Ht[i];
					}
				}

				int remIndex = 2;
				for (int i =0; i <3; ++i)
				{
					if (i != minHtIndex && i != nextMinHtIndex)
					{
						remIndex = i;
						break;
					}
				}
				
				Ht[minHtIndex] = H0;
				Ht[nextMinHtIndex] = fmax(H0, nextMinHt - 0.5*(Ht[minHtIndex] - minHt));
				Ht[remIndex] = Ht[remIndex] - (Ht[minHtIndex] - minHt) - (Ht[nextMinHtIndex] - nextMinHt); 

			}

			if (isnan(Ht[0]) || isnan(Ht[1]) || isnan(Ht[2]))
			{
				printf("water height negative during wetting and drying treatment for el %d with avgH = %e\n", el, avgH);
				exit(EXIT_FAILURE);
			}

			// Assign the nodal values
			junc->zeta[el][0] = Ht[0];
			junc->zeta[el][1] = Ht[1];
			junc->zeta[el][2] = Ht[2];
		}

		int NodalWDFlag[3];
		
		// Store the wet and dry nodes
		for (int i = 0; i < 3; ++i)
		{
			NodalWDFlag[i] = (Ht[i] > H0);

		}		
		
		// redistribute momentum if necessary

		
		if (NodalWDFlag[0] == 0 && NodalWDFlag[1] == 0 && NodalWDFlag[2] == 0)
		{
			for (int i = 0; i < 3; ++i)
			{
				junc->Qx[el][i] = 0;
				junc->Qy[el][i] = 0;
	
			}
		}

		else if (NodalWDFlag[0] == 0 || NodalWDFlag[1] == 0 || NodalWDFlag[2] == 0)
		{
			// Store the nodal values of the velocities in each direction
			double u[3];
			double v[3];

			double Qx[3] = {junc->Qx[el][0], junc->Qx[el][1], junc->Qx[el][2]};
			double Qy[3] = {junc->Qy[el][0], junc->Qy[el][1], junc->Qy[el][2]};

			for (int i = 0; i < 3; ++i)
			{
				u[i] = Qx[i]/Ht[i];
				v[i] = Qy[i]/Ht[i];
			}	


			int npos = 0;
			double delu = 0;
			double delv = 0;
			for (int i = 0; i < 3; ++i)
			{
				npos += NodalWDFlag[i];
				delu += u[i]*(1-NodalWDFlag[i])*(fabs(Qx[i]) > VELZERO);
				delv += v[i]*(1-NodalWDFlag[i])*(fabs(Qy[i]) > VELZERO);
				

			}			

			for (int i = 0; i < 3; ++i)
			{
				u[i] = NodalWDFlag[i]*(u[i]+delu/npos);
				v[i] = NodalWDFlag[i]*(v[i]+delv/npos);

				if (isnan(u[i]) || isnan(v[i]))
				{
					printf("WD-modified velocity NAN for element %d \n",el);
					exit(EXIT_FAILURE);
				}	
			}
			for (int i = 0; i < 3; ++i)
			{
				Qx[i] = Ht[i]*u[i];
				Qy[i] = Ht[i]*v[i];

			}
	
			for (int i = 0; i < 3; ++i)
			{
				junc->Qx[el][i] = Qx[i];
				junc->Qy[el][i] = Qy[i];

			}

		}

	
	}

}

#endif









