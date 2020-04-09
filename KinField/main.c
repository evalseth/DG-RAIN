/********************************************************************************//**
* @file main.c
*
* This file is the driver for ChanNet. When wetting and drying is on, the 
* minimum water threshold and the momentum threshold are also allocated in this file. 
* This needs to be changed later so that it will be a parameter than can be read 
* from a file.
************************************************************************************/
#include </usr/local/Cellar/python/3.7.7/Frameworks/Python.framework/Versions/3.7/include/python3.7m/Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SimulationSteps.h"
#include "Globals.h"

const double g = 9.810000000000;	/**< gravitational acceleration constant */
//const double g = 32.18; /***< gravitational acceleration constant in ft^2/s */

#ifdef WDON
const double H0 = 1e-3;     		/**< Minimum water height threshold that is maintained in all 
										 elements. Nodes whose water height is below this threshold 
										 are considered dry. */
const double VELZERO = 1e-2;		/**< a momentum threshold. If the magnitude of the momentum 
								  		at a dry node is below this threshold value, then this
								  		momentum will not be transferred to another wet node in
								  		the element.*/

#endif


int main(int argc, char **argv)
{

	setenv("PYTHONPATH", ".", 1);
	if (argc < 2)
	{
		printf("Need to provide at least the grid file\n");
		printf("Usage: ./KinField 'fort.14' 'ChannelCoordinates.in' (optional) 'JunctionMesh.in' (optional)\n");
		exit(EXIT_FAILURE);

	}
	
	char *fullMesh = argv[1];

	char *ChannelCoordinates;
	char *JunctionCoordinates;

	if (argc > 2)
		ChannelCoordinates = argv[2];
	else
		ChannelCoordinates = NULL;
	if (argc > 3)
		JunctionCoordinates = argv[3];
	else
		JunctionCoordinates = NULL;

	printf("Processing watershed file %s\n", fullMesh);
	printf("Channel Coordinates file = %s\n", ChannelCoordinates);
	printf("Junction Coordinates file = %s\n", JunctionCoordinates);
	// extract channel edges and create kinematic elements from the full mesh
	process_watershed(fullMesh, ChannelCoordinates, JunctionCoordinates);

	setup_1D_domains();
	setup_2D_domains();

	double FinalTime;
	printf("Enter the time you would like to run the simulation till:\n");
	scanf("%lf", &FinalTime);
	
	// initialize
	initialize_channels();
	initialize_junctions();	
	initialize_kinematicEls();
	//initialize_floodplains();

	// Step through time 

	time_evolution(FinalTime);

	return(0);

}


 
	
	
