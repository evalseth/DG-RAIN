#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include "Channels_KinRoutes.h"
#include "Simulation_Steps.h"

//const double g = 9.81;
const double g = 32.18;

int NumChannels;
int NumKinRoutes;

struct channel **ChannelList;
struct kinematicRoute **KinRouteList;

int main(int argc, char **argv)
{
	if (argc != 3)
	{
		printf("Usage: ./KinVenant 'ChannelNodes.in' 'KinRoutes.in'\n");
		exit(EXIT_FAILURE);
	}

	char *ChannelGrid = argv[1];
	char *KinGrid = argv[2];

	// create channels and kinematic flow routes from the  mesh
	create_flow_paths(ChannelGrid, KinGrid);
	preprocess_basis_functions();

	double FinalTime;
	printf("Enter the time you would like to run the simulation till: \n");
	scanf("%lf", &FinalTime);

	// initialize
	initialize_channels();
	initialize_kinRoutes();

	// Step through time
	time_evolution(FinalTime);



}
