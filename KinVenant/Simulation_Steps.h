#ifndef SIMULATION_STEPS

#define SIMULATION_STEPS

#include "Channels_KinRoutes.h"

extern void create_flow_paths(char* ChannelGrid, char* KinGrid);
extern void preprocess_basis_functions();
extern void initialize_channels();
extern void initialize_kinRoutes();
extern void time_evolution(double FinalTime);

#endif
