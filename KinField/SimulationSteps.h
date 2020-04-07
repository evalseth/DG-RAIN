/**************************************************************************//**
*
* @file SimulationSteps.h
*
* This file contains the function prototypes for functions that represent 
* different stages of the simulation, i.e. function to create the channels
* and junction structure, functions to initialize those structures and the
* function to evolve them in time.
*
* ***************************************************************************/

#ifndef SIMULATION_STEPS

#define SIMULATION_STEPS


/* @cond FUNCTION_PROTOTYPES */
extern void process_watershed(char* fullMesh, char* ChannelCoordinates, char* JunctionCoordinates);
extern void setup_1D_domains();
extern void setup_2D_domains();
extern void initialize_channels();
extern void initialize_junctions();
extern void initialize_floodplains();
extern void initialize_kinematicEls();
extern void time_evolution(double FinalTime);
/* @endcond */

#endif
