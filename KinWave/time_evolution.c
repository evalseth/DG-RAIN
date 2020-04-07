#include "globals.h"
#include "oneTimeStep.h"

double time;
void oneTimeStep(double time, double dt)
{
	// Intermediate RK step
	computeKinL(time);

	double oldH[TotalNumNodes];
	for (int i = 0; i < TotalNumNodes; i++)
	{
		oldH[i] = H[i];
		H[i] = H[i] + dt*KinRHS[i];
		//printf("KinRHS = %3.14f \t H = %3.14f\n", KinRHS[i], H[i]);
	}

	// apply slope limiter
	minmod();

	boundary_conditions();

	// Second stage RK 
	computeKinL(time);

	for (int i =0; i < TotalNumNodes; i++)
	{
		H[i] = 0.5*(oldH[i]+H[i]+dt*KinRHS[i]);
	}

	// slope limiter
	minmod();

	boundary_conditions();

}

void time_evolution(double FinalTime)
{
	time = 0.0; 
	dt = 1e-4;
	int Nstep = 0;

	while(time < FinalTime)
	{
		time = time + dt;
		oneTimeStep(time, dt);
		boundary_conditions();
	}

	//dt = 0.1*mindh/max_lambda;
	if ((time+dt) > FinalTime)
		dt = FinalTime - time;

}
