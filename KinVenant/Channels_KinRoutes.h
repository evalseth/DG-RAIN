#ifndef CHANNELS_KINROUTE

#define CHANNELS_KINROUTE

struct channel
{
	int NumEl;
	int NumEdges;
	int P;
	int Np;
	int NumNodes;

	double* x;
	double* y;
	double* z;
	double* b;
	double* m1;
	double* m2;

	double* dh;
	double* nFriction;

	double* NodalX;
	double* NodalY;
	double* NodalZ;
	double* NodalB;
	double* Nodalm1;
	double* Nodalm2;
	double* NodalnFriction;

	double* db;
	double* dm1;
	double* dm2;
	double* dz;

	double* A;
	double* Q;
	double* qL;
	
	double dt;
	double max_lambda;
	double mindh;

	int* GlobalNodeNum;

	double* beta;

	double collectedQL;

};

struct kinematicRoute
{
	
	int NumEl;
	int NumEdges;
	int NumNodes;
	int P;
	int Np;

	double* x;
	double* y;
	double* z;
	double* nFriction;
	double* dh;

	double* NodalX;
	double* NodalY;
	double* NodalZ;
	double* NodalnFriction;
	double* dz;

	double* H;
	double dt;
	double max_lambda;
	double mindh;

	int connectedChan;
	int connectedChanEl;

};


#endif
