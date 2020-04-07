#include <stdio.h>
#include <stdlib.h>
#include <math.h>	

#define PI 3.14159265359

struct node
{
	int num;
	struct node *next;
};


int NumEl;
int NumNodes;
int NumEdges;
double *x_coord;
double *y_coord;
double *z;
double *nFriction;
int *EltoVert;
int *EdgtoEls;
int *EdgtoVert;
int *EltoEdg;
int *FlowEdg;
int boundaryElNum;
struct node *boundaryEl; 

void *xcalloc(int items, int size)
{
	void *ptr = calloc(items, size);
	if (ptr == NULL)
	{
		printf("Unable to allocate memory\n");
		exit(EXIT_FAILURE);
	}
	return ptr;
}

void store_mesh(char *FullMesh)
{
	// read in the co-ordinates
	FILE *Coordinates = fopen(FullMesh,"r");
	
	// throw away the first line of the file
	char firstline[1000];
	fgets(firstline, sizeof(firstline), Coordinates);

	char secondline[1000];
	fgets(secondline, sizeof(secondline), Coordinates);
	sscanf(secondline, "%d  %d", &NumEl, &NumNodes);	
	
	// array to store the x, y and z values at the nodes
	x_coord = xcalloc(NumNodes,sizeof(double));
	y_coord = xcalloc(NumNodes,sizeof(double));
	z = xcalloc(NumNodes,sizeof(double));
	nFriction = xcalloc(NumNodes,sizeof(double));

	// variable to temporarily store the node numbers while reading the file
	int nodes;
	
	int j;
	for (j=0; j<NumNodes; ++j)
	{	
		char grid[1000];
		fgets(grid, sizeof(grid), Coordinates);
		sscanf(grid, "  %d %lf %lf %lf %lf", &nodes, &x_coord[j], &y_coord[j], &z[j], &nFriction[j]);
		z[j] = z[j];

	}

	EltoVert = xcalloc(3*NumEl,sizeof(int));
	int numNodesPerEl;
	int ElNum;
	for (j=0; j < NumEl; ++j)
	{
		char conn_table[1000];
		fgets(conn_table, sizeof(conn_table), Coordinates);
		int v1, v2, v3;
		sscanf(conn_table, "%d %d %d %d %d", &ElNum, &numNodesPerEl, &v1, &v2, &v3);
		// switch to indexing from 0
		EltoVert[j*3+0] = v1-1;
		EltoVert[j*3+1] = v2-1;
		EltoVert[j*3+2] = v3-1;
	}

	fclose(Coordinates);

}


void extract_mesh_attributes()
//void main()
{
	//int NumEl = 9;
	//int NumNodes = 9;

	//int EltoVert[] = {5,2,1,5,1,4,5,4,3,5,3,2};
	//int EltoVert[] = {0,6,7,1,0,7,8,1,7,2,1,8,3,2,8,3,8,7,3,7,4,4,7,5,5,7,6};

	// For each node, count the number of elements that share this node

		int **NtoNodes = xcalloc(NumNodes,sizeof(int*));
		for (int s = 0; s < NumNodes; s++)
		{
			NtoNodes[s] = xcalloc(10,sizeof(int));
		}
		int *NodesCount = xcalloc(NumNodes,sizeof(int));
	
{	
	int *ElCount = xcalloc(NumNodes,sizeof(int));
	for (int n = 0; n < NumNodes; ++n)
	{
		for (int j = 0; j < 10; ++j)
			NtoNodes[n][j] = -1;
	}
	for (int el = 0; el<NumEl; ++el)
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

	free(ElCount);

}	// Euler's relation: NumNodes - NumEdges + Numfaces = 1 -> NumEdges = NumNodes + NumEl - 1;
	NumEdges = NumNodes + NumEl - 1;
	EdgtoEls = xcalloc(2*NumEdges,sizeof(int));
	EdgtoVert = xcalloc(2*NumEdges,sizeof(int));
	EltoEdg = xcalloc(3*NumEl, sizeof(int));


	int EdgCount = 0;
	for (int n = 0; n < NumNodes; ++n)
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

	if (EdgCount != NumEdges)	
	{
		printf("EdgCount = %d\n", EdgCount);
		printf("Total number of edges is not correct\n");
		exit(EXIT_FAILURE);
	}	

	// Create EdgtoEls
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
		for (int el = 0; el < NumEl; ++el)
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
			if ((el == (NumEl-1)) && elfound == 1)
				EdgtoEls[e*2+1] = EdgtoEls[e*2];

		 }

	}

	// Create ElToEdg
	for (int el = 0; el < NumEl; el++)
	{
		EltoEdg[el*3]=-1;
		EltoEdg[el*3+1]=-1;
		EltoEdg[el*3+2]=-1;
	}

	boundaryEl = (struct node *) malloc(sizeof(struct node));
	struct node* boundaryElCounter = boundaryEl;

	for (int e = 0; e < NumEdges; e++)
	{
		static int bEl = 0;
		int el1 = EdgtoEls[e*2];
		int el2 = EdgtoEls[e*2+1];
		int edgNum = 0;
		while(EltoEdg[el1*3+edgNum] > -1  && edgNum < 2)
			edgNum++;

		EltoEdg[el1*3+edgNum] = e;

		if (el2 != el1)
		{
			int edgNum = 0;
			while(EltoEdg[el2*3+edgNum] >0 && edgNum < 2)
				edgNum++;
			EltoEdg[el2*3+edgNum] = e;			
		}
		else
		{
			boundaryElCounter->num = el1; 
			boundaryElCounter->next = (struct node *) malloc(sizeof(struct node));
			boundaryElCounter = boundaryElCounter->next;
			boundaryElNum++;
		}

	}
	boundaryElCounter->num = -1;

	free(NtoNodes);
	free(NodesCount);

}

void calculate_flow_edge()
{
	FlowEdg = xcalloc(NumEl, sizeof(int));
	for(int el = 0; el < NumEl; el++)
		FlowEdg[el] = -1;
	for (int el = 0; el < NumEl; el++)
	{
		// calculate the flow direction
		int begNode = el*3;
		int vrt1 = EltoVert[begNode];
		int vrt2 = EltoVert[begNode+1];
		int vrt3 = EltoVert[begNode+2];

		double x1 = x_coord[vrt1];
		double x2 = x_coord[vrt2];
		double x3 = x_coord[vrt3];
		double y1 = y_coord[vrt1];
		double y2 = y_coord[vrt2];
		double y3 = y_coord[vrt3];
		double z1 = -z[vrt1];
		double z2 = -z[vrt2];
		double z3 = -z[vrt3];

		double A = y1*(z2-z3) + y2*(z3-z1) + y3*(z1-z2);
		double B = z1*(x2-x3) + z2*(x3-x1) + z3*(x1-x2);
		double C = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);

		double flowVecx = B/C;
		double flowVecy = A/C;
		double flowDirAngle = atan2(flowVecy, flowVecx) + PI;

		double x_centroid = (x1+x2+x3)/3;
		double y_centroid = (y1+y2+y3)/3;

		int elVrt1 = EltoVert[el*3];
		int elVrt2 = EltoVert[el*3+1];
		int elVrt3 = EltoVert[el*3+2];

		// calculate the angles enclosed by each edge
		for(int ed = 0; ed< 3; ed++)
		{
			int edg = EltoEdg[el*3+ed];
			int vrt1 = EdgtoVert[edg*2];
			int vrt2 = EdgtoVert[edg*2+1];
				if((vrt1 == elVrt2 && vrt2 == elVrt1))
			{
				vrt1 = elVrt1;
				vrt2 = elVrt2;
			}

			else if ((vrt1 == elVrt3 && vrt2 == elVrt2))
			{
				vrt1 = elVrt2;
				vrt2 = elVrt3;
			}

			else if((vrt1 == elVrt1 && vrt2 == elVrt3))
			{
				vrt1 = elVrt3;
				vrt2 = elVrt1;
			}

			double x1 = x_coord[vrt1];
			double x2 = x_coord[vrt2];
			double y1 = y_coord[vrt1];
			double y2 = y_coord[vrt2];
			double ang0 = atan2(y1-y_centroid,x1-x_centroid) + PI;
			double ang1 = atan2(y2-y_centroid,x2-x_centroid) + PI;


			// order the angles in increasing order
	/*		if (ang0 > ang1)
			{
				double tmp = ang1;
				ang1 = ang0;
				ang0 = tmp;
			}

		*/

		/*	if((vrt1 == elVrt2 && vrt2 == elVrt1) || (vrt1 == elVrt3 && vrt2 == elVrt2))
			{
				
			}
	*/

			// if the flow direction is towards an edge, store the local edge number in
			// FlowEdg, otherwise store the global vertex number + 3 in FlowEdg

			double res1 = flowDirAngle - ang0;
			double res2 = ang1 - flowDirAngle;
			
			if (res1 > 0.01 && res2 > 0.01)
			{	
				FlowEdg[el] = edg;
				break;
			}
			else if(0 <= res1 && res1 <= 0.01)
			{
				FlowEdg[el] = vrt1+3;
				break;
			}
			else if(0 <= res2 && res2 <= 0.01)
			{
				FlowEdg[el] = vrt2+3;
				break;
			}
			else if (ang0 > ang1)
			{
				if( (0 < flowDirAngle && flowDirAngle < ang0) || (flowDirAngle > ang1))
				{
					FlowEdg[el] = edg;
					break;
				}
	
			}
		//	else if (el == 48 || el == 49 || el== 52)
		//		printf("el = %d flowDirAngle = %lf ang1 = %lf ang2 = %lf\n", el, flowDirAngle, ang0, ang1);
		
		}
		
		if(FlowEdg[el] == -1)
			printf("el = %d flowDirAngle = %lf\n", el, flowDirAngle);
		//printf("el = %d edg = %d\n", el, FlowEdg[el]);

	}

}

void add_to_node(struct node* root, int el)
{
	static int numadded = 0;
	struct node* ptr = root;
	int cont = 1;
	int add = 1;
	while (cont)
	{
		int num = ptr->num;
		if (el == num)
		{
			add = 0;
			cont = 0;
		}
		if (num == -1)
		{
			if (add)
			{
				ptr->num = el;
				ptr->next = (struct node*) malloc(sizeof(struct node));
				ptr = ptr->next;
				ptr->num = -1;

				numadded++;
			//	printf("added edge %d ", el);
			//	printf("numadded = %d\n", numadded);
				printf("%d \n%d \n", EdgtoVert[el*2]+1, EdgtoVert[el*2+1]+1);

			}
			cont = 0;
		}
		else
			ptr = ptr->next;
	}

	}

void calculate_channel_edges()
{
	struct node* ChannelEdges = (struct node*) malloc(sizeof(struct node));
	ChannelEdges->num = -1;
	struct node* ChannelEdgeCounter = ChannelEdges;

	for (int i = 0; i < NumEl; i++)
	{
		int flowEdg = FlowEdg[i];
		int edg1 = EltoEdg[i*3];
		int edg2 = EltoEdg[i*3+1];
		int edg3 = EltoEdg[i*3+2];
		
		int el11 = EdgtoEls[edg1*2];
		int el12 = EdgtoEls[edg1*2+1];
		
		int el21 = EdgtoEls[edg2*2];
		int el22 = EdgtoEls[edg2*2+1];

		int el31 = EdgtoEls[edg3*2];
		int el32 = EdgtoEls[edg3*2+1];

		int eln1, eln2, eln3;

		if (el11 == i && el12 !=i)
			eln1 = el12;
		else
			eln1 = el11;

		if (el21 == i && el22 != i)
			eln2 = el22;
		else
			eln2 = el21;

		if (el31 == i && el32 != i)
			eln3 = el32;
		else
			eln3 = el31;

		if (FlowEdg[eln1] == flowEdg)
		{
			if (eln1 != i)
				add_to_node(ChannelEdges, flowEdg);
		}

		else if (FlowEdg[eln2] == flowEdg)
		{
			if (eln2 != i)
				add_to_node(ChannelEdges, flowEdg);
		}

		else if (FlowEdg[eln3] == flowEdg)
		{
			if(eln3 != i)
				add_to_node(ChannelEdges, flowEdg);
		}
		
	
	}

}

void calculate_flow_path()
{
	struct node* FlowPaths = (struct node*) malloc(boundaryElNum*sizeof(struct node));
	struct node* boundaryElCounter = boundaryEl;
	
	struct node* ChannelEdges = (struct node*)malloc(sizeof(struct node));
	struct node* ChannelEdgeCounter = ChannelEdges;

	int i = 0;
	while (boundaryElCounter->num != -1)
	{
		//printf("i = %d\n", i);
		int el = boundaryElCounter->num;
		int cont = 1;

		struct node* FlowPathsTracker = &FlowPaths[i];
		while (cont)
		{
			FlowPathsTracker->num = el;
			printf("curr el = %d  ", el);
			FlowPathsTracker->next = (struct node*) malloc(sizeof(struct node));
			int edg1 = FlowEdg[el];
			int edg2;
			int el1 = EdgtoEls[2*edg1];
			int el2 = EdgtoEls[2*edg1+1];

			printf("el1 = %d el2 = %d ",el1, el2);
			if(el1 == el2)
			{
				cont = 0;
				printf("aborted\n");
			}
			else
			{

				if(el1==el)
				{
					edg2 = FlowEdg[el2];
					el = el2;
				}
				else
				{
					edg2 = FlowEdg[el1];
					el = el1;
				}
				
				//printf("edg1 = %d edge2 = %d\n", edg1, edg2);

				if (edg1 == edg2)
				{
					printf("channel edge\n");
					int vrt1 = EdgtoVert[2*edg1];
					int vrt2 = EdgtoVert[2*edg1+1];
					//printf("%d %d \n", vrt1, vrt2);

					ChannelEdgeCounter->num = edg1;
					ChannelEdgeCounter->next = (struct node*) malloc(sizeof(struct node));
					ChannelEdgeCounter = ChannelEdgeCounter->next;
					cont = 0;
				}
				else
				{
					printf("\n");
					FlowPathsTracker = FlowPathsTracker->next;
				}
			}
		}
		i++;
		boundaryElCounter = boundaryElCounter->next;
	}
	ChannelEdgeCounter->num = -1;
	printf("total num flowpaths = %d\n", i);

	ChannelEdgeCounter = ChannelEdges;
	int ch = 0;
	while(ChannelEdgeCounter->num != -1)
	{
		int edg = ChannelEdgeCounter->num;
		int vrt1 = EdgtoVert[2*edg];
		int vrt2 = EdgtoVert[2*edg+1];
		printf("%d %d \n", vrt1, vrt2);
		ch++;
		ChannelEdgeCounter = ChannelEdgeCounter->next;
	}
	printf("total num channel edges = %d\n", ch);

}

int main()
{
	store_mesh("fort.14");
	extract_mesh_attributes();
	calculate_flow_edge();
	//calculate_flow_path();
	calculate_channel_edges();
}
