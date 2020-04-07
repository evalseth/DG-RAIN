#include "globals.h"

void boundary_conditions()
{
	H[0] = 2;
	H[TotalNumNodes-1] = H[TotalNumNodes-2];

}
