#ifndef TWODPROCESSINGROUTINES

#define TWODPROCESSINGROUTINES

void read2DGridFile(char *gridFile, char* regionType);
int createEdgtoVertAndElCount(int* EltoVert, int currNumEl, int NumVerts, int* EdgtoVert, int* ElCount);
void createEdgtoEls(int *EltoVert, int* EdgtoVert, int currNumEl, int NumEdges, int* EdgtoEls);
void createGlobaltoLocalEdg(int* EltoVert, int* EdgtoVert, int *EdgtoEls, int NumEdges, int* GlobaltoLocalEdg);
void createVtoEl(int* EltoVert, int* ElCount, int NumEl, int NumVerts, int** VtoEl);

#endif
