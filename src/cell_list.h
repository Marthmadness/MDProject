#ifndef CELL_LIST
#define CELL_LIST
#include "md.h"
int** cell_list(double **state,int N, int nx, double L,double sigma);
int find_neighbor(int nx, int index, int delx, int dely, int delz);
#endif