#include "cell_list.h"
int** cell_list(double **state,int N,double L,double sigma){
	int nx = (int)L/(2.5*sigma);
	int** cells = malloc(nx*nx*nx*sizeof(int*));
	return cells;
}