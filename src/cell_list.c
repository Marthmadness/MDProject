#include "md.h"
#include "cell_list.h"
int** cell_list(double **state,int N, int nx, double L,double sigma){
	int total = nx*nx*nx;
	int** cells = malloc(total*sizeof(int*));
	int* counter = malloc(total*sizeof(int));
	for(int i = 0; i<total; i++){
		cells[i] = malloc(N*sizeof(int*));
		counter[i]=0;
		for(int j = 0; j<total; j++){
			cells[i][j] = -1;
		}
	}
	for(int i = 0; i<N; i++){
		int xindex = floor(state[i][0]*nx/L);
		int yindex = floor(state[i][1]*nx/L);
		int zindex = floor(state[i][2]*nx/L);
		int index = xindex + nx*yindex + nx*nx*zindex;
		cells[index][counter[index]]=i;
		counter[index]++;
	}
	return cells;
}
int find_neighbor(int nx, int index, int delx, int dely, int delz){
	int indexc[3];
	indexc[2] = index/(nx*nx);
	indexc[1] = (index%(nx*nx))/nx;
	indexc[0] = index%nx;
	indexc[0] += delx;
	indexc[1] += dely;
	indexc[2] += delz;
	for(int i =0; i<3; i++){
		if(indexc[i]<0){
			indexc[i] = nx-1;
		}
		if(indexc[i]==nx){
			indexc[i] =0;
		}
	}
	return indexc[0]+nx*indexc[1]+nx*nx*indexc[2];
}