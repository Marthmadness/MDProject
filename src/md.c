#include "md.h"
#include "verlet.h"

//initializing stuff
int main(){	
	//cube case
	double L = 100; //length of one side of the box
	int N_core = 3; // number of particles per side in a cubic lattice thing.
	int N = N_core*N_core*N_core;
	int testPart = 0;
	
	double **state;
	state = malloc(N*sizeof(double *));
	double **past;
	past = malloc(N*sizeof(double *));
	
	for(int i=0; i<N; i++){
		state[i] = malloc(6*sizeof(double));
		past[i] = malloc(3*sizeof(double));
		}
	
	int type[N];//for later stuff if I feel like it
	//initialize positions on grid

	//cube stuff
	for(int i = 0; i<N; i++){	
		double ran = 0;//for now keep velocity 0
		past[i][0] = (i%N_core)*L/N_core;
		past[i][1] = ((i/N_core)%N_core)*L/N_core;
		past[i][2] = ((i/(N_core*N_core))%N_core)*L/N_core;
		//in our case, start with slightly randomized velocity, nothing too fancy
		state[i][0] = past[i][0]+(ran*L*(rand()-0.5)/(double)RAND_MAX);
		state[i][1] = past[i][1]+(ran*L*(rand()-0.5)/(double)RAND_MAX);
		state[i][2] = past[i][2]+(ran*L*(rand()-0.5)/(double)RAND_MAX);
		}
	printf("%f %f %f %f \n",0.00,state[testPart][0],state[testPart][1],state[testPart][2]);


	/*
	//for pair interaction
	state[0][1] = L/32;
	state[0][2] = 0;
	state[0][3] = 0;
	state[0][4] = 0;
	state[0][5] = 0;
	state[0][0] = 0;
	state[1][1] = 31*L/32.0;
	state[1][2] = 0;
	state[1][3] = 0;
	state[1][4] = 0;
	state[1][5] = 0;
	state[1][0] = 0;
	past[0][1] = L/32;
	past[0][2] = 0;
	past[0][3] = 0;
	past[0][4] = 0;
	past[0][5] = 0;
	past[0][0] = 0;
	past[1][1] = 31*L/32.0;
	past[1][2] = 0;
	past[1][3] = 0;
	past[1][4] = 0;
	past[1][5] = 0;
	past[1][0] = 0;
	printf("%f %f %f %f \n",t,state[testPart][0],state[testPart][1],state[testPart][2]);
	*/
		
	verletengine(state, past, N, L);
	for(int i=0; i<N; i++){
		free(state[i]);
		free(past[i]);
		}
	free(state);
	free(past);
	return 0;
	}
