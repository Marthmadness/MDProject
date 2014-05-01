#include "md.h"
#include "verlet.h"

#define INPUT
#define VV

//initializing stuff
int main(){	
	#ifdef CUBE
	double L = 10; //length of one side of the box
	int N_core = 5; // number of particles per side in a cubic lattice thing.
	int N = N_core*N_core*N_core;
	int testPart = 2;
	#endif
	#ifdef INPUT
	FILE *fp;
	fp = fopen("coords_1_kj.dat","r");
	if(fp==NULL){
		perror("Failed to open file");
		return -1;
	}
	char buff[100];
	fgets(buff, 100, fp);
	fgets(buff, 100, fp);
	fgets(buff, 100, fp); 
	int N;
	fscanf(fp, "%d", &N);
	fgets(buff, 100, fp);
	fgets(buff, 100, fp);
	double xmin, xmax;
	fscanf(fp, "%lf %lf", &xmin, &xmax);
	fgets(buff, 100, fp);
	double L = xmax-xmin;
	printf("%d \n", N);
	printf("%lf \n", L);
	fgets(buff, 100, fp);
	fgets(buff, 100, fp);
	#endif
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
	#ifdef INPUT
	for(int i = 0; i<N; i++){
		int trash;
		double trash2;
		fgets(buff, 100, fp);//advance buffer
		fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf",&trash, &trash, &trash, &state[i][0], &state[i][1], &state[i][2],&trash2, &trash2, &trash2);
		for(int j = 0; j<3; j++){
			if(state[i][j]<0 || state[i][j]>=L){
				state[i][j] = state[i][j] - L*floor(state[i][j]/L); //enforce pbc
			}
		}
		past[i][0] = 0;
		past[i][1] = 0;
		past[i][2] = 0;
		state[i][3] = 0;
		state[i][4] = 0;
		state[i][5] = 0;		
	}
	int testPart = 0;
	fclose(fp);
	#endif
	
	#ifdef CUBE
	//cube stuff
	srand (time(NULL));//seed time
	#ifndef VV
	for(int i = 0; i<N; i++){	
		double ran = 0.0;//for now keep velocity 0
		past[i][0] = (i%N_core)*L/N_core;
		past[i][1] = ((i/N_core)%N_core)*L/N_core;
		past[i][2] = ((i/(N_core*N_core))%N_core)*L/N_core;
		//in our case, start with slightly randomized velocity, nothing too fancy
		state[i][0] = past[i][0]+(ran*L*((rand()/(double)RAND_MAX)-0.5));
		state[i][1] = past[i][1]+(ran*L*((rand()/(double)RAND_MAX)-0.5));
		state[i][2] = past[i][2]+(ran*L*((rand()/(double)RAND_MAX)-0.5));
		state[i][3] = 0;
		state[i][4] = 0;
		state[i][5] = 0;
		//enforce pbc on state vectors:		
		for(int j = 0; j<3; j++){
			if(state[i][j]<0 || state[i][j]>=L){
				state[i][j] = state[i][j] - L*floor(state[i][j]/L); //enforce pbc
			}
		}
		
	}
	#endif
	#ifdef VV
	for(int i = 0; i<N; i++){	
		double ran = 0.2;//for now keep velocity 0
		past[i][0] = 0;
		past[i][1] = 0;
		past[i][2] = 0;
		state[i][0] = (i%N_core)*L/N_core+(ran*L*((rand()/(double)RAND_MAX)-0.5));
		state[i][1] = ((i/N_core)%N_core)*L/N_core+(ran*L*((rand()/(double)RAND_MAX)-0.5));
		state[i][2] = ((i/(N_core*N_core))%N_core)*L/N_core+(ran*L*((rand()/(double)RAND_MAX)-0.5));
		state[i][3] = 0;
		state[i][4] = 0;
		state[i][5] = 0;
		//enforce pbc on state vectors:		
		for(int j = 0; j<3; j++){
			if(state[i][j]<0 || state[i][j]>=L){
				state[i][j] = state[i][j] - L*floor(state[i][j]/L); //enforce pbc
			}
		}
		//printf("%f %f %f \n", state[i][0], state[i][1], state[i][2]);
	}
	#endif
	#endif	
	printf("%lf %lf %lf %lf \n",0.00,state[testPart][0],state[testPart][1],state[testPart][2]);
	#ifndef VV	
	verletengine(state, past, N, L);
	#endif
	#ifdef VV
	vv(state, past, N,L);
	#endif
	for(int i=0; i<N; i++){
		free(state[i]);
		free(past[i]);
	}
	free(state);
	free(past);
	return 0;
	}
