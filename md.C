#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
//prototyping
void verletengine(double **state, double **past, int N, double L);

//initializing stuff
int main(void){
	
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
void verletengine(double **state, double **past, int N, double L){
// set conditions

	double t = 0; //initial time is 0
	double deltat = 0.001;// timestep
	int steps = 10000; //number of steps to simulate
	double mass = 1.0;
	srand(time(NULL));
	int testPart = 0;

	//begin the loop of the simulation
	for(int step=0; step<steps; step++){
		//zero out forces
		for(int i = 0; i<N; i++){	
			for(int j = 3; j<6; j++){
				state[i][j] = 0;
				}
			}
		//start calculating forces
		double epsilon = 1.0;
		double sigma = 1.0;
		for(int i = 0; i < N-1; i++){//loop over all pairs
			for(int j=i+1; j<N; j++){
					double rij[3];
					double rij2 = 0;
					for(int k=0; k<3; k++){
						rij[k] = state[j][k] - state[i][k]; 
						if(abs(rij[k]) > 0.5*L){
							rij[k] = (L - abs(rij[k]))*((rij[k] < 0) - (rij[k] > 0));//minimum image convention, requires that the particles are inside the cube
						}
						rij2 += rij[k]*rij[k];
					}
					double force = -24/(rij2) * (2*(pow(sigma,12)/(pow(rij2,6))) - (pow(sigma,6)/pow(rij2,3)));
					//now update each component of the force
					for(int k=3; k<6; k++){
						state[i][k] += force*rij[k-3];
						state[j][k] += -state[i][k];//newton's 3rd law
					}
				}
			}
		//now integrate with respect to time! Verlet style!
		double vTotal = 0.0;
		for(int i=0; i<N; i++){
				for(int j=0; j<3; j++){
					double next = 2*state[i][j] - past[i][j] + deltat*deltat*state[i][j+3]/mass;
					if(next<0 || next>L){
						next = next - L*floor(next/L); //enforce pbc
					}
					if(j==1){
						vTotal += next-state[i][j];
					}
					past[i][j] = state[i][j];
					state[i][j] = next;
				}
			}
		t += deltat;
		if(step%100 ==99)
		printf("%f %f %f %f %f \n",t,state[testPart][0],state[testPart][1],state[testPart][2],vTotal);
		
		
		}
}