#include "verlet.h"
#include "md.h"

void verletengine(double **state, double **past, int N, double L){
// set conditions
	
	double t = 0; //initial time is 0
	double deltat = 0.001;// timestep
	int steps = 10000; //number of steps to simulate
	double mass = 39.948;
	srand(time(NULL));
	int testPart = 0;
	double kb = 0.00831;
	double temp = 298;
	double epsilon = 120*kb*temp;
	double sigma = 0.34;

	//begin the loop of the simulation
	for(int step=0; step<steps; step++){
		//zero out forces
		for(int i = 0; i<N; i++){	
			for(int j = 3; j<6; j++){
				state[i][j] = 0;
			}
		}
		//start calculating forces
		double p = 0;
		lj(state, N, L, epsilon, sigma, &p);
		//now integrate with respect to time! Verlet style!
		p = p + (N/(L*L*L))*kb*temp;
		double temp = 0.0;
		for(int i=0; i<N; i++){
				for(int j=0; j<3; j++){
					double next = 2*state[i][j] - past[i][j] + deltat*deltat*state[i][j+3]/mass;
					if(next<0 || next>=L){
						next = next - L*floor(next/L); //enforce pbc
					}

					temp += 0.5*mass*((next-past[i][j])/(2*deltat))*((next-past[i][j])/(2*deltat));//calculate temperature
					temp/=N;
					past[i][j] = state[i][j];
					state[i][j] = next;
				}
			}
			
		t += deltat;
		if(step%100 ==99){
		printf("%f %f %f %f %f %f\n",t,state[testPart][0],state[testPart][1],state[testPart][2],temp,p);
		}
		
	}
}