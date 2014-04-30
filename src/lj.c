#include "md.h"
#include "lj.h"
void lj_full(double **state, int N, double L, double epsilon, double sigma, double *p){
	double pressure=0;
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
			double mag = 0;
			mag = force(rij2, epsilon, sigma);
			//now update each component of the force
			for(int k=3; k<6; k++){
				state[i][k] += mag*rij[k-3];
				state[j][k] += -mag*rij[k-3];//newton's 3rd law
				pressure += mag*rij[k-3]*rij[k-3];
			}
		}
	}
	*p = pressure/(3*L*L*L);
	
	//printf("%lf %lf %lf \n", state[1][3], state[1][4], state[1][5]);
}
double force(double rij2, double epsilon, double sigma){
return (-24*epsilon/(rij2) * (2*(pow(sigma,12)/(pow(rij2,6))) - (pow(sigma,6)/pow(rij2,3))));
}
