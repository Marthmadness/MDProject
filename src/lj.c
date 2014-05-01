#include "md.h"
#include "lj.h"
void lj_full(double **state, int N, double L, double epsilon, double sigma, double *p){
	double pressure=0;
	double test = 0;
	double test2 = 0;
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
			test = rij2;
			double mag = 0;
			mag = force(rij2, epsilon, sigma);
			test2 = mag;
			//now update each component of the force
			for(int k=3; k<6; k++){
				state[i][k] += mag*rij[k-3];
				state[j][k] += -mag*rij[k-3];//newton's 3rd law
				pressure += mag*rij[k-3]*rij[k-3];
			}
		}
	}
	*p = pressure/(3*L*L*L);
	int x = 0;
	printf("%lf %lf %lf %lf %lf \n", state[x][3], state[x][4], state[x][5], test, test2);
}
double force(double rij2, double epsilon, double sigma){
double r6 = rij2*rij2*rij2;
double r12 = r6*r6;
double sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
double sigma12 = sigma6*sigma6;
return (-24*epsilon/(rij2) * (2*(sigma12/r12) - (sigma6/r6)));
}
