#include "md.h"
#include "lj.h"
#include "cell_list.h"
void lj_full(double **state, int N, double L, double epsilon, double sigma, double *p){
	double pressure=0;

	for(int i = 0; i < N-1; i++){//loop over all pairs
		for(int j=i+1; j<N; j++){
			double rij[3];
			double rij2 = 0;
			for(int k=0; k<3; k++){
				rij[k] = state[j][k] - state[i][k];
				//printf("%lf AHHHHH  %lf\n", rij[k], fabs(rij[k]));
				if(fabs(rij[k]) > (0.5*L)){
					rij[k] = (L - fabs(rij[k]))*((rij[k] < 0) - (rij[k] > 0));//minimum image convention, requires that the particles are inside the cube
				}
				rij2 += rij[k]*rij[k];
			}
			double mag = 0;
			mag = force_cutoff(rij2, epsilon, sigma);
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
	//printf("%lf %lf %lf %lf %lf \n", state[x][3], state[x][4], state[x][5], test, test2);
}
void lj_cell(double **state, int** list, int N, int nx, double L, double epsilon, double sigma, double *p){
	double pressure=0;
	for(int i = 0; i < nx*nx*nx; i++){//loop each cell
		for(int x=0; x<=1; x++){//then loop over all possible neighboors
			for(int y= 0; y<=1; y++){
				for(int z=0; z<=1; z++){
					int j = find_neighbor(nx, i, x, y, z);
					for(int n1=0; n1<N; n1++){
						if(list[i][n1]==-1){
							break;
						}
						for(int n2=0; n2<N; n2++){
							if(list[j][n2]==-1){
								break;
							}
							double rij[3];
							double rij2 = 0;
							for(int k=0; k<3; k++){
								rij[k] = state[list[i][n1]][k] - state[list[j][n2]][k];
								//printf("%lf AHHHHH  %lf\n", rij[k], fabs(rij[k]));
								if(fabs(rij[k]) > (0.5*L)){
									rij[k] = (L - fabs(rij[k]))*((rij[k] < 0) - (rij[k] > 0));//minimum image convention, requires that the particles are inside the cube
								}
								rij2 += rij[k]*rij[k];
							}
							double mag = 0;
							mag = force_cutoff(rij2, epsilon, sigma);
							//now update each component of the force
							for(int k=3; k<6; k++){
								state[list[i][n1]][k] += mag*rij[k-3];
								state[list[j][n2]][k] += -mag*rij[k-3];//newton's 3rd law
								pressure += mag*rij[k-3]*rij[k-3];
							}
						}
					}
				}
			}
		}
	}
	*p = pressure/(3*L*L*L);
	//printf("%lf %lf %lf %lf %lf \n", state[x][3], state[x][4], state[x][5], test, test2);
}

void lj_cell_bit(double **state, int** list, int N, int nx, double L, double epsilon, double sigma, double *p){
	double pressure=0;
	for(int i = 0; i < nx*nx*nx; i++){//loop each cell
		for(char t=0; t<=8; t++){//then loop over all possible neighboors
			int x = t&0x01;
			int y = t&0x02;
			int z = t&0x04;
			int j = find_neighbor(nx, i, x, y, z);
			for(int n1=0; n1<N; n1++){
				if(list[i][n1]==-1){
					break;
				}
				for(int n2=0; n2<N; n2++){
					if(list[j][n2]==-1){
						break;
					}
					double rij[3];
					double rij2 = 0;
					for(int k=0; k<3; k++){
						rij[k] = state[list[i][n1]][k] - state[list[j][n2]][k];
						//printf("%lf AHHHHH  %lf\n", rij[k], fabs(rij[k]));
						if(fabs(rij[k]) > (0.5*L)){
							rij[k] = (L - fabs(rij[k]))*((rij[k] < 0) - (rij[k] > 0));//minimum image convention, requires that the particles are inside the cube
						}
						rij2 += rij[k]*rij[k];
					}
					double mag = 0;
					mag = force_cutoff(rij2, epsilon, sigma);
					//now update each component of the force
					for(int k=3; k<6; k++){
						state[list[i][n1]][k] += mag*rij[k-3];
						state[list[j][n2]][k] += -mag*rij[k-3];//newton's 3rd law
						pressure += mag*rij[k-3]*rij[k-3];
					}		
				}
			}
		}
	}
	*p = pressure/(3*L*L*L);
	//printf("%lf %lf %lf %lf %lf \n", state[x][3], state[x][4], state[x][5], test, test2);
}
double force(double rij2, double epsilon, double sigma){
double r6 = rij2*rij2*rij2;
double r12 = r6*r6;
double sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
double sigma12 = sigma6*sigma6;
return (-24*epsilon/(rij2) * (2*(sigma12/r12) - (sigma6/r6)));
}
double force_cutoff(double rij2, double epsilon, double sigma){
double r6 = rij2*rij2*rij2;
double r12 = r6*r6;
double sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
double sigma12 = sigma6*sigma6;
return ((-24*epsilon/(rij2) * (2*(sigma12/r12) - (sigma6/r6)))+force(2.5*sigma,epsilon,sigma));
}
