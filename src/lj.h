#ifndef LJ_H
#define LJ_H
#include "md.h"
void lj_full(double **state,int N, double L, double epsilon, double sigma, double *p);
double force(double rij2, double epsilon, double sigma);
double force_cutoff(double rij2, double epsilon, double sigma);
#endif