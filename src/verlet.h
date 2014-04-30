#ifndef VERLET_H
#define VERLET_H
void verletengine(double **state, double **past, int N, double L);
void vv(double **state, double **velocity, int N, double L);
#endif