#ifndef SPKMEANS_H_
#define SPKMEANS_H_

double** fit_wam(double** mat, int n, int d);
double** fit_ddg(double** mat, int n, int d);
double** fit_gl(double** mat, int n, int d);
double*** fit_jacobi(double** mat, int n);
double** fit_spk(int K, int maxiter, double EPS, int vectorLength,int N, double** clustersArr ,double** dataPointsArr); 

#endif
