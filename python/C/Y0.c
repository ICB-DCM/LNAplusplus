#include "Y0.h"
void Y0(const double *Theta, double varOut[2])
{
double k_m = Theta[0];
double k_p = Theta[1];
double g_m = Theta[2];
double g_p = Theta[3];
double phi0 = Theta[4];
double phi1 = Theta[5];
varOut[0] = k_m/g_m;
varOut[1] = k_m*k_p/(g_m*g_p);
}