#include "Channels_KinRoutes.h"

extern void RoeFluxChan(double *Fhat, double A_L, double A_R, double Q_L, double Q_R, double b, double m1val, double m2val, double g);
extern void chanLF(double *Fhat, double A_L, double A_R, double Q_L, double Q_R, double b, double g);
extern double upwindKinRoute(double H_L, double H_R, double S0_L, double S0_R, double n);
