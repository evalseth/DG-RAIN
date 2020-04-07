#ifndef NUMERICAL_FLUX

#define NUMERICAL_FLUX

extern void RoeFluxChan(double *Fhat, double A_L, double A_R, double Q_L, double Q_R, double b, double m1val, double m2val, double g);
extern void LFChan(double *Fhat, double A_L, double A_R, double Q_L, double Q_R, double b, double m1val, double m2val, double g);
extern double upwindKinRoute(double H_L, double H_R, double S0_L, double S0_R, double n);

extern double RoeFluxJunc(double zeta_in, double zeta_ex, double Qx_in, double Qx_ex,
	 double Qy_in, double Qy_ex, double z_edge, double nx, double ny, double localG, double *Fhatdotn);
extern void ComputeJuncInnerProducts(struct TwoDRegion *junc, double *Fhat1dotn, double *Fhat2dotn, 
							double *Fhat3dotn, int el, double *SI, double *VI);

#endif
