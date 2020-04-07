struct infil_func_params
{
	double K, psi, delta_theta, time; 
};

double infil_func(double x, void*params);

double infil_deriv(double x, void* params);

void infil_fdf(double x, void* params, double *y, double *dy);

double calculateInfiltration(double psi, double s_e, double theta_e, double K, double time);
