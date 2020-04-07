#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <math.h>

#include "Infiltration.h"

double infil_func(double x, void* params)
{
	struct infil_func_params *p = (struct infil_func_params *) params;
	double K = p->K;
	double psi = p->psi;
	double delta_theta = p->delta_theta;
	double t = p->time;

	double psi_delta_theta = psi*delta_theta;

	return (x - K*t - psi_delta_theta*log((1+x/(psi_delta_theta))));
}

double infil_deriv(double x, void* params)
{
	struct infil_func_params *p = (struct infil_func_params *) params;
	double psi = p->psi;
	double delta_theta = p->delta_theta;
	double psi_delta_theta = psi*delta_theta;

	double deriv = 1 - 1.0/(1+ x/psi_delta_theta);
	return deriv;
}

void infil_fdf(double x, void* params, double* y, double* dy)
{
	struct infil_func_params *p = (struct infil_func_params *) params;
	double K = p->K;
	double psi = p->psi;
	double delta_theta = p->delta_theta;
	double t = p->time;

	double psi_delta_theta = psi*delta_theta;

	*y = x - K*t - psi_delta_theta*log(1+x/psi_delta_theta);
	*dy = 1 - 1.0/(1+x/psi_delta_theta);

}

double calculateInfiltration(double psi, double s_e, double theta_e, double K, double time)
{
	double deltaTheta = (1-s_e)*theta_e;

	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fdfsolver_type *T;
	gsl_root_fdfsolver *s;

	gsl_function_fdf FDF;
	struct infil_func_params params = {K, psi, deltaTheta, time};
	FDF.f = &infil_func;
	FDF.df = &infil_deriv;
	FDF.fdf = &infil_fdf;
	FDF.params = &params;

	T = gsl_root_fdfsolver_newton;
	s = gsl_root_fdfsolver_alloc(T);

	// initial guess for F = Kt
	double x = K*time;
	double x0;
	gsl_root_fdfsolver_set(s, &FDF, x);

	do
	{
		iter++;
		status = gsl_root_fdfsolver_iterate(s);
		x0 = x;
		x = gsl_root_fdfsolver_root(s);
		status = gsl_root_test_delta(x, x0, 0, 1e-3);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fdfsolver_free(s);

	if (status != GSL_SUCCESS)
		printf("not converged\n");

	return x;
}
