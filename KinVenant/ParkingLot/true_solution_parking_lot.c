#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

// Limit of Iterative Times
#define MAXTIMES 1000

double func(double q, double* t)
{
	return(*t - 1800 - (600 - 21600*q)/(2.82*pow(q,2.0/5)));
	//return (2.802*pow(q,2.0/5)*(*t-30)+600-21600*q);
	//return(182.88-70921098*q-2.209636*pow(q,2.0/5)*(*t - 1800));
}

double dfunc(double q, double* t)
{
	return(600*2/(2.82*5)*pow(q,-7.0/5) - 21600*3/(5*2.82)*pow(q,-2.0/5));
	//return(-2.0/5*2.802/(pow(q,3.0/5))*(*t-30) - 21600);
	//return(-70921.98-2.2096361*2/5*(*t-1800)*pow(q,-3.0/5));
}

void fdfunc(double x, double* t, double *ret_f, double* ret_df)
{
	*ret_f = func(x,t);
	*ret_df = dfunc(x,t);
	return;
}

double root_finder(double t)
{
	int i, times, status;
	double tmp_t = t;
	gsl_function_fdf fdf;
	gsl_root_fdfsolver *workspace;
	double q;

	/* FDF Solver */
	//workspace = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);
	workspace = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_secant);
	fdf.f = &func;
	fdf.df = &dfunc;
	fdf.fdf = &fdfunc;
	fdf.params = &tmp_t;

	/* set initial value */
	q = 0.001;

	// set solver
	gsl_root_fdfsolver_set(workspace, &fdf, q);
	for (times = 0; times < MAXTIMES; times++)
	{
		status =gsl_root_fdfsolver_iterate(workspace);
		q = gsl_root_fdfsolver_root(workspace);
		status = gsl_root_test_residual(func(q, &tmp_t), 1.0e-5);

		if((status != GSL_CONTINUE) || (status == GSL_EZERODIV))
		{
			//printf("Status : %s\n", gsl_strerror(status));
			//printf("\n Root = %25.17e\n\n", q);
			break;
		}
	}

	// free
	gsl_root_fdfsolver_free(workspace);
	return q;
}

int main(int argc, char** argv)
{
	double t0 = 0;
	double tN = 3600;
	double dt = 30;

	int N = (tN)/dt;

	double q[N+1];
	double t[N+1];

	for (int i = 0; i <= N; i++)
	{
		if ( i > 0)
			t[i] = t[i-1] + dt;
		else
			t[i] = 0;

		if (t[i] < 1494.0)
		{
			double nval = 0.025*pow(3.28084,-1.0/3);
			double istar = 2.0/12/3600;
			q[i] = 1.0/nval*0.04*pow(istar*t[i],5.0/3);
			q[i] = q[i]*60;

		}
		else if (t[i] <= 1800)
			q[i] = 1.67;
		else
		{
			//double fval = func(0.002586,&t[i]);
			//printf("fval = %lf \n" ,fval);
			q[i] = root_finder(t[i]);
			q[i] = q[i]*60;

		}
		
		printf("%lf \t %lf \n", t[i]/60, q[i]);
		//printf("%lf \t %lf \n", t[i]/60, q[i]);
	}

	double integralTime;
	printf("Enter the time (in seconds) till when you would like to calculate the integral. It has to be greater than %lf:\n", dt);
	scanf("%lf", &integralTime);
	integralTime = integralTime/60;
	double integral = 0;
	double time = 0;
	int i = 1;
	dt = dt/60;

	FILE *IntegralFile = fopen("TrueIntegral.txt","w");
	while (time < integralTime)
	{
		if (time+dt > integralTime)
		{
			printf("time before = %lf\n", time);
			dt = integralTime - time;
		}
		time = time+dt;
		printf("time after = %lf\n", time);

		integral+= 0.5*(q[i] + q[i-1])*dt;
		i++;
		printf("time = %lf\n", time);
		fprintf(IntegralFile, "%lf \t %lf\n", time, integral);
	}
	printf("the integral is %lf \n", integral);

	fclose(IntegralFile);
}
