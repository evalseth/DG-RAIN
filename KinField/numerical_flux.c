#include <stdio.h>
#include <math.h>
#include "Math_Functions.h"
#include "Globals.h"
#include "Constitutive_Equations.h"
#include "Watershed.h"

void RoeFluxChan(double *Fhat, double A_L, double A_R, double Q_L, double Q_R, double b, 
	double m1, double m2, double localG)
{
	double I1_L, I1_R;
	double F1_L, F1_R;
	double F2_L, F2_R;
	double c_L, c_R;
	double u_L, u_R;
	double lambda1_L, lambda1_R, lambda2_L, lambda2_R;
	double uhat, chat;
	double lambdahat1, lambdahat2;
	double epsilon1, epsilon2;

	double h_L = getH(A_L,b,m1,m2);
	double h_R = getH(A_R, b, m1, m2);
	
	// width at the free surface
	double B_L = b+m1*h_L+m2*h_L;
	double B_R = b+m1*h_R+m2*h_R;

	I1_L = getI1(A_L, b, m1, m2);
	I1_R = getI1(A_R, b, m1, m2);

	F1_L = Q_L;
	F1_R = Q_R;

	F2_L = pow(Q_L,2)/A_L + localG*I1_L;
	F2_R = pow(Q_R,2)/A_R + localG*I1_R;
		
	c_L = sqrt(localG*A_L/B_L);
	c_R = sqrt(localG*A_R/B_R);
		
	u_L = Q_L/A_L;
	u_R = Q_R/A_R;
	
	lambda1_L = u_L + c_L;
	lambda1_R = u_R + c_R;
	lambda2_L = u_L - c_L;
	lambda2_R = u_R - c_R;
		
	uhat = (Q_L/sqrt(A_L) + Q_R/sqrt(A_R))/(sqrt(A_L) + sqrt(A_R));
	chat = sqrt(localG/2*(A_L/B_L + A_R/B_R));
		
	lambdahat1 = uhat + chat;
	lambdahat2 = uhat - chat;

	epsilon1 = max(0, (lambdahat1 - lambda1_L));
	epsilon1 = max(epsilon1, (lambda1_R - lambdahat1));
	lambdahat1 = max(epsilon1,fabs(lambdahat1));
		
	epsilon2 = max(0, (lambdahat2 - lambda2_L));
	epsilon2 = max(epsilon2, (lambda2_R - lambdahat2));
	lambdahat2 = max(epsilon2,fabs(lambdahat2));

	// Eigenvectors
	//double R1[2] = {1, uhat + chat};
	//double R2[2] = {1, uhat - chat};
	double Rhat[2][2] = {{1,1},{uhat+chat,uhat-chat}};
	double dtm = 1/(Rhat[0][0]*Rhat[1][1]-Rhat[0][1]*Rhat[1][0]);
	double Rhatinv[2][2] = {{dtm*Rhat[1][1], -dtm*Rhat[0][1]},{-dtm*Rhat[1][0], dtm*Rhat[0][0]}};

	double RoeJac[2][2];
	RoeJac[0][0] = lambdahat1*Rhat[0][0]*Rhatinv[0][0]+lambdahat2*Rhat[0][1]*Rhatinv[1][0]; 
	RoeJac[0][1] = lambdahat1*Rhat[0][0]*Rhatinv[0][1]+lambdahat2*Rhat[0][1]*Rhatinv[1][1]; 
	RoeJac[1][0] = lambdahat1*Rhat[1][0]*Rhatinv[0][0]+lambdahat2*Rhat[1][1]*Rhatinv[1][0]; 
	RoeJac[1][1] = lambdahat1*Rhat[1][0]*Rhatinv[0][1]+lambdahat2*Rhat[1][1]*Rhatinv[1][1]; 

	double jump[2] = {A_L - A_R, Q_L - Q_R};

	Fhat[0] = 0.5*(F1_L + F1_R) + 0.5*(RoeJac[0][0]*jump[0]+RoeJac[0][1]*jump[1]);
	Fhat[1] = 0.5*(F2_L + F2_R) + 0.5*(RoeJac[1][0]*jump[0]+RoeJac[1][1]*jump[1]);

}


void LFChan(double *Fhat, double A_L, double A_R, double Q_L, double Q_R, double b, double m1, double m2, double localG)
{
	double I1_L, I1_R;
	double F1_L, F1_R;
	double F2_L, F2_R;
	double c_L, c_R;
	double u_L, u_R;

	double h_L = getH(A_L,b, m1, m2);
	double h_R = getH(A_R,b, m1, m2);
	I1_L = getI1(A_L, b, m1, m2);
	I1_R = getI1(A_R, b, m1, m2);

	F1_L = Q_L;
	F1_R = Q_R;

	F2_L = pow(Q_L,2)/A_L + I1_L;
	F2_R = pow(Q_R,2)/A_R + I1_R;
	//F2_L = pow(Q_L,2)/A_L + localG*I1_L;
	//F2_R = pow(Q_R,2)/A_R + localG*I1_R;

	// free surface width
	double B_L = b+m1*h_L+m2*h_L;
	double B_R = b+m1*h_R + m2*h_R;

	c_L = sqrt(localG*A_L/B_L);
	c_R = sqrt(localG*A_R/B_R);
		
	u_L = Q_L/A_L;
	u_R = Q_R/A_R;

	double C = max(fabs(u_L + c_L), fabs(u_R+c_R));
	C = max(C, fabs(u_L - c_L));
	C = max(C, fabs(u_R - c_R));
	double jump[2] = {A_L - A_R, Q_L - Q_R};
	Fhat[0] = 0.5*(F1_L+F1_R) + 0.5*C*jump[0];
	Fhat[1] = 0.5*(F2_L+F2_R) + 0.5*C*jump[1];
	

}

/**************************************************************************************//**
*
* Function to calculate Roe's flux for the 2-D shallow water equations at an interface 
* (element edges).
* @param[in] zeta_in value of the water surface height at the edge coming from the interior
* element
* @param[in] zeta_ex value of the water surface height at the edge coming from the exterior
* element
* @param[in] Qx_in value of the momentum in the x-direction coming from the interior element
* @param[in] Qx_ex value of the momentum in the x-direction coming from the exterior element
* @param[in] Qy_in value of the momentum in the y-direction coming from the interior element
* @param[in] Qy_ex value of the momentum in the y-direction coming from the exterior element
* @param[in] z_edge value of the bathymetry at the edge. We assume that the bathymetry 
* information is provided at the edge and thus will be continuous at an edge. In the future
* we might need to revise this assumption and figure out how to handle it if the bathymetry
* information is provided at the nodes instead.
* @param[in] nx x-component of the outward unit normal vector to the interior element at the edge
* @param[in] ny y-component of the outward unit normal vector to the interior element at the edge
* @param[in] localG value of the gravitational acceleration constant. This will be the same as g,
* except when we set it to zero to handle wetting and drying.
* @param[out] Fhatdotn a pointer to an array of size 3 where the numerical flux computed
*  will be stored
*  @return current_max_lambda the maximum eigenvalue at the interface
* *****************************************************************************************/
double RoeFluxJunc(double zeta_in, double zeta_ex, double Qx_in, double Qx_ex,
	 double Qy_in, double Qy_ex, double z_edge, double nx, double ny, double localG, double *Fhatdotn)
{
	double H_in = zeta_in + z_edge;
	double H_ex = zeta_ex + z_edge;

	double u_in = Qx_in/H_in;
	double u_ex = Qx_ex/H_ex;

	double v_in = Qy_in/H_in;
	double v_ex = Qy_ex/H_ex;

	double jump[3] = {zeta_in - zeta_ex, Qx_in - Qx_ex, Qy_in - Qy_ex};

	double F1x_in = H_in*u_in;
	double F1x_ex = H_ex*u_ex;

	double F1y_in = H_in*v_in;
	double F1y_ex = H_ex*v_ex;
	 
	double F2x_in = H_in*u_in*u_in + 0.5*localG*(H_in*H_in-z_edge*z_edge);
	double F2x_ex = H_ex*u_ex*u_ex + 0.5*localG*(H_ex*H_ex-z_edge*z_edge);


	//printf("F2x_in = %lf localG = %lf z_edge = %lf\n", F2x_in, localG, z_edge);
	double F2y_in = H_in*u_in*v_in;
	double F2y_ex = H_ex*u_ex*v_ex;

	double F3x_in = H_in*u_in*v_in;
	double F3x_ex = H_ex*u_ex*v_ex;

	double F3y_in = H_in*v_in*v_in + 0.5*localG*(H_in*H_in-z_edge*z_edge);
 	double F3y_ex = H_ex*v_ex*v_ex + 0.5*localG*(H_ex*H_ex-z_edge*z_edge);

	double Fn_in[3] = {F1x_in*nx+F1y_in*ny, F2x_in*nx+F2y_in*ny, F3x_in*nx+F3y_in*ny};
		
	double Fn_ex[3] = {F1x_ex*nx+F1y_ex*ny, F2x_ex*nx+F2y_ex*ny, F3x_ex*nx+F3y_ex*ny};

		
	// Roe's averages
	double Hhat = 0.5*(H_in+H_ex);
	double chat = sqrt(g*Hhat);
	double denom =1.0/( sqrt(H_in) + sqrt(H_ex));
	double uhat = (u_in*sqrt(H_in)+u_ex*sqrt(H_ex))*denom;
	double vhat = (v_in*sqrt(H_in)+v_ex*sqrt(H_ex))*denom;

	//printf("H_in = %lf , H_ex = %lf, u_in = %lf, u_ex = %lf, v_in = %lf, v_ex = %lf\n", H_in, H_ex, u_in, u_ex, v_in, v_ex);
	//printf("Hhat = %lf, uhat = %lf, vhat = %lf\n", Hhat, uhat, vhat);

	double Rhat[3][3] = {{1, 0, 1},{uhat-chat*nx, -ny, uhat+chat*nx},{vhat-chat*ny, nx,vhat+chat*ny}};
	double Rhatinv[3][3];
	denom = -2*chat;
	Rhatinv[0][0] = (-chat-nx*uhat-ny*vhat)/denom;
	Rhatinv[0][1] = nx/denom;
	Rhatinv[0][2] = ny/denom;
	Rhatinv[1][0] = (-2*chat*ny*uhat+2*chat*nx*vhat)/denom;
	Rhatinv[1][1] = 2*chat*ny/denom;
	Rhatinv[1][2] = -2*chat*nx/denom;
	Rhatinv[2][0] = (-chat+nx*uhat+ny*vhat)/denom;
	Rhatinv[2][1] = -nx/denom;
	Rhatinv[2][2] = -ny/denom;

	double c_in = sqrt(g*H_in);
	double c_ex = sqrt(g*H_ex);
	double lambda1_in = u_in*nx+v_in*ny-c_in;
	double lambda2_in = u_in*nx+v_in*ny;
	double lambda3_in = u_in*nx+v_in*ny+c_in;

	double lambda1_ex = u_ex*nx+v_ex*ny-c_ex;
	double lambda2_ex = u_ex*nx+v_ex*ny;
	double lambda3_ex = u_ex*nx+v_ex*ny+c_ex;

	double lambdahat1 = uhat*nx+vhat*ny-chat;
	double lambdahat2 = uhat*nx+vhat*ny;
	double lambdahat3 = uhat*nx+vhat*ny+chat;

	double epsilon1 = fmax(0,(lambdahat1-lambda1_in));
	epsilon1 = fmax(epsilon1, (lambda1_ex-lambdahat1));
	lambdahat1 = fmax(epsilon1, fabs(lambdahat1));

	double epsilon2 = fmax(0,(lambdahat2-lambda2_in));
	epsilon2 = fmax(epsilon2, (lambda2_ex-lambdahat2));
	lambdahat2 = fmax(epsilon2, fabs(lambdahat2));

	double epsilon3 = fmax(0,(lambdahat3-lambda3_in));
	epsilon3 = fmax(epsilon3, (lambda3_ex-lambdahat3));
	lambdahat3 = fmax(epsilon3, fabs(lambdahat3));

	double AbsLambdaHat[3][3] = {{lambdahat1,0,0},{0,lambdahat2,0},{0,0,lambdahat3}};
	
	double prod1[3];
	double prod2[3];
	double R_Lam_Rinv_jump[3]; 

	MatrixVectorMultiply(&Rhatinv[0][0],jump,prod1,3,3);
	MatrixVectorMultiply(&AbsLambdaHat[0][0],prod1,prod2,3,3);
	MatrixVectorMultiply(&Rhat[0][0],prod2,R_Lam_Rinv_jump,3,3);

	for (int k=0; k < 3; ++k)
	{
		Fhatdotn[k] = 0.5*(Fn_in[k]+Fn_ex[k]) + 0.5*R_Lam_Rinv_jump[k];
		//printf("Fn_in = %lf Fn_ex = %lf\n", Fn_in[k], Fn_ex[k]);
		//printf("Fhat = %lf\n", Fhatdotn[k]);
	} 		

	double current_max_lam = max(fabs(lambda1_in), fabs(lambda1_ex));
	current_max_lam = max(current_max_lam, fabs(lambda2_in));
	current_max_lam = max(current_max_lam, fabs(lambda2_ex));
	current_max_lam = max(current_max_lam, fabs(lambda3_in));
	current_max_lam = max(current_max_lam, fabs(lambda3_ex));

	return current_max_lam;
}
