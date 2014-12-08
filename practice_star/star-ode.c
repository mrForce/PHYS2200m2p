/*
 * The function provides the right hand sides of the two differential equations in the system.
 * dp/dr = -3*(pi^2)*m(r)*p(r)/r^2
 * dm/dr = 4*pi*p(r)	
 * 
 */

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>


int func (double r, const double y[], double f[], 
	  void *params) 
{
  /*
   * y[0] - density, p
   * y[1] - mass, m
   */

 double gamma = pow(y[0], 2./3.)/(3*sqrt(1 + pow(y[0], 2./3.)));
 f[1] = pow(r, 2)*y[0];
 double p_central = *(double *) params;
 if(r < 1.e-8){
	f[0] = -1*y[0]*(4/3)*M_PI*r*p_central;
 }else{
 	f[0] = -1*y[1]*y[0]/(gamma*pow(r, 2));
 }
  return GSL_SUCCESS;
}

