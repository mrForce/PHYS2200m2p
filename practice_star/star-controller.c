/*
	This program solves for the internal structure of a white dwarf! It uses the following differential equations to model the density:

dm/dr = p(r)*r^2
dp/dr = -m(r)*p(r)/(g(p)*r^)

Where g(p) = (p^(2/3))/(3*sqrt(1 + p^(2/3)))
*/
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_const_mksa.h>

int func (double t, const double y[], double f[], void *params);

int main (void) 
{
  size_t neqs = 2;          /* number of equations */
  double eps_abs = 1.e-8, 
    eps_rel = 0.;           /* desired precision */
  double stepsize = 1e-6;   /* initial integration step */
  double t = 0.0, t1 = 2.; /* time interval */


  int status;
  /* 
   * Initial conditions 
   */
  double y[2] = {10.0, 0.0 };   /* for res1 */
  double p_central = y[0];

  /*
   * Explicit embedded Runge-Kutta-Fehlberg (4,5) method. 
   * This method is a good general-purpose integrator.
   */
  gsl_odeiv2_step    *s = gsl_odeiv2_step_alloc 
                            (gsl_odeiv2_step_rkf45, neqs);
  gsl_odeiv2_control *c = gsl_odeiv2_control_y_new (eps_abs, 
						    eps_rel);
  gsl_odeiv2_evolve  *e = gsl_odeiv2_evolve_alloc (neqs);
    
  gsl_odeiv2_system sys = {func, NULL, neqs, &p_central};
    
  /*
   * Evolution loop 
   */
  while ((y[0] >= 0.0) )
  {
    status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t,
                                      t1, &stepsize, y);

    if (status != GSL_SUCCESS) {
      printf ("Troubles: % .5e  % .5e  % .5e\n", 
              t, y[0], y[1]);
      break;
    }
    printf ("% .5e  % .5e  % .5e \n", 	      
            t, y[0], y[1]);
  }

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);

  return 0;
}
