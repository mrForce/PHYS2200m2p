
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_const_mksa.h>

struct doublePair{
	double first;
 	double second;
};
int func (double r, const double y[], double f[], void *params);

struct doublePair calculate_radius_mass (double p_central);


int main(void){

  double p;
  struct doublePair data;

  for(p = 0.1; p <= 1000000.0; p *= 1.1){

        data = calculate_radius_mass(p);
	printf ("% .5e  % .5e  % .5e\n", p, data.first, data.second);
 }
  return 0;

}

struct doublePair calculate_radius_mass (double p_central)
{
  size_t neqs = 2;          /* number of equations */
  double eps_abs = 1.e-8, 
    eps_rel = 0.;           /* desired precision */
  double stepsize = 1e-6;   /* initial integration step */
  double t = 0., t1 = 70000000.; /* time interval */
  int status;
  /* 
   * Initial conditions 
   */
  double y[2] = {p_central, 0.0 };   /* for res1 */


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
  double last_radius = 0;
  double last_mass = 0;
  while ((y[0] > 0.0) && t < t1 )
  {
    status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t,
                                      t1, &stepsize, y);

    if (status != GSL_SUCCESS) {
      printf ("Troubles: % .5e  % .5e  % .5e\n", 
              t, y[0], y[1]);
      break;
    }
    //printf ("% .5e  % .5e  % .5e\n", 
           //   t, y[0], y[1]);
    if(y[0] > 0){
	last_radius = t;
    	last_mass = y[1];
    }
    
  }
  struct doublePair data;
  data.first = last_radius;
  data.second = last_mass;

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);

  return data;
}

int func (double r, const double y[], double f[], 
	  void *params) 
{
  /*
   * y[0] - density, p
   * y[1] - mass, m
   */
 double gamma = pow(y[0], 2./3.)/(3*sqrt(1 + pow(y[0], 2./3.)));
 double p_central = *(double *) params;
 if(r < 1.e-8){
	f[0] = -1*y[0]*(4./3.)*M_PI*r*p_central;
 }else{
 	f[0] = -1*y[1]*y[0]/(gamma*pow(r, 2));
 }
f[1] = pow(r, 2)*y[0];
  return GSL_SUCCESS;
}

