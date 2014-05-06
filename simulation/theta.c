/***********************************/
/* Theta neuron simulation library */
/***********************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<ctype.h>
#include<limits.h>

/* Additional libraries (custom libraries) */
#include "utils.h"
#include "common.h"
#include "theta.h"
#include "intg.h"

/* === FUNCTION  th ====================
 * Description:  Theta neuron algorithm
 *   Variables:  Theta_i and parameters
 * ======================================= */

t_th THETA(t_data d,t_th th) {
  double (*fptr) (double, t_th);
  fptr = theta;
  if((th.th = rk4(th.th,th,d.dt,fptr)) >= PI) {
    th.spike2 = 1;
    th.th = -PI;
  }
  return th;
}



/* === FUNCTION  theta ====================
 * Description:  Theta neuron equation
 *   Variables:  Theta variable and parameters
 * ======================================= */

double theta(double x, t_th t) {
  double a, b, theta_punto;
  a = t.eta + t.J*t.r2 + t.g*t.V0*t.r2;
  b = t.g*t.r2;
  theta_punto = 1 - cos(x) + a*(1 + cos(x)) - b*sin(x);
  return theta_punto;
}


