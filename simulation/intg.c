/*********************************/
/* Numerical integration library */
/*********************************/

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
#include "intg.h"




/* === FUNCTION  Heun ====================
 * Description:  dif. eq. integration
 *   Variables:  dt,x
 * ======================================= */

double Heun(double x, t_qif prmts, double dt, double(* F)( double , t_qif)) {
  double xint = 0, x2 = 0;

  xint = x + dt*(*F)(x,prmts);
  sprintf(mesg,"xint = %lf",xint );
  DEBUG3(mesg);
  x2   = x + 0.2*dt*((*F)(xint,prmts) + (*F)(x,prmts));
  sprintf(mesg,"x2 = %lf",x2 );
  DEBUG3(mesg);
  return x2;
}

/* === FUNCTION  rk4 ====================
 * Description:  Runge-Kutta, order 4
 *   Variables:  dt, x, F
 * ======================================= */

double rk4(double x, t_qif prmts, double dt, double (*F)(double ,t_qif)) {
  double k1, k2, k3, k4;

  k1 = dt*(*F)(x,prmts);
  k2 = dt*(*F)(x+k1/2, prmts);
  k3 = dt*(*F)(x+k2/2,prmts);
  k4 = dt*(*F)(x+k3,prmts);
  return (x + k1/6 + k2/3 + k3/3 + k4/6);
}


/* === FUNCTION  rk42 ====================
 * Description:  Runge-Kutta, order 4
 *   Variables:  d, x, F
 * ======================================= */

double rk42(double x, t_qif prmts, t_data d, double (*F)(double ,t_qif,t_data)) {
  double k1, k2, k3, k4;

  k1 = d.dt*(*F)(x,prmts,d);
  k2 = d.dt*(*F)(x+k1/2, prmts,d);
  k3 = d.dt*(*F)(x+k2/2,prmts,d);
  k4 = d.dt*(*F)(x+k3,prmts,d);
  return (x + k1/6 + k2/3 + k3/3 + k4/6);
}


/* /\* === FUNCTION  rk4_void ==================== */
/*  * Description:  Runge-Kutta, order 4 */
/*  *   Variables:  d, x, F */
/*  * ======================================= *\/ */

/* double rk42(double x, void struct prmts, void struct d, double (*F)(double ,void struct,void struct)) { */
/*   double k1, k2, k3, k4; */

/*   k1 = d.dt*(*F)(x,prmts,d); */
/*   k2 = d.dt*(*F)(x+k1/2, prmts,d); */
/*   k3 = d.dt*(*F)(x+k2/2,prmts,d); */
/*   k4 = d.dt*(*F)(x+k3,prmts,d); */
/*   return (x + k1/6 + k2/3 + k3/3 + k4/6); */
/* } */
