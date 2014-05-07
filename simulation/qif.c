/***************************/
/*  QIF simulation library */
/***************************/

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
#include "qif.h"
#include "intg.h"


/* === FUNCTION  QIF ====================
 * Description:  QIF algorithm
 *   Variables:  Input data: vt,v0,N,dt,TT,th[]
 * ======================================= */

t_qif QIF(t_data d, t_qif th) {
  /* counters */

  double (* fptr) (double, t_qif);
  fptr = qif;

  /* Quadratic Integrate and Fire Algorithm */
  /* In each step, we must compute the qif eq. for each oscillator,
   * taking into account the S(t) function and I_i   */
  th.tr = 0;

  if((th.v = rk4(th.v,th,d.dt,fptr)) >= th.vp) {
    th.tr = 1;
    th.v = th.vr;
  }
  /* if((th.th = rk4(th.th,th,dt,fptr)) >= PI) */
  /*   th.th = -PI; */

  return th;
}


/*++++++++++++++++++++++++++++++++++++++
  Quadratic Integrate and Fire equation
  with an Input.

  double QIF: QIF model equation

  double I:   Constant input field (x[1])

  double J:   Mean-field input constant (J_i = J, for all i) (x[2])

  double r:  (...) Mean field... (x[3])

  double v:   Potential  (x[0])
  ++++++++++++++++++++++++++++++++++++++*/

double  qif(double x,t_qif p) {  
  double I;
  I = p.J*p.r + p.eta - p.g*p.r*(x - p.V0);
  return (x*x + I);
}
