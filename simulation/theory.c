/********************************/
/*  Rate eq. simulation library */
/********************************/

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
#include "theory.h"
#include "intg.h"


/* === FUNCTION  Theory ====================
 * Description:  Rate Eqs. algorithm
 *   Variables:  Input data: vt,v0,N,dt,TT,th[]
 * ======================================= */

t_qif Theory(t_data d, t_qif th) {
  /* counters */

  double (* fptr) (double, t_qif, t_data);
  fptr = theory1;
  th.rh2 = th.rh;
  th.vh2 = th.vh;
  th.rh = rk42(th.rh2,th,d,fptr);

  fptr = theory2;
  th.vh = rk42(th.vh2,th,d,fptr);
  return th;
}


/*++++++++++++++++++++++++++++++++++++++

  ++++++++++++++++++++++++++++++++++++++*/

double  theory1(double x,t_qif p, t_data d) {  
  double I;
  I = d.Deta/(PI) + 2*p.rh2*p.vh2 + d.DJ*p.rh2 - d.g*p.rh2*(p.rh2*(PI) - d.DE);
  return (I);
}

double  theory2(double x,t_qif p, t_data d) {  
  double I;
  I = d.eta + p.vh2*p.vh2 - p.rh2*p.rh2*(PI)*(PI) + d.J*p.rh2*(PI) - d.g*p.rh2*(PI)*(p.vh2 - d.E);
  return (I);
}
