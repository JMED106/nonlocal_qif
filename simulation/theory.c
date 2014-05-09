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

T_FR Theory(t_data *d, T_FR fr) {
  /* counters */

  d->rx = fr.r;
  d->vx = fr.v;

  /* Calculamos la integral mediante una suma de Riemann */
  fr.S = Coupling(*d,fr);
  d->S = fr.S;

  fr.r = rk4_void(d->rx,d->dt,rdot,d);
  fr.v = rk4_void(d->vx,d->dt,vdot,d);

  return fr;
}

double Coupling(t_data d, T_FR fr) {
  int i;
  double s = 0;
  double norm = 2*PI;
  for(i = 0; i < d.l; i++) 
    s += J_x(d,fr.x,(*(d.FR))[i].x)*(*(d.FR))[i].r;
  /* sprintf(mesg,"La suma para fr.x: %lf da %lf ",fr.x,s ); */
  /* DEBUG(mesg); */
  return s/norm;		/* Falta la normalización... */
}

double J_x(t_data d,double x, double x_prima) {
  return (d.J0 + d.J1*cos(x_prima-x));
}
  
/*++++++++++++++++++++++++++++++++++++++

  ++++++++++++++++++++++++++++++++++++++*/

double rdot(double x, t_data *d) {
  return d->Deta/PI + 2.0*d->vx*x;
}

double vdot(double x, t_data *d) {
  return d->eta + x*x-d->rx*d->rx*PI*PI + d->S;
}

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
