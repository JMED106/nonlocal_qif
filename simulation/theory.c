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

T_FR Theory(t_data *d, T_FR *fr) {
  /* counters */

  d->rx = fr->rp[d->t];
  d->vx = fr->v;

  /* Calculamos la integral mediante una suma de Riemann */
  fr->S = Coupling(*d,*fr);
  d->S = fr->S;
  /* printf("\nSuma: %lf",d->S); */
  /* ESPERA */
  /* fr->rp[d->t+1] = rk4_void(d->rx,d->dt,rdot,d); */
  /* fr->v        = rk4_void(d->vx,d->dt,vdot,d); */

  /* fr->rp[d->t+1] = Heun_void(d->rx,d->dt,rdot,d); */
  /* fr->v        = Heun_void(d->vx,d->dt,vdot,d); */

  fr->rp[d->t+1] = Euler_void(d->rx,d->dt,rdot,d);
  fr->v        =   Euler_void(d->vx,d->dt,vdot,d);



}

double Coupling(t_data d, T_FR fr) {
  int i;
  float s = 0.;
  double norm = d.l;
  
  for(i = 1; i <= d.l; i++)  {
    s += J_x(d,fr.x,i)*((*(d.FR))[i].rp[d.t]);
    /* if(d.t%10 == 0) { */
  }

    /* if(d.t < 5 || d.t > 200) { */
    /*   sprintf(mesg,"La suma para fr.x %d en t: %d: da %lf ",fr.x,d.t,s ); */
    /*   DEBUG(mesg); */
    /* } */

    /* if(d.t == 5 && fr.x == d.l) ESPERA; */
    /* } */
  s = s/norm;
  return s;		/* Falta la normalización... */
}

double J_x(t_data d,int x, int x_prima) {
  return (d.J0 +2.0*d.J1*cos((x_prima-x)*d.dx));
}
  
/*++++++++++++++++++++++++++++++++++++++

  ++++++++++++++++++++++++++++++++++++++*/

double rdot(double x, t_data *d) {
  return d->Deta/M_PI + 2.0*d->vx*x;
}

double vdot(double x, t_data *d) {
  float v;
  v = d->eta + x*x-d->rx*d->rx*M_PI*M_PI + d->S;

  return v;
}

double  theory1(double x,t_qif p, t_data d) {  
  double I;
  I = d.Deta/(M_PI) + 2.*p.rh2*p.vh2 + d.DJ*p.rh2 - d.g*p.rh2*(p.rh2*(M_PI) - d.DE);
  return (I);
}

double  theory2(double x,t_qif p, t_data d) {  
  double I;
  I = d.eta + p.vh2*p.vh2 - p.rh2*p.rh2*(M_PI)*(M_PI) + d.J*p.rh2*(M_PI) - d.g*p.rh2*(M_PI)*(p.vh2 - d.E);
  return (I);
}
