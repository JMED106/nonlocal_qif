#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "qif.h"
#include "nrutil.h"
#include "nr.h"

void F_QIF(double *QIF, double *qif_pr,double Phi, int pop, t_data d) {
  int j;
  double QIFJ;

  double tau_v;
  
  double (*fptr)(double , double, double,  t_data);
		/*  v      Phi     eta    parameters  */
  fptr = Vdot;
    
  for(j=1 ;j<=d.N ;j++ ) {	/* For each neuron */
    if(Spikes[pop][j] == 0UL) {
      QIFJ = QIF[j];
      if((QIF[j] = Euler(QIFJ,qif_pr[j],Phi, d,fptr)) >= d.vp ) { 
	tau_v = 1.0/(QIF[j]);
	Spikes[pop][j] = 1UL<<(d.wait - (int)ceil((d.tau_p-tau_v)/d.dt));
      }
    }
    else {      
      Spikes[pop][j] = Spikes[pop][j]>>1;
      if(Spikes[pop][j] == 1UL<<((int)(d.wait/2.0))) 
	QIF[j] = d.vr;
    }
  }

}


double Vdot(double v, double Phi, double eta,t_data d) {
  return (pow(v,2.0) + Phi + eta);
}


double Euler(double x, double eta, double Phi, t_data d, double (*F)(double v, double Phi, double eta, t_data d)) {
  return x + d.dt*(F(x,Phi,eta,d));
}
	   

