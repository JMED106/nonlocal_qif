#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include "field.h"

void N_F(double *J,double *rx,int i,  t_data d) {
  (void) signal(SIGFPE,WarningFPE);
  int j;
  
  J[i] = 0.0;
  for(j=1 ;j<=d.l ;j++ ) 
    J[i] += J_F(d,i,j,rx)/d.l;
}


double J_F(t_data d, int i, int j, double *rx) { 
  (void) signal(SIGFPE,WarningFPE); 
  return (d.J0 + 2.0*d.J1*cos((i-j)*d.dx))*rx[j];
}
