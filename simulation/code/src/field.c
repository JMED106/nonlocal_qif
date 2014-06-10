#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "field.h"

void N_F(double *J,double *rx,int i,  t_data d) {
  int j;
  J[i] = 0.0;
/* #pragma omp parallel pivate(j) */
  for(j=1 ;j<=d.l ;j++ ) 
    J[i] += J_F(d,i,j,rx);
  J[i] /= d.l;
}


double J_F(t_data d, int i, int j, double *rx) {  
  return (d.J0 + 2.0*d.J1*cos((i-j)*d.dx))*rx[j];
}
