#ifndef __INTG__H__
#define __INTG__H__

#include "common.h"

double Heun(double , t_th, double, double (* F)( double , t_th) );
double rk4(double , t_th, double ,double (*F)(double ,t_th));
double rk42(double , t_th, t_data ,double (*F)(double ,t_th,t_data));
#endif
