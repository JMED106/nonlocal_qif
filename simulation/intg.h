#ifndef __INTG__H__
#define __INTG__H__

#include "common.h"

double Heun(double , t_qif, double, double (* F)( double , t_qif) );
double rk4(double , t_qif, double ,double (*F)(double ,t_qif));
double rk42(double , t_qif, t_data ,double (*F)(double ,t_qif,t_data));
 /* double rk4_void(double , void struct, void struct ,double (*F)(double ,void struct,void struct)); */

#endif
