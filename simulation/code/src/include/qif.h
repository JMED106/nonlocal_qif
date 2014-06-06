#ifndef __QIF__H__
#define __QIF__H__

#include "common.h"


#ifndef debug
extern int debug;
extern int d_level;
#endif
#ifndef GNUPLOT_PATH
#define GNUPLOT_PATH "/usr/bin/gnuplot"		/* Ruta del binario de gnuplot */
#endif


#ifndef PI
#define PI (4*atan(1.0))
#endif
extern char DATA_FILE[50];
extern FILE *gp;
extern unsigned long int **Spikes;
extern int *spike2;
extern int **spike;


/* Functions declaration. */
double Euler(double x, double eta, double Phi, t_data d, double (*F)(double v, double Phi, double eta, t_data d));
double Vdot(double v, double Phi, double eta,t_data d);
void F_QIF(double *QIF, double *qif_pr,double Phi, int pop, t_data d);

#endif
