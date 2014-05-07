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

/* Functions declaration. */
t_qif QIF(t_data, t_qif );
double qif (double ,t_qif );

#endif

