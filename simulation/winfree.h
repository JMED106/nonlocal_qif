#ifndef __WINFREE__H__
#define __WINFREE__H__

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


#endif

