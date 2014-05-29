#ifndef __UTILS__H__
#define __UTILS__H__

#include <stdio.h>
#include "common.h"

#ifndef DEBUG
extern int debug;
extern int d_level;
extern char mesg[1024];
#endif

#define ARG_ERROR {fprintf(stderr,"ERROR: Bad format for arguments.\nTry '--help' option for more information\n");exit(1);}
#define READ_ERROR {fprintf(stderr,"ERROR: The file is not readable or it does not exist.\nDo you have reading permissions?\n");exit(1);}
#define TRUE 1
#define FALSE 0

#ifndef DEBUG
#define DEBUG(msg) if (debug == 1) {fprintf(stderr, "\n%s:%3d %s", __FILE__, __LINE__, (msg)); fflush(stderr);}
#define DEBUG2(msg) if (debug == 1 && d_level > 1) {fprintf(stderr, "\n%s:%3d %s", __FILE__, __LINE__, (msg)); fflush(stderr);}
#define DEBUG3(msg) if (debug == 1 && d_level > 2) {fprintf(stderr, "\n%s:%3d %s", __FILE__, __LINE__, (msg)); fflush(stderr);}
#endif
#define ESPERA while(getchar() != '\n');


#ifndef GNUPLOT_PATH
#define GNUPLOT_PATH "/usr/bin/gnuplot"		/* Ruta del binario de gnuplot */
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define NP_CHECK(ptr) \
    { \
        if (NULL == (ptr)) { \
            fprintf(stderr, "%s:%d NULL POINTER: %s n\n", \
                __FILE__, __LINE__, #ptr); \
            exit(-1); \
        } \
    } \


extern char DATA_FILE[50];
int Help(void);
t_data Arg(char *, t_data );
t_data Variables(char , char * , t_data d);
char Archivo(char *, int , FILE **);
t_data Scan_Data(char *, t_data );
void Gnuplot_Init(t_data d, char *day, char *hour);
int Create_Dir(char dir[]);
#endif
