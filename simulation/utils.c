/** Various generic utilities */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


#include "utils.h"
#include "gnuplot_i.h"



/* Help del programa */
int Help() {
  exit(1);
}

/* Función para gestionar los archivos de resultados */
char Archivo(char archivo[], int s, FILE **fin)
{ 
  char archivo0[20];
  sprintf(archivo0,"%s_%d.txt",archivo,s);

  printf("\nOpening results file    [%15s] .............",archivo0);
  if((*fin = fopen(archivo0, "r")) == NULL) {		/* Intentamos leer el archivo */
    printf("FAILED");
    printf("\nCreating results file   [%15s] .............",archivo0);
    if ((*fin = fopen (archivo0, "w")) == NULL) {	/* Creamos el archivo de resultados */
	printf("FAILED");
	fprintf (stderr, "\nERROR: Impossible to create file '%s' ......... TERMINATED\n\n", archivo0);
	exit(1);
    } else printf("%6s","OK"); 
  } else {
    printf("%6s","OK"); 
    printf("\nWriting at results file [%15s] .............",archivo0);
    if((*fin = fopen(archivo0, "a")) == NULL) {
      printf(" FAILED");
      fprintf (stderr, "\nERROR: Impossible to create file '%s' ......... TERMINATED\n\n", archivo0);
      exit(1);
    } else printf("%6s","OK");
  }
  return 0;
}


/* === FUNCTION  Gnuplot ====================
 * Description:  Plotting function
 *   Variables:  ???????
 * ======================================= */

void Gnuplot_Init(t_data d, char *day, char *hour) {
  
  /* Gnuplot Plotting (Experimental) */
  gnuplot_ctrl *g1, *g2, *g3, *g4;
  
  g1 = gnuplot_init();
  g2 = gnuplot_init();
  g3 = gnuplot_init();
  g4 = gnuplot_init();
    
  gnuplot_close(g1);
  gnuplot_close(g2);
  gnuplot_close(g3);
  gnuplot_close(g4);
}

/* === FUNCTION  Create ====================
 * Description:  Creates a new folder
 *   Variables:  Name of the folder
 * ======================================= */

int Create_Dir(char dir[]) {
  char command[100];
  sprintf(command,"cd %s 1>/dev/null 2>&1",dir);
  fflush(stdout); fflush(stdin);
  if(system(command) != 0) {
    sprintf(mesg,"Creating new directory ./%s", dir);
    DEBUG(mesg);
    sprintf(command,"mkdir %s",dir);
    system(command);
  }
  return 0;
}



/* === FUNCTION  Arg ====================
 * Description:  Handles arguments
 *   Variables:  argument, data structure
 * ======================================= */

t_data Arg(char *argv, t_data d) {
  int i = 0;
  char valor[50];
  t_data data;
  data = d;
  if(isdigit(argv[0]))
    return d;
  if(argv[0] == 'S');
  else {
    while(argv[i+2] != '\0') {
      valor[i] = argv[i+2];
      i++;
    }
    valor[i] = '\0';
  }

  if(argv[0] == '-' && argv[1] == 'd') { 
    debug = 1;
    d_level = atoi(valor);
    if(d_level == 0) d_level =1;
    return d;
  }  
  data= Variables(argv[0], valor,d); 
  return data;
}


/* === FUNCTION  Scan Data ===============
 * Description:  Scans data from file
 *   Variables:  File name, data structure
 * ======================================= */

t_data Scan_Data(char *file, t_data d) {
  char c, mesg[100],arg[10], valor[10],aux[100];
  FILE *fin;
  int lines = 0, *linea_vacia,  i = 0, j = 0, error = 0;
  t_data data;  
  data = d;
  linea_vacia = (int*) calloc (1, sizeof(int));
  
  printf("\nOpening data file %s ...............",file);
  if ((fin = fopen (file, "r")) == NULL) {			/* Abrimos el archivo de datos */
    printf(" FAILED\n");
    fprintf (stderr, "\nERROR: File ’%s’ does not exist or it can't be read......... TERMINATED\n\n", file);
    exit(1);
  }
  else
    printf(" OK\n");

  printf("Reading data file %s ...............",file);
  fflush(stdin);
  fflush(stdout);

  while((c = getc(fin)) != EOF) { /* Leemos el archivo de datos */
    if(c == '#'){ /* Si se detecta una # como primer carácter de línea... */
      linea_vacia[lines] = 1;
      DEBUG3("Comentario Encontrado");
    }
    if(c == '\n') {
      DEBUG3("Línea encontrada");
      lines++;
      linea_vacia = (int*) realloc (linea_vacia,(lines+1)*sizeof(int));
      linea_vacia[lines] = 0;
    }
  }
  rewind(fin);
  for(j = 0; j < lines; j++) {
    sprintf(mesg,"línea %d: vacía: %d ",j,linea_vacia[j] );
    DEBUG3(mesg);
    if(linea_vacia[j] == 0)  {
      DEBUG("Escaneando línea");
      fscanf(fin,"%s\n",arg);
      if(isdigit(arg[0])) {
	printf("\nERROR de lectura: la línea %d tiene formato incorrecto\n",lines);
	error = 1;
      }
      else {
	i = 0;
	while(arg[i+2] != '\0') {
	  valor[i] = arg[i+2];
	  i++;
	}
	valor[i] = '\0';
	sprintf(mesg,"Valor de %c = %s",arg[0],valor);
	DEBUG(mesg);
	sprintf(mesg,"Valor de %c = %lf (convertido)\n",arg[0],atof(valor));
	DEBUG(mesg);
      }
      data = Variables(arg[0],valor, d);
      sprintf(mesg,"data.init_dist = %d\t data.vr = %lf ",data.init_dist,data.vr );
  DEBUG(mesg);
      d = data;
sprintf(mesg,"data.init_dist = %d\t data.vr = %lf ",data.init_dist,data.vr );
 DEBUG(mesg);
    } else {
       fscanf(fin," %[^\n]s",aux);
    }
  }
  
  fclose(fin);
  if(error == 0) printf(" OK\n");
  else		 printf(" ERROR\n");

  sprintf(mesg,"data.vr = %lf ",data.vr );
  DEBUG(mesg);

  return data;
}


/* === FUNCTION  Variables================
 * Description:  Assigns value to variables
 *   Variables:  argument type, data st., value
 * ======================================= */

t_data Variables(char c, char *v, t_data d) {
  t_data data;
  data = d;

  switch(c) {
  case 'a':			/* (Provisional) J0 */
    data.J0 = atof(v);
    return data;
    break;
  case 'b':
    data.J1 = atof(v);		/* (Provisional) J1 */
    return data;
    break;
  case 'c':
    data.J2 = atof(v);		/* (Provisional) J2 */
    return data;
    break;
  case 'd':
    /* assignment */
    return data;
    break;
  case 'e':
    /* assignment */
    return data;
    break;
  case 'f':
    /* assignment */
    return data;
    break;
  case 'g':
    /* assignment */
    return data;
    break;
  case 'h':			/* External current (eta) */
    data.eta = atof(v);
    return data;
    break;
  case 'i':			/* Initial condition (distribution) */
    data.init_dist = atoi(v);
    return data;
    break;
  case 'j':
    /* assignment */
    return data;
    break;
  case 'k':
    /* assignment */
    return data;
    break;
  case 'l':			/* Number of clusters (size) */
    data.l = atoi(v);
    return data;
    break;
  case 'm':			/* Minimum value of the variable to scan */
    data.min = atof(v);
    return data;
    break;
  case 'n':
    /* assignment */
    return data;
    break;
  case 'o':
    /* assignment */
    return data;
    break;
  case 'p':			/* Peak potential QIF */
    data.vp = atof(v);
    return data;
    break;
  case 'q':
    /* assignment */
    return data;
    break;
  case 'r':			/* Reset potential QIF */
    data.vr = atof(v);
    return data;
    break;
  case 's':			/* Step of the scanned variable */
    data.step = atof(v);
    return data;
    break;
  case 't':			/* Time step of the simulation */
    data.dt = atof(v);
    return data;
    break;
  case 'u':		       
    /* assignment */
    return data;
    break;
  case 'v':
    /* assignment */
    return data;
    break;
  case 'w':
    /* assignment */
    return data;
    break;
  case 'x':
    /* assignment */
    return data;
    break;
  case 'y':
    /* assignment */
    return data;
    break;
  case 'z':
    /* assignment */
    return data;
    break;
  case 'A':
    /* assignment */
    return data;
    break;
  case 'B':
    /* assignment */
    return data;
    break;
  case 'C':
    /* assignment */
    return data;
    break;
  case 'D':
    /* assignment */
    return data;
    break;
  case 'E':
    /* assignment */
    return data;
    break;
  case 'F':
    /* assignment */
    return data;
    break;
  case 'G':
    /* assignment */
    return data;
    break;
  case 'H':			/* Delta-Eta (if H ≠ 0 distribution is enabled) */
    data.Deta = atof(v);
    return data;
    break;
  case 'I':
    /* assignment */
    return data;
    break;
  case 'J':
    /* assignment */
    return data;
    break;
  case 'K':
    /* assignment */
    return data;
    break;
  case 'L':
    /* assignment */
    return data;
    break;
  case 'M':			/* Maximum value of the scnanned variable */
    data.max = atof(v);
    return data;
    break;
  case 'N':			/* Number of neurons QIF */
    data.N = atoi(v);
    return data;
    break;
  case 'O':
    /* assignment */
    return data;
    break;
  case 'P':			/* Amplitude of the perturbation QIF */
    data.pert_amplitude = atof(v);
    return data;
    break;
  case 'Q':
    /* assignment */
    return data;
    break;
  case 'R':			/* Boolean: periodic boundaries */
    data.ring = atoi(v);
    return data;
    break;
  case 'S':			/* Scan mode boolean */
    data.scan_mode = atoi(v);
    return data;
    break;
  case 'T':			/* Total time of the simulation */
    data.TT = atof(v);
    return data;
    break;
  case 'U':
    /* assignment */
    return data;
    break;
  case 'V':
    /* assignment */
    return data;
    break;
  case 'W':
    /* assignment */
    return data;
    break;
  case 'X':
    /* assignment */
    return data;
    break;
  case 'Y':
    /* assignment */
    return data;
    break;
  case 'Z':
    /* assignment */
    return data;
    break;
  default:
    ARG_ERROR;
  }
}
