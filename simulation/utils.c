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
  int option = 0;
  char answ[100];
  char pdf1[100], pdf2[100],pdf3[100];
  double x, y;
  x = sqrt(d.eta/2.0 + 0.5*sqrt(d.eta*d.eta + d.eta_sigma*d.eta_sigma));
  y = d.eta_sigma/(2.0*x);
  
  /* Gnuplot Plotting (Experimental) */
  gnuplot_ctrl *g1, *g2, *g3, *g4;
  
  g1 = gnuplot_init();
  g2 = gnuplot_init();
  g3 = gnuplot_init();
  g4 = gnuplot_init();

  gnuplot_cmd(g2,"set title \"{/Symbol h} = %.2lf, {/Symbol s} = %.2lf\"",d.eta,d.eta_sigma);
  gnuplot_cmd(g1, "plot \"./results/%s_%s_raster.dat\" w p lw 0 title \"{/Symbol h} = %.2lf\\ {/Symbol s} = %.2lf\"",day,hour,d.eta,d.eta_sigma);
  gnuplot_cmd(g2, "plot \"./results/%s_%s_kuramoto.dat\" u 1:($1<7?1/0:$3) w l title \"r(t)\", \"\" u 1:($1<7?1/0:$2) w l title \"order\", \"\" u 1:($1<7?1/0:$4/pi) w l title \"{/Symbol Y}(t)/{/Symbol p} \" ",day,hour);
  gnuplot_cmd(g3,"set title \"{/Symbol h} = %.2lf, {/Symbol s} = %.2lf\"",d.eta,d.eta_sigma);
  gnuplot_cmd(g3,"set samples 10000");
  gnuplot_cmd(g3,"set arrow from %lf,%lf to %lf,%lf nohead lc rgb \"blue\"\n",-y,0.0,-y,1.0/(PI*x));
  gnuplot_cmd(g3,"plot ((%lf/pi)/(((%lf)+(x))**2+%lf**2)), \"./results/hist_simu.dat\";",x,y,x);

  do{
    printf("\nDo you want to export these plots to a PDF file?\n(Y/N)"); scanf("%s",answ);
    fflush(stdout); fflush(stdin);
    option = (((strcasecmp(answ,"y") == 0) || (strcasecmp(answ,"n") == 0))?1:0);
  } while (option != 1);

  printf("\n");
  if((option = ((strcasecmp(answ,"y") == 0)?1:0) ) == 1) {
    Create_Dir("plots");

    sprintf(pdf1,"./plots/%s_%s_raster.pdf",day,hour);
    sprintf(pdf2,"./plots/%s_%s_order.pdf",day,hour);
    sprintf(pdf3,"./plots/%s_%s_hist.pdf",day,hour);

    gnuplot_cmd(g1, "set encoding iso_8859_15");
    gnuplot_cmd(g1, "set grid back");
    gnuplot_cmd(g1, "set terminal postscript eps enhanced color");
    gnuplot_cmd(g1, "set output \"| epstopdf --filter --outfile=%s\"",pdf1);

    gnuplot_cmd(g2, "set encoding iso_8859_15");
    gnuplot_cmd(g2, "set grid back");
    gnuplot_cmd(g2, "set terminal postscript eps enhanced color");
    gnuplot_cmd(g2, "set output \"| epstopdf --filter --outfile=%s\"",pdf2);

    gnuplot_cmd(g3, "set encoding iso_8859_15");
    gnuplot_cmd(g3, "set grid back");
    gnuplot_cmd(g3, "set terminal postscript eps enhanced color");
    gnuplot_cmd(g3, "set output \"| epstopdf --filter --outfile=%s\"",pdf3);

    gnuplot_cmd(g1, "plot \"./results/%s_%s_raster.dat\" w p pt 0 title \"{/Symbol h} = %.2lf\\ {/Symbol s} = %.2lf\"",day,hour,d.eta,d.eta_sigma);
    gnuplot_cmd(g2, "plot \"./results/%s_%s_kuramoto.dat\" u 1:($1<7?1/0:$3) w l title \"r(t)\", \"\" u 1:($1<7?1/0:$2) w l title \"order\", \"\" u 1:($1<7?1/0:$4/pi) w l title \"{/Symbol Y}(t)/{/Symbol p}\" ",day,hour);

    gnuplot_cmd(g3,"plot ((%lf/pi)/((%lf+x)**2+%lf**2)) title \"Theory\", \"./results/hist_simu.dat\" title \"Simulation\";",x,y,x);
  }
  fflush(stdin); fflush(stdout);

  sprintf(answ,"cp ./results/hist_simu.dat ./results/%s_%s_hist_simu.dat",day,hour);
  system(answ);
  sprintf(answ,"cp ./results/volt.dat ./results/%s_%s_voltdist.dat",day,hour);
  system(answ);
    
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
    case 'i':			/* Initial configuration state */
      data.init_dist = atoi(v);
      return data;
      break;
    case 'r':			/* Reseting Potential V_res */
      data.vr = atof(v);
      return data;
      break;
    case 'p':			/* Peak Potential V_peak */
      data.vp = atof(v);
      return data;
      break;
    case 'j':			/* Value of J0 */
      data.J = atof(v);
      return data;
      break;
    case 'J':			/* J_dist {0,1} */
      data.J_dist = atoi(v);
      return data;
      break;
    case 'm':			/* J_gamma */
      data.J_sigma = atof(v);
      return data;
      break;
    case 'N':			/* Total number of neurons */
      data.N = atoi(v);
      return data;
      break;
    case 'T':			/* Simulation time */
      data.TT = atof(v);
      return data;
      break;
    case 't':			/* Time step */
      data.dt = atof(v);
      return data;
      break;
    case 'e':			/* Value of Eta0 */
      data.eta = atof(v);
      return data;
      break;
    case 'E':			/* eta_dist */
      data.eta_dist = atoi(v);
      return data;
      break;
    case 'd':			/* eta_gamma */
      data.eta_sigma = atof(v);
      return data;
      break;
    case 'g':			/* Value of g0 */
      data.g = atof(v);
      return data;
      break;
    case 'v':			/* Value of V0 */
      data.v0 = atof(v);
      return data;
      break;
    case 'V':			/* V0_dist */
      data.v0_dist = atoi(v);
      return data;
      break;
    case 'B':			/* V0_gamma */
      data.v0_sigma = atof(v);
      return data;
      break;
    case 's':			/* Scan mode: ... */
      data.scan_mode = atoi(v);
      return data;
      break;
    case 'z':
      data.variable = atoi(v);	/* Variable to scan */
      return data;
      break;
    case 'x':			/* Minimum value */
      data.min_x = atof(v);
      return data;
      break;
    case 'Z':			/* Maximum value */
      data.max_x = atof(v);
      return data;
      break;
    case 'X':			/* Increment */
      data.dx = atof(v);
      return data;
      break;	
    case 'f':			/* Amplitude of the perturbation */
      data.pert_amplitude = atof(v);
      return data;
      break;
    case 'P':			/* FR *** Amplitude of the perturbation of FR */
      data.perturbation_FR = atoi(v);
      return data;
      break;
    case 'R':			/* OPTIONAL 1 */
      /* Assing something */
      return data;
      break;
    case 'l':			/* OPTIONAL 2 */
      /* Assing something */
      return data;
      break;
    case 'G':			/* OPTIONAL 3 */
      /* Assing something */
      return data;
      break;
    case 'b':			/* OPTIONAL 4 */
      /* Assing something */
      return data;
      break;
    default:
      ARG_ERROR;
  }
}
