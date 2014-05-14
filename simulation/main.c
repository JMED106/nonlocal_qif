/*********************************************/
/* Simulation of the non-localized FR model. */
/*********************************************/

/*******************/
/* Version alfa 1. */
/*******************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<gsl/gsl_rng.h>		/* GNU Scientific library */
#include<gsl/gsl_randist.h>
#include<ctype.h>
#include<limits.h>
#include<omp.h>		        /* Parallel OMP library */
#include"gnuplot_i.h"		/* Gnuplot pipeline suit */
#include<unistd.h>

#define PI 3.14159265358979323846

/* Global variables declarationx */
int debug = 0,			/* Boolean: debugging activator */
  d_level = 0;			/* Debugging level */
char DATA_FILE[50] = "input_data.conf";
char mesg[1024];
char mesg2[1024];

/* Additional libraries (custom libraries) */
#include "utils.h"
#include "file.h"
#include "qif.h"
#include "theory.h"
#include "common.h"


/* Functions declaration. */
void    Intro(t_data *d, int *Nscan, int *tr_TT);
double *InitCond(double *p , int N, int dist_type,  double center, double gamma);
void    InitialState(t_qif * th, t_data , int type);
void InitialState_FR(T_FR *fr,t_data d, int type);
char   *DataDebug(t_data d, FILE **file);
double  MeanField(t_qif *th, t_data , int type);
t_data *Var_update(t_data *d);
void    Data_Files(t_data *d);
void    R_calc(t_data d,double *fr, double *voltage);
void    R_script(t_data d, double x, double y);
void    Perturbation(t_qif **th, t_data d,int type);

/* === FUNCTION  main ====================
 * Description:  Program Body
 *   Variables:  It accepts external arguments
 * ======================================= */

int main(int argc, char **argv) {

#ifdef _OPENMP    		/* Compilation with OMP */
  int numthreads;
  numthreads = omp_get_max_threads();
  omp_set_num_threads(numthreads);
#endif  

  /* +++++++++++++++++ External Things ++++++++++++++++++++++ */

  time_t ti;   struct tm *tm;  ti = time(NULL); tm = localtime(&ti);
  if((argc >1) && (argv[1][1] == 'd')) {debug = 1; d_level = 1;} /* Initial debugging */

  srand(time(NULL));
  gsl_rng_env_setup();		/* Random generator initialization */
  gsl_rng_default_seed = rand()*RAND_MAX;
  sprintf(mesg,"Seed: %ld ",gsl_rng_default_seed);
  DEBUG(mesg);
  
  /* Date variables */
  char day[100], hour[100];
  strftime(day, 100, "%m-%d-%Y",tm);
  sprintf(hour,"%d.%02d",tm->tm_hour,tm->tm_min);

  /* Data store */
  Create_Dir("results");
  FILE *file[6];
  FILE *fileS[4];
  FILE *volt;

  t_file FileT,FileS;

  /* +++++++++++++++++ Simulation variables ++++++++++++++++ */
  int i,j,t, total_t;
  int def;
  int time_correction;

  /* Parameters */
  t_data *d,*data;
  d = malloc(sizeof(t_data));
  data = malloc(sizeof(t_data));

  /* Results Variables */
  int Nscan;

  /* Dynamic variables (system variables) */
  t_qif **neur;
  T_FR *FR;
  double time, x;
  
  /* ++++++++++++++++++++++++++++++++++++++++++++++++++ */
  
  /* Creating data files */
  sprintf(d->file,"%s_%s",day,hour);
  sprintf(mesg2,"results/%s",d->file);
  Create_Dir(mesg2);
  Create_Dir("temp");
  sprintf(mesg2,"cp volt_avg.R ./temp"); system(mesg2);

  /*********************************/
  /* Program arguments assignation */
  /*********************************/
  *data = Scan_Data(DATA_FILE,*d);
  d = data;

  while((argc--) != 1) {	/* Terminal arguments can be handled */
    *data = Arg(argv[argc], *d);
    d = data;
    if(argv[1][0] != '-') def = 0;
  }

  FR = (T_FR*) calloc (d->l,sizeof(T_FR)); /* We create our spatially extended systems (number of columns) */
  /* Each FR[i] represents a cluster of neurons, i.e. the columns */

  neur = (t_qif**) calloc(d->l,sizeof(t_qif*));
  for(j=0 ;j<d->l ;j++ ) {
    neur[j] = (t_qif*) calloc (d->N,sizeof(t_qif));
    FR[j].rp = (double*) calloc ((int)((double)d->TT/d->dt) + 2, sizeof(double));
  }

  /* Initial condition of the firing rate is stored in the 0th position of the vector */
  d->FR = &FR;
  d->QIF = &neur;

#ifdef _OPENMP			/* Work is divided between the cores */
  int chunksize = d->l/numthreads;
  if(chunksize > 5) chunksize = 5;
#endif

  /* Initial tunning: step size, and scanning issues */
  Intro(d,&Nscan,&time_correction);
  total_t = (int)((float)d->TT/d->dt);
  d->t = 1;
  d->DX = 2.0*PI;
  d->dx = d->DX/(d->l*1.0);


  /**********************************/
  /* New simulations can start here */
  /**********************************/
  do {    
    if(d->scan_mode >= 1) 
      d = Var_update(d);
    Data_Files(d);		/* Paths to the results' files */
    if(def == 0)		/* Debug message: data display */
      DEBUG2(DataDebug(*d,&(file[0])));

    FileT = LoadFileLibrary(d->ftime, d->tmodes);
    FileT.multiopen(&file,&FileT);
    /* InitialState_FR(FR,*d,4); */
    for(i=0 ;i<d->l ;i++ ) {
      FR[i].rp[0] = (d->J0 + sqrt(pow(d->J0,2) + 4.0*PI*PI*d->eta))/(2*PI*PI)+ 0.0001;
      FR[i].v = 0.0;
    }
    FR[50].rp[0]+=0.0000;

    /* Reset counters */
    t = 0;
    d->t = 0;
    time = t*d->dt;
    /* Run */
    do {			/* Tiempo */
      if(t%(total_t/10) == 0) { /* Control point */
	sprintf(mesg,"%d%% ",(int)(t*100.0/total_t));
	DEBUG(mesg);
      }

#pragma omp parallel for schedule(dynamic,chunksize)
      for(i=0 ;i<d->l ;i++ ) {	/* Espacio */
	/* Asignamos la posición en la que nos encontramos x_i = -PI + i*dx, dónde dx = 2PI/l */
	FR[i].x =  -PI + i*d->dx;
	/* Ejecutamos la función de las eqs. FR */
	FR[i] = Theory(d,FR[i]);
      }
      for(i=0 ;i<d->l ;i++ ) {
	fprintf(file[1],"%lf\t",FR[i].rp[d->t]);
	fprintf(file[2],"%lf\t",FR[i].v);
      }
      fprintf(file[1],"\n");
      fprintf(file[2],"\n");

      t++; d->t++;
      time = t*d->dt;
    } while(time < d->TT);	/* Paso temporal (se puede hacer de la misma manera que en el QIF */
    sprintf(mesg,"%d%% ",(int)(t*100.0/total_t));
    DEBUG(mesg);

    for(i=0 ;i<d->l ;i++ ) {
      x = -PI + i*d->dx;
      fprintf(file[0],"%lf\t%lf\t%lf\t%lf\n",x,FR[i].rp[d->t],FR[i].v,J_x(*d,x,0));
    }
	
    FileT.closeall(&file,&FileT);
    d->scan++;
  } while (d->scan < Nscan);  /* Simulation ends here */

  system("rm -r ./temp");
  sprintf(mesg,"cp ./results/%s/*.txt ~/Escritorio/gp/",d->file);
  system(mesg);
  printf("\n");
  return 0;  
}


/* === FUNCTION  InitCond ====================
 * Description:	 Sets initial condition of the 
 *               neuron population given a 
 *               distribution.
 *   Variables:  Neuron vector, distr ID.
 *               min, max, wide(gamma)
 * ======================================= */

double *InitCond(double *p, int N, int distr_type, double center, double gamma) {
  int i;			/* Counters */
  double k;
  gsl_rng *r;

  gsl_rng_default_seed = rand()*RAND_MAX;

  const gsl_rng_type *T;
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  
  NP_CHECK(p);
  switch(distr_type) {
  case 0:			/* Cauchy ??? */
    for(i=0 ;i<N ;i++ ) 
      p[i] = center + gsl_ran_cauchy(r, gamma);
    break;
  case 1:			/* Gaussian (random) */
    for(i=0 ;i<N ;i++ ) 
      p[i] = center + gsl_ran_gaussian(r,gamma);
    break;
  case 2:			/* Cauchy ordered */
    for(i=0 ;i<N ;i++ ) {
      k = (2.0*(i+1) - N -1.0)/(N+1.0);
      p[i] = center + gamma*tan((PI/2.0)*k);
    }
    break;
  default:
    ARG_ERROR;
    break;
  }

  return p;
}


/* === FUNCTION  InitialState ====================
 * Description:  Oscillators initial state setup
 *   Variables:  none
 * ======================================= */

void InitialState(t_qif *th,t_data d, int type) {
  int i;
  double k = 0, h, v;
  gsl_rng *r;
  gsl_rng_env_setup();		/* Random generator initialization */
  gsl_rng_default_seed = rand()*RAND_MAX;

  const gsl_rng_type *T;
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  /* We must define the relationship between v-th
   * v = tg(th/2). 
   * But th goes from 0 to PI 
   * v goes from v_reset to v_peak (finites) */

  switch(type) {
  case 0:			/* Uniform distribution between reset and peak values */
    for(i=0 ;i < d.N ;i++ ) {
      th[i].v = d.vr + (d.vp - d.vr)*gsl_rng_uniform(r);
    }
    break;
  case 2:			/* Null distribution: constant 0 valued. */
    for(i=0 ;i < d.N ;i++ ) {
      th[i].v = 0.0;
    }
    break;
  case 4:			/* Lorentzian distribution centered at 10 and width of 5 */
    h = 5;			/* (take into account that the distribution is cut at vr */
    v= 10;			/*  and vp) */
    for(i=0 ;i < d.N ;i++ ) {
      k = (2.0*(i+1) -d.N -1.0)/(d.N+1.0);
      th[i].v = v + h*tan((PI/2.0)*k);
      if(fabs(th[i].v) > d.vp) th[i].v = v;
    }
    break;
  default:
    ARG_ERROR;
    break;
  }
}

void InitialState_FR(T_FR *fr,t_data d, int type) {
  int i;
  double k = 0, h, v;
  gsl_rng *r;
  gsl_rng_env_setup();		/* Random generator initialization */
  gsl_rng_default_seed = rand()*RAND_MAX;

  const gsl_rng_type *T;
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  /* We must define the relationship between v-th
   * v = tg(th/2). 
   * But th goes from 0 to PI 
   * v goes from v_reset to v_peak (finites) */

  switch(type) {
  case 0:			/* Uniform distribution between reset and peak values */
    for(i=0 ;i < d.l ;i++ ) {
      fr[i].rp[0] = 1.0 + 1.0*gsl_rng_uniform(r);
    }
    break;
  case 2:			/* Null distribution: constant 0 valued. */
    for(i=0 ;i < d.l ;i++ ) {
      fr[i].rp[0] = 0.0;
    }
    break;
  case 4:			/* Lorentzian distribution centered at 10 and width of 5 */
    h = 1;			/* (take into account that the distribution is cut at vr */
    v= 1;			/*  and vp) */
    for(i=0 ;i < d.l ;i++ ) {
      k = (2.0*(i+1) -d.l -1.0)/(d.l+1.0);
      fr[i].rp[0] = v + h*tan((PI/2.0)*k);
      if(fabs(fr[i].rp[0]) > 1) fr[i].rp[0] = v;
    }
    break;
  default:
    ARG_ERROR;
    break;
  }
}


/* === FUNCTION  MeanField ====================
 * Description:  Computes S
 *   Variables:  all th_i
 * ======================================= */

double MeanField(t_qif *th, t_data d,int type) {
  int i;
  double s = 0;
  double delta = 1.0/d.dt;
  double norm = (1.0*PI/d.N)*delta;

  if(type == 1) {
    th[0].global_s1 = 0;
    for(i=0 ;i<d.N ;i++ ) {
      if(th[i].spike == 1)
	s+= 1;
    }
    th[0].global_s1 = s;
  } 

  s = s*norm;		/* Something wrong with couplinh (MF is not OK) */
  return s;
}

/* === FUNCTION  data ====================
 * Description:  Updates target variable in scan mode
 *   Variables:  Data structure
 * ======================================= */

t_data *Var_update(t_data *d) {
  t_data *d2;
  d2 = malloc(sizeof(*d2));
  char *var_name[4] = {"J","eta","g","E"};
  *d2 = *d;
  switch(d->variable) {
  case 1:			/* J */
    if(d->scan == 0)
      d2->J = d->min;
    else
      d2->J = d2->J + (d->step);
    sprintf(mesg,"Variable [%s] changed from: %lf to %lf.",
	    var_name[d->variable-1],d->J,d2->J);
    d2->var_value = d2->J;
    break;
  case 2:			/* eta */
    if(d->scan == 0)
      d2->eta = d->min;
    else
      d2->eta = d2->eta + (d->step);
    sprintf(mesg,"Variable [%s] changed from: %lf to %lf.",
	    var_name[d->variable-1],d->eta,d2->eta);
    d2->var_value = d2->eta;
    break;
  case 3:			/* g */
    if(d->scan == 0)
      d2->g = d->min;
    else
      d2->g += d->step;
    sprintf(mesg,"Variable [%s] changed from: %lf to %lf.",
	    var_name[d->variable-1],d->g,d2->g);
    d2->var_value = d2->g;
    break;
  case 4:			/* E */
    if(d->scan == 0)
      d2->E = d->min;
    else
      d2->E += d->step;
    sprintf(mesg,"Variable [%s] changed from: %lf to %lf.",
	    var_name[d->variable-1],d->E,d2->E);
    d2->var_value = d2->E;
    break;
  default:
    d2->scan_mode = 0;
    break;
  }

  DEBUG(mesg);
  return d2;
}


/* === FUNCTION  Intro ====================
 * Description:  Intro to the program
 *   Variables:  Data structure d->
 * ======================================= */

void Intro(t_data *d, int *Nscan, int *tr_TT) {

  int max_bits = 8*sizeof(unsigned long int); /* Architecture of the machine */
  char *var_name[4] = {"J","eta","g","E"};
  double tr, tr1 = 0, tr2 = 0;	/* Refractory period */
  double dt, scanstep = 0;


  /** DEBUG ***********************************************/
  sprintf(mesg,"Entering debug mode [%d] ...",d_level);
  DEBUG(mesg);
  sprintf(mesg,"We are running on a %d bits system\n",max_bits);
  DEBUG(mesg);

  if(d->scan_mode >= 1) {
    sprintf(mesg,"Scanning mode: ON. Variable: [%s]", var_name[d->variable-1]);
    DEBUG(mesg);
  } else {
    sprintf(mesg,"Scanning mode: OFF. ");
    DEBUG(mesg);
  }

  scanstep = d->step;
  *Nscan = (int)ceil(1+(fabs(d->max - d->min)/((float)d->step)));
  if(d->max < d->min) d->step = (-1.0)*scanstep;
  if(d->scan_mode == 0) *Nscan = 0;

  /* /\* Time step setup *\/ */
  /* dt = d->dt; */
  /* /\* dt(tr) adjustment *\/ */
  /* tr1 = 1.0/fabs(d->vp); */
  /* tr2 = 1.0/fabs(d->vr); */
  /* tr = (tr1 < tr2)? tr1:tr2; */
  
  /* dt = tr; */
  /* if(dt != d->dt) { */
  /*   d->dt = dt; */
  /*   sprintf(mesg,"Time step dt adjusted to: dt = %lf ",dt ); */
  /*   DEBUG(mesg); */
  /* } */

  /* *tr_TT = (int)(tr/dt); */
  /* sprintf(mesg,"Refractory period (in steps): %d ",*tr_TT); */
  /* DEBUG(mesg); */

  if(d->N*d->TT > d->MaxDim)
    d->disable_raster = 1;
  else
    d->disable_raster = 0;
}


/* === FUNCTION  DataDebug ====================
 * Description:  Displays selected data
 *   Variables:  data structure
 * ======================================= */

char *DataDebug(t_data d,FILE **file) {
  OpenFile(file,d.file_parameters,"w");  
  sprintf(mesg,"\nCustom parameters\n"
	  "-----------------\n"
	  "\t Reseting potential (r ) = %5.2lf\n"
	  "\t Peak potential (p)      = %5.2lf\n"
	  "\t Value of J (J)          = %5.2lf\n"
	  "\t Value of DJ (m)         = %5.2lf\n"
	  "\t Value of eta (e)        = %5.2lf\n"
	  "\t Value of Deta (d)       = %5.2lf\n"
	  "\t Value of E (v)          = %5.2lf\n"
	  "\t Value of DE (B)         = %5.2lf\n"
	  "\t Value of g (g)          = %5.2lf\n"
	  "\t Total number of Neurons = %5d\n"
	  "\t Simulation Time         = %5.2lf\n"
	  "\t Time step               = %5.4lf\n", d.vr, d.vp, d.J,d.DJ,d.eta,d.Deta,d.E,d.DE,d.g, d.N, d.TT,d.dt);
  fprintf(*file,mesg);
  CloseFile(file);
  return mesg;
}

void Data_Files(t_data *d) {
  int i;
  /* File paths */
  d->ftime = (char**) malloc (7*sizeof(char*));
  d->tmodes = (char**) malloc (6*sizeof(char*));
  d->fscan = (char**) malloc (4*sizeof(char*));
  d->smodes = (char**) malloc (4*sizeof(char*));
  for(i = 0; i < 7; i++) 
    (d->ftime)[i] = (char*) malloc(256*sizeof(char));
  (d->ftime)[6] = NULL;

  sprintf(d->file_parameters,"./results/%s/parameters.txt",d->file);

  sprintf(d->file_rvJ_FR,"./results/%s/rvJ_FR_%.lf.txt",d->file,d->var_value);
  (d->ftime)[0] = d->file_rvJ_FR;
  (d->tmodes)[0] = "w";
  sprintf(d->file_rt_FR,"./results/%s/rt_FR_%.lf.txt",d->file,d->var_value);
  (d->ftime)[1] = d->file_rt_FR;
  (d->tmodes)[1] = "w";
  sprintf(d->file_vt_FR,"./results/%s/vt_FR_%.lf.txt",d->file,d->var_value);
  (d->ftime)[2] = d->file_vt_FR;
  (d->tmodes)[2] = "w";
  sprintf(d->file_rvJ_QIF,"./results/%s/rvJ_QIF_%.lf.txt",d->file,d->var_value);
  (d->ftime)[3] = d->file_rvJ_QIF;
  (d->tmodes)[3] = "w";
  sprintf(d->file_rt_QIF,"./results/%s/rt_QIF_%.lf.txt",d->file,d->var_value);
  (d->ftime)[4] = d->file_rt_QIF;
  (d->tmodes)[4] = "w";
  sprintf(d->file_vt_QIF,"./results/%s/vt_QIF_%.lf.txt",d->file,d->var_value);
  (d->ftime)[5] = d->file_vt_QIF;
  (d->tmodes)[5] = "w";

  if(d->scan_mode > 0) {    
    sprintf(d->file_rvp_FR,"./results/%s/rvp_FR_%.lf.txt",d->file,d->var_value);
    sprintf(d->file_rvp_QIF,"./results/%s/rvp_QIF_%.lf.txt",d->file,d->var_value); 
  }
}


/* === FUNCTION  R ====================
 * Description:  Calls R script and store variables x, y
 *   Variables:  
 * ======================================= */

void  R_calc(t_data d, double *fr, double *voltage) {
  FILE *Rcmd;
  char cmd[30], cmd2[30], cmd3[30], cmd4[30];
  Rcmd = popen("R -e \"source('./volt_avg.R')\" -q -s","r");
  fscanf(Rcmd,"%s %s\n%s %s",cmd, cmd2,cmd3,cmd4);
  pclose(Rcmd);
  *fr =atof(cmd2);
  *voltage = atof(cmd4);
}


/* === FUNCTION  R ====================
 * Description:  Creates R script for R_calc
 *   Variables:  vpeak, aprox x and y (theta)
 * ======================================= */

void R_script(t_data d, double x, double y) {
  FILE *Rscript;
  Rscript = fopen("volt_avg.R","w");
  
  fprintf(Rscript,"setwd('./temp');\n"
	  "volt_dist_files = list.files(pattern= 'volt_dist*');\n"
	  "voltages <- numeric(length(volt_dist_files));\n"
	  "firing_rates <- numeric(length(volt_dist_files));\n"
	  "for (i in volt_dist_files) \{\n"
	  "    data<-read.table(i);\n"
	  "    data2 <- data[abs(data[,1])<=%d,];\n"
	  "    data3 <- hist(data2,breaks=seq(-%d,%d,by=0.1),freq=FALSE);\n"
	  "    x <-data3$breaks[1:length(data3$counts)]\n"
	  "    datax <- data.frame(x,data3$density)\n"
	  "    lorentz <- nls(data3$density ~(1/pi)*a/((x-b)*(x-b)+a*a),data=datax,start=list(a=%lf,b=%lf),algorithm=\"plinear\", nls.control(maxiter = %d , tol = %e, minFactor = 1/%d, printEval = FALSE, warnOnly = FALSE))\n"
	  "    voltages[i] <- coef(lorentz)[2];\n"
	  "    firing_rates[i] <- coef(lorentz)[1];\n"
	  "}\n"
	  "avg_volt <- mean(voltages);\n"
	  "avg_fr <- mean(firing_rates);\n"
	  "print(avg_fr);\n"
	  "print(avg_volt);\n",(int)d.vp,(int)d.vp,(int)d.vp,x*(PI),y,50,1e-5,2048);
  fclose(Rscript);
}


/* === FUNCTION  Perturbation ====================
 * Description:  Perturbation of the system
 *   Variables:  Voltages, amplitud
 * ======================================= */

void Perturbation(t_qif **th, t_data d,int type)  {
  int i;
  if(type == 0) 
    for(i=0 ;i<d.N ;i++ ) {
      ((*th)[i].v) += d.pert_amplitude;
    }
  else if(type == 1)
    (* th)[0].FR = d.pert_amplitude;
}
