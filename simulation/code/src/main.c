/*********************************************/
/* Simulation of the non-localized FR model. */
/*********************************************/

/*******************/
/* Version beta 0. */
/*******************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<gsl/gsl_rng.h>		/* GNU Scientific library */
#include<gsl/gsl_randist.h>
#include<omp.h>		        /* Parallel OMP library */



/* Global variables declaration */
int debug = 0,			/* Boolean: debugging activator */
  d_level = 0;			/* Debugging level */
char DATA_FILE[50] = "input_data.conf";
char mesg[1024];
char mesg2[1024];

/* Additional libraries (custom libraries) */
#include "common.h"		/* This must be first */
#include "utils.h"
#include "file.h"
#include "nrutil.h"
#include "nr.h"
#include "field.h"

/* Functions declaration. */
double *InitCond(double *p, int N, int distr_type, double center, double gamma);
void    InitialState(t_qif *th,t_data d, int type);
double  MeanField(t_qif *th, t_data d,int type);
t_data *Var_update(t_data *d);
void    Intro(t_data *d, int *Nscan, int *tr_TT);
char   *DataDebug(t_data d,FILE **file);
void    Data_Files(t_data *d);

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
  sprintf(mesg,"Seed: %ld ",gsl_rng_default_seed);  DEBUG(mesg);
  
  /* Date variables */
  char day[100], hour[100];  strftime(day, 100, "%m-%d-%Y",tm);  sprintf(hour,"%d.%02d",tm->tm_hour,tm->tm_min);

  /* Data store */
  Create_Dir("results");

  /* +++++++++++++++++ Simulation variables ++++++++++++++++ */
  int i,j,t, t_max;
  int def = 0;
  int time_correction = 0;
  int TIME;

  /* Parameters (all the parameters are stored in a structure) */
  t_data *d,*data;
  d = malloc(sizeof(t_data));
  data = malloc(sizeof(t_data));

  /* Results Variables */
  int Nscan;

  FILE *file[6];
  for(i=0 ;i<6 ;i++ ) {
    file[i] = malloc(sizeof(FILE));
  }
  t_file FileT;


  /* Dynamic variables (system variables) */
  double *R,*R2,*V;
  double *J;

  double **qif;
  
  
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

  d->ftime = (char**) malloc (6*sizeof(char*));
  d->tmodes = (char**) malloc (6*sizeof(char*));
  d->fscan = (char**) malloc (4*sizeof(char*));
  d->smodes = (char**) malloc (4*sizeof(char*));
  for(i = 0; i < 5; i++) 
    (d->ftime)[i] = (char*) malloc(256*sizeof(char));
  (d->ftime)[5] = NULL;

  while((argc--) != 1) {	/* Terminal arguments can be handled */
    *data = Arg(argv[argc], *d);
    d = data;
    if(argv[1][0] != '-') def = 0;
  }

  R = dvector(1,d->l);
  R2 = dvector(1,d->l);
  V = dvector(1,d->l);
  J = dvector(1,d->l);

  qif = dmatrix(1,d->l,1,d->N);

  d->R = &R;
  d->V = &V;


#ifdef _OPENMP			/* Work is divided between the cores */
  int chunksize = d->l/numthreads;
  if(chunksize > 10) chunksize = 10;
#endif

  /* Initial tunning: step size, and scanning issues */
  Intro(d,&Nscan,&time_correction);
  t_max = (int)((float)d->TT/d->dt);

  /**********************************/
  /* New simulations can start here */
  /**********************************/
  do {    
    /* InitialState_FR(FR,*d,4); */
#pragma omp parallel for schedule(dynamic,chunksize)
    for(i=1 ;i<=d->l ;i++ ) {
      R[i]=(d->J0+sqrt(d->J0*d->J0+4.*M_PI*M_PI*d->eta))/(2.*M_PI*M_PI)+0*1E-6;
      R[i]+=-1E-4*cos(i*d->dx);
      R2[i] = R[i];
      V[i] = 0.0;
    }
    /* printf("\nR[1] = %lf\td->R[0] = %lf\n",R[50],(*d->R)[50]); */
    /* File names */
    Data_Files(d);
    
    /* Debbug  and openning files*/
    if(def == 0) DEBUG(DataDebug(*d,&file[0]));
    FileT = LoadFileLibrary(d->ftime, d->tmodes);
    FileT.multiopen(&file,&FileT);

    /* Reset counters */
    t = 0;
    TIME = 0.0;
    /* Run */
    do {			/* Tiempo */
      if(t%(t_max/10) == 0) { /* Control point */
	sprintf(mesg,"%d%% ",(int)(t*100.0/t_max));
	DEBUG(mesg);
      }
#pragma omp parallel for private(j) schedule(dynamic,chunksize)
      for(i=1 ;i<=d->l ;i++ ) { 	/* Espacio */
	N_F(J,R2,i,*d);
	/* We integrate the ODEs */
	R[i] += d->dt*(d->Deta/M_PI + 2.*R[i]*V[i]);
	V[i] += d->dt*(d->eta + pow(V[i],2) - pow(R[i]*M_PI,2) + J[i]);
      }
      /* Results are stored */
      for(i=1 ;i<d->l ;i++ ) {
	R2[i] = R[i];
	fprintf(file[1] ,"%lf ",R[i]);
	fprintf(file[2] ,"%lf ",V[i]);
      }
      fprintf(file[1] ,"\n");
      fprintf(file[2] ,"\n");
      
      t++;
      TIME = t*d->dt;
    } while(TIME < d->TT);	/* Paso temporal (se puede hacer de la misma manera que en el QIF */

    /* Final shape is stored */
    for(i=1 ;i<=d->l ;i++ ) 
      fprintf(file[0] ,"%lf\t%lf\t%lf\t%lf\n",-M_PI + i*d->dx,R[i],V[i],J[i]);
    
    sprintf(mesg,"%d%% ",(int)(t*100.0/t_max));
    DEBUG(mesg);
    FileT.closeall(&file,&FileT);
  } while (d->scan < Nscan);  /* Simulation ends here */


  free_dvector(R,1,d->l);
  free_dvector(R2,1,d->l);
  free_dvector(V,1,d->l);
  free_dvector(J,1,d->l);


  free(d);
  system("rm -r ./temp");
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
      p[i] = center + gamma*tan((M_PI/2.0)*k);
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
   * But th goes from 0 to M_PI 
   * v goes from v_reset to v_peak (finites) */

  switch(type) {
  case 0:			/* Uniform distribution between reset and peak values */
    for(i=0 ;i < d.N ;i++ ) 
      th[i].v = d.vr + (d.vp - d.vr)*gsl_rng_uniform(r);    
    break;
  case 1:			/* Null distribution: constant 0 valued. */
    for(i=0 ;i < d.N ;i++ ) 
      th[i].v = 0.0;
    break;
  case 2:			/* Lorentzian distribution centered at 10 and width of 5 */
    h = 5;			/* (take into account that the distribution is cut at vr */
    v= 10;			/*  and vp) */
    for(i=0 ;i < d.N ;i++ ) {
      k = (2.0*(i+1) -d.N -1.0)/(d.N+1.0);
      th[i].v = v + h*tan((M_PI/2.0)*k);
      /* if(fabs(th[i].v) > d.vp) th[i].v = d.vr + (d.vp - d.vr)*gsl_rng_uniform(r); /\* If we want to have ou values between vr and vp *\/ */
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
  double norm = (1.0*M_PI/d.N)*delta;

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

/* === FUNCTION  Var_Update ====================
 * Description:  Updates target variable in scan mode
 *   Variables:  Data structure
 * ======================================= */

t_data *Var_update(t_data *d) {
  t_data *d2;
  d2 = malloc(sizeof(*d2));
  char *var_name[4] = {"J","eta"};
  *d2 = *d;
  switch(d->variable) {
  case 1:			/* J */
    if(d->scan == 0)
      d2->J0 = d->min;
    else
      d2->J0 = d2->J0 + (d->step);
    sprintf(mesg,"Variable [%s] changed from: %lf to %lf.",
	    var_name[d->variable-1],d->J0,d2->J0);
    d2->var_value = d2->J0;
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
  char *var_name[4] = {"J","eta"};
  double scanstep = 0;


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

  if(d->N*d->TT > d->MaxDim)
    d->disable_raster = 1;
  else
    d->disable_raster = 0;

  /* Space step size (dx) */
  d->dx = 2.0*M_PI/d->l;
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
	  "\t Value of eta (e)        = %5.2lf\n"
	  "\t Value of Deta (d)       = %5.2lf\n"
	  "\t Value of J0 (a)         = %5.2lf\n"
	  "\t Value of J1 (b)         = %5.2lf\n"
	  "\t Value of  l (l)         = %5d\n"
	  "\t Total number of Neurons = %5d\n"
	  "\t Simulation Time         = %5.2lf\n"
	  "\t Time step               = %5.4lf\n", d.vr, d.vp,
	  d.eta,d.Deta,d.J0,d.J1,d.l, d.N, d.TT,d.dt);
  fprintf(*file,mesg);
  CloseFile(file);
  return mesg;
}

void Data_Files(t_data *d) {
  /* File paths */
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
  /* (d->ftime)[5] = d->file_vt_QIF; */
  /* (d->tmodes)[5] = "w"; */

  if(d->scan_mode > 0) {    
    sprintf(d->file_rvp_FR,"./results/%s/rvp_FR_%.lf.txt",d->file,d->var_value);
    sprintf(d->file_rvp_QIF,"./results/%s/rvp_QIF_%.lf.txt",d->file,d->var_value); 
  }
}

