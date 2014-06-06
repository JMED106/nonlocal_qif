/*********************************************/
/* Simulation of the non-localized FR model. */
/*********************************************/

/*************************/
/* Version 2.0-qif-beta. */
/*************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <limits.h>
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


int **spike;			/* Vector that contains number of spikes of each pop */
int *spike2;
unsigned long int  **Spikes;    /* Matrix that contains spikes of each neuron in each pop */

/* Additional libraries (custom libraries) */
#include "common.h"		/* This must be first */
#include "utils.h"
#include "file.h"
#include "nrutil.h"
#include "nr.h"
#include "field.h"
#include "qif.h"

/* Functions declaration. */
double *InitCond(double *p, int N, int distr_type, double center, double gamma);
double *InitialState(double *p, t_data d, int distr_type, double center, double gamma);
double  MeanField(t_data d,int pop);
double at(double t, t_data d);
t_data *Var_update(t_data *d);
void    Intro(t_data *d, int *Nscan, int *tr_TT);
char   *DataDebug(t_data d,FILE **file);
void    Data_Files(t_data *d);
void FiringRate(int *count, int *Tspikes, t_data d, FILE **file);

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

  double *rqif, *rqif2, *vqif;
  double *FR;
  double **QIF;
  double **qif_eta, *Jqif;
    
  double *pt,*pt2;	      /* Auxiliar pointer for initial dist. */

  int *count;
  
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

  QIF = dmatrix(1,d->l,1,d->N);	    /* QIF neurons matrix */
  qif_eta = dmatrix(1,d->l,1,d->N); /* QIF parameter matrix */
  rqif = dvector(1,d->l);
  FR = dvector(1,d->l);
  rqif2 = dvector(1,d->l);
  vqif = dvector(1,d->l);
  Jqif = dvector(1,d->l);

  pt = dvector(1,d->N);		/* Allocation of the auxiliary pointer */
  pt2 = dvector(1,d->N);	/* Allocation of the auxiliary pointer */

  count = ivector(1,d->l);  

  d->R = &R;
  d->V = &V;

#ifdef _OPENMP			/* Work is divided between the cores */
  int chunksize = d->l/numthreads;
  if(chunksize > 5) chunksize = 5;
#endif

  /* Initial tunning: step size, and scanning issues */
  Intro(d,&Nscan,&time_correction);
  t_max = (int)((float)d->TT/d->dt);

  /* Matrix for time correction */
  Spikes = ulimatrix(1,d->l,1,d->N);
  spike = imatrix(0,t_max,1,d->l);
  spike2 = ivector(1,d->l);
  /**********************************/
  /* New simulations can start here */
  /**********************************/
  do {        
#pragma omp parallel for schedule(dynamic,chunksize)
    for(i=1 ;i<=d->l ;i++ ) {
      R[i]=(d->J0+sqrt(d->J0*d->J0+4.*M_PI*M_PI*d->eta))/(2.*M_PI*M_PI)+0*1E-6;
      R[i]+=-1E-4*cos(i*d->dx);
      rqif[i] = R[i];
      rqif2[i] = R[i];
      R2[i] = R[i];
      V[i] = 0.0;
      vqif[i] = 0.0;
      count[i] = 0;
      spike2[i] = 0;
      for(j=1 ;j<=d->N ;j++ ) 
	Spikes[i][j] = 0;
    }
    pt2 = InitCond(pt2,d->N,2,d->eta,d->Deta);      
    /* Initial state of the QIF neurons */
    for(i=1 ;i<=d->l ;i++ ) {
      /* Each position has a given R, which is the width of .. */
      /* .. the Lorentzian distribution */
      pt = InitialState(pt,*d,2,V[i],R[i]);      
      for(j=1 ;j<=d->N ;j++ ) {
	QIF[i][j] = pt[j];
	qif_eta[i][j] = pt2[j];
      }
    }
    FILE *hola;
    hola = fopen("inicial.txt","w");
    for(i=1 ;i<=d->N ;i++ ) {
      fprintf(hola ,"%lf\t%lf\n ",QIF[1][i],qif_eta[1][i]);
    }
    fclose(hola);
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
    d->t = 0;
    /* Run */
    do {			/* Tiempo */
      if(t%(t_max/10) == 0) { /* Control point */
	sprintf(mesg,"%d%% ",(int)(t*100.0/t_max));
	DEBUG(mesg);
      }
#pragma omp parallel for private(j) schedule(dynamic,chunksize)
      for(i=1 ;i<=d->l ;i++ ) { 	/* Espacio */
	spike[t][i] = 0;
	N_F(J,R2,i,*d);
	N_F(Jqif,rqif2,i,*d);
	/* We integrate the ODEs */
	R[i] += d->dt*(d->Deta/M_PI + 2.*R[i]*V[i]);
	V[i] += d->dt*(d->eta + pow(V[i],2) - pow(R[i]*M_PI,2) + J[i]);
	/* QIF integration */
	F_QIF(QIF[i],qif_eta[i],Jqif[i],i,*d);
	/* r and v are computed */
	rqif[i] = MeanField(*d,i);
	count[i]++;
	/* v here */
      }
      /* Results are stored */
      for(i=1 ;i<=d->l ;i++ ) {
	R2[i] = R[i];
	rqif2[i] = rqif[i];
	fprintf(file[1] ,"%lf ",R[i]);
	fprintf(file[2] ,"%lf ",V[i]);

	/* fprintf(file[4] ,"%lf ",rqif[i]); */

	/* Firing rate of QIF */
	if(t%10 == 0) 
	  FiringRate(&count[i], &spike2[i], *d, &file[4]);
      }
      fprintf(file[1] ,"\n");
      fprintf(file[2] ,"\n");
      if(t%10 == 0) 
	fprintf(file[4] ,"\n");

      d->t++;
      t++;
      TIME = t*d->dt;
    } while(TIME < d->TT);	/* Paso temporal (se puede hacer de la misma manera que en el QIF */

    /* Final shape is stored */
    for(i=1 ;i<=d->l ;i++ ) 
      fprintf(file[0] ,"%lf\t%lf\t%lf\t%lf\n",-M_PI + i*d->dx,R[i],V[i],J[i]);

    /* Final shape is stored */
    for(i=1 ;i<=d->l ;i++ ) 
      fprintf(file[3] ,"%lf\t%lf\t%lf\t%lf\n",-M_PI + i*d->dx,R[i],V[i],Jqif[i]);
    
    sprintf(mesg,"%d%% ",(int)(t*100.0/t_max));
    DEBUG(mesg);
    FileT.closeall(&file,&FileT);
    free_dvector(pt,1,d->N);
    free_dvector(pt2,1,d->N);
  } while (d->scan < Nscan);  /* Simulation ends here */


  free_dvector(R,1,d->l);
  free_dvector(R2,1,d->l);
  free_dvector(V,1,d->l);
  free_dvector(J,1,d->l);
  free_dmatrix(qif_eta,1,d->l,1,d->N);
  free_dmatrix(QIF,1,d->l,1,d->N);

  free_dvector(rqif,1,d->l);
  free_dvector(FR,1,d->l);
  free_dvector(rqif2,1,d->l);
  free_dvector(vqif,1,d->l);
  free_dvector(Jqif,1,d->l);

  free_ivector(count,1,d->l);
  free_ulimatrix(Spikes,1,d->l,1,d->N);
  free_imatrix(spike,0,t_max,1,d->l);
  free_ivector(spike2,1,d->l);


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
    for(i=1 ;i<=N ;i++ ) 
      p[i] = center + gsl_ran_cauchy(r, gamma);
    break;
  case 1:			/* Gaussian (random) */
    for(i=1 ;i<=N ;i++ ) 
      p[i] = center + gsl_ran_gaussian(r,gamma);
    break;
  case 2:			/* Cauchy ordered */
    for(i=1 ;i<=N ;i++ ) {
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

double *InitialState(double *p, t_data d, int distr_type, double center, double gamma) {
  int i;
  double k = 0;
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

  switch(distr_type) {
  case 0:			/* Uniform distribution between reset and peak values */
    for(i=1 ;i <= d.N ;i++ ) 
      p[i] = d.vr + (d.vp - d.vr)*gsl_rng_uniform(r);    
    break;
  case 1:			/* Null distribution: constant 0 valued. */
    for(i=1 ;i <=d.N ;i++ ) 
      p[i] = 0.0;
    break;
  case 2:			/* Lorentzian distribution*/
    for(i=1 ;i <= d.N ;i++ ) {
      k = (2.0*(i+1) -d.N -1.0)/(d.N+1.0);
      p[i] = center + gamma*tan((M_PI/2.0)*k);
      if(fabs(p[i]) > d.vp) p[i] = 0.0; /* If we want to have ou values between vr and vp */
    }
    break;
  default:
    ARG_ERROR;
    break;
  }
  return p;
}


/* === FUNCTION  MeanField ====================
 * Description:  Computes S
 *   Variables:  all th_i
 * ======================================= */

double MeanField(t_data d,int pop) {
  int i;
  double s = 0;
  double delta = 1.0/d.dt;
  double norm = (1.0/d.N)*delta;
  double norm2 = (1.0*M_PI/d.N)*delta;
  double t_min = 0;
  spike[d.t][pop] = 0;

  for(i=1 ;i<=d.N ;i++ ) { 
    if(Spikes[pop][i] == d.t_spike) 
      s+= 1;
  }
  spike2[pop] += s;
  spike[d.t][pop] = s;

  /* s = 0; */
  /* if(d.t_min >= d.t) t_min = d.t; */
  /* else t_min = d.t_min; */
  /* for(i=d.t-t_min ;i<=d.t ;i++ )  */
  /*   s+=spike[i][pop]; */

  /* s+=at((d.t-i)*d.dt,d)*spike[i-1][pop]; */
  s = s*norm;	
  return s;
}

/* === FUNCTION  at ====================
 * Description:  Exponential function for s(t)
 *   Variables:  time and tau
 * ======================================= */

double at(double t, t_data d) {
  return (1.0/d.tau)*exp(-t/d.tau);
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
  double min_dt = 0.0;


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
  /** ***** ***********************************************/

  scanstep = d->step;
  *Nscan = (int)ceil(1+(fabs(d->max - d->min)/((float)d->step)));
  if(d->max < d->min) d->step = (-1.0)*scanstep;
  if(d->scan_mode == 0) *Nscan = 0;

  if(d->N*d->TT > d->MaxDim)
    d->disable_raster = 1;
  else
    d->disable_raster = 0;

  /* Time step, t_spike ... */
  d->tau_p = 1.0/d->vp;
  d->tau_r = 1.0/d->vr;
  d->wait = 2*(int)(ceil(d->tau_p/d->dt));
  sprintf(mesg,"Refractory stepts adjusted to: %d",d->wait);
  DEBUG(mesg);
  d->t_spike = 1UL<<(d->wait/2);
  /* Space step size (dx) */
  d->dx = 2.0*M_PI/d->l;

  /* d.tau from at */
  d->tau = 10*d->dt;
  d->t_min = 100;
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

void FiringRate(int *count, int *Tspikes, t_data d, FILE **file) {
  double inst_fr_avg = 0;

  inst_fr_avg  = (double)*Tspikes/(double)*count;
  inst_fr_avg  = inst_fr_avg/(d.N*d.dt);
  *count = 0;
  fprintf(*file,"%lf ",inst_fr_avg);
  inst_fr_avg = 0;
  *Tspikes = 0;
}
