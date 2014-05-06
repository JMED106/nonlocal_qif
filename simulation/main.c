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

#define PI 4*atan(1)

/* Global variables declarationx */
int debug = 0,			/* Boolean: debugging activator */
  d_level = 0;			/* Debugging level */
char DATA_FILE[50] = "input_data.conf";
char mesg[1024];
char mesg2[1024];

/* Additional libraries (custom libraries) */
#include "utils.h"
#include "qif.h"
#include "theory.h"
#include "common.h"


/* Functions declaration. */
void    Intro(t_data *d, int *Nscan, int *tr_TT);
double *InitCond(double *p , int N, int dist_type,  double center, double gamma);
void    InitialState(t_th * th, int type);
char   *DataDebug(t_data d, FILE *(*files)[]);
double  MeanField(t_th *th, double dt, int type);
t_data *Var_update(t_data *d);
char   *File_exists(char file_name[], double prmts);
void    Data_Files(FILE *(*files)[], t_data d, int action);
void    R_calc(t_data d,double *fr, double *voltage);
void    R_script(t_data d, double x, double y);
void    Perturbation(t_th **th, t_data d,int type);

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

  time_t t;   struct tm *tm;  t = time(NULL); tm = localtime(&t);
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
  FILE *file[8];

  /* +++++++++++++++++ Simulation variables ++++++++++++++++ */

  /* Parameters */
  t_data *d;

  /* Results Variables */
  int Nscan;

  /* Dynamic variables (system variables) */
  
  /* ++++++++++++++++++++++++++++++++++++++++++++++++++ */
  
  /* Creating data files */
  sprintf(d->file,"%s_%s",day,hour);
  sprintf(mesg2,"results/%s",d->file);
  Create_Dir(mesg2);
  Create_Dir("temp");
  sprintf(mesg2,"cp volt_avg.R ./temp");
  system(mesg2);

  /*********************************/
  /* Program arguments assignation */
  /*********************************/
  /* *data = Scan_Data(DATA_FILE,*d); */
  /* d = data; */

  /* while((argc--) != 1) {	/\* Terminal arguments can be handled *\/ */
  /*   *data = Arg(argv[argc], *d); */
  /*   d = data; */
  /*   if(argv[1][0] != '-') def = 0; */
  /* } */

#ifdef _OPENMP			/* Work is divided between the cores */
  int chunksize = d->N/numthreads;
  if(chunksize > 10) chunksize = 10;
#endif

  /* Initial tunning: step size, and scanning issues */


  /**********************************/
  /* New simulations can start here */
  /**********************************/
  do {    

  } while (d->scan < Nscan);  /* Simulation ends here */

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

void InitialState(t_th *th, int type) {
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
  sprintf(mesg,"th[0].vp = %lf\tth[0].vr = %lf ",th[0].vp,th[0].vr );
  DEBUG(mesg);

  switch(type) {
  case 0:			/* Uniform distribution between reset and peak values */
    for(i=0 ;i < th[0].N ;i++ ) {
      th[i].v = th[0].vr + (th[0].vp - th[0].vr)*gsl_rng_uniform(r);
    }
    break;
  case 2:			/* Null distribution: constant 0 valued. */
    for(i=0 ;i < th[0].N ;i++ ) {
      th[i].v = 0.0;
    }
    break;
  case 4:			/* Lorentzian distribution centered at 10 and width of 5 */
    h = 5;			/* (take into account that the distribution is cut at vr */
    v= 10;			/*  and vp) */
    for(i=0 ;i < th[0].N ;i++ ) {
      k = (2.0*(i+1) -th[0].N -1.0)/(th[0].N+1.0);
      th[i].v = v + h*tan((PI/2.0)*k);
      if(fabs(th[i].v) > th[0].vp) th[i].v = v;
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

double MeanField(t_th *th, double dt,int type) {
  int i;
  double s = 0;
  double delta = 1.0/dt;
  double norm = (1.0*PI/th[0].N)*delta;

  if(type == 1) {
    th[0].global_s1 = 0;
    for(i=0 ;i<th[0].N ;i++ ) {
      if(th[i].spike == 1)
	s+= 1;
    }
    th[0].global_s1 = s;
  } else if (type == 2) {
    th[0].global_s2 = 0;
    for(i=0 ;i<th[0].N ;i++ ) {
      if(th[i].spike2 == 1)
	s+= 1;
    }
    th[0].global_s2 = s;
  }

  s = s*norm;		/* Something wrong with couplinh (MF is not OK) */
  return s;
}

/* === FUNCTION  File ====================
 * Description:  Check whether a file exists
 *   Variables:  file_name, prmts
 * ======================================= */

 char *File_exists(char file_name[], double prmts) {
   int exist = 0;
   int count = 0;
   char *cmd;
   cmd = (char*) malloc(100*sizeof(char));
   do {
     count++;
     sprintf(cmd,"ls ./%s/%s_sigma-%.4lf_%d.dat",file_name,file_name,prmts,count);
     exist = system(cmd);
   } while(exist == 0);
   sprintf(cmd, "./%s/%s_sigma-%.4lf_%d.dat",file_name,file_name,prmts,count);
   return cmd;
 }


/* === FUNCTION  data ====================
 * Description:  Updates target variable in scan mode
 *   Variables:  Data structure
 * ======================================= */

t_data *Var_update(t_data *d) {
  t_data *d2;
  d2 = malloc(sizeof(*d2));
  char *var_name[4] = {"J","eta","g","v0"};
  *d2 = *d;
  switch(d->variable) {
  case 1:			/* J */
    if(d->scan == 0)
      d2->J = d->min_x;
    else
      d2->J = d2->J + (d->dx);
    sprintf(mesg,"Variable [%s] changed from: %lf to %lf.",
	    var_name[d->variable-1],d->J,d2->J);
    d2->var_value = d2->J;
    break;
  case 2:			/* eta */
    if(d->scan == 0)
      d2->eta = d->min_x;
    else
      d2->eta = d2->eta + (d->dx);
    sprintf(mesg,"Variable [%s] changed from: %lf to %lf.",
	    var_name[d->variable-1],d->eta,d2->eta);
    d2->var_value = d2->eta;
    break;
  case 3:			/* g */
    if(d->scan == 0)
      d2->g = d->min_x;
    else
      d2->g += d->dx;
    sprintf(mesg,"Variable [%s] changed from: %lf to %lf.",
	    var_name[d->variable-1],d->g,d2->g);
    d2->var_value = d2->g;
    break;
  case 4:			/* v0 */
    if(d->scan == 0)
      d2->v0 = d->min_x;
    else
      d2->v0 += d->dx;
    sprintf(mesg,"Variable [%s] changed from: %lf to %lf.",
	    var_name[d->variable-1],d->v0,d2->v0);
    d2->var_value = d2->v0;
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
  char *var_name[4] = {"J","eta","g","v0"};
  double tr, tr1 = 0, tr2 = 0;	/* Refractory period */
  double dt, scandx = 0;


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

  scandx = d->dx;
  *Nscan = (int)ceil(1+(fabs(d->max_x - d->min_x)/((float)d->dx)));
  if(d->max_x < d->min_x) d->dx = (-1.0)*scandx;

  /* Time step setup */
  dt = d->dt;
  /* dt(tr) adjustment */
  tr1 = 1.0/fabs(d->vp);
  tr2 = 1.0/fabs(d->vr);
  tr = (tr1 < tr2)? tr1:tr2;
  
  dt = tr;
  if(dt != d->dt) {
    d->dt = dt;
    sprintf(mesg,"Time step dt adjusted to: dt = %lf ",dt );
    DEBUG(mesg);
  }

  *tr_TT = (int)(tr/dt);
  sprintf(mesg,"Refractory period (in steps): %d ",*tr_TT);
  DEBUG(mesg);

  if(d->N*d->TT > d->MaxDim)
    d->disable_raster = 1;
  else
    d->disable_raster = 0;
}


/* === FUNCTION  DataDebug ====================
 * Description:  Displays selected data
 *   Variables:  data structure
 * ======================================= */

char *DataDebug(t_data d,FILE *(*files)[]) {
  sprintf(mesg,"\nCustom parameters\n"
	  "-----------------\n"
	  "\t Reseting potential (r ) = %5.2lf\n"
	  "\t Peak potential (p)      = %5.2lf\n"
	  "\t Value of J (J)          = %5.2lf\n"
	  "\t Value of DJ (m)         = %5.2lf\n"
	  "\t Value of eta (e)        = %5.2lf\n"
	  "\t Value of Deta (d)       = %5.2lf\n"
	  "\t Value of V0 (v)         = %5.2lf\n"
	  "\t Value of DV0 (B)        = %5.2lf\n"
	  "\t Value of g (g)          = %5.2lf\n"
	  "\t Total number of Neurons = %5d\n"
	  "\t Simulation Time         = %5.2lf\n"
	  "\t Time step               = %5.4lf\n", d.vr, d.vp, d.J,d.J_sigma,d.eta,d.eta_sigma,d.v0,d.v0_sigma,d.g, d.N, d.TT,d.dt);
  fprintf((*files)[5] ,"%s",mesg);
  return mesg;
}

void Data_Files(FILE *(*files)[], t_data d, int action) {
  char cmd[200];
  if(action == 0) {		/* Creating flies */
    /* Dynamics of a given neuron (test) */
    sprintf(mesg,"./results/%s/neuron-%d.dat",d.file,d.N/4);
    (*files)[0] = fopen(mesg,"w");
    fprintf((*files)[0],"# Time (ms)\tV(t) (mV) (QIF)\tV(t) (mV) (Theta)\n");

    /* Final voltage distribution */
    sprintf(mesg,"./results/%s/final_volt_dist_%.2lf.dat",d.file,d.var_value);
    (*files)[1] = fopen(mesg,"w");
    fprintf((*files)[1],"# Neuron (i)\tV_f (QIF)\tV_f (Theta)\n");

    /* Raster plot */
    if(d.disable_raster == 0){     
      sprintf(mesg,"./results/%s/raster_%.2lf.dat",d.file,d.var_value);
      (*files)[2] = fopen(mesg,"w");
      fprintf((*files)[2],"# Time (ms)\tNeuron (i,qif)\tNeuron (i,Theta)\n");
    }

    /* Firing rate (and other magnitudes) */
    sprintf(mesg,"./results/%s/firing_rate_%.2lf.dat",d.file,d.var_value);
    (*files)[3] = fopen(mesg,"w");
    fprintf((*files)[3],"# Time (ms)\tFR (QIF)\tR (QIF)\tPhi (QIF)\tR (Theta)\tR (Theta)\tPhi (Theta)\n");    

    /* Parameters */
    sprintf(mesg,"./results/%s/parameters_%.2lf.txt",d.file,d.var_value);
    (*files)[5] = fopen(mesg,"w");

    /* Scan mode */
    if(d.scan_mode >= 1 && d.scan == 0) {
      switch(d.variable) {
      case 1:
	sprintf(mesg,"./results/%s/fr_vs_J.dat",d.file);
	sprintf(cmd,"# J\t");
	break;
      case 2:
	sprintf(mesg,"./results/%s/fr_vs_eta.dat",d.file);
	sprintf(cmd,"# Eta\t");
	break;
      case 3:
	sprintf(mesg,"./results/%s/fr_vs_g.dat",d.file);
	sprintf(cmd,"# g\t");
	break;
      case 4:
	sprintf(mesg,"./results/%s/fr_vs_V0.dat",d.file);
	sprintf(cmd,"# V0\t");
	break;
      default:
	d.scan_mode = 0;
	break;
      }
      (*files)[4] = fopen(mesg,"w");
      fprintf((*files)[4] ,"%sr (QIF)\tr (theta)\n",cmd);
      sprintf(cmd,"Scan data will be saved in %s", mesg);
      DEBUG(cmd);
    }
  } else if(action == 1) {
    fclose((*files)[0]);
    fclose((*files)[1]);
    fclose((*files)[2]);
    fclose((*files)[3]);
    fclose((*files)[5]);
    if(d.scan_mode >= 1 && d.scan >= d.scan_max)
      fclose((*files)[4]);
  } else if(action == 2) {
    sprintf(mesg,"./temp/volt_dist%d.dat",d.voltdist);
    (*files)[6] = fopen(mesg,"w");
  } else if(action == 3) {
    fclose((*files)[6]);
  } else if(action == 4) {
    sprintf(mesg,"./results/%s/inst_fr_%.2lf.dat",d.file,d.var_value);
    (*files)[7] = fopen(mesg,"w");
  } else if(action == 5) {
    fclose((*files)[7]);
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

void Perturbation(t_th **th, t_data d,int type)  {
  int i;
  if(type == 0) 
    for(i=0 ;i<d.N ;i++ ) {
      ((*th)[i].v) += d.pert_amplitude;
    }
  else if(type == 1)
    (* th)[0].FR = d.pert_amplitude;
}
