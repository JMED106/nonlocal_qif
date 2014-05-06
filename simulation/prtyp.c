/**********************************/
/* Simulation of N neuron system. */
/**********************************/

/*********************************/
/* Features:			 */
/* 				 */
/* - Simulation of QIF		 */
/* - Simulation of Theta neurons */
/* - Simulation of FR equations	 */
/*********************************/


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
#include "winfree.h"
#include "theta.h"
#include "common.h"


/* Functions declaration. */
void    Intro(t_data *d, int *Nscan, int *tr_TT);
double *InitCond(double *p , int N, int dist_type,  double center, double gamma);
void    InitialState(t_th * th, int type);
char   *DataDebug(t_data d, FILE *(*files)[]);
double  MeanField(t_th *th, double dt, int type);
double  OrderParameter(t_th *th,double *phi);
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

  time_t t;   struct tm *tm;  t = time(NULL); tm = localtime(&t);
  if((argc >1) && (argv[1][1] == 'd')) {debug = 1; d_level = 1;} /* Initial debugging */

  srand(time(NULL));
  gsl_rng_env_setup();		/* Random generator initialization */
  gsl_rng_default_seed = rand()*RAND_MAX;
  sprintf(mesg,"Seed: %ld ",gsl_rng_default_seed);
  DEBUG(mesg);
  
  /* Declaration of variables */
  /* Contadores */
  int i, time, tr_TT = 0;
  int Nscan = 0;
  /* Others */
  int def = 1;			/* Default = 1 */
  int temporal_correction = 1;
  int total_spikes1 = 0, total_spikes2 = 0;

  double var_value;
  double *final_conf,  *final_conf2;

  /* Date variables */
  char day[100], hour[100];
  strftime(day, 100, "%m-%d-%Y",tm);
  sprintf(hour,"%d.%02d",tm->tm_hour,tm->tm_min);

  /* Data store */
  Create_Dir("results");
  FILE *file[8];
  
  /* System variables (magnitudes) */
  t_th *th;			/* neuron vector */
  double R = 0;			/* Kuramoto order parameter */

  double fr_volt[2];		/* Stores width and center of the Lorentzian distr. (R) */
  double fr_x;			/* Center of the Lorentzian distr. */
  double fr_y;			/* Width  " " " " */
  double perturbation = 0;

  double fr, fr2;		/* Mean field (firing rate) */
  double avg_fr1,avg2_fr1;
  int medida;
  double inst_fr_avg;
  int intervalo;
  int total_spikes_inst;

  double avg_v1 = 0, avg_v2 = 0;
  long double v1_center = 0, v2_center = 0;

  double dt, *p, phi = 0;

  double rh_final = 1, vh_final = 1; /* FR *** initial conds. for r and v ( FR equations) */

  t_data *d, *data;
  d = malloc(sizeof(*d));
  data = malloc(sizeof(*data));

  d->scan = 0;			/* Exploring is unabled by default */
  d->MaxDim = 5000000;		/* Raster plot will be disabled for higher dim */
  
  sprintf(d->file,"%s_%s",day,hour);
  sprintf(mesg2,"results/%s",d->file);
  Create_Dir(mesg2);
  Create_Dir("temp");
  sprintf(mesg2,"cp volt_avg.R ./temp");
  system(mesg2);

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

#ifdef _OPENMP			/* Work is divided between the cores */
  int chunksize = d->N/numthreads;
  if(chunksize > 10) chunksize = 10;
#endif

  /* Initial tunning: step size, and scanning issues */
  Intro(d,&Nscan,&tr_TT);
  if(d_level == 4){ temporal_correction = 0; DEBUG("Temporal correction is OFF");}
  d->var_value = d->eta;
  d->scan_max = Nscan;
  d->pert = 0;

  if(d->scan_mode == 3) {
    d->pert = 1;
    d->scan_mode = 2;
  }
  double dTT = d->TT;
  double dTT2 = d->TT*2.0;

  /* Configuration at the end of simulation */
  final_conf  = (double*) malloc (d->N*sizeof(double));
  final_conf2 = (double*) malloc (d->N*sizeof(double));

  /**********************************/
  /* New simulations can start here */
  /**********************************/
  do {    
    if(d->scan_mode >= 1) 
      d = Var_update(d);
    Data_Files(&file,*d,0);
    Data_Files(&file,*d,4);
    total_spikes1 = 0;
    total_spikes2 = 0;
    
    /********************/
    /* Oscillator setup */
    /********************/
    th  = (t_th*) calloc (d->N,sizeof(t_th)); 
    for(i=0 ;i<d->N ;i++ ) {		/* The total number (N) is stored in each package */
      th[i].N = d->N;
    }
    th[0].pert = 0;
    th[0].rh = 1;
    th[0].vh = 1;                       /* FR **** Initial condition for the FR equations */

    p = (double*) malloc ((d->N+1)*sizeof(double));

    /* QIF: vr and vp and g */
    for(i=0 ;i<d->N ;i++ ) {
      th[i].vr = d->vr;
      th[i].vp = d->vp;
      th[i].g = d->g;
    }
    
    
    if(d->scan_mode == 2 && d->scan > 0) {
      for(i=0 ;i<d->N ;i++ ) {
	th[i].v = final_conf[i];                /* We mantain the last configuration of voltages */
	th[i].th = final_conf2[i];
      }
    } else {
      InitialState(th,d->init_dist);		/* Initial state configuration */
      /* InitialState(th,d->init_dist+1); */
    }

    for(i=0 ;i<d->N ;i++ ) {
      th[i].spike = 0;		                /* Spike */
      th[i].tr = 0;
      th[i].tr2 = 0;
      th[i].ph = VarChange(th[i]);
    }
    
    dt = d->dt;

    if(def == 0)		/* Debug message: data display */
      DEBUG2(DataDebug(*d,&file));

    if(d->scan_mode == 2 && d->scan > 0) {
      th[0].rh = rh_final;	/* FR **** final configuration of the pevious ... */
      th[0].vh = vh_final;	/* ... simulation is taken as the initial configuration */
    }
    
    /** Distributions *******************/   

    /* QIF: J */
    sprintf(mesg,"Building distribution for J_i. J0 = %lf ",d->J );
    if(d->J_dist == 1)  DEBUG(mesg);
    p = InitCond(p,d->N,2,d->J,d->J_sigma);
    if(d->J_dist == 0) {
      for(i=0 ;i<d->N ;i++ ) 
	th[i].J = d->J;
    } else
      for(i=0 ;i<d->N ;i++ ) 
	th[i].J = p[i];
  
    /* QIF: eta */
    sprintf(mesg,"Building distribution for eta_i. eta0 = %lf ",d->eta );
    if(d->eta_dist == 1)  DEBUG(mesg);
    p = InitCond(p,d->N,2,d->eta,d->eta_sigma);
    if(d->eta_dist == 0) {
      for(i=0 ;i<d->N ;i++ ) 
	th[i].eta = d->eta;
    } else {
      for(i=0 ;i<d->N ;i++ ) 
	th[i].eta = p[i];
    }

    /* QIF: V0 */
    sprintf(mesg,"Building distribution for V0_i. V0 = %lf ",d->v0 );
    if(d->v0_dist == 1)  DEBUG(mesg);
    p = InitCond(p,d->N,2,d->v0,d->v0_sigma);
    if(d->v0_dist == 0) {
      for(i=0 ;i<d->N ;i++ ) 
	th[i].V0 = d->v0;
    } else
      for(i=0 ;i<d->N ;i++ ) 
	th[i].V0 = p[i];


    /* Simulation */
    sprintf(mesg,"Simulating dynamics of: v' = v² + I ");    DEBUG(mesg);
    sprintf(mesg,"I = J*r + n - g*r(v - v0)");    DEBUG(mesg);
    sprintf(mesg,"I = %.3lf*r + %.3lf - %.3lf*r(v - %.3lf) ",d->J,d->eta,d->g,d->v0);    DEBUG(mesg);

    time = 0; fr = 0; fr2 = 0;  avg_fr1 = 0; 
    avg2_fr1 = 0; medida = 0;
    v1_center = 0; v2_center = 0; 
    d->voltdist = 0;
    intervalo = 0;
    inst_fr_avg = 0;
    total_spikes_inst = 0;

    do {
      /* In each step we must compute the mean field potential (r) */
      if(time%((int)((float)d->TT/(10.0*dt))) == 0) { /* Control point */
	sprintf(mesg,"%d%% ",(int)(((time*dt)/(float)d->TT)*100) );
	DEBUG(mesg);
      }
      /* Parallelizable, watch out with pointers and shared variables */
#pragma omp parallel for schedule(dynamic,chunksize)
      for(i=0 ;i<d->N ;i++ ) {	/* i represents each oscillator */
	th[i].r = fr;	
	if(th[0].pert == 1) 
	  th[i].r = th[0].FR;
	th[i].r2 = fr2;
	th[i].spike = 0;
	th[i].spike2 = 0;

	if(temporal_correction == 1) {
	  if(th[i].tr == 1) {
	    th[i].tr = 2;
	    th[i].spike = 1;
	  } else if(th[i].tr == 2) {
	    th[i].tr = 0;
	  } else 
	    th[i] = QIF(*d,th[i]);
	} else {
	  th[i] = QIF(*d, th[i]);
	  if(th[i].tr == 1) th[i].spike = 1;
	  th[i].tr = 0;
	}

	th[i] = THETA(*d,th[i]);

	if(time*dt >= (d->TT/4.0)*3.0) 
	  th[i].total_spikes += th[i].spike;
      }
      /* Simulation of rate equations */
      th[0] = Theory(*d,th[0]);
      /* Ends here */

      fr = MeanField(th,d->dt,1);   fr2 = MeanField(th,d->dt,2);
      if((time*d->dt)/d->TT > 0.9) {
	total_spikes1 += th[0].global_s1;
	total_spikes2 += th[0].global_s2;
	for(i=0 ;i<d->N ;i++ ) {
	  v1_center += th[i].v;
	  v2_center += th[i].th;
	}

	v1_center /= d->N;
	v2_center /= d->N;
	medida++;
	avg_v1 += v1_center;
	avg_v2 += v2_center;
      }
      /* Inst FR */
      total_spikes_inst += th[0].global_s1;
      intervalo++;

      if(abs(time - (int)(d->TT/d->dt)) <= 50) {
	d->voltdist++;
	Data_Files(&file,*d,2);
	for(i=0 ;i<d->N ;i++ ) 
	  fprintf(file[6],"%lf\n",th[i].v);
	Data_Files(&file,*d,3);
      }

      if(d->disable_raster == 0)
	for(i=0 ;i<d->N ;i++ ) {
	  if(th[i].spike == 1)
	    fprintf(file[2] ,"%lf\t%d\t-1\n",time*d->dt,i);
	  if(th[i].spike2 == 1)
	    fprintf(file[2] ,"%lf\t-1\t%d\n",time*d->dt,i);
	}
      R = OrderParameter(th,&phi);
      th[0].FR = 0;
      th[0].pert = 0;
      /* Perturbation introduced here */
      if((time > (d->TT/d->dt)-100) && d->pert == 1 && d->dx >= 0) {
	if(time == (d->TT/d->dt)-100 +1) {
	  avg_fr1 = (double)total_spikes1/((double)medida*d->dt);
	  avg_fr1 /= d->N;
	  printf("\n Average FR 0: %.8lf",avg_fr1);
	  fflush(stdout);	  
	}
	if(abs(time-(int)((d->TT/d->dt)-100))%5 == 0) {
	  Perturbation(&th,*d,1);
	  th[0].pert = 1;
	  printf("\nPerturbation: %lf!!",d->pert_amplitude);
	  perturbation = d->pert_amplitude;
	}
	if(time == (d->TT/d->dt)-1) {
	  d->TT = dTT2;
	  d->pert = 2;
	}
      } else if((time == (d->TT/d->dt)-1) && d->pert == 1) {
	Perturbation(&th,*d,0);
	th[0].vh += d->perturbation_FR; /* FR *** Perturbation for FR equations */
	th[0].pert = 1;
	printf("\nPerturbation: %lf!!",d->pert_amplitude);
	perturbation = d->pert_amplitude;
	d->TT = dTT2;
	avg_fr1 = (double)total_spikes1/((double)medida*d->dt);
	avg_fr1 /= d->N;
	printf("\n Average FR 0: %.8lf",avg_fr1);
	fflush(stdout);
	d->pert = 2;
      }
      fprintf(file[3] ,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",time*dt,fr,R,phi,fr2,R,phi,th[0].rh,th[0].vh);
      fprintf(file[0] ,"%lf\t%lf\t%lf\t%lf\n",time*d->dt,th[d->N/4].v,th[d->N/4].th,perturbation);

      if(time%50 == 0) {
	inst_fr_avg  = (double)total_spikes_inst/(double)intervalo;
	inst_fr_avg  = inst_fr_avg/(d->N*d->dt);
	intervalo = 0;
	fprintf(file[7],"%lf\t%lf\n",time*d->dt,inst_fr_avg);
	inst_fr_avg = 0;
	total_spikes_inst = 0;
      }
      perturbation = 0;

      time++;    
    } while(time*d->dt< d->TT);

    if(d->pert == 2) {
      d->TT = dTT;
      d->pert = 1;
    }
    avg_fr1 = (double)total_spikes1/((double)medida*d->dt);
    avg_fr1 /= d->N;
    printf("\n Average FR 1: %.8lf",avg_fr1);

    avg2_fr1 = (double)total_spikes2/((double)medida*d->dt);
    avg2_fr1 /= d->N;
    printf("\n Average FR Theta: %.8lf\n",avg2_fr1);

    avg_v1 /= medida;
    avg_v2 /= medida;
    avg_v2 = tan(avg_v2*0.5);

    R_script(*d,avg2_fr1,avg_v2);
    R_calc(*d,&fr_x,&fr_y); /* fr: 0  voltage: 1 */
    fr_volt[0] = fr_x;
    fr_volt[1] = fr_y;
    fr_x = 2*fr_x/(PI);
    fr_y = 2*fr_y;

    printf("\n Average Voltage (v) : %lf\n Average Voltage (Theta) : %lf",avg_v1, avg_v2);
    printf("\n Average Voltage (v) using v distribution: %lf\n Average FR using v dist : %lf\n",fr_volt[1]*2,fr_volt[0]*(2.0/(PI)));

    
    if(d->scan_mode >= 1) {
      switch(d->variable) {
      case 1:
	var_value = d->J;
	break;
      case 2:
	var_value = d->eta;
	break;
      case 3:
	var_value = d->g;
	break;
      case 4:
	var_value = d->v0;
	break;
      }
      fprintf(file[4] ,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",var_value,avg_fr1,avg2_fr1,avg_v1, avg_v2,fr_x,fr_y);
    }
    
    printf("\n");
    free(p);
    for(i=0 ;i<d->N ;i++ ) {
      fprintf(file[1] ,"%d\t%lf\t%lf\n",i,th[i].v,th[i].th);
      final_conf[i]  = th[i].v;
      final_conf2[i] = th[i].th;
    }

    free(th);
    d->scan++;
    Data_Files(&file,*d,1);
    Data_Files(&file,*d,5);
    if(d->scan_mode == 0)
      break;

  } while (d->scan < Nscan);
  /* Simulation ends here */


  free(final_conf);
  free(final_conf2);

  /* system("R -f histogram.R --quiet --slave"); */
  /* Gnuplot_Init(*d,day,hour); */
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
      th[i].th = VarChange( th[i]);
      th[i].ph = th[i].th;
    }
    break;
  case 1:			/* Uniform distribution for theta neurons */
    for(i=0 ;i < th[0].N ;i++ ) {
      th[i].th = -PI + 2*PI*gsl_rng_uniform(r);
    }
    break;
  case 2:			/* Null distribution: constant 0 valued. */
    for(i=0 ;i < th[0].N ;i++ ) {
      th[i].v = 0.0;
      th[i].ph = th[i].th;
    }
    break;
  case 3:			/* Null distribution for theta neurons. Constant 0 valued. */
    for(i=0 ;i < th[0].N ;i++ ) {
      th[i].th = 0.0;
    }
    break;
  case 4:			/* Lorentzian distribution centered at 10 and width of 5 */
    h = 5;			/* (take into account that the distribution is cut at vr */
    v= 10;			/*  and vp) */
    for(i=0 ;i < th[0].N ;i++ ) {
      k = (2.0*(i+1) -th[0].N -1.0)/(th[0].N+1.0);
      th[i].v = v + h*tan((PI/2.0)*k);
      if(fabs(th[i].v) > th[0].vp) th[i].v = v;

      th[i].ph = th[i].th;
    }
    break;
  case 5:			/* Lorentzian distribution fot theta neurons. */
    h = 0.1;
    v= 0;
    for(i=0 ;i < th[0].N ;i++ ) {
      k = (2.0*(i+1) -th[0].N -1.0)/(th[0].N+1.0);
      th[i].th = v + h*tan((PI/2.0)*k);
      if(fabs(th[i].th) > PI) th[i].th = 0;
      th[i].ph = th[i].th;
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

/* === FUNCTION  VarChange ====================
 * Description:  Variable change: v = tg(th/2)
 *   Variables:  th structure
 * ======================================= */

double VarChange(t_th th) {
  return( 2*atan(th.v));
}

/* === FUNCTION  OrderParameter ====================
 * Description:  Calculates r (r*exp(i\phi))
 *   Variables:  oscillator phases
 * ======================================= */

double OrderParameter(t_th *th,double *phi)  {
  int i;
  double order_cos = 0, order_sin = 0;
  for(i=0 ;i<th[0].N ;i++ ) {
    order_cos += cos(th[i].ph);
    order_sin += sin(th[i].ph);
  }
    
  *phi = atan(order_sin/order_cos);
  return sqrt((1.0/th[0].N)*(1.0/th[0].N)*((order_cos*order_cos) + (order_sin*order_sin)));
  /* modulo del número complejo */
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
