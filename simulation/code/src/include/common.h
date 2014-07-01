#ifndef __COMMON__
#define __COMMON__


void WarningFPE(int);


typedef struct {
  int init_dist;		/* i: Type of initial distribution */
  double vr;			/* r: Reseting potential */
  double vp;			/* p: Peak potential */

  double eta;			/* h: External current */
  double Deta;		        /* H: Width of the external current */
  int eta_dist;			/* H ≠ 0 Type of distribution of eta */

  int N;			/* N: Number of neurons */
  double TT;			/* T: Total time of the simulation */
  double dt;			/* t: Time step */
  double dt2;			/* f(t) Auxiliary time step (not used) */

  int scan_mode;		/* S: Scan mode boolean */
  int scan;			/* counter */
  int scan_max;
  int variable;
  double var_value;

  double min;	                /* m: Minimum value of the variable to scan */
  double max;			/* M: Maximum value of the variable to scan */
  double step;			/* s: Increment of the value in each simulation */
  int MaxDim;			/* Maximum dimension of the raster plot */
  int disable_raster;
  char file[200];
  int voltdist;			/* Voltage distribution integer */

  /* Perturbation and external inputs */
  int perturbation_ON;  	/* Pertubation is ON */
  int perturbation_time;	/* Counter */
  int perturbation_max_time;	/* Duration of the perturbation */
  int pert_start;

  double R_perturbation;	/* Amplitude for perturbations at R */
  int pert_typeR;
  double V_perturbation;	/* Amplitude for perturbations at V */
  int pert_typeV;

  /* Related to the non-localized system */
  int l;			/* l: Number of clusters (columns) */
  double J0,J1,J2;		/* a,b,c: Amplitudes of the modes of the coupling function */
  double dx;			/* Space interval */
  double DX;			/* Size of the system */

  /* File names */
  char file_parameters[256];
  char file_rvJ_FR[256];
  char file_rt_FR[256];
  char file_vt_FR[256];
  char file_rvJ_QIF[256];
  char file_rt_QIF[256];
  char file_vt_QIF[256];
  char file_rvp_FR[256];
  char file_rvp_QIF[256];

  char **ftime; char **tmodes;
  char **fscan; char **smodes;
  
  double rx;
  double vx;

  double **R;
  double **V;
  double S;
  int t;			/* Time counter */
  double tau;			/* characteristic time in at */

  double t_min;

  unsigned long int t_spike;	/* Time when the spike occurs */
  double tau_p;			/* Time to reach infinity from vpeak (t_spike) */
  double tau_r;			/* Time to reach vreset from -infinity */
  int wait;			/* Refractory Time steps from vp to vr */

  double sample;
} t_data;

  
#ifndef GSL_SIGN
#define GSL_SIGN(x) ((x) >= 0 ? 1: -1)
#endif
#endif
