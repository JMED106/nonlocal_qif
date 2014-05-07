#ifndef __COMMON__
#define __COMMON__

typedef struct {
  double v;

  double eta;

  double r;			/* Firing Rate saved in each oscillator */
  double spike;			/* Has the oscillator spiked? */
  int tr;			/* How many steps ago has the oscillator spiked? */
  int total_spikes;	        /* Total number of spike of this oscillator */
  int global_s1;

  int N;			/* Total number of oscillators (redundant) */
  double FR;
  int pert;
} t_qif;

typedef struct {
  int init_dist;		/* i: Type of initial distribution */
  double vr;			/* r: Reseting potential */
  double vp;			/* p: Peak potential */

  double J;			/*  */
  double DJ;			/* Dispersion of the J */
  int J_dist;			/* Type of distribution of J */

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

  double min;	                /* m: Minimum value of the variable to scan */
  double max;			/* M: Maximum value of the variable to scan */
  double step;			/* s: Increment of the value in each simulation */
  int MaxDim;			/* Maximum dimension of the raster plot */
  int disable_raster;

  char file[200];

  int voltdist;			/* Voltage distribution integer */

  int pert;
  double pert_amplitude;        /* P: Amplitude of the perturbation QIF */
  double perturbation_FR;	/* FR *** Amplitude of the perturbation FR*/

  /* Related to the non-localized system */
  int l;			/* l: Number of clusters (columns) */
  double J0,J1,J2;		/* a,b,c: Amplitudes of the modes of the coupling function */
  int ring;			/* R: Boolean to determine periodic boundaries */
} t_data;


typedef struct {
  double r, r2;
  double v, v2;
} t_FR;


#ifndef GSL_SIGN
#define GSL_SIGN(x) ((x) >= 0 ? 1: -1)
#endif

double VarChange(t_th th);

#endif
