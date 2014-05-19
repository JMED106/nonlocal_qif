#ifndef __COMMON__
#define __COMMON__

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

  double g;
  double E;		        /*  */
  double DE;		        /*  */
  int E_dist;			/*  */

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

  int pert;
  double pert_amplitude;        /* P: Amplitude of the perturbation QIF */
  double perturbation_FR;	/* FR *** Amplitude of the perturbation FR*/

  /* Related to the non-localized system */
  int l;			/* l: Number of clusters (columns) */
  double J0,J1,J2;		/* a,b,c: Amplitudes of the modes of the coupling function */
  int ring;			/* R: Boolean to determine periodic boundaries */
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
  struct  t_FR **FR;		/* We will use this object to carry info about r and v */
  struct T_qif ***QIF;
  
  double rx;
  double vx;
  
  double S;
  int t;			/* Time counter */
} t_data;

typedef struct {
  double r;
  struct  t_FR **FR;
} str_pr;

struct T_qif {
  double v;
  double x;

  double eta;

  double r;			/* Firing Rate saved in each oscillator */
  double spike;			/* Has the oscillator spiked? */
  int tr;			/* How many steps ago has the oscillator spiked? */
  int total_spikes;	        /* Total number of spike of this oscillator */
  int global_s1;

  int N;			/* Total number of oscillators (redundant) */
  double FR;
  int pert;

  /* For g */
  double E;
  double g;
  double J;

  /* Provisionales (a borrar) */
  double rh,rh2,vh,vh2;
};

struct t_FR {
  double r, r2,*rp;
  double v, v2;
  double S;			/* Suma de Riemann (integral) */
  int x;
} ;

typedef struct t_FR T_FR;
typedef struct T_qif t_qif;
  


#ifndef GSL_SIGN
#define GSL_SIGN(x) ((x) >= 0 ? 1: -1)
#endif
#endif
