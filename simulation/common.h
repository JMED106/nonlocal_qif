#ifndef __COMMON__
#define __COMMON__

typedef struct {
  int id;			/* Id of the oscillator */
  int id_rand;
  double v;
  double ph;
  double th;			/* Phase of the osc. */
  double w;			/* w of the osc */
  double vp;			/* peak potential i */
  double vr;			/* reset potential i */
  double J;
  double eta;
  double g;
  double V0;
  double r;			/* Firing Rate saved in each oscillator */
  double r2;
  double spike;			/* Has the oscillator spiked? */
  double spike2;
  int tr;			/* How many steps ago has the oscillator spiked? */
  int tr2;
  int total_spikes;			/* Total number of spike of this oscillator */
  int global_s1;
  int global_s2;
  int N;			/* Total number of oscillators (redundant) */
  double FR;
  int pert;

  double rh;			/* r and v of the firing rate equations. */
  double vh;
  double rh2, vh2; 		/* Copy of the previous ones */

} t_th;

typedef struct {
  int init_dist;		/* Type of initial distribution */
  double vr;			/* Reseting potential */
  double vp;			/* Peak potential */

  double J;			/* J */
  double J_sigma;		/* Dispersion of the J */
  int J_dist;			/* Type of distribution of J */

  double eta;			/* External current */
  double eta_sigma;		/* Dispersion of the external current */
  int eta_dist;			/* Type of distribution of eta */

  double g;			/* g */
  double g_sigma;		/* Dispersion of the g */
  int g_dist;			/* Type of distribution of g */

  double v0;			/* v0 */
  double v0_sigma;		/* Dispersion of the v0 */
  int v0_dist;			/* Type of distribution of v0 */

  int N;			/* Number of neurons */
  double TT;			/* Total time of the simulation */
  double dt;			/* Time step */
  double dt2;			/* Auxiliary time step (not used) */

  int scan_mode;		/* Scan mode boolean */
  int scan;			/* counter */
  int scan_max;

  int variable;			/* Variable to scan (...?) */
  double var_value;
  double min_x;			/* Minimum value of the variable to scan */
  double max_x;			/* Maximum value of the variable to scan */
  double dx;			/* Increment of the value in each simulation */
  int MaxDim;			/* Maximum dimension of the raster plot */
  int disable_raster;

  char file[200];

  int voltdist;			/* Voltage distribution integer */

  int pert;
  double pert_amplitude;        /* Amplitude of the perturbation QIF */
  double perturbation_FR;	/* FR *** Amplitude of the perturbation FR*/
} t_data;


typedef struct {
  double r;
  double v;
} t_FR;

typedef struct {
  int l;			/* Number of functional columns */
} t_system;

#ifndef GSL_SIGN
#define GSL_SIGN(x) ((x) >= 0 ? 1: -1)
#endif

double VarChange(t_th th);

#endif
