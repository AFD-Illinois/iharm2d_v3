
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <string.h>

/*************************************************************************
      COMPILE-TIME PARAMETERS : 
*************************************************************************/
/** SWITCHES **/

/* include wind source terms */
/* turning wind off in global BH models has bad effect on funnel */
#define WIND        (1)

/* use optically thin cooling */
/* if cooling is on then munit has to be set */
#define COOLING     (0)

/* select reconstruction routine; options are
	reconstruct_lr_lin  - piecewise linear (choice of slope limiters)
	reconstruct_lr_par  - piecewise parabolic (recommended)
	reconstruct_lr_weno - weno5
	reconstruct_lr_mp5  - mp5
*/
#define RECONSTRUCT_LR reconstruct_lr_par

/** END SWITCHES **/

/** dimensions of problem **/
#define N1       (32)		/* number of physical zones in X1-direction */
#define N2       (32)		/* number of physical zones in X2-direction */
#define NMAX     (N1 > N2 ? N1 : N2) /* this sizes 1D slices */

#define NPR        (8) 		/* number of primitive variables */
#define NDIM       (4)		/* number of total dimensions.  Never changes */
#define NPG        (4)		/* number of positions on grid for grid functions */
#define NG         (3)		/* number of ghost zones */

/* A numerical convenience to represent a small non-zero quantity compared to unity:*/
#define SMALL	(1.e-20)

/* Max. value of gamma, the lorentz factor */
#define GAMMAMAX (50.)

/* maximum fractional increase in timestep per timestep */
#define SAFE	(1.3)

/*************************************************************************
    MNEMONICS SECTION 
*************************************************************************/
/* mnemonics for primitive vars; conserved vars */
#define RHO	(0)
#define UU	(1)
#define U1	(2)
#define U2	(3)
#define U3	(4)
#define B1	(5)
#define B2	(6)
#define B3	(7)

/* mnemonics for dimensional indices */
#define TT	(0)
#define RR	(1)
#define TH	(2)
#define PH	(3)

/* mnemonics for centering of grid functions */
#define FACE1	(0)
#define FACE2	(1)
#define CORN	(2)
#define CENT	(3)

/* mnemonics for failure modes */
#define FAIL_UTOPRIM        (1)
#define FAIL_VCHAR_DISCR    (2)
#define FAIL_COEFF_NEG	    (3)
#define FAIL_COEFF_SUP	    (4)
#define FAIL_GAMMA          (5)
#define FAIL_METRIC         (6)

/*************************************************************************
    GLOBAL ARRAY SECTION 
*************************************************************************/
extern double p[N1 + 2*NG][N2 + 2*NG][NPR];	/* space for primitive vars */
extern double dq[N1 + 2*NG][N2 + 2*NG][NPR];	/* slopes */
extern double F1[N1 + 2*NG][N2 + 2*NG][NPR];	/* fluxes */
extern double F2[N1 + 2*NG][N2 + 2*NG][NPR];	/* fluxes */
extern double ph[N1 + 2*NG][N2 + 2*NG][NPR];	/* half-step primitives */
extern int    pflag[N1 + 2*NG][N2 + 2*NG];	/* identifies failure points */

/* for debug & diagnostics */
extern double sum_boundary_flux[2*NDIM][NPR];	/* fluxes of each variable on each boundary */
extern double psave[N1+2*NG][N2+2*NG][NPR];   /* stores old data for time derivatives */
extern double dtsave ;
extern double Jcon[N1+2*NG][N2+2*NG][NDIM];

/* particle-related global variables */
#define NPTOT	0
extern double xp[NPTOT][NDIM];

/* fluid physics parameters */
extern double gam;

/* numerical parameters, counter */
extern double cour;
extern double dx[NPR], startx[NPR];
extern double dt;
extern double t, tf;
extern int nstep;

/* output parameters */
extern double DTd;
extern double DTl;
extern double DTi;
extern double DTp;
extern int DTr;

extern int dump_cnt;
extern int pdump_cnt;
extern int image_cnt;
extern int rdump_cnt;

extern int nstroke;

/* global variables that indicate local position */
extern int icurr, jcurr;

/* geometry data structure */
struct of_geom {
	double gcon[NDIM][NDIM];	/* contravariant (index up) metric */
	double gcov[NDIM][NDIM];	/* covariant (index dn) metric */
	double g;			/* sqrt(-det(g_{\mu\nu})) */
	double alpha;			/* lapse function = (-g_{tt})^{-1/2} */
};

/* state of fluid */
struct of_state {
	double ucon[NDIM];		/* contravariant (index up) four-velocity */
	double ucov[NDIM];		/* covariant (index dn) four-velocity */
	double bcon[NDIM];		/* contravariant (index up) B field four-vector */
	double bcov[NDIM];		/* covariant (index dn) B field four-vector */
};

/* grid functions */
extern double conn[N1 + 2*NG][N2 + 2*NG][NDIM][NDIM][NDIM];
extern struct of_geom ggeom[N1+2*NG][N2+2*NG][NPG] ;

/*************************************************************************
    MACROS
*************************************************************************/
#define START1	(0+NG)
#define START2	(0+NG)

/* loop over all active zones */
#define ZLOOP for(i=0+START1;i<N1+START1;i++)for(j=0+START2;j<N2+START2;j++)

/* specialty loop */
extern int istart, istop, jstart, jstop;
#define ISLOOP(istart,istop) \
        for(i=istart+START1;i<=istop+START1;i++)
#define JSLOOP(jstart,jstop) \
	for(j=jstart+START2;j<=jstop+START2;j++)
#define ZSLOOP(istart,istop,jstart,jstop) \
        for(i=istart+START1;i<=istop+START1;i++)\
	for(j=jstart+START2;j<=jstop+START2;j++)

/* loop over Primitive variables */
#define PLOOP  for(k=0;k<NPR;k++)
/* loop over all Dimensions; second rank loop */
#define DLOOP  for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP  for(j=1;j<NDIM;j++) for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA for(j=1;j<NDIM;j++)

#define delta(i,j) ( (i == j) ? 1. : 0.)
#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])


/*************************************************************************
    FUNCTION DECLARATIONS 
*************************************************************************/
double bsq_calc(double *pr, struct of_geom *geom);
double divb_calc(int i, int j);
double gdet_func(double lgcov[][NDIM]);
double get_spin();
double ran_uniform() ;
double slope_lim(double y1, double y2, double y3);
double synchrotron_cooling_func(double rho, double u, double bsq, double r);

int mhd_gamma_calc(double *pr, struct of_geom *geom, double *gamma);
int restart_init(void);

struct of_geom *get_geometry(int i, int j, int loc);

void advance_particles(double prim[][N2 + 2*NG][NPR], double Dt);
void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon);
void BL_coord(double *X, double *r, double *th);
void bound_prim(double pr[][N2 + 2*NG][NPR]);
void current_calc() ;
void conn_func(double *X, struct of_geom *geom,
	       double lconn[][NDIM][NDIM]);
void cool_down(double prim[][N2 + 2*NG][NPR], double Dt);
void coord(int i, int j, int loc, double *X);
void diag_dump();
void diag_image();
void diag_log();
void diag_pdump();
void diag_flux(double F1[][N2 + 2*NG][NPR], double F2[][N2 + 2*NG][NPR]);
void dump(FILE * fp);
void fail(int fail_type);
void find_model_limits(int i, int j, double *rhoflr, double *uuflr);
void fixup(double (*pv)[N2 + 2*NG][NPR]);
void fixup_utoprim(double (*pv)[N2 + 2*NG][NPR]);
void gcon_func(double lgcov[][NDIM], double lgcon[][NDIM]);
void gcov_func(double *X, double lgcov[][NDIM]);
void get_state(double *pr, struct of_geom *geom, struct of_state *q);
void image_ppm(double *f, char *image_descr);
void init(void);
void init_particles(void);
void init_rng(int seed) ;
void linear_mc(double x1, double x2, double x3, double *lout, double *rout) ;
void lower(double *a, struct of_geom *geom, double *b);
void lr_to_flux(double p_l[NPR], double p_r[NPR], struct of_geom *geom, int dir, 
	double Flux[NPR], double *max_speed) ;
void ludcmp(double **a, int n, int *indx, double *d);
void mhd_calc(double *pr, int dir, struct of_state *q, double *mhd);
void para(double x1, double x2, double x3, double x4, 
	double x5, double *lout, double *rout) ;
void pdump(FILE * fp) ;
void primtoflux(double *pa, struct of_state *q, int dir,
		struct of_geom *geom, double *fl);
void primtoU(double *p, struct of_state *q, struct of_geom *geom,
	     double *U);
void raise(double *v1, struct of_geom *geom, double *v2);
void reconstruct_lr_lin(double ptmp[][NPR], int N, double p_l[][NPR], double p_r[][NPR]);
void reconstruct_lr_mp5(double ptmp[][NPR], int N, double p_l[][NPR], double p_r[][NPR]);
void reconstruct_lr_par(double ptmp[][NPR], int N, double p_l[][NPR], double p_r[][NPR]);
void reconstruct_lr_weno(double ptmp[][NPR], int N, double p_l[][NPR], double p_r[][NPR]);
void restart_write(void);
void restart_read(FILE * fp);
void zero_arrays(void);
void set_coordinates(void);
void set_geometry(void);
void step_ch(void);
void source(double *pa, struct of_geom *geom, int ii, int jj, double *Ua, double Dt);
void timestep(void);
void u_to_v(double *pr, int i, int j);
void ucon_calc(double *pr, struct of_geom *geom, double *ucon);
void uRcov_calc(double *pr, struct of_geom *geom, double *uRcov);
void usrfun(double *pr, int n, double *beta, double **alpha);
void weno(double x1, double x2, double x3, double x4, 
	double x5, double *lout, double *rout) ;
int Utoprim(double *Ua, struct of_geom *geom, double *pa);

void mhd_vchar(double *pr, struct of_state *q, struct of_geom *geom,
	   int dir, double *cmax, double *cmin);
void rad_vchar(double *pr, struct of_state *q, struct of_geom *geom,
	   int dir, double *cmax, double *cmin);
void wind_source(double *ph, struct of_geom *geom, struct of_state *q,
	int ii, int jj, double *dU);

