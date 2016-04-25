
/* global dependent variable arrays */
double p[N1 + 2*NG][N2 + 2*NG][NPR];	/* space for primitive vars */
double dq[N1 + 2*NG][N2 + 2*NG][NPR];	/* slopes */
double F1[N1 + 2*NG][N2 + 2*NG][NPR];	/* fluxes */
double F2[N1 + 2*NG][N2 + 2*NG][NPR];	/* fluxes */
double ph[N1 + 2*NG][N2 + 2*NG][NPR];	/* half-step primitives */
int pflag[N1 + 2*NG][N2 + 2*NG];	/* identifies failure points */

/* for debug & diagnostics */
double psave[N1+2*NG][N2+2*NG][NPR];
double Jcon[N1+2*NG][N2+2*NG][NDIM];
double dtsave ;

/* grid functions */
double conn[N1 + 2*NG][N2 + 2*NG][NDIM][NDIM][NDIM];
struct of_geom ggeom[N1+2*NG][N2+2*NG][NPG] ;

/* particles */
double xp[NPTOT][NDIM];

/*************************************************************************
    GLOBAL VARIABLES SECTION 
*************************************************************************/
/* physics parameters */
double gam;

/* numerical parameters */
double cour;
double dx[NPR], startx[NPR];
double dt;
double t, tf;
int istart, istop, jstart, jstop;
int nstep;

/* output parameters */
double DTd;
double DTp;
double DTl;
double DTi;
int DTr;
int dump_cnt;
int pdump_cnt;
int image_cnt;
int rdump_cnt;
int nstroke;

/* diagnostics */
double sum_boundary_flux[2*NDIM][NPR];	/* fluxes of each variable on each boundary */

/* current local position */
int icurr, jcurr;

