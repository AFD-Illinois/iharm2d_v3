
/* 
   include automatic measurement of zone cycles/second 
   
   cfg 12.24.14 
   
*/

#include "decs.h"
#include "defs.h"	
#include <time.h>

/* this indirection is used to select the inverter
   and the reconstruction routine with the preprocessor */
#define QUOTE(macro) #macro
#define STR(macro) QUOTE(macro)
#define RECONSTR STR(RECONSTRUCT_LR)

#define  TIMER_NSTEP	(256)

/*****************************************************************/
/*****************************************************************
   main():
   ------

     -- Initializes, time-steps, and concludes the simulation. 
     -- Handles timing of output routines;
     -- Main is main, what more can you say.  

-*****************************************************************/
int main(int argc, char *argv[])
{
	double tdump, tpdump, timage, tlog, zcps;
	time_t ti_t,t0_t,t1_t,tf_t ;

	/* initialize code timer */
	ti_t = time(NULL) ;	
	t0_t = ti_t ;

        /** begin report on options, to stderr **/
        fprintf(stderr,"interpolation: %s\n", RECONSTR) ;
	if(!strncmp(RECONSTR ,"reconstruct_lr_lin",18) && NG == 2) {
		fprintf(stderr,"not enough ghost zones!\n") ;
		fprintf(stderr,"WENO or PARA + NG=2\n") ;
		exit(1);
	}

        fprintf(stderr,"Wind source: ") ;
        if(WIND) fprintf(stderr,"ENABLED\n") ;
        else fprintf(stderr,"DISABLED\n") ;

        fprintf(stderr,"Cooling: ") ;
        if(COOLING) fprintf(stderr,"ENABLED\n") ;
        else fprintf(stderr,"DISABLED\n") ;

        fprintf(stderr,"Lagrangian tracers: ") ;
        if(NPTOT <= 0) fprintf(stderr,"DISABLED\n") ;
        else fprintf(stderr,"ENABLED\n") ;
	/** end report on options **/

	/* initialize, either directly from init.c or from restart (checkpoint) dump */
	system("mkdir -p dumps images");
	zero_arrays();
	if (!restart_init()) {
		init_rng(1) ;	/* standard seed for reproducibility */
		init();
		if(NPTOT > 0) init_particles();
	}

	/* emit initial diagnostics */
	diag_log();
	diag_dump();
	diag_image();
	if(NPTOT > 0) diag_pdump();

	/* initialize diagnostic timers */
	tdump = t + DTd;
	if(NPTOT > 0) tpdump = t + DTp;
	timage = t + DTi;
	tlog = t + DTl;

	/** begin evolution **/
	fprintf(stderr,"t, tf: %g %g\n",t,tf) ;
	while (t < tf) {

		/* this pushes variables forward in time */
		step_ch();
		nstep++;

		/* uncomment to get diagnostics on every step */
		//diag_dump();
		//diag_image();

		fprintf(stderr, "%10.5g %10.5g %8d\n", t, dt, nstep) ;

		/* begin diagnostics */
		if (t >= tdump) {
			diag_dump();
			tdump += DTd;
		}
		if (t >= timage) {
			diag_image();
			timage += DTi;
		}
		if (t >= tlog) {
			diag_log();
			tlog += DTl;
		}
		if (NPTOT > 0 && t >= tpdump) {
			diag_pdump();
			tpdump += DTp;
		}
		/* end diagnostics */

		/* restart dump at fixed timestep interval */
		if (nstep % DTr == 0)
			restart_write();

		/* report on speed */
		//if(nstep == 1024) exit(0) ; 
		if(nstep%TIMER_NSTEP == 0) {
			t1_t = time(NULL) ;
			zcps = TIMER_NSTEP*N1*N2/((double)(t1_t-t0_t)) ;
			fprintf(stderr,"Zone cycles/second: %g\n", zcps) ;
			t0_t = t1_t ;
		}
	}
	/** end evolution **/

	/* final speed report */
	tf_t = time(NULL) ;
	zcps = nstep*N1*N2/((double)(tf_t-ti_t)) ;
	fprintf(stderr,"Run average zone cycles/second: %g\n", zcps) ;

	/* emit final diagnostics */
	diag_log();
	diag_dump();
	diag_image();
	if(NPTOT > 0) diag_pdump();

	/* done! */
	return (0);
}

/* initialize primitive variable arrays, slopes, fluxes, flags */
void zero_arrays()
{
        int i, j, k;

        /* everything must be initialized to zero */
        ZSLOOP(-NG, N1-1 + NG, -NG, N2-1 + NG) {
                PLOOP {
                        p[i][j][k] = 0.;
                        ph[i][j][k] = 0.;
                        dq[i][j][k] = 0.;
                        F1[i][j][k] = 0.;
                        F2[i][j][k] = 0.;
                }
                pflag[i][j] = 0;
        }
}

