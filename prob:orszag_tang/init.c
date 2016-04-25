
/*

initial conditions for relativistic orszag-tang vortex.

*/

#include "decs.h"

void init()
{
	int i, j;
	double X[NDIM];

	/* physical parameters */
	gam = 4./3.;

	/* set up grid functions */
	set_geometry();

	/* end of simulation.  time unit = G M/c^3 */
	t = 0.;
	tf = 4.0*M_PI ;
	dtsave = dt = 1.e-4;
	cour = 0.9;

	/* output frequencies */
	DTd = tf/10.;		/* dumping frequency, in units of M */
	DTl = tf/100.;		/* logfile frequency, in units of M */
	DTi = tf/100.;		/* image file frequ., in units of M */
	DTr = 512;		/* restart file frequ., in timesteps */

	/* start diagnostic counters */
	nstep = 0;
	dump_cnt = 0;
	image_cnt = 0;
	rdump_cnt = 0;

	/* relativistic Orszag-Tang vortex */
	ZLOOP {
		coord(i, j, CENT, X);

                p[i][j][U1] = -sin(X[2] + M_PI);
                p[i][j][U2] = sin(X[1] + M_PI);
                p[i][j][U3] = 0. ;
                p[i][j][B1] = -sin(X[2] + M_PI) ;
                p[i][j][B2] = sin(2.*(X[1] + M_PI)) ;
                p[i][j][B3] = 0. ;

                p[i][j][RHO] = 25./9. ;
                p[i][j][UU] = (5./3.)/(gam - 1.) ;
	}

	/* enforce boundary conditions */
	fixup(p);
	bound_prim(p);

}

