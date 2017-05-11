
/*

initial conditions for relativistic orszag-tang vortex.

description follows initial conditions for the athena
	test problem

*/

#include "decs.h"

void init()
{
	int i, j;
	double X[NDIM];

	/* physical parameters */
	gam = 5./3.;

	/* set up grid functions */
	set_geometry();

	/* problem is relativistic at tscale = 1; 
		nonrelativistic for tscale << 1. */
	double tscale = 0.05;

	/* end of simulation.  time unit = G M/c^3 */
	t = 0.;
	tf = 1.5*M_PI/tscale;
	dtsave = dt = 1.e-4;
	cour = 0.9;

	/* output frequencies */
	DTd = tf/10.;		/* dumping frequency, in units of M */
	DTl = tf/200.;		/* logfile frequency, in units of M */
	DTi = tf/200.;		/* image file frequ., in units of M */
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
                p[i][j][U2] =  sin(X[1] + M_PI) ;
                p[i][j][U3] = 0.;

                p[i][j][B1] = -sin(X[2] + M_PI);
                p[i][j][B2] =  sin(2.*(X[1] + M_PI));
                p[i][j][B3] = 0. ;

                p[i][j][RHO] = 25./9.;
                p[i][j][UU] = 5./(3.*(gam - 1.));

                p[i][j][U1] *= tscale ;
                p[i][j][U2] *= tscale ;
                p[i][j][U3] *= tscale ;
                p[i][j][B1] *= tscale ;
                p[i][j][B2] *= tscale ;
                p[i][j][B3] *= tscale ;
                p[i][j][UU] *= tscale*tscale ;
	}

	/* enforce boundary conditions */
	fixup(p);
	bound_prim(p);

}

