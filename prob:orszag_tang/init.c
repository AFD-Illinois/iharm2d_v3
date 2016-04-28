
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

	/* end of simulation.  time unit = G M/c^3 */
	t = 0.;
	tf = 2.;
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

		double v1 = -0.5*sin(2.*M_PI*X[2]);
                double v2 =  0.5*sin(2.*M_PI*X[1]);
                double lorentzFactor = 1./sqrt(1. - v1*v1 - v2*v2);

                p[i][j][U1] = lorentzFactor*v1;
                p[i][j][U2] = lorentzFactor*v2;
                p[i][j][U3] = 0.;

                p[i][j][B1] = -1./sqrt(4.*M_PI) * sin(2.*M_PI*X[2]);
                p[i][j][B2] =  1./sqrt(4.*M_PI) * sin(4.*M_PI*X[1]);
                p[i][j][B3] = 0. ;

                p[i][j][RHO] = 25./(36.*M_PI);
                p[i][j][UU] = 5./(12.*M_PI*(gam - 1.));
	}

	/* enforce boundary conditions */
	fixup(p);
	bound_prim(p);

}

