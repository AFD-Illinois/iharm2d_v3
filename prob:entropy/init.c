
/*

test problem:

advection of a magnetized entropy wave.

coordinates are x1,x2

initial conditions:
x1 < 0: B2 = B, rho = 0.01,
x1 > 0: B2 = 0, rho = 1, 
Ptot = B^2/2 + Pgas = const.

u1 = 0.3

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
	tf = 10.;
	dtsave = dt = 1.e-4;
	cour = 0.8;

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

	/* entropy wave */
	double BB,PB;
	BB = 0.5;
	PB = BB*BB*0.5;
	double v1 = 2./tf ;
	double u1 = v1/sqrt(1. - v1*v1);
	double rhomin = 0.5;
	double rhomax = 1.;
	double width = 0.02;
	double ptot = PB*1.5;
	ZLOOP {
		coord(i, j, CENT, X);

		double f = tanh(sin(2.*M_PI*X[1])/(width*2.*M_PI));

                p[i][j][U1] = u1;
                p[i][j][U2] = 0.;
                p[i][j][U3] = 0.;

                p[i][j][B1] = 0.;
                p[i][j][B2] = (1. - f)*0.5*BB ;
                p[i][j][B3] = 0. ;

                p[i][j][RHO] = (1. + f)*0.5*rhomax + rhomin ;
		double pgas = ptot - p[i][j][B2]*p[i][j][B2]/2.;
                p[i][j][UU] = pgas/(gam - 1.);
	}

	/* enforce boundary conditions */
	fixup(p);
	bound_prim(p);

}

