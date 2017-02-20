
/*

Kelvin-Helmholtz instability problem

follows initial conditions from Lecoanet et al. 2015,
	MNRAS 455, 4274.

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

	/* scale toward nonrelativistic limit */
	double tscale;
	tscale = 1. ; // relativistic
	tscale = 0.01 ; // nonrelativistic

	/* end of simulation */
	t = 0.;
	tf = 6.0/tscale;
	dtsave = dt = 1.e-4;
	cour = 0.95;

	/* output frequencies */
	DTd = tf/10.;		/* dumping frequency, in units of M */
	DTl = tf/100.;		/* logfile frequency, in units of M */
	DTi = tf/400.;		/* image file frequ., in units of M */
	DTr = 512;		/* restart file frequ., in timesteps */

	/* start diagnostic counters */
	nstep = 0;
	dump_cnt = 0;
	image_cnt = 0;
	rdump_cnt = 0;

	/* follows notation of Lecoanet et al. eq. 8 et seq. */
	double rho0 = 1.;	
	double Drho = 0.1;
	double P0 = 10.;
	double uflow = 1.;
	double a = 0.05;
	double sigma = 0.2;
	double A = 0.01;
	double z1 = 0.5;
	double z2 = 1.5;
	double x,z;

	ZLOOP {
		coord(i, j, CENT, X);

		/* Lecoanet's x <-> x1; z <-> x2 */
		x = X[1];
		z = X[2];

                p[i][j][RHO] = 1. + (Drho/rho0)*0.5*(
			tanh((z - z1)/a) - tanh((z - z2)/a));
                p[i][j][UU] = P0/(gam - 1.);

		p[i][j][U1] = uflow*(tanh((z - z1)/a) -
			tanh((z - z2)/a) - 1.);
		p[i][j][U2] = A*sin(2.*M_PI*x)*(
			exp(-(z-z1)*(z-z1)/(sigma*sigma)) +
			exp(-(z-z2)*(z-z2)/(sigma*sigma)));
                p[i][j][U3] = 0.;

                p[i][j][B1] = 0.;
                p[i][j][B2] = 0.;
                p[i][j][B3] = 0.;


		/* rescale */
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

