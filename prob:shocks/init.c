
/*

test problem:

shock tubes

coordinates are x1,x2

implemented shocks:

- sod shock tube (SOD)
	as specified in Sod 1978, JCP 27, 1
- Brio & Wu shock tube (BRIOWU)

- TODO: komissarov shocks

*/

#include "decs.h"

/* pick your poison */
#define SOD 	1
#define BRIOWU 	0

void init()
{
	int i, j;
	double X[NDIM];

	/** generic bits **/
	
	/* set up grid functions */
	set_geometry();

	/* start diagnostic counters */
	nstep = 0;
	dump_cnt = 0;
	image_cnt = 0;
	rdump_cnt = 0;

#if SOD
	/* physical parameters */
	gam = 1.4;

	/* rescale for relativistic/nonrelativistic version */
	tf = 0.2;
	double tscale = 0.01; 		// nonrelativistic

	/* entropy wave */
	ZLOOP {
		coord(i, j, CENT, X);

                p[i][j][U1] = 0.;
                p[i][j][U2] = 0.;
                p[i][j][U3] = 0.;

                p[i][j][B1] = 0.;
                p[i][j][B2] = 0.;
                p[i][j][B3] = 0.;

                p[i][j][RHO] = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.125;

		double pgas = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.1;
                p[i][j][UU] = pgas/(gam - 1.);
	}
#endif

#if BRIOWU
	/* physical parameters */
	gam = 2.0;

	/* rescale for relativistic/nonrelativistic version */
	tf = 0.1;
	double tscale = 1.e-2; 	// nonrelativistic

	/* entropy wave */
	ZLOOP {
		coord(i, j, CENT, X);

                p[i][j][U1] = 0.;
                p[i][j][U2] = 0.;
                p[i][j][U3] = 0.;

                p[i][j][B1] = 0.75;
                p[i][j][B2] = 0.;
                p[i][j][B3] = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : -1.0;

                p[i][j][RHO] = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.125;
		double pgas =  (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.1;
                p[i][j][UU] = pgas/(gam - 1.);
	}
#endif
	/** more generic bits **/

	t = 0.;
	tf /= tscale;
	dtsave = dt = tf/1.e6;
	cour = 0.5;

	/* output frequencies */
	DTd = tf/10.;		/* dumping frequency, in units of M */
	DTl = tf/100.;		/* logfile frequency, in units of M */
	DTi = tf/100.;		/* image file frequ., in units of M */
	DTr = 512;		/* restart file frequ., in timesteps */

	/* rescale */
	ZLOOP {
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

