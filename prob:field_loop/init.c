
/*

smooth field loop advection test problem

description follows initial conditions for the athena
	test problem

alternative description uses smooth current configuration.

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
	double eps = 0.05;
	t = 0.;
	tf = 2./sin(M_PI/3.) * 1./eps;
	dtsave = dt = 0.01/dx[1];
	cour = 0.95;

	/* output frequencies */
	DTd = tf/10.;		/* dumping frequency, in units of M */
	DTl = tf/100.;		/* logfile frequency, in units of M */
	DTi = tf/200.;		/* image file frequ., in units of M */
	DTr = 512;		/* restart file frequ., in timesteps */

	/* start diagnostic counters */
	nstep = 0;
	dump_cnt = 0;
	image_cnt = 0;
	rdump_cnt = 0;

	/* field loop */
	double A[N1+2*NG][N2+2*NG];
	double r;
	for(i=0;i<N1+2*NG;i++)
	for(j=0;j<N2+2*NG;j++) 
		A[i][j] = 0.;
	for(i=0;i<N1+2*NG;i++)
	for(j=0;j<N2+2*NG;j++) {
		coord(i, j, CORN, X);

		r = sqrt(X[1]*X[1] + X[2]*X[2]);

		/* original gardner & stone version */
		/* current is discontinuous */
		A[i][j] = eps*1.e-4*fmax(0.3-r, 0.);

		/* smooth field loop */
		/* current is continuous */
		A[i][j] = eps*1.e-3*exp(-r*r/(2.*0.2*0.2));

		fprintf(stderr,"%d %d %g %g\n",i,j,r,A[i][j]) ;

	}
	ZLOOP {
		coord(i, j, CENT, X);

                p[i][j][U1] = eps*sin(M_PI/3.);
                p[i][j][U2] = eps*cos(M_PI/3.);
                p[i][j][U3] = 0.;

                p[i][j][B1] = -(A[i][j] - A[i][j + 1]
                                + A[i + 1][j] - A[i + 1][j + 1]) / 
                                (2. * dx[2]);
                p[i][j][B2] = (A[i][j] + A[i][j + 1]
                                - A[i + 1][j] - A[i + 1][j + 1]) / 
                                (2. * dx[1]) ;
                p[i][j][B3] = 0. ;

                p[i][j][RHO] = 1.;
                p[i][j][UU] = eps*eps/(gam - 1.);

	}

	/* enforce boundary conditions */
	fixup(p);
	bound_prim(p);

}

