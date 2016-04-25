
/*

initial conditions for Komissarov cylindrical
explosion problem.

see Komissarov 1999, MNRAS 303, 343

explosion problem is described on p. 360.

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
	tf = 4. ;
	dtsave = dt = 1.e-4*dx[1]/1.;
	cour = 0.1;

	/* output frequencies */
	DTd = tf/10.;		/* dumping frequency, in units of M */
	DTl = tf/10.;		/* logfile frequency, in units of M */
	DTi = tf/50.;		/* image file frequ., in units of M */
	DTr = 512;		/* restart file frequ., in timesteps */

	/* start diagnostic counters */
	nstep = 0;
	dump_cnt = 0;
	image_cnt = 0;
	rdump_cnt = 0;

	double r;
	double fr(double R,double min,double max) ;
	ZLOOP {
		coord(i, j, CENT, X);
                r = sqrt(X[1]*X[1] + X[2]*X[2]) ;

                p[i][j][U1] = 0. ;
                p[i][j][U2] = 0. ;
                p[i][j][U3] = 0. ;

		/* strong: 1; relatively weak: 0.1; weak: 0.01 */
                p[i][j][B1] = 0.01 ;
                p[i][j][B2] = 0. ;
                p[i][j][B3] = 0. ;

                p[i][j][RHO] = fr(r,1.e-4,1.e-2) ;
                p[i][j][UU] = fr(r,3.e-5,1.)/(gam - 1.) ;
	}

	/* enforce boundary conditions */
	fixup(p);
	bound_prim(p);

}

/* this is the initial density, internal energy profile
   described on Komissarov p. 360 */
double fr(double R,double min,double max)
{
        double Rout,Rin,dR ;

        Rin = 0.8 ;
        Rout = 1.0 ;
        dR = (Rout - Rin) ;

        if(R > Rout) return(min) ;
        if(R < Rin) return(max) ;
        else {
                return( exp(
                        log(max)*(Rout - R)/dR +
                        log(min)*(R - Rin)/dR)) ;
        }
}

