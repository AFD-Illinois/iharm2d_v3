
#include "decs.h"

/* 

   source term for plasma in 
   continuity, momentum, and energy equations.

   goal is to minimize invocation of floors

   problem-dependent; designed for torus evolution.

*/

void wind_source(double *prim, struct of_geom *geom, struct of_state *q,
                        int ii, int jj, double *dU)
{
        double X[NDIM], drhopdt, Tp ;
        double r, th ;
        double dp[NPR],U[NPR] ;
	struct of_state qdp ;
	double cth ;
        int k ;

        /* need coordinates to evaluate particle addtn rate */
        coord(ii, jj, CENT, X);
        BL_coord(X, &r, &th);
	cth = cos(th) ;

	/* here is the rate at which we're adding particles */
	/* this function is designed to concentrate effect in the
	   funnel in black hole evolutions */
        drhopdt = 2.e-4*cth*cth*cth*cth/pow(1. + r*r,2) ;  

        dp[RHO] = drhopdt ;

        Tp = 10. ;  /* temp, in units of c^2, of new plasma */
        dp[UU] = drhopdt*Tp*3. ;

	/* add in particles in normal observer frame */
        dp[U1] = 0. ;
        dp[U2] = 0. ;
        dp[U3] = 0. ;
	/* don't change the magnetic field */
        dp[B1] = 0. ;
        dp[B2] = 0. ;
        dp[B3] = 0. ;

        /* add in plasma to the T^t_a component of the stress-energy tensor */
	/* notice that U already contains a factor of sqrt{-g} */
	get_state(dp, geom, &qdp) ;	/* can't use q(unperturbed fluid) */
	primtoU(dp, &qdp, geom, U);
        for(k=0;k<NPR;k++) dU[k] += U[k] ;

}
