
#include "decs.h"

/*
 
 this file contains all the coordinate dependent
 parts of the code, except the initial and boundary
 conditions 

 gcov_func: return line element
 set_coordinates: fixes grid position
 get_spin: return black hole spin a
 BL_coord: return Boyer-Lindquist coordinates from x1,x2
 model_limits: set floors on U and rho in a coordinate-dependent way
 
*/

/* coordinate-dependent parameters: 
   not accessible elsewhere in the code */

/* mandatory routines */
/* insert covariant (line element) form of metric here */

void gcov_func(double *X, double gcov[][NDIM])
{
	int j,k;

	/* Minkowski metric */
	DLOOP gcov[j][k] = 0.;
	gcov[0][0] = -1.;
	gcov[1][1] = 1.;
	gcov[2][2] = 1.;
	gcov[3][3] = 1.;
}

/* some coordinate parameters, grid location, dx */
void set_coordinates()
{

	startx[1] = -6.;
	startx[2] = -6.;
	startx[3] = 0.;

	dx[1] = 12. / N1;
	dx[2] = 12. / N2;
	dx[3] = 1.;	/* used to evaluate integrals */

}

#define RHOMIN  (1.e-5)
#define UUMIN   (1.e-7)

void find_model_limits(int i, int j, double *rhoflr, double *uuflr)
{

        /* set floors */
        *rhoflr = RHOMIN ;
        *uuflr = UUMIN ;
        /* end set floors */
}

#undef UUMIN
#undef RHOMIN

