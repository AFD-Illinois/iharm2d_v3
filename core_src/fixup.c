
#include "decs.h"

/* 

   fixup:
   check a set of primitive variables.

   set minimum density (set "floor" on density)
   set minimum internal energy 
   set maximum lorentz factor wrt normal observer
   	(normal observer has u_{\mu} = {n, 0, 0, 0})


   fixup_utoprim:
   on zones where utoprim fails, replace primitives
   other than B with averages over nearby non-failed
   zones, if any exist. if none exist, then set
   velocity to normal observer velocity and 
   leave density and internal energy alone.
   at the end, run fixup.

*/

void fixup(double (*pv)[N2 + 2*NG][NPR])
{
	int i, j;
	int nfix = 0;
	int fixup1zone(int i, int j, double pv[NPR]) ;

	/* loop over zones */
	ZLOOP nfix += fixup1zone(i, j, pv[i][j]);

	/* report how many zones were fixed on each pass */
	if(nfix > 0) fprintf(stderr,"[fix: %d]",nfix) ;	
}

int fixup1zone(int i, int j, double pv[NPR])
{
	double rhoflr, uuflr;
	double f, gamma;
	struct of_geom *geom;

	/* set density, internal energy floors */
	find_model_limits(i, j, &rhoflr, &uuflr);

	/* implement floor on rho, u (momentum *not* conserved) */
	int fix_flag = 0;
	if (pv[RHO] < rhoflr || isnan(pv[RHO])) {
		pv[RHO] = rhoflr;
		fix_flag++;
	}
	if (pv[UU] < uuflr || isnan(pv[UU])) {
		pv[UU] = uuflr;
		fix_flag++;
	}
	/* end rho, u floors */

	/* limit lorentz factor of fluid wrt normal observer */
	geom = get_geometry(i, j, CENT) ;
	if (mhd_gamma_calc(pv, geom, &gamma)) {
		/* Treat gamma failure here as "fixable" for fixup_utoprim() */
		pflag[i][j] = -333;
	} else {
		if (gamma > GAMMAMAX) {
			f = sqrt((GAMMAMAX * GAMMAMAX - 1.) / (gamma * gamma - 1.));
			pv[U1] *= f;
			pv[U2] *= f;
			pv[U3] *= f;

			fix_flag++;
		}
	}
	/* end lorentz factor limit */

	return(fix_flag);
}


/**************************************************************************************
 INTERPOLATION STENCILS:  
 ------------------------
   -- let the stencils be characterized by the following numbering convention:

           1 2 3 
           8 x 4      where x is the point at which we are interpolating 
           7 6 5
*******************************************************************************************/

/* 2468  */
#define AVG4_1(pr,i,j,k) (0.25*(pr[i][j+1][k]+pr[i][j-1][k]+pr[i-1][j][k]+pr[i+1][j][k]))

/* 1357  */
#define AVG4_2(pr,i,j,k) (0.25*(pr[i+1][j+1][k]+pr[i+1][j-1][k]+pr[i-1][j+1][k]+pr[i-1][j-1][k]))

/*******************************************************************************************
  fixup_utoprim(): 

    -- figures out (w/ pflag[]) which stencil to use to interpolate bad point from neighbors;

    -- here we use the following numbering scheme for the neighboring cells to i,j:  

                      1  2  3 
                      8  x  4        where "x" is the (i,j) cell or the cell to be interpolated
                      7  6  5

 *******************************************************************************************/

void fixup_utoprim(double (*pv)[N2 + 2*NG][NPR])
{
	int i, j, k;
	static int pf[9];
	int fcount ;

	/* count the number of failed zones */
	fcount = 0 ;
	ZSLOOP(-NG, (N1-1 + NG), -NG, (N2-1 + NG)) {
		if(pflag[i][j]) fcount++ ;
	}
	if(fcount > 0) fprintf(stderr,"[fail: %d]",fcount) ;

	/* Fix interior points only */
	ZSLOOP(0, (N1 - 1), 0, (N2 - 1)) {
		if (pflag[i][j]) {
			pf[1] = !pflag[i - 1][j + 1];
			pf[2] = !pflag[i][j + 1];
			pf[3] = !pflag[i + 1][j + 1];
			pf[8] = !pflag[i - 1][j];
			pf[4] = !pflag[i + 1][j];
			pf[7] = !pflag[i - 1][j - 1];
			pf[6] = !pflag[i][j - 1];
			pf[5] = !pflag[i + 1][j - 1];

			/* pf's  are true if they represent good points */
			if (pf[2] && pf[4] && pf[6] && pf[8]) {
				/* leave B alone */
				for(k=0;k<B1;k++) pv[i][j][k] = AVG4_1(pv, i, j, k);
			} else if (pf[1] && pf[3] && pf[5] && pf[7]) {
				/* leave B alone */
				for(k=0;k<B1;k++) pv[i][j][k] = AVG4_2(pv, i, j, k);
			} else {
				/* all is lost.
				   average rho, u over all nearby zones
				   set v = normal observer values
				   leave B alone. */
				for (k = RHO; k <= UU; k++) {
					pv[i][j][k] = 0.5 * (AVG4_1(pv, i, j, k) +
						   AVG4_2(pv, i, j, k));
				}
				pv[i][j][U1] = pv[i][j][U2] = pv[i][j][U3] = 0.;
			}
			pflag[i][j] = 0; /* cell fixed, can use it for interpolation 
					   elsewhere */

			/* finally: apply limits */
			fixup1zone(i, j, pv[i][j]);  
		}
	}

	return;
}

