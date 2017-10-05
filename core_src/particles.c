
#include "decs.h"

/*

 advance particle positions using fluid half-step primitives 

 revised cfg 11 apr 2016 to improve interpolation scheme 

*/

void advance_particles(double ph[][N2 + 2*NG][NPR], double Dt)
{
	int k, l, i, j;
	double ucon[NDIM], vel[NDIM];
	double stopx[NDIM];
	double f1,f2;
	struct of_geom *geom;

	stopx[1] = startx[1] + N1*dx[1] ;
	stopx[2] = startx[2] + N2*dx[2] ;

	for (l = 0; l < NPTOT; l++) {

		/* don't update particles that are off-grid */
		if(xp[l][1] > startx[1] && xp[l][2] > startx[2] &&
		   xp[l][1] < stopx[1]  && xp[l][2] < stopx[2]) {

			/* the four-velocities are zone-centered */
			f1 = (xp[l][1] - startx[1] + 0.5*dx[1]) / dx[1];
			f2 = (xp[l][2] - startx[2] + 0.5*dx[2]) / dx[2];

		   	/* find nearest zone center */
			i = lround( f1 ) ;
			j = lround( f2 ) ;

                        geom = get_geometry(i, j, CENT);
                        ucon_calc(ph[i][j], geom, ucon);
                        for (k = 1; k < NDIM; k++) vel[k] = ucon[k]/ucon[0];

			/* push particle forward */
			for (k = 1; k < NDIM; k++)
				xp[l][k] += Dt * vel[k];

		}

	}
	/* done! */
}

/* 
	output coordinates and velocities of particles 
*/

void pdump(FILE * fp)
{
	int l;

	for (l = 0; l < NPTOT; l++) {
		fprintf(fp, "%g %g %g\n", xp[l][1], xp[l][2], xp[l][3]);
	}
}

/*

 initialize Lagrangian tracer particles
 cfg 4 feb 09

 simplified 10 apr 2016 cfg

*/

void init_particles()
{
	int i, j, Np;
	double dmass, mass, Nexp, sample_factor, X[NDIM];
	struct of_geom *geom ;
	struct of_state q ;
	double U[NPR];

	/* global variables */
	pdump_cnt = 0;
	DTp = 20.;

	/* assign particles according to restmass density */
	/* first find total rest-mass on grid */
	mass = 0.;
	ZLOOP {
		geom = get_geometry(i, j, CENT) ;
		get_state(p[i][j], geom, &q);
		primtoU(p[i][j], &q, geom, U) ;

		dmass = U[RHO]*dx[1]*dx[2]*dx[3] ;
		mass += dmass;
	}

	/* normalization factor so that we get the # of
	   particles we want */
	sample_factor = NPTOT / mass;

	/* now sample, keeping running total of particles
	   created */
	Np = 0;
	Nexp = 0.;
	ZLOOP {
		geom = get_geometry(i, j, CENT) ;
		get_state(p[i][j], geom, &q);
		primtoU(p[i][j], &q, geom, U) ;

		dmass = U[RHO]*dx[1]*dx[2]*dx[3] ;

		Nexp += dmass * sample_factor;

		/* assign particle to random position in cell */
		while(Nexp >= 0.5) {

			/* lower left corner of zone is here */
			coord(i, j, CORN, X);

			xp[Np][0] = 0.;
			xp[Np][1] = X[1] + ran_uniform() * dx[1];
			xp[Np][2] = X[2] + ran_uniform() * dx[2];
			xp[Np][3] = startx[3] + ran_uniform() * dx[3];

			Np++;
			Nexp -= 1.;
		}
	}

	/* report back! */
	fprintf(stderr, "made %d particles of %d (%g)\n", Np, NPTOT, Nexp);

	/* done! */
}
