
/*
 generic code for advancing primitive variables 
*/

#include "decs.h"

/* local functions for step_ch */
double advance(double pi[][N2 + 2*NG][NPR], double pb[][N2 + 2*NG][NPR],
	       double Dt, double pf[][N2 + 2*NG][NPR]);
double fluxcalc(double pr[][N2 + 2*NG][NPR], double F[][N2 + 2*NG][NPR],
		int dir);
void   flux_ct(double F1[][N2 + 2*NG][NPR], double F2[][N2 + 2*NG][NPR]);

/* maximum fractional increase in timestep per timestep */
#define SAFE    (1.3)


/*
  step_ch():
  ---------
     -- handles the sequence of making the time step, the fixup of unphysical values, 
        and the setting of boundary conditions;

     -- also sets the dynamically changing time step size;

*/

void step_ch()
{
	double ndt;
	int i, j, k;

	/* step primitive variables to the half step */
	fprintf(stderr, "h");
	ndt = advance(p, p, 0.5 * dt, ph);	

	fixup(ph);		/* Set updated densities to floor, set limit for gamma */
	bound_prim(ph);		/* Set boundary conditions for primitive variables, flag bad ghost zones */
	fixup_utoprim(ph);	/* Fix the failure points using interpolation */
	bound_prim(ph);		/* Reset boundary conditions with fixed up points */

	/* step Lagrangian tracer particles forward */
	if(NPTOT > 0) advance_particles(ph, dt);

	/* take second-order-accurate full time step */
	fprintf(stderr, "f");
	ZSLOOP(-NG,N1-1+NG,-NG,N2-1+NG) PLOOP psave[i][j][k] = p[i][j][k];
	ndt = advance(p, ph, dt, p);

	dtsave = dt ;

#if ( COOLING )
	/* operator-split the cooling. 
	   integrate internal energy equation in fluid frame */
	cool_down(p, dt);
#endif

	fixup(p);
	bound_prim(p);
	fixup_utoprim(p);
	bound_prim(p);

	/* check timestep */
	if (dt < 1.e-15*tf) {
		fprintf(stderr, "timestep too small\n");
		exit(11);
	}

	/* increment time */
	t += dt;

	/* set next timestep */
	if (ndt > SAFE * dt) ndt = SAFE * dt;	/* allow only moderate growth in timestep */
	if (t + ndt > tf) ndt = tf - t;	        /* don't step beyond end of run */
	dt = ndt;

	/* done! */
}

/***********************************************************************************************
  advance():
  ---------
     -- responsible for what happens during a time step update, including the flux calculation, 
         the constrained transport calculation (aka flux_ct()), the finite difference 
         form of the time integral, and the calculation of the primitive variables from the 
         update conserved variables;
     -- also handles the "fix_flux()" call that sets the boundary condition on the fluxes;

***********************************************************************************************/
double advance(double pi[][N2 + 2*NG][NPR],
	       double pb[][N2 + 2*NG][NPR],
	       double Dt, 
	       double pf[][N2 + 2*NG][NPR])
{
	int i, j, k;
	double ndt, ndt1, ndt2, U[NPR], dU[NPR];
	struct of_state q;

	ZLOOP PLOOP pf[i][j][k] = pi[i][j][k];	/* needed for Utoprim */

	fprintf(stderr, "0");

	ndt1 = fluxcalc(pb, F1, 1);
	ndt2 = fluxcalc(pb, F2, 2);

	/* implement Toth's flux_ct constrained transport */
	flux_ct(F1, F2);

	/* evaluate diagnostics based on fluxes */
	diag_flux(F1, F2);

	fprintf(stderr, "1");
#pragma omp parallel \
 shared ( pi, pb, pf, F1, F2, ggeom, pflag, dx, Dt ) \
 private ( i, j, k, dU, q, U )
{
#pragma omp for
	/** now update pi to pf **/
	ZLOOP {

		source(pb[i][j], &(ggeom[i][j][CENT]), i, j, dU, Dt);
		get_state(pi[i][j], &(ggeom[i][j][CENT]), &q);
		primtoU(pi[i][j], &q, &(ggeom[i][j][CENT]), U);

		PLOOP {
			U[k] +=
			    Dt * (-(F1[i + 1][j][k] - F1[i][j][k]) / dx[1]
				  - (F2[i][j + 1][k] - F2[i][j][k]) / dx[2]
				  + dU[k]
			    );
		}

		pflag[i][j] = Utoprim(U, &(ggeom[i][j][CENT]), pf[i][j]);
	}
}

	ndt = 1. / (1. / ndt1 + 1. / ndt2);
	//fprintf(stderr,"\nndt: %g %g %g\n",ndt,ndt1,ndt2) ;
	fprintf(stderr, "2");

	return (ndt);
}


/*
  fluxcalc():
  ---------
     -- sets the numerical fluxes, avaluated at the cell boundaries using the slope limiter
        slope_lim();

     -- only has HLL and Lax-Friedrichs  approximate Riemann solvers implemented;
        
*/

#define OFFSET	1
double fluxcalc(double pr[][N2 + 2*NG][NPR], double F[][N2 + 2*NG][NPR], int dir)
{
	int i, j, k ;
	double p_l[NMAX+2*NG][NPR], p_r[NMAX+2*NG][NPR] ;
	double ndt,dtij;
	static double cmax[N1+2*NG][N2+2*NG] ;
	double ptmp[NMAX+2*NG][NPR] ;

	ndt = 1.e3 ;

	if(dir == 1) {
		/* loop over other direction */
#pragma omp parallel \
 shared ( pr, ggeom, F, dir, cmax ) \
 private ( ptmp, p_l, p_r, i, j, k )
{
#pragma omp for
		JSLOOP(-1,N2) {

			/* copy out variables and operate on those */
			ISLOOP(-NG,N1-1+NG) PLOOP ptmp[i][k] = pr[i][j][k] ;

			/* reconstruct left & right states */
			RECONSTRUCT_LR(ptmp, N1, p_l, p_r) ;

			ISLOOP(0,N1) {
				lr_to_flux(p_r[i-1],p_l[i], &(ggeom[i][j][FACE1]), dir, F[i][j], &cmax[i][j]) ;
			}
		}
}
	}
	else if(dir == 2) {
		/* loop over other direction */
#pragma omp parallel \
 shared ( pr, ggeom, F, dir, cmax ) \
 private ( ptmp, p_l, p_r, i, j, k )
{
#pragma omp for
		ISLOOP(-1,N1) {

			/* copy out variables and operate on those */
			JSLOOP(-NG,N2-1+NG) PLOOP ptmp[j][k] = pr[i][j][k] ;

			/* reconstruct left & right states */
			RECONSTRUCT_LR(ptmp, N2, p_l, p_r) ;

			JSLOOP(0,N2) 
				lr_to_flux(p_r[j-1],p_l[j], &(ggeom[i][j][FACE2]), dir, F[i][j], &cmax[i][j]) ;

		}
}
	}
	ZLOOP {
		dtij = cour * dx[dir] / cmax[i][j];
		if (dtij < ndt) ndt = dtij;
	}

	return (ndt);

}

/***********************************************************************************************/
/***********************************************************************************************
  flux_ct():
  ---------
     -- performs the flux-averaging used to preserve the del.B = 0 constraint (see Toth 2000);
        
***********************************************************************************************/
void flux_ct(double F1[][N2 + 2*NG][NPR], double F2[][N2 + 2*NG][NPR])
{
	int i, j;
	static double emf[N1 + 2*NG][N2 + 2*NG];

	/* calculate EMFs */
	/* Toth approach: just average to corners */
	/* best done in scalar mode! */
	ZSLOOP(0, N1, 0, N2) {
		emf[i][j] = 0.25 * (F1[i][j][B2] + F1[i][j - 1][B2]
		                  - F2[i][j][B1] - F2[i - 1][j][B1]);
	}

	/* rewrite EMFs as fluxes, after Toth */
	ZSLOOP(0, N1, 0, N2 - 1) {
		F1[i][j][B1] = 0.;
		F1[i][j][B2] = 0.5 * (emf[i][j] + emf[i][j + 1]);
	}
	ZSLOOP(0, N1 - 1, 0, N2) {
		F2[i][j][B1] = -0.5 * (emf[i][j] + emf[i + 1][j]);
		F2[i][j][B2] = 0.;
	}

}

void lr_to_flux(double p_l[NPR], double p_r[NPR], struct of_geom *geom, 
		int dir, double Flux[NPR], double *max_speed)
{
	struct of_state state_l, state_r ;
	double F_l[NPR],F_r[NPR],U_l[NPR],U_r[NPR] ;
	double cmax_l,cmax_r,cmin_l,cmin_r,ctop ;
	int k ;

	get_state(p_l, geom, &state_l);
	get_state(p_r, geom, &state_r);

	primtoflux(p_l, &state_l, dir, geom, F_l);
	primtoflux(p_r, &state_r, dir, geom, F_r);

	primtoflux(p_l, &state_l, 0, geom, U_l);
	primtoflux(p_r, &state_r, 0, geom, U_r);

	mhd_vchar(p_l, &state_l, geom, dir, &cmax_l, &cmin_l);
	mhd_vchar(p_r, &state_r, geom, dir, &cmax_r, &cmin_r);

	ctop = fmax(fabs(cmax_l), 
		fmax( fabs(cmax_r),
		 fmax( fabs(cmin_r), fabs(cmin_l) ))) ;

	/* local Lax-Friedrichs flux */
	PLOOP Flux[k] = 0.5*(F_l[k] + F_r[k] - ctop*(U_r[k] - U_l[k]));

	*max_speed = ctop ;

	return ;
}

/* sum fluxes over boundaries.  coordinate-independent */
void diag_flux(double F1[][N2 + 2*NG][NPR], double F2[][N2 + 2*NG][NPR])
{
        int i,j,k,l;

	for(l=0;l<NDIM;l++)
	for(k=0;k<NPR;k++)
		sum_boundary_flux[l][k] = 0. ;

        for (j = START2; j < N2+START2; j++) {
		for(k=0;k<NPR;k++) sum_boundary_flux[2][k] += F1[0+START1][j][k]    * dx[2] * dx[3];
		for(k=0;k<NPR;k++) sum_boundary_flux[3][k] += F1[0+START1+N1][j][k] * dx[2] * dx[3];
        }
        for (i = START1; i < N1+START1; i++) {
		for(k=0;k<NPR;k++) sum_boundary_flux[4][k] += F2[0+START2][i][k]    * dx[1] * dx[3];
		for(k=0;k<NPR;k++) sum_boundary_flux[5][k] += F2[0+START2+N2][i][k] * dx[1] * dx[3];
	}
}

