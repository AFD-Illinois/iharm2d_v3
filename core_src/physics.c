
#include "decs.h"

/*

  primtoflux():
  ---------
   --  calculate fluxes in direction dir, 
        
*/

void primtoflux(double *pr, struct of_state *q, int dir,
		struct of_geom *geom, double *flux)
{
	int k;
	double mhd[NDIM];

	/* particle number flux */
	flux[RHO] = pr[RHO] * q->ucon[dir];

	mhd_calc(pr, dir, q, mhd);

	/* MHD stress-energy tensor w/ first index up, 
	 * second index down. */
	flux[UU] = mhd[0];
	flux[U1] = mhd[1];
	flux[U2] = mhd[2];
	flux[U3] = mhd[3];

	/* dual of Maxwell tensor */
	flux[B1] = q->bcon[1] * q->ucon[dir] - q->bcon[dir] * q->ucon[1];
	flux[B2] = q->bcon[2] * q->ucon[dir] - q->bcon[dir] * q->ucon[2];
	flux[B3] = q->bcon[3] * q->ucon[dir] - q->bcon[dir] * q->ucon[3];

	PLOOP flux[k] *= geom->g;

}

/* calculate "conserved" quantities; provided strictly for
 * historical reasons */
void primtoU(double *pr, struct of_state *q, struct of_geom *geom,
	     double *U)
{

	primtoflux(pr, q, 0, geom, U);
	return;
}

/* calculate magnetic field four-vector */
void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon)
{
	int j;

	bcon[0] = pr[B1] * ucov[1] + pr[B2] * ucov[2] + pr[B3] * ucov[3];
	for (j = 1; j < 4; j++)
		bcon[j] = (pr[B1 - 1 + j] + bcon[0] * ucon[j]) / ucon[0];

	return;
}

/* MHD stress tensor, with first index up, second index down */
/* notice that factor of sqrt(4 pi) is absorbed into defn of b, which
   implies that it must be reinserted to get B in Gauss -- this
   is origin of the 4 pi in Bunit = c sqrt[4 Pi rhounit] */
void mhd_calc(double *pr, int dir, struct of_state *q, double *mhd)
{
	int j;
	double r, u, P, w, bsq, eta, ptot;

	r = pr[RHO];
	u = pr[UU];
	P = (gam - 1.) * u;
	w = P + r + u;
	bsq = dot(q->bcon, q->bcov);
	eta = w + bsq;
	ptot = P + 0.5 * bsq;

	/* single row of mhd stress tensor, 
	 * first index up, second index down */
	DLOOPA mhd[j] = eta * q->ucon[dir] * q->ucov[j]
	    + ptot * delta(dir, j) - q->bcon[dir] * q->bcov[j];

}

/* add in source terms to equations of motion */
void source(double *ph, struct of_geom *geom, int ii, int jj, double *dU,
	    double Dt)
{
	double mhd[NDIM][NDIM];
	int j, k;
	struct of_state q;

	get_state(ph, geom, &q);
	mhd_calc(ph, 0, &q, mhd[0]);
	mhd_calc(ph, 1, &q, mhd[1]);
	mhd_calc(ph, 2, &q, mhd[2]);
	mhd_calc(ph, 3, &q, mhd[3]);

	/* contract mhd stress tensor with connection */
	PLOOP dU[k] = 0.;
	DLOOP {
		dU[UU] += mhd[j][k] * conn[ii][jj][k][0][j];
		dU[U1] += mhd[j][k] * conn[ii][jj][k][1][j];
		dU[U2] += mhd[j][k] * conn[ii][jj][k][2][j];
		dU[U3] += mhd[j][k] * conn[ii][jj][k][3][j];
	}
	
	PLOOP dU[k] *= geom->g;

#if ( WIND )
	/* uncomment to put new source terms in */
	wind_source(ph, geom, &q, ii, jj, dU);
#endif

	/* done! */
}

/* returns b^2 (i.e., twice magnetic pressure) */
double bsq_calc(double *pr, struct of_geom *geom)
{
	struct of_state q;

	get_state(pr, geom, &q);
	return (dot(q.bcon, q.bcov));
}

/* find ucon, ucov, bcon, bcov from primitive variables */
void get_state(double *pr, struct of_geom *geom, struct of_state *q)
{

	/* get ucon */
	ucon_calc(pr, geom, q->ucon);
	lower(q->ucon, geom, q->ucov);
	bcon_calc(pr, q->ucon, q->ucov, q->bcon);
	lower(q->bcon, geom, q->bcov);

	return;
}

/* find contravariant four-velocity */
void ucon_calc(double *pr, struct of_geom *geom, double *ucon)
{
	double alpha, gamma;
	double beta[NDIM];
	int j;

	alpha = geom->alpha ;
	SLOOPA beta[j] = geom->gcon[0][j] * alpha * alpha;

	/*
	if (mhd_gamma_calc(pr, geom, &gamma)) {
		fflush(stderr);
		fprintf(stderr, "\nucon_calc(): gamma failure \n");
		fflush(stderr);
		fail(FAIL_GAMMA);
	}
	*/

	mhd_gamma_calc(pr, geom, &gamma) ;
	ucon[0] = gamma / alpha;
	SLOOPA ucon[j] = pr[U1 + j - 1] - gamma * beta[j] / alpha;

	return;
}

/* find gamma-factor wrt normal observer */
int mhd_gamma_calc(double *pr, struct of_geom *geom, double *gamma)
{
	double qsq;

	qsq = geom->gcov[1][1] * pr[U1] * pr[U1]
	    + geom->gcov[2][2] * pr[U2] * pr[U2]
	    + geom->gcov[3][3] * pr[U3] * pr[U3]
	    + 2. * (geom->gcov[1][2] * pr[U1] * pr[U2]
		    + geom->gcov[1][3] * pr[U1] * pr[U3]
		    + geom->gcov[2][3] * pr[U2] * pr[U3]);

	/*
	if (qsq < 0.) {
		if (fabs(qsq) > 1.E-10) {	// then assume not just machine precision
			fprintf(stderr,
				"gamma_calc():  failed: i,j,qsq = %d %d %28.18e \n",
				icurr, jcurr, qsq);
			fprintf(stderr,
				"v[1-3] = %28.18e %28.18e %28.18e  \n",
				pr[U1], pr[U2], pr[U3]);
			*gamma = 1.;
			return (1);
		} else
			qsq = 1.E-10;	// set floor
	}
	*/

	*gamma = sqrt(1. + qsq);

	return (0);
}

/*  
 * VCHAR():
 * 
 * calculate components of magnetosonic velocity 
 * corresponding to primitive variables p 
 *
 * cfg 7-10-01
 * 
 */

void mhd_vchar(double *pr, struct of_state *q, struct of_geom *geom, int js,
	   double *vmax, double *vmin)
{
	double discr, vp, vm, bsq, EE, EF, va2, cs2, cms2, rho, u;
	double Acov[NDIM], Bcov[NDIM], Acon[NDIM], Bcon[NDIM];
	double Asq, Bsq, Au, Bu, AB, Au2, Bu2, AuBu, A, B, C;
	int j;


	DLOOPA Acov[j] = 0.;
	Acov[js] = 1.;
	raise(Acov, geom, Acon);

	DLOOPA Bcov[j] = 0.;
	Bcov[0] = 1.;
	raise(Bcov, geom, Bcon);

	/* find fast magnetosonic speed */
	bsq = dot(q->bcon, q->bcov);
	rho = fabs(pr[RHO]);
	u = fabs(pr[UU]);
	EF = rho + gam * u;
	EE = bsq + EF;
	va2 = bsq / EE;
	cs2 = gam * (gam - 1.) * u / EF;
	cms2 = cs2 + va2 - cs2 * va2;	/* and there it is... */

	//cms2 *= 1.1 ;

	/* check on it! */
	/*
	if (cms2 < 0.) {
		fprintf(stderr,"\n\ncms2: %g %g %g\n\n",gam,u,EF) ;
		fail(FAIL_COEFF_NEG);
		cms2 = SMALL;
	}
	if (cms2 > 1.) {
		fail(FAIL_COEFF_SUP);
		cms2 = 1.;
	}
	*/

	/* now require that speed of wave measured by observer 
	   q->ucon is cms2 */
	Asq = dot(Acon, Acov);
	Bsq = dot(Bcon, Bcov);
	Au = dot(Acov, q->ucon);
	Bu = dot(Bcov, q->ucon);
	AB = dot(Acon, Bcov);
	Au2 = Au * Au;
	Bu2 = Bu * Bu;
	AuBu = Au * Bu;

	A = Bu2 - (Bsq + Bu2) * cms2;
	B = 2. * (AuBu - (AB + AuBu) * cms2);
	C = Au2 - (Asq + Au2) * cms2;

	discr = B * B - 4. * A * C;
	/*
	if ((discr < 0.0) && (discr > -1.e-10))
		discr = 0.0;
	else if (discr < -1.e-10) {
		fprintf(stderr, "\n\t %g %g %g %g %g\n", A, B, C, discr,
			cms2);
		fprintf(stderr, "\n\t q->ucon: %g %g %g %g\n", q->ucon[0],
			q->ucon[1], q->ucon[2], q->ucon[3]);
		fprintf(stderr, "\n\t q->bcon: %g %g %g %g\n", q->bcon[0],
			q->bcon[1], q->bcon[2], q->bcon[3]);
		fprintf(stderr, "\n\t Acon: %g %g %g %g\n", Acon[0],
			Acon[1], Acon[2], Acon[3]);
		fprintf(stderr, "\n\t Bcon: %g %g %g %g\n", Bcon[0],
			Bcon[1], Bcon[2], Bcon[3]);
		fail(FAIL_VCHAR_DISCR);
		discr = 0.;
	}
	*/

	discr = sqrt(discr);
	vp = -(-B + discr) / (2. * A);
	vm = -(-B - discr) / (2. * A);

	*vmax = fmax(vp,vm) ;
	*vmin = fmin(vp,vm) ;

	/*
	if (vp > vm) {
		*vmax = vp;
		*vmin = vm;
	} else {
		*vmax = vm;
		*vmin = vp;
	}
	*/

	/*
	fprintf(stderr,"vmax, vmin: %g %g\n",*vmax,*vmin) ;
	*/

	return;
}

/*  
   estimate div . b, using stencil that is consistent
   with emf 
*/
double divb_calc(int i, int j)
{
	double divb ;

	/* divb flux-ct defn; corner-centered */
	divb = 0.5*(
		  p[i][j][B1]    *ggeom[i][j][CENT].g
		+ p[i][j-1][B1]  *ggeom[i][j-1][CENT].g
		- p[i-1][j][B1]  *ggeom[i-1][j][CENT].g
		- p[i-1][j-1][B1]*ggeom[i-1][j-1][CENT].g
		)/dx[1] +
		0.5*(
		  p[i][j][B2]    *ggeom[i][j][CENT].g
		+ p[i-1][j][B2]  *ggeom[i-1][j][CENT].g
		- p[i][j-1][B2]  *ggeom[i][j-1][CENT].g
		- p[i-1][j-1][B2]*ggeom[i-1][j-1][CENT].g
		)/dx[2];

	return(divb);
}

/* Add any additional source terms (e.g. cooling functions) */

