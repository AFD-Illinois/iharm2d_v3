
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
double a, R0, Rin, Rout, hslope ;

/* mandatory routines */
/* insert covariant (line element) form of metric here */

void gcov_func(double *X, double gcov[][NDIM])
{

	/* transform x1, x2 to boyer-lindquist coordinates */
	double r,th;
	BL_coord(X, &r, &th);

	/** begin set metric in BL coordinates **/
	double sth,cth,s2,rho2 ;
	cth = cos(th);
	sth = sin(th);
	s2 = sth * sth;
	rho2 = r * r + a * a * cth * cth;

	int j, k;
	DLOOP gcov[j][k] = 0.;

	gcov[TT][TT] = (-1. + 2. * r / rho2) ;
	gcov[TT][1] = (2. * r / rho2) ;
	gcov[TT][3] = (-2. * a * r * s2 / rho2) ;

	gcov[1][TT] = gcov[TT][1];
	gcov[1][1] = (1. + 2. * r / rho2) ;
	gcov[1][3] = (-a * s2 * (1. + 2. * r / rho2)) ;

	gcov[2][2] = rho2 ;

	gcov[3][TT] = gcov[TT][3];
	gcov[3][1] = gcov[1][3];
	gcov[3][3] = s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2)) ;
	/** end set metric in BL coordinates **/

	/** transform to modified Kerr-Schild coordinates **/
	/* because each coordinate is independently stretched, 
	   transformation matrix is diagonal **/
	double dXpdX[NDIM];
        dXpdX[0] = 1. ;
	dXpdX[1] = r - R0;
	dXpdX[2] = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
	dXpdX[3] = 1.;

	for(j=0;j<NDIM;j++) 
	for(k=0;k<NDIM;k++) 
		gcov[j][k] *= dXpdX[j]*dXpdX[k] ;

}

/* some coordinate parameters, grid location, dx */
void set_coordinates()
{

	/* coordinate parameters */
	a = 0.9375;
	R0 = 0.0;
	Rin = 0.98 * (1. + sqrt(1. - a*a)) ;
	Rout = 40.;
	hslope = 0.3;

	startx[1] = log(Rin - R0);
	startx[2] = 0.;
	startx[3] = 0.;

	dx[1] = log((Rout - R0) / (Rin - R0)) / N1;
	dx[2] = 1. / N2;
	dx[3] = 2.*M_PI;	/* used to evaluate integrals */

	/* report! */
	fprintf(stderr,"set_coordinates: \n") ;
	fprintf(stderr,"a: %g, R0: %g, Rin: %g, Rout: %g, hslope: %g \n",
			a, R0, Rin, Rout, hslope) ;

}

/** end mandatory routines **/

double get_spin()
{
	return(a);
}

double BL_gdet_func(double r, double th);
void BL_gcov_func(double r, double th, double gcov[][NDIM]) ;
void BL_gcon_func(double r, double th, double gcon[][NDIM]) ;

/* x1,x2 -> boyer-lindquist r,theta */
void BL_coord(double *X, double *r, double *th)
{

	*r = exp(X[1]) + R0;
	*th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);

	// avoid singularity at polar axis
	if (fabs(*th) < SMALL) {
		if ((*th) >= 0)
			*th = SMALL;
		if ((*th) < 0)
			*th = -SMALL;
	}
	if (fabs(M_PI - (*th)) < SMALL) {
		if ((*th) >= M_PI)
			*th = M_PI + SMALL;
		if ((*th) < M_PI)
			*th = M_PI - SMALL;
	}

	return;
}

#define RHOMINLIMIT (1.e-20)
#define UUMINLIMIT  (1.e-20)
#define RHOMIN  (1.e-5)
#define UUMIN   (1.e-7)

void find_model_limits(int i, int j, double *rhoflr, double *uuflr)
{
	double X[NDIM];
	double r,th,rhoscal,uuscal ;

        coord(i, j, CENT, X);
        BL_coord(X, &r, &th);

        /* this is a "natural" scaling for inflow problems */
        rhoscal = pow(r, -1.5);
        uuscal = rhoscal/r ;

        /* set floors */
        *rhoflr = RHOMIN * rhoscal;
        *uuflr = UUMIN * uuscal;
        if (*rhoflr < RHOMINLIMIT) *rhoflr = RHOMINLIMIT;
        if (*uuflr < UUMINLIMIT) *uuflr = UUMINLIMIT;
        /* end set floors */

}

#undef UUMIN
#undef RHOMIN
#undef UUMINLIMIT
#undef RHOMINLIMIT

/* 
 this starts w/ BL 4-velocity and
 converts to 3-velocities in modified
 Kerr-Schild coordinates 
*/

void BL_4vel_to_prim(double ur, double uh, double up, double *X, double *pr)
{
	double r, th, ucon[NDIM], tmp[NDIM];
	double AA, BB, CC, discr;
	double alpha, gamma, beta[NDIM];
	int j, k;

	struct of_geom BL_geom;
	BL_coord(X, &r, &th);
	BL_geom.g = BL_gdet_func(r, th);
	BL_gcov_func(r, th, BL_geom.gcov);
	BL_gcon_func(r, th, BL_geom.gcon);

	ucon[1] = ur;
	ucon[2] = uh;
	ucon[3] = up;

	/* find ucon[TT] */
	AA = BL_geom.gcov[TT][TT];
	BB = 2. * (BL_geom.gcov[TT][1] * ucon[1] +
		   BL_geom.gcov[TT][2] * ucon[2] +
		   BL_geom.gcov[TT][3] * ucon[3]);
	CC = 1. +
	    BL_geom.gcov[1][1] * ucon[1] * ucon[1] +
	    BL_geom.gcov[2][2] * ucon[2] * ucon[2] +
	    BL_geom.gcov[3][3] * ucon[3] * ucon[3] +
	    2. * (BL_geom.gcov[1][2] * ucon[1] * ucon[2] +
		  BL_geom.gcov[1][3] * ucon[1] * ucon[3] +
		  BL_geom.gcov[2][3] * ucon[2] * ucon[3]);

	discr = BB * BB - 4. * AA * CC;
	ucon[TT] = (-BB - sqrt(discr)) / (2. * AA);
	/* now we've got ucon in BL coords */

	/* transform to Kerr-Schild */
	/* make transform matrix */
	double BL_to_KS[NDIM][NDIM];
	DLOOP BL_to_KS[j][k] = 0.;
	DLOOPA BL_to_KS[j][j] = 1.;
	BL_to_KS[0][1] = 2. * r / (r * r - 2. * r + a * a);
	BL_to_KS[3][1] = a / (r * r - 2. * r + a * a);

	/* transform ucon */
	DLOOPA tmp[j] = 0.;
	DLOOP tmp[j] += BL_to_KS[j][k] * ucon[k];
	DLOOPA ucon[j] = tmp[j];
	/* now we've got ucon in KS coords */

	/* transform to KS' coords */
	double KS_to_KSp[NDIM];
	KS_to_KSp[0] = 1. ;
	KS_to_KSp[1] = (1. / (r - R0));
	KS_to_KSp[2] = (1. / (M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2])));
	KS_to_KSp[3] = 1.;
	DLOOPA ucon[j] *= KS_to_KSp[j];

	/* now solve for v-- we can use the same u^t because
	 * it didn't change under KS -> KS' */
	struct of_geom KSp_geom;
	gcov_func(X, KSp_geom.gcov) ;
	gcon_func(KSp_geom.gcov, KSp_geom.gcon) ;
	alpha = 1.0 / sqrt( -KSp_geom.gcon[0][0] ) ;
	gamma = ucon[TT] * alpha;

	/* get shift = -beta */
	beta[1] = alpha * alpha * KSp_geom.gcon[0][1];
	beta[2] = alpha * alpha * KSp_geom.gcon[0][2];
	beta[3] = alpha * alpha * KSp_geom.gcon[0][3];

	/* subtract off shift */
	pr[U1] = ucon[1] + beta[1] * gamma / alpha;
	pr[U2] = ucon[2] + beta[2] * gamma / alpha;
	pr[U3] = ucon[3] + beta[3] * gamma / alpha;

	/* done! */
}

/* some functions to set BL geometry.  Needed only
   in the initial conditions */

double BL_gdet_func(double r, double th)
{
	double a2, r2;

	a2 = a * a;
	r2 = r * r;
	return (r * r * fabs(sin(th)) *
		(1. + 0.5 * (a2 / r2) * (1. + cos(2. * th)))
	    );
}

void BL_gcov_func(double r, double th, double gcov[][NDIM])
{
	int j, k;
	double sth, cth, s2, a2, r2, DD, mu;

	DLOOP gcov[j][k] = 0.;

	sth = fabs(sin(th));
	s2 = sth * sth;
	cth = cos(th);
	a2 = a * a;
	r2 = r * r;
	DD = 1. - 2. / r + a2 / r2;
	mu = 1. + a2 * cth * cth / r2;

	gcov[TT][TT] = -(1. - 2. / (r * mu));
	gcov[TT][3] = -2. * a * s2 / (r * mu);
	gcov[3][TT] = gcov[TT][3];
	gcov[1][1] = mu / DD;
	gcov[2][2] = r2 * mu;
	gcov[3][3] =
	    r2 * sth * sth * (1. + a2 / r2 + 2. * a2 * s2 / (r2 * r * mu));

}

void BL_gcon_func(double r, double th, double gcon[][NDIM])
{
	int j, k;
	double sth, cth, a2, r2, r3, DD, mu;

	DLOOP gcon[j][k] = 0.;

	sth = sin(th);
	cth = cos(th);

	a2 = a * a;
	r2 = r * r;
	r3 = r2 * r;
	DD = 1. - 2. / r + a2 / r2;
	mu = 1. + a2 * cth * cth / r2;

	gcon[TT][TT] = -1. - 2. * (1. + a2 / r2) / (r * DD * mu);
	gcon[TT][3] = -2. * a / (r3 * DD * mu);
	gcon[3][TT] = gcon[TT][3];
	gcon[1][1] = DD / mu;
	gcon[2][2] = 1. / (r2 * mu);
	gcon[3][3] = (1. - 2. / (r * mu)) / (r2 * sth * sth * DD);

}
