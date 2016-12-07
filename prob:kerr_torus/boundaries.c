
#include "decs.h"

void bound_prim(double prim[][N2 + 2*NG][NPR])
{
	int i, j, k, iz, jrefl;
	double rescale_fac;
	void inflow_check(double *pr, int ii, int jj, int type);

	/* inner r boundary condition: copy last value off grid */
	JSLOOP(0,N2-1) {
		iz = 0+START1 ;
		ISLOOP(-NG,-1) {
			PLOOP prim[i][j][k] = prim[iz][j][k];
			pflag[i][j] = pflag[iz][j];

			rescale_fac = ggeom[iz][j][CENT].g/ggeom[i][j][CENT].g ;
			prim[i][j][B1] *= rescale_fac ;
			prim[i][j][B2] *= rescale_fac ;
			prim[i][j][B3] *= rescale_fac ;
		}
	}

	/* outer r BC: copy last value off grid */
	JSLOOP(0,N2-1) {
		iz = N1-1+START1 ;
		ISLOOP(N1,N1-1+NG) {
			PLOOP prim[i][j][k] = prim[iz][j][k];
			pflag[i][j] = pflag[iz][j];

			rescale_fac = ggeom[iz][j][CENT].g/ggeom[i][j][CENT].g ;
			prim[i][j][B1] *= rescale_fac ;
			prim[i][j][B2] *= rescale_fac ;
			prim[i][j][B3] *= rescale_fac ;
		}
	}

	/* loop over ghost zones to zero inflow at the inner boundary */
	ISLOOP(-NG,-1) 
		JSLOOP(-NG,N2-1+NG) 
			inflow_check(prim[i][j], i, j, 0);

	/* loop over ghost zones to zero inflow at the outer boundary */
	ISLOOP(N1,N1-1+NG)
		JSLOOP(-NG,N2-1+NG) 
			inflow_check(prim[i][j], i, j, 1);

	/* polar BCs */
	ISLOOP(-NG,N1-1+NG) {
		JSLOOP(-NG,-1) {
			jrefl = -j+2*NG-1 ;
			PLOOP prim[i][j][k] = prim[i][jrefl][k];
			pflag[i][j] = pflag[i][jrefl];
		}
		JSLOOP(N2,N2-1+NG) {
			jrefl = -j+2*(N2+START2)-1 ;
			PLOOP prim[i][j][k] = prim[i][jrefl][k];
			pflag[i][j] = pflag[i][jrefl];
		}
	}

	/* enforce b and u antisymmetry at the poles */
	/* notice that 3-components do not change sign at poles */
	ISLOOP(-NG,N1-1+NG) {
		JSLOOP(-NG,-1) {
			prim[i][j][U2] *= -1.;
			prim[i][j][B2] *= -1.;
		}
		JSLOOP(N2,N2-1+NG) {
			prim[i][j][U2] *= -1.;
			prim[i][j][B2] *= -1.;
		}
	}

}

void inflow_check(double *pr, int ii, int jj, int type)
{
	struct of_geom *geom;
	double ucon[NDIM];
	int j, k;
	double alpha, beta1, gamma, vsq;

	geom = get_geometry(ii, jj, CENT);
	ucon_calc(pr, geom, ucon);

	/* inflow condition is u^1 directed toward the grid */
	if (((ucon[1] > 0.) && (type == 0))
	    || ((ucon[1] < 0.) && (type == 1))) {

	    	/* 
		   BC idea: 
		   if u^1 is onto the grid;
		   zero radial component of four-velocity
		   without changing any other "component" of the
		   velocity with respect to the normal observer
		*/

		/* find gamma and remove it from primitives */
		if (mhd_gamma_calc(pr, geom, &gamma)) {
			fflush(stderr);
			fprintf(stderr,
				"\ninflow_check(): gamma failure \n");
			fflush(stderr);
			fail(FAIL_GAMMA);
		}
		pr[U1] /= gamma;
		pr[U2] /= gamma;
		pr[U3] /= gamma;
		alpha = geom->alpha ;
		beta1 = geom->gcon[0][1] * alpha * alpha;

		/* reset radial velocity so radial 4-velocity
		 * is zero */
		pr[U1] = beta1 / alpha;

		/* now find new gamma and put it back in */
		vsq = 0.;
		SLOOP vsq +=
		    geom->gcov[j][k] * pr[U1 + j - 1] * pr[U1 + k - 1];
		if (fabs(vsq) < 1.e-13)
			vsq = 1.e-13;
		if (vsq >= 1.) {
			vsq = 1. - 1. / (GAMMAMAX * GAMMAMAX);
		}
		gamma = 1. / sqrt(1. - vsq);
		pr[U1] *= gamma;
		pr[U2] *= gamma;
		pr[U3] *= gamma;

		/* done */
	} else
		return;

}
