
/*

Kelvin-Helmholtz instability problem

follows initial conditions from Lecoanet et al. 2015,
	MNRAS 455, 4274.

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

	/* end of simulation */
	t = 0.;
	tf = 8.0;
	dtsave = dt = 1.e-4;
	cour = 0.9;

	/* output frequencies */
	DTd = tf/10.;		/* dumping frequency, in units of M */
	DTl = tf/400.;		/* logfile frequency, in units of M */
	DTi = tf/400.;		/* image file frequ., in units of M */
	DTr = 4096;		/* restart file frequ., in timesteps */

	/* start diagnostic counters */
	nstep = 0;
	dump_cnt = 0;
	image_cnt = 0;
	rdump_cnt = 0;

	/* here are the dimensionless parameters we want to control */
	/* zones are designated H (high density) and L (low density) */
	/* first specify high side parameters */
	double betaH = 1.0;
	double csH = 0.1;
	double vxH = 0.;
	double rhoH = 1.0 ;

	/* which implies... */
	double csqH = csH*csH ;          /* aux quants */ 
	double fH = gam - 1. - csqH ; 
	double uH = rhoH*csqH/(gam*fH);
	double bsqH = 2.*rhoH*csqH*(gam - 1.)/(betaH*gam*fH) ;

	/* which in turn implies... */
	double pH = (gam - 1.)*uH ;
	double pHtot = 0.5*bsqH + pH ;

	/* equilibrium condition */
	double pLtot = pHtot ;

	/* then low side */
	double csL = 0.5;
	double sigmaL = 1.0;
	double vxL = 0.5;
	double csqL = csL*csL ;
	double fL = gam - 1. - csqL ; /* aux */
	double rhoL = rhoH*csqH*(gam - 1.)*fL*(1. + betaH)/
		(betaH*fH*(sigmaL*gam*(gam - 1.) + csqL*(gam - 1. - gam*sigmaL))) ;
	double bsqL = -2.*rhoH*csqH*(1. + betaH)*fL*(gam - 1.)*sigmaL/
		(betaH*fH*(csqL*(gam*(sigmaL - 1.) + 1.) - sigmaL*gam*(gam - 1.))) ;
	double uL = -rhoH*csqH*csqL*(1. + betaH)*(gam - 1.)/
		(betaH*fH*gam*(csqL*(1. + gam*(sigmaL - 1.)) - (gam - 1.)*gam*sigmaL)) ;

	/* report! */
	fprintf(stderr,"high side rho, u, bsq: %g %g %g\n",rhoH, uH, bsqH) ;
	fprintf(stderr,"low side rho, u, bsq: %g %g %g\n",rhoL, uL, bsqL) ;

	/* convert velocity to four-velocity */
	double uxL = vxL/sqrt(1. - vxL*vxL) ;
	double uxH = vxH/sqrt(1. - vxH*vxH) ;

	/* location of discontinuities */
	double y1 = 0.5;
	double y2 = 1.5;
	double lL = y2 - y1 ;
	double lH = 2. - lL ;

	/* amplitude of noise */
	double eps = 0.05;  

	double x,y,phase;
	double dphase = 0.3*M_PI/2.;
	double q ;
	double bH = sqrt(bsqH);
	double bL = sqrt(bsqL);
	double sp = 0.;

	init_rng(1);

	ZLOOP {
		coord(i, j, CENT, X);

		x = X[1];
		y = X[2];

		if (y > y1 && y < y2) {  // zone L
			p[i][j][RHO] = rhoL ;
			p[i][j][UU] = uL * (1. + eps*(ran_uniform()-0.5));
			p[i][j][U1] = uxL;
			p[i][j][U2] = 0.;
			p[i][j][U3] = 0.;
			phase = M_PI*(y - y1)/lL + dphase;
			phase = 0.;
			sp = sin(phase);
			q = 1./sqrt(1. - vxL*vxL*sp*sp);
			p[i][j][B1] = bL*sp*q ;
			p[i][j][B2] = 0.;
			p[i][j][B3] = bL*cos(phase)*q;
		}
		else {  // zone H
			p[i][j][RHO] = rhoH ;
			p[i][j][UU] = uH * (1. + eps*(ran_uniform()-0.5));
			p[i][j][U1] = uxH;
			p[i][j][U2] = 0.;
			p[i][j][U3] = 0.;
			phase = (y > y2) ? M_PI*(y - y2)/lH : M_PI*(y + 1. - y2)/lH + dphase ;
			phase = 0. ;
			sp = sin(phase);
			q = 1./sqrt(1. - vxH*vxH*sp*sp);
			p[i][j][B1] = bH*sp*q ;
			p[i][j][B2] = 0.;
			p[i][j][B3] = bH*cos(phase)*q;
		}
	}

	/* enforce boundary conditions */
	fixup(p);
	bound_prim(p);

}

