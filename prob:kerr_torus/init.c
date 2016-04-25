
/*
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with physical but low density atmosphere
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"

/* boyer-lindquist related functions; needed to
   set up disk in BL coord then transform to Kerr-Schild prime coordinates */
void   BL_4vel_to_prim(double ur,double uh,double up, double *X, double *pr);

/* returns specific angular momentum of fishbone-moncrief torus */
double lfish_calc(double rmax) ;

double a ;

void init()
{
	int i, j;
	double r, th, sth, cth;
	double ur, uh, up, u, rho;
	double X[NDIM];
	struct of_geom *geom;

	/* for disk interior */
	double l, rin, lnh, expm2chi, up1;
	double DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
	double hm1;

	/* for magnetic field */
	static double A[N1 + 2*NG][N2 + 2*NG];
	double rho_av, rhomax, umax, beta, bsq, bsq_max, norm, q,
	    beta_act;
	double rmax ;

	/* physical parameters */
	gam = 13. / 9.;  /* adiabatic index: relativistic electrons; nonrel. ions */

	/* set up grid functions */
	set_geometry();

	/* disk parameters (use fishbone.m to select new solutions) */
	a = get_spin();  /* set_geometry must be called first */
	rin = 6.;        /* inner edge of disk */
	rmax = 12.;      /* pressure maximum */
	l = lfish_calc(rmax);   /* solve for angular momentum of fishbone-moncreif torus */
	beta = 1.e2;     /* plasma beta for field strength */

	/* end of simulation.  time unit = G M/c^3 */
	t = 0.;
	tf = 2000.0;
	dtsave = dt = 1.e-8;
	cour = 0.4;

	/* output frequencies */
	DTd = 50.;		/* dumping frequency, in units of M */
	DTl = 2.0;		/* logfile frequency, in units of M */
	DTi = 2.0;		/* image file frequ., in units of M */
	DTr = 512;		/* restart file frequ., in timesteps */

	/* start diagnostic counters */
	nstep = 0;
	dump_cnt = 0;
	image_cnt = 0;
	rdump_cnt = 0;

	rhomax = 0.;
	umax = 0.;
	ZLOOP {
		coord(i, j, CENT, X);
		BL_coord(X, &r, &th);

		sth = sin(th);
		cth = cos(th);

		/* calculate lnh */
		DD = r * r - 2. * r + a * a;
		AA = (r * r + a * a) * (r * r + a * a) -
		    DD * a * a * sth * sth;
		SS = r * r + a * a * cth * cth;

		thin = M_PI / 2.;
		sthin = sin(thin);
		cthin = cos(thin);
		DDin = rin * rin - 2. * rin + a * a;
		AAin = (rin * rin + a * a) * (rin * rin + a * a)
		    - DDin * a * a * sthin * sthin;
		SSin = rin * rin + a * a * cthin * cthin;

		if (r >= rin) {
			lnh =
			    0.5 *
			    log((1. +
				 sqrt(1. +
				      4. * (l * l * SS * SS) * DD / (AA *
								     sth *
								     AA *
								     sth)))
				/ (SS * DD / AA))
			    - 0.5 * sqrt(1. +
					 4. * (l * l * SS * SS) * DD /
					 (AA * AA * sth * sth))
			    - 2. * a * r * l / AA -
			    (0.5 *
			     log((1. +
				  sqrt(1. +
				       4. * (l * l * SSin * SSin) * DDin /
				       (AAin * AAin * sthin * sthin))) /
				 (SSin * DDin / AAin))
			     - 0.5 * sqrt(1. +
					  4. * (l * l * SSin * SSin) *
					  DDin / (AAin * AAin * sthin *
						  sthin))
			     - 2. * a * rin * l / AAin);
		} else
			lnh = 1.;


		/* regions outside torus */
		if (lnh < 0. || r < rin) {
			/* nominal values; real value set by fixup */
			rho = SMALL;
			u = SMALL;

			/* these normal observer values are demonstrably physical
			   for all points on grid */
			ur = 0.;
			uh = 0.;
			up = 0.;

			p[i][j][RHO] = rho;
			p[i][j][UU] = u;
			p[i][j][U1] = ur;
			p[i][j][U2] = uh;
			p[i][j][U3] = up;
		}
		/* region inside magnetized torus; u^i is calculated in
		 * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
		 * so it needs to be transformed at the end */
		else {
			hm1 = exp(lnh) - 1.;
			rho = pow(hm1 * (gam - 1.) / gam, 1. / (gam - 1.));
			u = pow(rho, gam) / (gam - 1.);
			ur = 0.;
			uh = 0.;

			/* calculate u^phi */
			expm2chi = SS * SS * DD / (AA * AA * sth * sth);
			up1 =
			    sqrt((-1. +
				  sqrt(1. + 4. * l * l * expm2chi)) / 2.);
			up = 2. * a * r * sqrt(1. +
					       up1 * up1) / sqrt(AA * SS *
								 DD) +
			    sqrt(SS / AA) * up1 / sth;


			p[i][j][RHO] = rho;
			if (rho > rhomax) rhomax = rho;
			p[i][j][UU] = u * (1. + 4.e-2 * (ran_uniform() - 0.5));
			if (u > umax && r > rin) umax = u;
			p[i][j][U1] = ur;
			p[i][j][U2] = uh;
			p[i][j][U3] = up;

			/* convert from 4-vel to 3-vel */
			BL_4vel_to_prim(ur,uh,up, X, p[i][j]);
		}

		p[i][j][B1] = 0.;
		p[i][j][B2] = 0.;
		p[i][j][B3] = 0.;
	}

	/* Normalize the densities so that max(rho) = 1 */
	fprintf(stderr, "rhomax: %g\n", rhomax);
	ZSLOOP(0, N1 - 1, 0, N2 - 1) {
		p[i][j][RHO] /= rhomax;
		p[i][j][UU] /= rhomax;
	}
	umax /= rhomax;
	rhomax = 1.;
	fixup(p);
	bound_prim(p);

	/* first find corner-centered 3-component of vector potential */
	ZSLOOP(0, N1, 0, N2) A[i][j] = 0.;
	ZSLOOP(0, N1, 0, N2) {
		/* vertical field through disk */
		/*
		   coord(i,j,CORN,X) ;
		   BL_coord(X,&r,&th) ;

		   A[i][j] = 0.5*r*sin(th) ;
		 */

		/* field follows isodensity contours */
		rho_av = 0.25 * (p[i][j][RHO] + p[i - 1][j][RHO] +
				 p[i][j - 1][RHO] + p[i - 1][j - 1][RHO])
				* (1. + 0.0*(ran_uniform()-0.5)) 
				;

		q = rho_av / rhomax - 0.2;
		if (q > 0.) A[i][j] = q;
	}

	/* now differentiate to find cell-centered B,
	   and begin normalization */
	bsq_max = 0.;
	ZLOOP {
		geom = get_geometry(i, j, CENT) ;

		/* flux-ct */
		p[i][j][B1] = -(A[i][j] - A[i][j + 1]
				+ A[i + 1][j] - A[i + 1][j + 1]) / 
				(2. * dx[2] * geom->g);
		p[i][j][B2] = (A[i][j] + A[i][j + 1]
			       - A[i + 1][j] - A[i + 1][j + 1]) / 
			       (2. * dx[1] * geom->g);

		p[i][j][B3] = 0.;

		bsq = bsq_calc(p[i][j], geom);
		if (bsq > bsq_max) bsq_max = bsq;
	}
	fprintf(stderr, "initial bsq_max: %g\n", bsq_max);

	/* finally, normalize to set field strength */
	beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
	fprintf(stderr, "initial beta: %g (should be %g)\n", beta_act,
		beta);
	norm = sqrt(beta_act / beta);
	bsq_max = 0.;
	ZLOOP {
		p[i][j][B1] *= norm;
		p[i][j][B2] *= norm;

		geom = get_geometry(i, j, CENT) ;
		bsq = bsq_calc(p[i][j], geom);
		if (bsq > bsq_max)
			bsq_max = bsq;
	}
	beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
	fprintf(stderr, "final beta: %g (should be %g)\n", beta_act, beta);

	/* enforce boundary conditions */
	fixup(p);
	bound_prim(p);

}

double lfish_calc(double r)
{
	return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
		 ((-2. * a * r *
		   (pow(a, 2) - 2. * a * sqrt(r) +
		    pow(r,
			2))) / sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
		  ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) + pow(a, 2) * 
		  (2. + r))) / sqrt(1 + (2. * a) / pow (r, 1.5) - 3. / r)))
		/ (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
		   (pow(a, 2) + (-2. + r) * r))
	    );
}

