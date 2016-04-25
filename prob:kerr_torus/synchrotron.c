
#include "decs.h"

/* evaluate synchrotron cooling rate in fluid frame */

/* input and output are in code units */

double synchrotron_cooling_func(double rho, double u, double bsq, double r)
{

	double GNEWT, MSUN, CL, ME, MP, M_SGRA, L_UNIT, M_UNIT;
	double RHO_UNIT, BSQ_UNIT ;
	double Theta_e, ne, lambda_cgs, lambda_sim;
	double expr, em;
	double Y_crit, f, lgt, c0, EE, nuc, nus, L ;
	static int firstc = 1;
	double find_lgY_crit(double lgt) ;

	/* physical constants */
	GNEWT = 6.67407e-8;
	MSUN = 1.989e33;
	CL = 2.997924e10;
	MP = 1.6726231e-24;
	ME = 9.10956e-28;
	EE = 4.80320680e-10;

	/** begin set units **/
	/* characteristic values for Sgra* */
	M_SGRA = 4.0e6 * MSUN;

	/* fundamental units */
	L_UNIT = GNEWT * M_SGRA / (CL * CL);
	M_UNIT = 1.e19;

	if (firstc) {
		fprintf(stderr, "M_UNIT: %g\n", M_UNIT);
		firstc = 0;
	}

	/* derived units */
	RHO_UNIT = M_UNIT / (L_UNIT * L_UNIT * L_UNIT);
	BSQ_UNIT = RHO_UNIT * CL * CL * 4. * M_PI;
	/* 4 pi necessary to get B in gauss */
	/** end set units **/

	/** begin nonself-absorbed cooling estimate **/
	/* assumes pure hydrogen plasma */
	Theta_e = (gam - 1.) * u / rho * (0.5 * MP / ME);
	ne = rho * RHO_UNIT / MP;
	bsq *= BSQ_UNIT;

	/* coefficient here is 4 e^4/(3 me^2 c^3) */
	em = 4. / 3.;
	/* this smoothly interpolates between Theta_e << 1 and 
	   Theta_e >> 1 */
	expr = pow(Theta_e, em) + pow(2. * Theta_e, 2. * em);
	lambda_cgs = 3.173475e-15 * ne * bsq * pow(expr, 1. / em);
	lambda_sim = lambda_cgs / (RHO_UNIT * CL * CL * CL / L_UNIT);
	/** end nonself-absorbed cooling estimate **/

	/** begin self-absorption correction **/
	/* don't count cooling if emitted photon mean free path is
	   smaller than L; coefficient is a parameter of the model */
	L = 0.1 * r * L_UNIT;
	nuc = EE * sqrt(bsq) / 2. / M_PI / ME / CL;
	nus = 2. / 9. * nuc * Theta_e * Theta_e;
	c0 = EE * EE * EE * sqrt(bsq) * ne / 27. / sqrt(2.) / ME / CL / CL;

	if (nus < SMALL || nuc < SMALL || c0 < SMALL * SMALL) {
		Y_crit = 0.0;
	} else {
		lgt = log(L * c0 / 2. / Theta_e / nus / nus / ME);

		/* use lookup table to find characteristic frequency */
		Y_crit = exp(find_lgY_crit(lgt));
	}

	/* find correction using analytic fitting formula */
	f = 0.5 * (exp(-Y_crit / 82.) + exp(-Y_crit / 360));

	/* reduce cooling rate */
	lambda_sim *= f;
	/** done with self-absorption correction **/

	return (lambda_sim);
}

/* advance internal energy over dt using cooling function */

void cool_down(double pi[][N2 + 2*NG][NPR], double dt)
{
	int i, j;
	double bsq, lambda, tau_cool, p_half;
	double r, th;
	double X[NDIM], ucon[NDIM];
	struct of_geom *geom;

	ZLOOP {
		geom = get_geometry(i, j, CENT) ;
		bsq = bsq_calc(pi[i][j], geom);
		ucon_calc(pi[i][j], geom, ucon);

		coord(i, j, CENT, X);
		BL_coord(X, &r, &th);

		/** take half-step **/
		lambda = synchrotron_cooling_func(pi[i][j][RHO], pi[i][j][UU], bsq, r);
		lambda /= ucon[0];	/* cooling rate wrt t rather than tau */
		tau_cool = pi[i][j][UU] / lambda;

		p_half = pi[i][j][UU] * exp(-dt * 0.5 / tau_cool);

		/** take full step **/
		lambda = synchrotron_cooling_func(pi[i][j][RHO], p_half, bsq, r);
		lambda /= ucon[0];	/* cooling rate wrt t rather than tau */
		tau_cool = p_half / lambda;

		pi[i][j][UU] *= exp(-dt / tau_cool);

		/* done! */
	}

}

/* lookup table for frequency below which cooling is
   reduced due to finite optical depth */
double find_lgY_crit(double lgt)
{

	int i;
	double x1, x2, f1, f2;
	double lgY;
	double lgY_crit[] = {
		0.626088, 0.899897, 1.16985, 1.43559, 1.6968,
		1.95314, 2.20434, 2.45013, 2.69029, 2.92463,
		3.15302, 3.37535, 3.59157, 3.80167, 4.00566,
		4.20359, 4.39556, 4.58167, 4.76205, 4.93684,
		5.10622
	};

	if (lgt < 0) {

		x1 = 0.0;
		x2 = 0.5;
		f1 = lgY_crit[1];
		f2 = lgY_crit[2];


	} else if (lgt > 10) {

		x1 = 9.5;
		x2 = 10.;
		f1 = lgY_crit[19];
		f2 = lgY_crit[20];

	} else {

		i = (int) floor((lgt) / 0.5 + 1.);	// larger that x[i]
		x1 = (i) * 0.5;
		x2 = (i + 1) * 0.5;
		f1 = lgY_crit[i];
		f2 = lgY_crit[i + 1];

	}

	lgY =
	    (f2 - f1) / (x2 - x1) * lgt + (f1 * x2 - x1 * f2) / (x2 - x1);

	return lgY;

}
