
#include "decs.h"

/* formatting strings that give the necessary precision */
#define FMT_DBL_OUT "%28.18e"
#define FMT_INT_OUT "%10d"

void dump(FILE * fp)
{
	int i, j, k;
	double divb;
	double X[NDIM];
	double vmin, vmax;
	struct of_geom *geom;
	struct of_state q;

	/***************************************************************
	  Write header information : 
	***************************************************************/

	fprintf(fp, FMT_DBL_OUT, t);
	fprintf(fp, FMT_INT_OUT, N1);
	fprintf(fp, FMT_INT_OUT, N2);
	fprintf(fp, FMT_DBL_OUT, startx[1]);
	fprintf(fp, FMT_DBL_OUT, startx[2]);
	fprintf(fp, FMT_DBL_OUT, dx[1]);
	fprintf(fp, FMT_DBL_OUT, dx[2]);
	fprintf(fp, FMT_DBL_OUT, tf);
	fprintf(fp, FMT_INT_OUT, nstep);
	fprintf(fp, FMT_DBL_OUT, gam);
	fprintf(fp, FMT_DBL_OUT, cour);
	fprintf(fp, FMT_DBL_OUT, DTd);
	fprintf(fp, FMT_DBL_OUT, DTl);
	fprintf(fp, FMT_DBL_OUT, DTi);
	fprintf(fp, FMT_INT_OUT, DTr);
	fprintf(fp, FMT_INT_OUT, dump_cnt);
	fprintf(fp, FMT_INT_OUT, image_cnt);
	fprintf(fp, FMT_INT_OUT, rdump_cnt);
	fprintf(fp, FMT_DBL_OUT, dt);

	fprintf(fp, "\n");

	/* calculate charge density */
	current_calc() ;

	/***************************************************************
	  Write grid data:
	***************************************************************/

	ZSLOOP(0, N1 - 1, 0, N2 - 1) {
		coord(i, j, CENT, X);

		fprintf(fp, FMT_DBL_OUT, X[1]);
		fprintf(fp, FMT_DBL_OUT, X[2]);

		/* emit full set of primitives */
		PLOOP fprintf(fp, FMT_DBL_OUT, p[i][j][k]);

		/* divb flux-ct defn; corner-centered.  Use
		   only interior corners */
		if (i > 0+NG && j > 0+NG && i < N1+NG && j < N2+NG) {
			divb = fabs(divb_calc(i,j));
		} else
			divb = 0.;
		
		fprintf(fp, FMT_DBL_OUT, divb);

		geom = get_geometry(i, j, CENT) ;
		get_state(p[i][j], geom, &q);

		/* emit 4-velocity and b-field in covariant & contravariant form */
		for (k = 0; k < NDIM; k++)
			fprintf(fp, FMT_DBL_OUT, q.ucon[k]);
		for (k = 0; k < NDIM; k++)
			fprintf(fp, FMT_DBL_OUT, q.ucov[k]);
		for (k = 0; k < NDIM; k++)
			fprintf(fp, FMT_DBL_OUT, q.bcon[k]);
		for (k = 0; k < NDIM; k++)
			fprintf(fp, FMT_DBL_OUT, q.bcov[k]);

		/* begin characteristic speed calculation */
		mhd_vchar(p[i][j], &q, geom, 1, &vmax, &vmin);
		fprintf(fp, FMT_DBL_OUT, vmin);
		fprintf(fp, FMT_DBL_OUT, vmax);

		mhd_vchar(p[i][j], &q, geom, 2, &vmax, &vmin);
		fprintf(fp, FMT_DBL_OUT, vmin);
		fprintf(fp, FMT_DBL_OUT, vmax);
		/* end characteristic speed calculation */

		/* useful to have \sqrt{-g} in output */
		fprintf(fp, FMT_DBL_OUT, geom->g);

		/* begin current density calculation */
		for (k = 0; k < NDIM; k++)
			fprintf(fp, FMT_DBL_OUT, Jcon[i][j][k]) ;

		double jcov[NDIM] ;
		lower(Jcon[i][j], geom, jcov) ;
		for (k = 0; k < NDIM; k++)
			fprintf(fp, FMT_DBL_OUT, jcov[k]) ;
		/* end current density calculation */

		fprintf(fp, "\n");
	}
}

#undef FMT_DBL_OUT 
#undef FMT_INT_OUT 
