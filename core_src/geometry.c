/*

	contains coordinate-independent, geometrical evaluations

	problem-independent

	removed gsl dependency
	7.19.21

*/

#include "decs.h"



double MINOR(double m[16], int r0, int r1, int r2, int c0, int c1, int c2);
void adjoint(double m[16], double adjOut[16]);
double determinant(double m[16]);

/*
  set_geometry():
  ----------

       -- calculates all grid functions that remain constant 
          over time, such as the metric (gcov), inverse metric 
          (gcon), connection coefficients (conn), and sqrt of 
          the metric's determinant (gdet).

*/
void set_geometry()
{
	int i, j ;
	double X[NDIM] = {0.} ;
	int k, l ;

	/* for zone averaging connection */
	double Xs[NDIM] = {0.};
	struct of_geom ggeoms ;
	double conns[NDIM][NDIM][NDIM] ;
	int m, n, o ;
	int NZA = 5 ;	/* number of points to average connection over */

	/* set up coordinate grid */
	set_coordinates();

	/* loop over grid and set geometry on grid */
	ZSLOOP(-NG, N1-1 + NG, -NG, N2-1 + NG) {

		/* zone-centered */
		coord(i, j, CENT, X);
		gcov_func(X, ggeom[i][j][CENT].gcov);
		ggeom[i][j][CENT].g = gdet_func(ggeom[i][j][CENT].gcov);
		gcon_func(ggeom[i][j][CENT].gcov, ggeom[i][j][CENT].gcon);
		ggeom[i][j][CENT].alpha = 1.0 / sqrt(-ggeom[i][j][CENT].gcon[0][0]);

		/* connection at zone center */
		//conn_func(X, &ggeom[i][j][CENT], conn[i][j]);

		/* begin calculating connection averaged over zone volume */
		/* calculate <\Gamma \sqrt{-g}>/\sqrt{-g} */
		for(m=0;m<NDIM;m++)
		for(n=0;n<NDIM;n++)
		for(o=0;o<NDIM;o++) conn[i][j][m][n][o] = 0. ;

		coord(i,j,CORN,X) ;
		for(k=0;k<NZA;k++) 
		for(l=0;l<NZA;l++) {
			Xs[1] = X[1] + (k+0.5)*dx[1]/NZA ;
			Xs[2] = X[2] + (l+0.5)*dx[2]/NZA ;
			gcov_func(Xs, ggeoms.gcov);
			ggeoms.g = gdet_func(ggeoms.gcov) ;
			gcon_func(ggeoms.gcov, ggeoms.gcon) ;
			conn_func(Xs, &ggeoms, conns) ;
			
			for(m=0;m<NDIM;m++)
			for(n=0;n<NDIM;n++)
			for(o=0;o<NDIM;o++)
				conn[i][j][m][n][o] += ggeoms.g*conns[m][n][o] ;
		}
		for(m=0;m<NDIM;m++)
		for(n=0;n<NDIM;n++)
		for(o=0;o<NDIM;o++)
			conn[i][j][m][n][o] /= (ggeom[i][j][CENT].g*NZA*NZA) ;
		/* end calculating connection averaged over zone volume */

		/* corner-centered */
		coord(i, j, CORN, X);
		gcov_func(X, ggeom[i][j][CORN].gcov);
		ggeom[i][j][CORN].g = gdet_func(ggeom[i][j][CORN].gcov);
		gcon_func(ggeom[i][j][CORN].gcov, ggeom[i][j][CORN].gcon);
		ggeom[i][j][CORN].alpha = 1.0 / sqrt(-ggeom[i][j][CORN].gcon[0][0]);

		/* r-face-centered */
		coord(i, j, FACE1, X);
		gcov_func(X, ggeom[i][j][FACE1].gcov);
		ggeom[i][j][FACE1].g = gdet_func(ggeom[i][j][FACE1].gcov);
		gcon_func(ggeom[i][j][FACE1].gcov, ggeom[i][j][FACE1].gcon);
		ggeom[i][j][FACE1].alpha = 1.0 / sqrt(-ggeom[i][j][FACE1].gcon[0][0]);

		/* theta-face-centered */
		coord(i, j, FACE2, X);
		gcov_func(X, ggeom[i][j][FACE2].gcov);
		ggeom[i][j][FACE2].g = gdet_func(ggeom[i][j][FACE2].gcov);
		gcon_func(ggeom[i][j][FACE2].gcov, ggeom[i][j][FACE2].gcon);
		ggeom[i][j][FACE2].alpha = 1.0 / sqrt(-ggeom[i][j][FACE2].gcon[0][0]);
	}

	/* done! */
}


/*
    coord():
    -------
       -- given the indices i,j and location in the cell, return with 
          the values of X1,X2 there;  
       -- the locations are defined by : 
           -----------------------
           |                     |
           |                     |
           |FACE1   CENT         |
           |                     |
           |CORN    FACE2        |
           ----------------------

*/
void coord(int i, int j, int loc, double *X)
{
	if (loc == FACE1) {
		X[1] = startx[1] + (i - START1) * dx[1];
		X[2] = startx[2] + (j + 0.5 - START2) * dx[2];
	} else if (loc == FACE2) {
		X[1] = startx[1] + (i + 0.5 - START1) * dx[1];
		X[2] = startx[2] + (j - START2) * dx[2];
	} else if (loc == CENT) {
		X[1] = startx[1] + (i + 0.5 - START1) * dx[1];
		X[2] = startx[2] + (j + 0.5 - START2) * dx[2];
	} else {
		X[1] = startx[1] + (i - START1) * dx[1];
		X[2] = startx[2] + (j - START2) * dx[2];
	}

	return;
}

/* assumes gcov has been set first; returns determinant */
/* 
	removed gsl dependency
	7.19.21
*/
double gdet_func(double gcov[][NDIM])
{
	double det = determinant(&gcov[0][0]);
	if(det == 0) {
		fprintf(stderr,
			"gdet_func(): singular matrix encountered! \n");
		fail(FAIL_METRIC);
	}
	return sqrt(fabs(det));
}


/* invert gcov to get gcon */
/*

	removed gsl dependency
	7.19.21

*/

void gcon_func(double gcov[][NDIM], double gcon[][NDIM])
{
  double *m = &gcov[0][0];
  double *invOut = &gcon[0][0];
  adjoint(m, invOut);

  double det = determinant(m);
  if(det == 0.) {
	fprintf(stderr,
		"gcon_func(): singular matrix encountered! \n");
	fail(FAIL_METRIC);
  }
  double inv_det = 1. / det;
  for (int i = 0; i < 16; ++i) {
    invOut[i] = invOut[i]*inv_det;
  }
}



/*
  conn_func():
  -----------

   -- this gives the connection coefficient
	\Gamma^{i}_{j,k} = conn[..][i][j][k]
   --  where i = {1,2,3,4} corresponds to {t,r,theta,phi}


  notice: this version uses *numerical differentiation* to
  obtain the connection from the line element

*/

/* Sets the spatial discretization in numerical derivatives : */
#define DELTA 1.e-5

/* NOTE: parameter hides global variable */
void conn_func(double *X, struct of_geom *geom, double conn[][NDIM][NDIM])
{
	int i, j, k, l;
	double tmp[NDIM][NDIM][NDIM];
	double Xh[NDIM], Xl[NDIM];
	double gh[NDIM][NDIM];
	double gl[NDIM][NDIM];

	for (k = 0; k < NDIM; k++) {
		for (l = 0; l < NDIM; l++)
			Xh[l] = X[l];
		for (l = 0; l < NDIM; l++)
			Xl[l] = X[l];
		Xh[k] += DELTA;
		Xl[k] -= DELTA;
		gcov_func(Xh, gh);
		gcov_func(Xl, gl);

		for (i = 0; i < NDIM; i++)
			for (j = 0; j < NDIM; j++)
				conn[i][j][k] =
				    (gh[i][j] - gl[i][j]) / (Xh[k] -
							     Xl[k]);
	}

	/* now rearrange to find \Gamma_{ijk} */
	for (i = 0; i < NDIM; i++)
		for (j = 0; j < NDIM; j++)
			for (k = 0; k < NDIM; k++)
				tmp[i][j][k] =
				    0.5 * (conn[j][i][k] + conn[k][i][j] -
					   conn[k][j][i]);

	/* finally, raise index */
	for (i = 0; i < NDIM; i++)
		for (j = 0; j < NDIM; j++)
			for (k = 0; k < NDIM; k++) {
				conn[i][j][k] = 0.;
				for (l = 0; l < NDIM; l++)
					conn[i][j][k] +=
					    geom->gcon[i][l] *
					    tmp[l][j][k];
			}

	/* done! */
}

#undef DELTA

/* Lowers a contravariant rank-1 tensor to a covariant one */
void lower(double *ucon, struct of_geom *geom, double *ucov)
{

	ucov[0] = geom->gcov[0][0] * ucon[0]
	    + geom->gcov[0][1] * ucon[1]
	    + geom->gcov[0][2] * ucon[2]
	    + geom->gcov[0][3] * ucon[3];
	ucov[1] = geom->gcov[1][0] * ucon[0]
	    + geom->gcov[1][1] * ucon[1]
	    + geom->gcov[1][2] * ucon[2]
	    + geom->gcov[1][3] * ucon[3];
	ucov[2] = geom->gcov[2][0] * ucon[0]
	    + geom->gcov[2][1] * ucon[1]
	    + geom->gcov[2][2] * ucon[2]
	    + geom->gcov[2][3] * ucon[3];
	ucov[3] = geom->gcov[3][0] * ucon[0]
	    + geom->gcov[3][1] * ucon[1]
	    + geom->gcov[3][2] * ucon[2]
	    + geom->gcov[3][3] * ucon[3];

	return;
}

/* Raises a covariant rank-1 tensor to a contravariant one */
void raise(double *ucov, struct of_geom *geom, double *ucon)
{

	ucon[0] = geom->gcon[0][0] * ucov[0]
	    + geom->gcon[0][1] * ucov[1]
	    + geom->gcon[0][2] * ucov[2]
	    + geom->gcon[0][3] * ucov[3];
	ucon[1] = geom->gcon[1][0] * ucov[0]
	    + geom->gcon[1][1] * ucov[1]
	    + geom->gcon[1][2] * ucov[2]
	    + geom->gcon[1][3] * ucov[3];
	ucon[2] = geom->gcon[2][0] * ucov[0]
	    + geom->gcon[2][1] * ucov[1]
	    + geom->gcon[2][2] * ucov[2]
	    + geom->gcon[2][3] * ucov[3];
	ucon[3] = geom->gcon[3][0] * ucov[0]
	    + geom->gcon[3][1] * ucov[1]
	    + geom->gcon[3][2] * ucov[2]
	    + geom->gcon[3][3] * ucov[3];

	return;
}

/* return pointer to local geometry */
struct of_geom *get_geometry(int ii, int jj, int kk)
{
	icurr = ii;
	jcurr = jj;

	return( &(ggeom[ii][jj][kk]) ) ;
}




/* cofactor helper function */
inline double MINOR(double m[16], int r0, int r1, int r2, int c0, int c1, int c2)
{
  return m[4*r0+c0]*(m[4*r1+c1]*m[4*r2+c2] - m[4*r2+c1]*m[4*r1+c2]) -
         m[4*r0+c1]*(m[4*r1+c0]*m[4*r2+c2] - m[4*r2+c0]*m[4*r1+c2]) +
         m[4*r0+c2]*(m[4*r1+c0]*m[4*r2+c1] - m[4*r2+c0]*m[4*r1+c1]);
}
/* returns determinant of 4 by 4 matrix */
inline double determinant(double m[16])
{
  return m[0]*MINOR(m,1,2,3,1,2,3) -
         m[1]*MINOR(m,1,2,3,0,2,3) +
         m[2]*MINOR(m,1,2,3,0,1,3) -
         m[3]*MINOR(m,1,2,3,0,1,2);
}
/* sets adjOut as adjoint of m */
inline void adjoint(double m[16], double adjOut[16])
{
  adjOut[ 0] =  MINOR(m,1,2,3,1,2,3);
  adjOut[ 1] = -MINOR(m,0,2,3,1,2,3);
  adjOut[ 2] =  MINOR(m,0,1,3,1,2,3);
  adjOut[ 3] = -MINOR(m,0,1,2,1,2,3);

  adjOut[ 4] = -MINOR(m,1,2,3,0,2,3);
  adjOut[ 5] =  MINOR(m,0,2,3,0,2,3);
  adjOut[ 6] = -MINOR(m,0,1,3,0,2,3);
  adjOut[ 7] =  MINOR(m,0,1,2,0,2,3);

  adjOut[ 8] =  MINOR(m,1,2,3,0,1,3);
  adjOut[ 9] = -MINOR(m,0,2,3,0,1,3);
  adjOut[10] =  MINOR(m,0,1,3,0,1,3);
  adjOut[11] = -MINOR(m,0,1,2,0,1,3);

  adjOut[12] = -MINOR(m,1,2,3,0,1,2);
  adjOut[13] =  MINOR(m,0,2,3,0,1,2);
  adjOut[14] = -MINOR(m,0,1,3,0,1,2);
  adjOut[15] =  MINOR(m,0,1,2,0,1,2);
}

