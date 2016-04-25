
#include "decs.h"

/*

subroutines related to 
special diagnostic routine to obtain current
and charge density in fluid frame

ASSUMES: p, psave, and dtsave set.

RETURNS: 4-current values through global array Jcon 

Does not require BCs to be set correctly; 
current calculated only at interior points.

Written to be readable rather than efficient.

cfg 25 june 10

- fixed bug in calculation of DF, cfg 18 june 12

*/

void current_calc()
{
        static double pa[N1+2*NG][N2+2*NG][NPR] ;
	double Fcon_calc(double *prim, int i1, int i2, int i, int j) ;
	double gF0p[NDIM],gF0m[NDIM],gF1p[NDIM],gF1m[NDIM],gF2p[NDIM],gF2m[NDIM] ;
	struct of_geom *geom ;
	int i,j,k ;
	
	/* calculate a time-centered p */
	ZSLOOP(-1,N1-1+NG,-1,N2-1+NG) PLOOP pa[i][j][k] = 0.5*(p[i][j][k] + psave[i][j][k]) ;

	/* calculate J using centered differences; interior zones only */
	ZLOOP for(k=0;k<NDIM;k++) Jcon[i][j][k] = 0. ;
	ZLOOP {

		/* get gdet * Fmunu at neighboring points */
		/* first time direction */
		for(k = 0 ; k < NDIM ; k++) gF0p[k] = Fcon_calc(p[i][j],  0, k, i, j) ;
		for(k = 0 ; k < NDIM ; k++) gF0m[k] = Fcon_calc(psave[i][j], 0, k, i, j) ;

		/* x1 direction */
		for(k = 0 ; k < NDIM ; k++) gF1p[k] = Fcon_calc(pa[i+1][j], 1, k, i+1, j) ;
		for(k = 0 ; k < NDIM ; k++) gF1m[k] = Fcon_calc(pa[i-1][j], 1, k, i-1, j) ;

		/* x2 direction */
		for(k = 0 ; k < NDIM ; k++) gF2p[k] = Fcon_calc(pa[i][j+1], 2, k, i, j+1) ;
		for(k = 0 ; k < NDIM ; k++) gF2m[k] = Fcon_calc(pa[i][j-1], 2, k, i, j-1) ;

		/* get gdet at point */
		geom = get_geometry(i, j, CENT) ;

		/* difference ; use Maxwell in the form D_b F^{ab} = 4\pi J^a,
		   assuming symmetry along the 3-axis */
		for(k = 0 ; k < NDIM ; k++) {
			Jcon[i][j][k] = (1./(4.*M_PI*geom->g))*(
				(gF0p[k] - gF0m[k])/dtsave +
				(gF1p[k] - gF1m[k])/(2.*dx[1]) +
				(gF2p[k] - gF2m[k])/(2.*dx[2]) 
				) ;
		}
	}

	return ;
}

/* return single component of the contravariant maxwell tensor at position i,j
	component i1,i2, constructed from primitives prim */
double Fcon_calc(double *prim, int i1, int i2, int i, int j) 
{
	struct of_geom *geom ;
	double ucon[NDIM],ucov[NDIM],bcon[NDIM],bcov[NDIM] ;
	double Fcon,gFcon,dFcon ;
	int k,l ;
	int antisym(int a, int b, int c, int d) ;

	if(i1 == i2) return(0.) ;

	geom = get_geometry(i,j, CENT) ;
        ucon_calc(prim, geom, ucon);
        lower(ucon, geom, ucov);
        bcon_calc(prim, ucon, ucov, bcon);
        lower(bcon, geom, bcov);

	Fcon = 0. ;
	for(k=0;k<4;k++)
	for(l=0;l<4;l++) {
		dFcon = (-1./geom->g)*antisym(i1,i2,k,l)*ucov[k]*bcov[l] ;
		Fcon += dFcon ;
	}

	gFcon = Fcon * geom->g ;

	return(gFcon) ;
}


/* completely antisymmetric symbol in 4D.
   verified against mathematica */
int antisym(int a, int b, int c, int d)
{
	int pp(int n, int *P) ;

        /** check for a valid permutation **/
        /* range? */
        if(a < 0 || a > 3) return(100) ;
        if(b < 0 || b > 3) return(100) ;
        if(c < 0 || c > 3) return(100) ;
        if(d < 0 || d > 3) return(100) ;

        /* entries different? */
        if(a == b) return(0) ;
        if(a == c) return(0) ;
        if(a == d) return(0) ;
        if(b == c) return(0) ;
        if(b == d) return(0) ;
        if(c == d) return(0) ;

        /* determine parity of permutation */
        int p[4] = {a,b,c,d} ;
        return(pp(4,p)) ;
}

/* algorithm tracks back to Norm Hardy; good for general n */
int pp(int n, int P[n])
{
        int j,x ;
        int p = 0;
        int v[n] ;

        for(j=0;j<n;j++) v[j] = 0 ;

        for(j=0;j<n;j++) {
                if (v[j]) p++;
                else {
                        x = j;
                        do {
                                x = P[x];
                                v[x] = 1;
                        } while (x != j);
                }
        }

        if(p%2 == 0) return(1) ;
        else return(-1) ;
}

