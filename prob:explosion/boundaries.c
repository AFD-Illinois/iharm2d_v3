
#include "decs.h"

void bound_prim(double prim[][N2 + 2*NG][NPR])
{
	int i, j, k;

	/* periodic boundary conditions */
	JSLOOP(0,N2-1) {
		ISLOOP(-NG,-1) {
			PLOOP prim[i][j][k] = prim[i+N1][j][k];
			pflag[i][j] = pflag[i+N1][j];
		}
		ISLOOP(N1,N1-1+NG) {
			PLOOP prim[i][j][k] = prim[i-N1][j][k];
			pflag[i][j] = pflag[i-N1][j];
		}
	}

	ISLOOP(-NG,N1-1+NG) {
		JSLOOP(-NG,-1) {
			PLOOP prim[i][j][k] = prim[i][j+N2][k];
			pflag[i][j] = pflag[i][j+N2];
		}
		JSLOOP(N2,N2-1+NG) {
			PLOOP prim[i][j][k] = prim[i][j-N2][k];
			pflag[i][j] = pflag[i][j-N2];
		}
	}
}
