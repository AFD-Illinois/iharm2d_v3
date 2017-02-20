
#include "decs.h"

void bound_prim(double prim[][N2 + 2*NG][NPR])
{
	int i, j, k;

	static double primL[NPR];
	static double primR[NPR];
	static int firstc = 1;

#if 0
	/* just copy out the values in the last zones */
	if(nstep<2) {
		PLOOP primL[k] = prim[0][0][k];
		PLOOP primR[k] = prim[N1-1][0][k];
		firstc = 0;
	}

	/* constant boundary conditions */
	JSLOOP(0,N2-1) {
		ISLOOP(-NG,-1) {
			PLOOP prim[i][j][k] = primL[k];
			pflag[i][j] = pflag[0][j];
		}
		ISLOOP(N1,N1-1+NG) {
			PLOOP prim[i][j][k] = primR[k];
			pflag[i][j] = pflag[N1-1][j];
		}
	}
#endif
#if 1
	JSLOOP(0,N2-1) {
		ISLOOP(-NG,-1) {
			PLOOP prim[i][j][k] = prim[i+N1][j][k];
			pflag[i][j] = pflag[0][j];
		}
		ISLOOP(N1,N1-1+NG) {
			PLOOP prim[i][j][k] = prim[i-N1][j][k];
			pflag[i][j] = pflag[N1-1][j];
		}
	}
#endif

	/* periodic boundary conditions */
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
