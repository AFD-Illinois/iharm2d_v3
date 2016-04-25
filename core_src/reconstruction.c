
/*

includes all interpolation routines 

- linear (MC or other limiter)
- parabolic (from collela and woodward)
- weno5
- mp5

*/

#include "decs.h"

/* mnemonics for slope limiter */
#define MC	(0)
#define VANL	(1)
#define MINM	(2)

/* set the limiter */
#define LIMITER  MC

/* performs the slope-limiting for the numerical flux calculation */

double slope_lim(double y1, double y2, double y3)
{
	double Dqm, Dqp, Dqc, s;

	/* woodward, or monotonized central, slope limiter */
	if (LIMITER == MC) {
		Dqm = 2. * (y2 - y1);
		Dqp = 2. * (y3 - y2);
		Dqc = 0.5 * (y3 - y1);
		s = Dqm * Dqp;
		if (s <= 0.)
			return 0.;
		else {
			if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
				return (Dqm);
			else if (fabs(Dqp) < fabs(Dqc))
				return (Dqp);
			else
				return (Dqc);
		}
	}
	/* van leer slope limiter */
	else if (LIMITER == VANL) {
		Dqm = (y2 - y1);
		Dqp = (y3 - y2);
		s = Dqm * Dqp;
		if (s <= 0.)
			return 0.;
		else
			return (2. * s / (Dqm + Dqp));
	}

	/* minmod slope limiter (crude but robust) */
	else if (LIMITER == MINM) {
		Dqm = (y2 - y1);
		Dqp = (y3 - y2);
		s = Dqm * Dqp;
		if (s <= 0.)
			return 0.;
		else if (fabs(Dqm) < fabs(Dqp))
			return Dqm;
		else
			return Dqp;
	}

	fprintf(stderr, "unknown slope limiter\n");
	exit(10);

	return (0.);


}

void linear_mc(double x1, double x2, double x3, double *lout, double *rout) 
{
	double Dqm,Dqp,Dqc,s;

	Dqm = 2. * (x2 - x1);
	Dqp = 2. * (x3 - x2);
	Dqc = 0.5 * (x3 - x1);

	s = Dqm * Dqp;

	if (s <= 0.)
		s = 0.;
	else {
		if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
			s = Dqm;
		else if (fabs(Dqp) < fabs(Dqc))
			s = Dqp;
		else
			s = Dqc;
	}

	/* reconstruct left, right */
	*lout = x2 - 0.5*s;
	*rout = x2 + 0.5*s;
}

/*
 * parabolic interpolation subroutin  
 * ref. Colella && Woodward's PPM paper
 *
 * using zone-centered value of 5 continuous zones 
 * to get left and right value of the middle zone.
 *  
 * 
 */
/* author: Xiaoyue Guan */

void para(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{
         int i ;
         double y[5], dq[5];
         double Dqm, Dqc, Dqp, aDqm,aDqp,aDqc,s,l,r,qa, qd, qe;

         y[0]=x1;
         y[1]=x2;
         y[2]=x3;
         y[3]=x4;
         y[4]=x5;

         /* CW 1.7 */
         for(i=1 ; i<=3 ; i++) {
               Dqm = 2. *(y[i]-y[i-1]);
               Dqp = 2. *(y[i+1]-y[i]);
               Dqc = 0.5 *(y[i+1]-y[i-1]);
               aDqm = fabs(Dqm) ;
               aDqp = fabs(Dqp) ;
               aDqc = fabs(Dqc) ;
               s = Dqm*Dqp;

               if (s <=0.) dq[i]=0.;       //CW1.8
               else dq[i]=copysign(fmin(aDqc, fmin(aDqm, aDqp)), Dqc);
	       //else dq[i]=MY_MIN(aDqc,MY_MIN(aDqm,aDqp))*MY_SIGN(Dqc);
         }

         /* CW 1.6 */
         l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/6.0;
         r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/6.0;

         qa=(r-y[2])*(y[2]-l);
         qd=(r-l);
         qe=6.0*(y[2]-0.5*(l+r));

         if (qa <=0. ) {
                l=y[2];
                r=y[2];
         }

         if (qd*(qd-qe)<0.0) l=3.0*y[2]-2.0*r;
         else if (qd*(qd+qe)<0.0) r=3.0*y[2]-2.0*l;

         *lout=l;   //a_L,j
         *rout=r;
}

/*
 * left and right state reconsitruction using WENO-5,
 * all numbers from WHAM paper
 *
 * using zone-centered value of n=5 continuous zones 
 * to get left and right value of the middle zone.
 *  
 */

/* author: Monika Moscibrodzka */
/* modified and sign error corrected by CFG 12.28.13 */

void weno(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{

    /* based on shu scholarpedia article, eqs 1,2,3 */
    /* 3rd order interpolations */
    double vr[3],vl[3];
    vr[0] =  (3./8.)*x1-(5./4.)*x2+(15./8.)*x3;
    vr[1] = (-1./8.)*x2+(3./4.)*x3+(3./8.)*x4;
    vr[2] =  (3./8.)*x3+(3./4.)*x4-(1./8.)*x5;

    vl[0] =  (3./8.)*x5-(5./4.)*x4+(15./8.)*x3;
    vl[1] = (-1./8.)*x4+(3./4.)*x3+(3./8.)*x2;
    vl[2] =  (3./8.)*x3+(3./4.)*x2-(1./8.)*x1;

    /* smoothness indicators, from Tchekh et al. A18, equiv to Shu eq. 8 */
    double beta[3];
    beta[0]=(13./12.)*pow(x1-2.*x2+x3,2)+(1./4.)*pow(x1-4.*x2+3.*x3,2);
    beta[1]=(13./12.)*pow(x2-2.*x3+x4,2)+(1./4.)*pow(x4-x2,2);
    beta[2]=(13./12.)*pow(x3-2.*x4+x5,2)+(1./4.)*pow(x5-4.*x4+3.*x3,2);
    
    /* nonlinear weights, after shu eq. 9 */
    double den,wtr[3],Wr,wr[3],wtl[3],Wl,wl[3],eps;
    eps=1.e-26;

    den = eps+beta[0]; den *= den; wtr[0] = (1./16.)/den;
    den = eps+beta[1]; den *= den; wtr[1] = (5./8.)/den;
    den = eps+beta[2]; den *= den; wtr[2] = (5./16.)/den;
    Wr = wtr[0]+wtr[1]+wtr[2];
    wr[0] = wtr[0]/Wr ;
    wr[1] = wtr[1]/Wr ;
    wr[2] = wtr[2]/Wr ;

    den = eps+beta[2]; den *= den; wtl[0] = (1./16.)/den;
    den = eps+beta[1]; den *= den; wtl[1] = (5./8.)/den;
    den = eps+beta[0]; den *= den; wtl[2] = (5./16.)/den;
    Wl = wtl[0]+wtl[1]+wtl[2];
    wl[0] = wtl[0]/Wl ;
    wl[1] = wtl[1]/Wl ;
    wl[2] = wtl[2]/Wl ;

    *lout = vl[0]*wl[0]+vl[1]*wl[1]+vl[2]*wl[2];
    *rout = vr[0]*wr[0]+vr[1]*wr[1]+vr[2]*wr[2];
}

/* left and right state reconstruction using MP5
   from pluto

   via Mani Chandra,
   CFG 11-25-14
*/

/* The following code taken from PLUTO. Needs to be cleaned up. */
#define MINMOD(a,b)  ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0)
double Median (double a, double b, double c)
{
  return (a + MINMOD(b-a,c-a));
}

void mp5(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{
        double MP5_Reconstruct(double Fjm2, double Fjm1, double Fj, 
                double Fjp1, double Fjp2) ;

        *rout = MP5_Reconstruct(x1,x2,x3,x4,x5) ;
        *lout = MP5_Reconstruct(x5,x4,x3,x2,x1) ;

        return;
}

double MP5_Reconstruct(double Fjm2, double Fjm1, double Fj, double Fjp1, double Fjp2)
{
  double f, d2, d2p, d2m; 
  double dMMm, dMMp;
  double scrh1,scrh2, Fmin, Fmax; 
  double fAV, fMD, fLC, fUL, fMP;
  static double alpha = 4.0, epsm = 1.e-12;

  f  = 2.0*Fjm2 - 13.0*Fjm1 + 47.0*Fj + 27.0*Fjp1 - 3.0*Fjp2;
  f /= 60.0;   

  fMP = Fj + MINMOD(Fjp1-Fj,alpha*(Fj-Fjm1));

  if ((f - Fj)*(f - fMP) <= epsm) return f;

  d2m = Fjm2 + Fj - 2.0*Fjm1;    /* -- Eq. (2.19) -- */
  d2  = Fjm1 + Fjp1 - 2.0*Fj;
  d2p = Fj + Fjp2 - 2.0*Fjp1;    /* -- Eq. (2.19) -- */

  scrh1 = MINMOD(4.0*d2 - d2p, 4.0*d2p - d2);
  scrh2 = MINMOD(d2, d2p);
  dMMp  = MINMOD(scrh1,scrh2);   /* -- Eq. (2.27) -- */

  scrh1 = MINMOD(4.0*d2m - d2, 4.0*d2 - d2m);
  scrh2 = MINMOD(d2, d2m);
  dMMm  = MINMOD(scrh1,scrh2);   /* -- Eq. (2.27) -- */

  fUL = Fj + alpha*(Fj - Fjm1);   /* -- Eq. (2.8) -- */
  fAV = 0.5*(Fj + Fjp1);        /* -- Eq. (2.16) -- */
  fMD = fAV - 0.5*dMMp; /* -- Eq. (2.28) -- */
  fLC = 0.5*(3.0*Fj - Fjm1) + 4.0/3.0*dMMm;  /* -- Eq. (2.29) -- */

  scrh1 = fmin(Fj, Fjp1); scrh1 = fmin(scrh1, fMD);
  scrh2 = fmin(Fj, fUL);    scrh2 = fmin(scrh2, fLC);
  Fmin  = fmax(scrh1, scrh2);  /* -- Eq. (2.24a) -- */

  scrh1 = fmax(Fj, Fjp1); scrh1 = fmax(scrh1, fMD);
  scrh2 = fmax(Fj, fUL);    scrh2 = fmax(scrh2, fLC);
  Fmax  = fmin(scrh1, scrh2);  /* -- Eq. (2.24b) -- */

  f = Median(f, Fmin, Fmax); /* -- Eq. (2.26) -- */
  return f;

}

#undef MINMOD

/***

Reconstruction routines for linear, parabolic, weno, mp5

***/

void reconstruct_lr_lin(double ptmp[NMAX+2*NG][NPR], int N, 
	double p_l[NMAX+2*NG][NPR], double p_r[NMAX+2*NG][NPR]) 
{
        double dqtmp[NMAX+2*NG][NPR] ;
        int i,k ;

	ISLOOP(-1,N) PLOOP dqtmp[i][k] = 0. ;

        /* calculate slopes */
	ISLOOP(-1,N) PLOOP dqtmp[i][k] = slope_lim(ptmp[i-1][k],ptmp[i][k],ptmp[i+1][k]) ;

        /* reconstruct left */
	ISLOOP(0,N) PLOOP p_l[i][k] = ptmp[i][k] - 0.5*dqtmp[i][k];

        /* reconstruct right */
	ISLOOP(-1,N-1) PLOOP p_r[i][k] = ptmp[i][k] + 0.5*dqtmp[i][k];
}


void reconstruct_lr_par(double ptmp[NMAX+2*NG][NPR], int N, 
	double p_l[NMAX+2*NG][NPR], double p_r[NMAX+2*NG][NPR]) 
{
	int i,k ;

	ISLOOP(-1,N) {
		PLOOP {
			para(   ptmp[i-2][k], 
				ptmp[i-1][k], 
				ptmp[i][k], 
				ptmp[i+1][k], 
				ptmp[i+2][k],
				&p_l[i][k], 
				&p_r[i][k]) ;
		}
	}
}

void reconstruct_lr_weno(double ptmp[NMAX+2*NG][NPR], int N, 
	double p_l[NMAX+2*NG][NPR], double p_r[NMAX+2*NG][NPR]) 
{
	int i,k ;

	ISLOOP(-1,N) {
		PLOOP {
			weno(   ptmp[i-2][k], 
				ptmp[i-1][k], 
				ptmp[i][k], 
				ptmp[i+1][k], 
				ptmp[i+2][k],
				&p_l[i][k], 
				&p_r[i][k]) ;
		}
	}
}

void reconstruct_lr_mp5(double ptmp[NMAX+2*NG][NPR], int N, 
        double p_l[NMAX+2*NG][NPR], double p_r[NMAX+2*NG][NPR]) 
{
        int i,k ;

        ISLOOP(-1,N) {
                PLOOP {
                        mp5(   ptmp[i-2][k], 
                                ptmp[i-1][k], 
                                ptmp[i][k], 
                                ptmp[i+1][k], 
                                ptmp[i+2][k],
                                &p_l[i][k], 
                                &p_r[i][k]) ;
                }
        }
}

