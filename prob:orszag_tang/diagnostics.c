
/* 
   diagnostics routines.

   problem-independent *except* for choice
   	of information output to log file
	and dump files 
   
*/

#include "decs.h"

/* output full dump of simulation state */
void diag_dump()
{
	char dump_file_name[100];
	FILE *dump_file;

	/* emit dump file */
	sprintf(dump_file_name, "dumps/dump%03d", dump_cnt);
	fprintf(stderr, "DUMP     file=%s\n", dump_file_name);
	dump_file = fopen(dump_file_name, "w");

	if (dump_file == NULL) {
		fprintf(stderr, "error opening dump file\n");
		exit(2);
	}

	dump(dump_file);
	fclose(dump_file);

	/* increment global dump counter */
	dump_cnt++;
}

/* output to logfile */
void diag_log()
{
	double U[NPR],sum_U[NPR],divb;
	double divbmax = 0. ;
	int imax = 0;
	int jmax = 0;
	int i,j,k;
	struct of_geom *geom;
	struct of_state q;
	static FILE *log_file;
	static int firstc=1;

	PLOOP sum_U[k] = 0. ;
	
	/* open diagnostics output file */
	if(firstc) {
		firstc = 0;

		log_file = fopen("ener.out","a") ;
		if (log_file == NULL) {
			fprintf(stderr,"error opening energy output file\n"); 
			exit(1);
		}
	}

	/* evaluate diagnostics for log file */
	ZSLOOP(0, N1 - 1, 0, N2 - 1) {
		geom = get_geometry(i, j, CENT) ;
		get_state(p[i][j], geom, &q);
		primtoU(p[i][j], &q, geom, U);

		/* diagnostic 1: sum mass, angular momentum, etc., over grid */
		PLOOP sum_U[k] += U[k] * dx[1]*dx[2]*dx[3];

		/* diagnostic 2: evaluate maximum divb */
		/* Use only interior corners */
		if(i > 0+NG && j > 0+NG && i < N1+NG && j < N2+NG) {
			divb = fabs(divb_calc(i,j));

			if (divb > divbmax) {
				imax = i;
				jmax = j;
				divbmax = divb;
			}
		}
	}

	/* alert user in stderr */
	fprintf(stderr, "LOG      t=%g \t divbmax: %d %d %g\n", t,imax,jmax,divbmax);

	/* begin log file output */
	fprintf(log_file, "%15.8g ",t) ;

	/* output integrated quantities to log file */
	fprintf(log_file, "%15.8g %15.8g %15.8g %15.8g %15.8g ",
			sum_U[RHO],
			sum_U[UU],
			sum_U[U1],
			sum_U[U2],
			sum_U[U3]
			) ;
	/* output boundary fluxes of interest to file */

	fprintf(log_file, "\n");
	fflush(log_file);
	/* end log file output */

	/* done! */
}


/* output lagrangian tracer particle dump */
void diag_pdump()
{
	char pdump_file_name[100];
	FILE *pdump_file;

	sprintf(pdump_file_name, "dumps/pdump%03d", pdump_cnt);
	fprintf(stderr, "PDUMP     file=%s\n", pdump_file_name);
	pdump_file = fopen(pdump_file_name, "w");

	if (pdump_file == NULL) {
		fprintf(stderr, "error opening pdump file\n");
		exit(2);
	}

	pdump(pdump_file);
	fclose(pdump_file);

	pdump_cnt++;
}

/* output desired images */
void diag_image()
{
	int i,j,k;
	double dimage[N1*N2];

	/** for each type of image file, 
		- assign floating point values to image array dimage
		- assign descriptor to image file, e.g. rho
		- emit image **/

	/* density */
	k = 0;
	for(j=0+START2;j<N2+START2;j++)for(i=0+START1;i<N1+START1;i++) {
		dimage[k] = p[i][j][RHO];
		k++;
	}
	image_ppm(dimage, "rho");

	/* log density */
	k = 0;
	for(j=0+START2;j<N2+START2;j++)for(i=0+START1;i<N1+START1;i++) {
		dimage[k] = log(p[i][j][RHO]);
		k++;
	}
	image_ppm(dimage, "logrho");

	/* b . b */
	k = 0;
	struct of_geom *geom;
	double bsq;
	for(j=0+START2;j<N2+START2;j++)for(i=0+START1;i<N1+START1;i++) {
		geom = get_geometry(i, j, CENT);
		bsq = bsq_calc(p[i][j],geom);
		dimage[k] = bsq;
		k++;
	}
	image_ppm(dimage, "bsq");
	
	/* log( b . b ) */
	k = 0;
	for(j=0+START2;j<N2+START2;j++)for(i=0+START1;i<N1+START1;i++) {
		geom = get_geometry(i, j, CENT);
		bsq = bsq_calc(p[i][j],geom);
		dimage[k] = log( fabs( bsq ) + 1.e-15 );
		k++;
	}
	image_ppm(dimage, "logbsq");

	/* temperature */
	k = 0;
	double T;
	double Pressure_rho0_u(double rho, double u) ;
	for(j=0+START2;j<N2+START2;j++)for(i=0+START1;i<N1+START1;i++) {
		T = Pressure_rho0_u(p[i][j][RHO], p[i][j][UU]) / p[i][j][RHO];
		dimage[k] = T;
		k++;
	}
	image_ppm(dimage, "T");

	
	/* log temperature */
	k = 0;
	for(j=0+START2;j<N2+START2;j++)for(i=0+START1;i<N1+START1;i++) {
		T = Pressure_rho0_u(p[i][j][RHO], p[i][j][UU]) / p[i][j][RHO];
		dimage[k] = log( fabs(T) + 1.e-15 );
		k++;
	}
	image_ppm(dimage, "logT");
	
	/* current in frame of fluid, squared */
	current_calc();
	double Jsq,Jdu,gJsq,jcov[NDIM];
	struct of_state q;
	k = 0;
	for(j=0+START2;j<N2+START2;j++)for(i=0+START1;i<N1+START1;i++) {
		geom = get_geometry(i, j, CENT);
        	get_state(p[i][j], geom, &q);
                lower(Jcon[i][j], geom, jcov) ;
                Jsq = jcov[0]*Jcon[i][j][0] + jcov[1]*Jcon[i][j][1] + 
                        jcov[2]*Jcon[i][j][2] + jcov[3]*Jcon[i][j][3] ;
                Jdu = jcov[0]*q.ucon[0] + jcov[1]*q.ucon[1] + 
                        jcov[2]*q.ucon[2] + jcov[3]*q.ucon[3] ;
                gJsq = geom->g * (Jsq + Jdu*Jdu) ;
		dimage[k] = gJsq;

		k++;
	}
	image_ppm(dimage, "gJsq");

	/* log current */
	k = 0;
	for(j=0+START2;j<N2+START2;j++)for(i=0+START1;i<N1+START1;i++) {
		geom = get_geometry(i, j, CENT);
        	get_state(p[i][j], geom, &q);
                lower(Jcon[i][j], geom, jcov) ;
                Jsq = jcov[0]*Jcon[i][j][0] + jcov[1]*Jcon[i][j][1] + 
                        jcov[2]*Jcon[i][j][2] + jcov[3]*Jcon[i][j][3] ;
                Jdu = jcov[0]*q.ucon[0] + jcov[1]*q.ucon[1] + 
                        jcov[2]*q.ucon[2] + jcov[3]*q.ucon[3] ;
                gJsq = geom->g * (Jsq + Jdu*Jdu) ;
		dimage[k] = log(fabs(gJsq) + 1.e-15);

		k++;
	}
	image_ppm(dimage, "loggJsq");

	/* lorentz factor wrt normal observer */
	k = 0;
	double gamma;
	for(j=0+START2;j<N2+START2;j++)for(i=0+START1;i<N1+START1;i++) {
		geom = get_geometry(i, j, CENT) ;
                
                if (mhd_gamma_calc(p[i][j], geom, &gamma)) {
                        gamma = 1.;
                }
		dimage[k] = gamma;

		k++;
	}
	image_ppm(dimage, "gamma");

	/* log lorentz factor wrt normal observer */
	k = 0;
	for(j=0+START2;j<N2+START2;j++)for(i=0+START1;i<N1+START1;i++) {
		geom = get_geometry(i, j, CENT) ;
                
                if (mhd_gamma_calc(p[i][j], geom, &gamma)) {
                        gamma = 1.;
                }
		dimage[k] = log(fabs(gamma) + 1.e-15);

		k++;
	}
	image_ppm(dimage, "loggamma");

	/* done! increment image counter */
	image_cnt++;
}

/* error handling: run has failed, emit diagnostics and exit */
void fail(int fail_type)
{
	void map_failed_zone(int i, int j);

	fprintf(stderr, "\n\nfail: %d %d %d\n", icurr, jcurr, fail_type);

	map_failed_zone(icurr, jcurr);

	fprintf(stderr, "fail\n");

	diag_log();
	diag_image();
	diag_dump();
	diag_pdump();

	/* for diagnostic purposes */
	exit(0);
}

/* map out region around failure point */
void map_failed_zone(int i, int j)
{
	int k;

	fprintf(stderr, "area map\n");

	PLOOP {
		fprintf(stderr, "variable %d \n", k);
		fprintf(stderr, "i = \t %12d %12d %12d\n", i - 1, i,
			i + 1);
		fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j + 1,
			p[i - 1][j + 1][k], p[i][j + 1][k],
			p[i + 1][j + 1][k]);
		fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j,
			p[i - 1][j][k], p[i][j][k],
			p[i + 1][j][k]);
		fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j - 1,
			p[i - 1][j - 1][k], p[i][j - 1][k],
			p[i + 1][j - 1][k]);
	}

	/* other diagnostics here: */
}

