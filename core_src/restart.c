
/* restart functions: restart_init and restart_dump */

#include "decs.h"

/* formatting strings that give the necessary precision */
#define FMT_DBL_OUT "%28.18e"
#define FMT_INT_OUT "%10d"

/* local function decs */
FILE *get_restart_file(char *restart_filename);

/***********************************************************************/
/***********************************************************************
  restart_write():
     -- writes current state of primitive variables to the 
        checkpointing/restart file. 
     -- uses ASCII text format ;
     -- when changing this routine, be sure to make analogous changes 
        in restart_read();
************************************************************************/
void restart_write()
{
	FILE *fp;
	int i, j, k, l;


	if (rdump_cnt % 2 == 0) {
		
		fp = fopen("dumps/rdumpA", "wt");
		fprintf(stderr, "RESTART  file=dumps/rdumpA\n");
	} else {
		fp = fopen("dumps/rdumpB", "wt");
		fprintf(stderr, "RESTART file=dumps/rdumpA\n");
	}

	if (fp == NULL) {
		fprintf(stderr, "Cannot open restart file\n");
		exit(2);
	}

  /*************************************************************
	  Write the header of the restart file: 
  *************************************************************/
	fprintf(fp, FMT_INT_OUT, N1);
	fprintf(fp, FMT_INT_OUT, N2);
	fprintf(fp, FMT_DBL_OUT, t);
	fprintf(fp, FMT_DBL_OUT, tf);
	fprintf(fp, FMT_INT_OUT, nstep);
	fprintf(fp, FMT_DBL_OUT, gam);
	fprintf(fp, FMT_DBL_OUT, cour);
	fprintf(fp, FMT_DBL_OUT, DTd);
	fprintf(fp, FMT_DBL_OUT, DTl);
	fprintf(fp, FMT_DBL_OUT, DTi);
	fprintf(fp, FMT_DBL_OUT, DTp);
	fprintf(fp, FMT_INT_OUT, DTr);
	fprintf(fp, FMT_INT_OUT, dump_cnt);
	fprintf(fp, FMT_INT_OUT, image_cnt);
	fprintf(fp, FMT_INT_OUT, rdump_cnt);
	fprintf(fp, FMT_DBL_OUT, dt);

	fprintf(fp, "\n");

  /*************************************************************
	  Write the body of the restart file: 
  *************************************************************/
	ZSLOOP(-NG, N1-1 + NG, -NG, N2-1 + NG) {
		PLOOP fprintf(fp, FMT_DBL_OUT, p[i][j][k]);
		fprintf(fp, "\n");
	}

	for(l=0;l<NPTOT;l++) {
		for(k=1;k<NDIM;k++) fprintf(fp, FMT_DBL_OUT, xp[l][k]) ;
	}

	fclose(fp);

	rdump_cnt++;

	return;
}

/***********************************************************************/
/***********************************************************************
  restart_init():
     -- main driver for setting initial conditions from a checkpoint 
        or restart file. 
     -- determines if there are any restart files to use and then 
         lets the  user choose if there are more than one file. 
     -- then calls initializes run with restart data;
************************************************************************/
int restart_init()
{
	FILE *fp;
	FILE *fpA = NULL;
	FILE *fpB = NULL;
	char ans[100];
	int i, j, k;

	/* initialize dependent variable arrays */
	zero_arrays();

	/* set up grid functions */
	set_geometry();


  /********************************************************************
   Check to see which restart files exist. 
   Use the only one that exists, else prompt user to decide 
     which one to use if we have a choice : 
  ********************************************************************/

        fpA = get_restart_file("dumps/rdumpA");
        fpB = get_restart_file("dumps/rdumpB");

	if ((fpA == NULL) && (fpB == NULL)) {
		fprintf(stderr, "No restart file\n");
		return (0);
	} else {
		fprintf(stderr, "\nRestart file exists! \n");
		if (fpA == NULL) {
			fprintf(stderr, "Using dumps/rdumpB ... \n");
			fp = fpB;
		} else if (fpB == NULL) {
			fprintf(stderr, "Using dumps/rdumpA ... \n");
			fp = fpA;
		} else {
			/* need to figure out which dump is which */

			fprintf(stderr,
				"Use dumps/rdumpA (0) or dumps/rdumpB (1)?   [A|B]  \n");
			fscanf(stdin, "%s", ans);
			if (strncmp(ans, "A", 1) == 0) {
				fp = fpA;
			} else {
				fp = fpB;
			}
		}
	}

  /********************************************************************
   Now that we know we are restarting from a checkpoint file, then 
   we need to read in data, assign grid functions and define the grid: 
  ********************************************************************/
	restart_read(fp);
	fclose(fp);

	/* set half-step primitives everywhere */
	ZSLOOP(-NG, N1-1 + NG, -NG, N2-1 + NG) PLOOP ph[i][j][k] = p[i][j][k];

	/* set metric functions */
	set_coordinates();

	/* bound */
	bound_prim(p);
	bound_prim(ph);


  /***********************************************************************
    Make any changes to parameters in restart file  here: 
      e.g., cour = 0.4 , tf ...
  ************************************************************************/
	//cour = 0.9 ;
	//tf = 4000. ;

	fprintf(stderr, "done with restart init.\n");

	/* done! */
	return (1);

}

/***********************************************************************/
/***********************************************************************
  restart_read():
     -- reads in data from the restart file, which is specified in 
         restart_init() but is usually named "dumps/rdump[0,1]" 
************************************************************************/
void restart_read(FILE * fp)
{
	int idum, i, j, k, l;

  /*************************************************************
	  READ the header of the restart file: 
  *************************************************************/
	fscanf(fp, "%d", &idum);
	if (idum != N1) {
		fprintf(stderr,
			"error reading restart file; N1 differs\n");
		exit(3);
	}
	fscanf(fp, "%d", &idum);
	if (idum != N2) {
		fprintf(stderr,
			"error reading restart file; N2 differs\n");
		exit(4);
	}

	fscanf(fp, "%lf", &t);
	fscanf(fp, "%lf", &tf);
	fscanf(fp, "%d",  &nstep);
	fscanf(fp, "%lf", &gam);
	fscanf(fp, "%lf", &cour);
	fscanf(fp, "%lf", &DTd);
	fscanf(fp, "%lf", &DTl);
	fscanf(fp, "%lf", &DTi);
	fscanf(fp, "%lf", &DTp);
	fscanf(fp, "%d",  &DTr);
	fscanf(fp, "%d",  &dump_cnt);
	fscanf(fp, "%d",  &image_cnt);
	fscanf(fp, "%d",  &rdump_cnt);
	fscanf(fp, "%lf", &dt);

	fprintf(stderr,"%g %g %d %g %g\n", t,tf,nstep,gam,cour) ;
	fprintf(stderr,"%g %g %g %g %d\n", DTd,DTl,DTi,DTp,DTr) ;
	fprintf(stderr,"%d %d %d %g\n", dump_cnt,image_cnt,rdump_cnt,dt);


  /*************************************************************
	  READ the body of the restart file: 
  *************************************************************/
	ZSLOOP(-NG, N1-1 + NG, -NG, N2-1 + NG) 
		PLOOP fscanf(fp, "%lf", &(p[i][j][k]));

	for(l=0;l<NPTOT;l++)
	for(k=1;k<NDIM;k++) fscanf(fp, "%lf", &(xp[l][k])) ;

	return;
}

#undef FMT_DBL_OUT
#undef FMT_INT_OUT

FILE *get_restart_file(char *restart_filename)
{
	char *header_line = NULL;
	size_t linecap = 0;
	int idum;
	double fdum;
	FILE *fp;

	fp = fopen(restart_filename, "r") ;

	if(fp != NULL) {
		getline(&header_line, &linecap, fp);

		/* 	assumes 
			3rd header field is time
			5th header field is nstep 
		*/
		sscanf(header_line,"%d %d %lf %lf %d",
			&idum, &idum, &t, &fdum, &nstep) ;
		fprintf(stderr,"%s: t = %g, nstep = %d\n", 
			restart_filename, t, nstep) ;

		rewind(fp);
	}

	return(fp);
}
