

/*

modified 17 June 2012 CFG
        - eliminate r8 option (never used).
        - eliminate read-in of color map (map.ppm no longer required)
        - insert analytic function for john.pal color map
        - removed failimage capability (probably unwise, but eliminates
                global arrays for failimage) 

*/

#include "decs.h"
#include <ctype.h>

/******************************************************************************
  image_ppm(): 
  -----------
     -- generates a color mapped "image" or pixelated output file following 
        the "raw" PPM  format ("man ppm" for more details). 

     -- color map determined by get_color_map()

    CFG 14 Sept 07: modified to go from almost-min to almost-max

******************************************************************************/
void image_ppm(double *f, char *image_descr)
{
	int i ;
	double max, min ;
	static double q[N1 * N2];
        int red, green, blue ;
	char image_file_name[100];
	FILE *image_file;
	int compare_doubles(const void *a, const void *b);
        void john_pal(double data, double min, double max, 
		int *pRed, int *pGreen, int *pBlue) ;


	/* open file; check */
	sprintf(image_file_name,"images/im_%s_%04d.ppm", image_descr, image_cnt);

	if ((image_file = fopen(image_file_name, "w")) == NULL) {
		fflush(stderr);
		fprintf(stderr, "image(): Cannot open %s !! \n", image_file_name);
		fflush(stderr);
		return;
	}
	/* end open file */

	/*  mapping is in 255 steps from max to min */
	/*  here, set min and max */
	for (i = 0; i < N1 * N2; i++) q[i] = f[i];
                // this bit sorts the image values and lops off the top and bottom
                // end of the distribution to create max and min.
		// qsort is a standard unix utility now
	qsort(q, N1 * N2, sizeof(double), compare_doubles);
	min = q[N1 * N2 / 128];
	max = q[N1 * N2 - N2 * N2 / 128];
	/*  end set min and max */

	/* emit ppm header */
	fprintf(image_file, "P6\n#  min=%g  , max=%g \n%d %d\n%d\n", min, max, N1, N2, 255);
	fflush(image_file);

	for (i = 0; i < N1 * N2; i++) {
                john_pal(f[i], min, max ,&red,&green,&blue) ;
		fputc((char) red, image_file);
		fputc((char) green, image_file);
		fputc((char) blue, image_file);
	}

	/* tidy up */
	fclose(image_file);

	return;
}

/* this is needed for qsort, which is used to set colormap min and max */
int compare_doubles(const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;

	return (*da > *db) - (*da < *db);
}

/* 
color palette, based on John Hawley's john.pal

        input: integer 0-255 
        output: red, green, blue integers 

        author: Bryan M. Johnson
*/

void john_pal(double data, double min, double max, int *pRed, int *pGreen, int *pBlue)
{
  double a, b, c, d, e, f;
  double x, y;
  double max_min = max - min ;

  if(max_min > 0.0) { // trust no one
    x = (data - min)/(max_min) ;
    
    /* ========== Red ============ */
    a = 4.0*x - 1.52549019607844;
    b = 4.52941176470589 - 4.0*x;
    y = a < b ? a : b;
    *pRed = (int)(255.0*y);
    *pRed = *pRed >   0 ? *pRed :   0;
    *pRed = *pRed < 255 ? *pRed : 255;

    /* ========== Green ========== */
    a = 4.0*x - 0.521568627450979;
    b = 2.52549019607844 - 4.0*x;
    c = a < b ? a : b;
    d = 4.0*x - 1.53725490196073;
    e = 3.52941176470581 - 4.0*x;
    f = d < e ? d : e;
    y = c > f ? c : f;
    *pGreen = (int)(255.0*y);
    *pGreen = *pGreen >   0 ? *pGreen :   0;
    *pGreen = *pGreen < 255 ? *pGreen : 255;

    /* ========== Blue =========== */
    a = 4.0*x + 0.498039215686276;
    b = 2.50980392156862 - 4.0*x;
    y = a < b ? a : b;
    *pBlue = (int)(255.0*y);
    *pBlue = *pBlue >   0 ? *pBlue :   0;
    *pBlue = *pBlue < 255 ? *pBlue : 255;

  }
  else {
    *pRed = *pGreen = *pBlue = (data > max ? 255: 0) ;
  }

  return;
}

#undef IRHO
#undef ITT
#undef IBSQ
#undef IGAM
#undef GJSQ
