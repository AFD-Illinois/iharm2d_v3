
/*

code to 
1. read in ppm image in x1,x2 (native code) coordinates
2. remap to r,z (cylindrical) coordinates
3. emit a ppm image in r,z coordinates

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SMALL 	1.e-14 ;
#define RED	0
#define GREEN 	1
#define BLUE	2

int oN1, oN2, nN1, nN2;
double dx1, dx2, rin, rout, dx, dy, hslope;

int main(int argc, char *argv[])
{
	int i, j, iold, jold, index;
	double xmax, ymax, dx2, x1max;
	double r, th, x1, x2;
	void ij_to_rth(int i, int j, double *r, double *th);
	void rth_to_x1x2(double r, double th, double *x1, double *x2);
	void get_ppm(char *comment_line, int *N1, int *N2,
		     char orig_image[][3]);
	static char oldimage[4096 * 4096][3];
	char newimage[3];
	char comment_line[100];

	/* get size of image from arguments */
	if (argc < 7) {
		fprintf(stderr, "Usage: \n");
		fprintf(stderr,
			"\t rz_image nN1 nN2 Rin Rout xmax ymax hslope < infile.ppm > outfile.ppm\n");
		fprintf(stderr, "  where \n");
		fprintf(stderr, "\t nN1     = # of new image pixels in x-dir \n");
		fprintf(stderr, "\t nN2     = # of new image pixels in y-dir \n");
		fprintf(stderr, "\t Rin     = Rin used in simulation \n");
		fprintf(stderr, "\t Rout    = Rout used in simulation\n");
		fprintf(stderr,
			"\t xmax    = interpolated image will span x=[0,xmax] \n");
		fprintf(stderr,
			"\t ymax    = interpolated image will span y=[-ymax,ymax] \n");
		fprintf(stderr, "\t hslope  = hslope value from run\n");
		fprintf(stderr, "\n\n ");
		exit(0);
	}
	sscanf(argv[1], "%d", &nN1);
	sscanf(argv[2], "%d", &nN2);
	sscanf(argv[3], "%lf", &rin);
	sscanf(argv[4], "%lf", &rout);
	sscanf(argv[5], "%lf", &xmax);
	sscanf(argv[6], "%lf", &ymax);
	sscanf(argv[7], "%lf", &hslope);

	get_ppm(comment_line, &oN1, &oN2, oldimage);

	x1max = log(rout / rin);
	dx1 = x1max / oN1;
	dx2 = 1. / (double) oN2;

	dx = xmax / (double) nN1;
	dy = 2. * ymax / (double) nN2;

	/* emit PPM header. simple! */
	fprintf(stdout, "P6\n%s\n%d %d\n%d\n", comment_line, nN1, nN2,
		255);

	/* interpolate to new image */
	for (j = nN2 - 1; j >= 0; j--)
		for (i = 0; i < nN1; i++) {
			ij_to_rth(i, j, &r, &th);
			rth_to_x1x2(r, th, &x1, &x2);

			/* set to black outside domain */
			if (x1 < 0. || x1 >= x1max || x2 < 0. || x2 >= 1.)
				newimage[RED] = newimage[GREEN] =
				    newimage[BLUE] = 0;

			/* set to old image pixel at center of zone */
			else {

				iold = (int) (x1 / dx1 - 1.e-20);
				jold = (int) (x2 / dx2 - 1.e-20);

				index = iold + oN1 * jold;

				newimage[RED] = oldimage[index][RED];
				newimage[GREEN] = oldimage[index][GREEN];
				newimage[BLUE] = oldimage[index][BLUE];
			}

			/* and dump it */
			fputc(newimage[RED], stdout);
			fputc(newimage[GREEN], stdout);
			fputc(newimage[BLUE], stdout);
		}

	return (0) ;
}

/* go from pixel on new image to r,th */
void ij_to_rth(int i, int j, double *r, double *th)
{
	double x, y;

	x = i * dx;
	y = (j - nN2 / 2) * dy;

	*r = sqrt(x * x + y * y);
	*th = atan2(x, y);	/* deliberately interchanged args */

	return;
}

/* go from r,th to x1,x2 coordinate */
void rth_to_x1x2(double r, double th, double *x1, double *x2)
{
	double root_find(double th);

	*x1 = log(r / rin);
	*x2 = root_find(th);

	return;
}

/* nonlinear root finding needed to convert x2 to theta */
double root_find(double th)
{
	int i;
	double X2a, X2b, X2c, tha, thb, thc;
	double dthdX2, dtheta_func(double y), theta_func(double y);

	if (th < M_PI / 2.) {
		X2a = 0. - SMALL;
		X2b = 0.5 + SMALL;
	} else {
		X2a = 0.5 - SMALL;
		X2b = 1. + SMALL;
	}
	tha = theta_func(X2a);
	thb = theta_func(X2b);

	/* bisect for a bit */
	for (i = 0; i < 8; i++) {
		X2c = 0.5 * (X2a + X2b);
		thc = theta_func(X2c);

		if ((thc - th) * (thb - th) < 0.)
			X2a = X2c;
		else
			X2b = X2c;

	}

	/* now do a couple of newton-raphson strokes */
	tha = theta_func(X2a);
	for (i = 0; i < 2; i++) {
		dthdX2 = dtheta_func(X2a);
		X2a -= (tha - th) / dthdX2;
		tha = theta_func(X2a);
	}

	return (X2a);
}

/* describes x2 -> theta transformation */
double theta_func(double x2)
{
	return (M_PI * x2 + 0.5 * (1. - hslope) * sin(2. * M_PI * x2));
}
double dtheta_func(double x2)
{
	return (M_PI * (1. + (1. - hslope) * cos(2. * M_PI * x2)));
}

/* read in the old (x1x2) ppm image */
void get_ppm(char *comment_line, int *N1, int *N2, char orig_image[][3])
{
	int i, ppm_ncolors;
	FILE *fp;
	char file_type[100];
	char line[100];

	/* get from stdin */
	fp = stdin;

  	/** get header **/
	/* file type */
	fgets(line, 90, fp);
	sscanf(line, "%s", file_type);
	if (strncmp(file_type, "P6", 2) != 0) {
		fprintf(stderr, "wrong input file type, P6 != %s\n",
			file_type);
		exit(1);
	}

	/* min,max comment line */
	fgets(comment_line, 90, fp);
	fprintf(stderr, "comment_line = %s\n", comment_line);

	/* get N1, N2 */
	fgets(line, 90, fp);
	sscanf(line, "%d %d", N1, N2);
	fprintf(stderr, "N1 = %d, N2 = %d\n", *N1, *N2);

	/* get ppm_ncolors */
	fgets(line, 90, fp);
	sscanf(line, "%d", &ppm_ncolors);
	if (ppm_ncolors != 255) {
		fprintf(stderr, "expected ppm_ncolors = 255, but %d\n",
			ppm_ncolors);
		exit(1);
	}

	for (i = 0; i < *N1 * (*N2); i++) {
		orig_image[i][RED] = fgetc(fp);
		orig_image[i][GREEN] = fgetc(fp);
		orig_image[i][BLUE] = fgetc(fp);
	}

	return;
}

