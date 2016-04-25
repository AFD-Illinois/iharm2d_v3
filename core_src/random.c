
/* 
	tools for returning random variables
	interface to gsl's random number utilities 
*/

#include "decs.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng *r ;

void init_rng(int seed)
{
        r = gsl_rng_alloc(gsl_rng_mt19937);     /* use Mersenne twister */
        gsl_rng_set(r, seed);
}

/* return pseudo-random value between 0 and 1 */
double ran_uniform()
{
        return (gsl_rng_uniform(r));
}

