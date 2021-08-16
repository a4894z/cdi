
/* 
 * to compile, I use: 
 * mex -O CFLAGS="\$CFLAGS -std=c99" random_poisson.c -lgsl -lgslcblas -lm  
 *
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "mex.h"

/* Input Arguments */

#define	LAMBDA_IN   prhs[0]

/* Output Arguments */

#define	PN_OUT    plhs[0]

/************************************************************************************************/

static void generate_poisson_noise (float noisy_signal_out[], const gsl_rng_type * gsl_type, unsigned int num_amt, float llambda[] ) {

    gsl_rng * rng;
      
    int i;
    unsigned int seed;
    //int * lam_poiss = malloc( num_amt  * sizeof(int));

    srand(time(NULL));                  // create random seed based on time
    seed  = rand(); 
  
    rng = gsl_rng_alloc (gsl_type);     // initialize RNG
    gsl_rng_set (rng,seed);   
    
    for (i = 0; i < num_amt; i++)
    {
      noisy_signal_out[i] = (float) gsl_ran_poisson( rng, (double) llambda[i] );  
    }

    return;
} 

/************************************************************************************************/

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    float *noisy_signal_out; 
    float *llambda;  
    unsigned int num_amt; 

    mwSize mLam, nLam;
    mLam = mxGetM(LAMBDA_IN); // get rows 
    nLam = mxGetN(LAMBDA_IN); // get cols

    if (nLam != 1)
	mexErrMsgTxt("Input Poisson means must be in a column vector!");
  
    llambda = mxGetData(LAMBDA_IN);

    mwSignedIndex dims[2];
    dims[0] = mLam;
    dims[1] = nLam;
    PN_OUT = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    noisy_signal_out = mxGetData(PN_OUT);

    num_amt = (unsigned int) mLam;  // size of the lambda input vector 
 
    generate_poisson_noise( noisy_signal_out, gsl_rng_taus, num_amt, llambda );
    //generate_poisson_noise( noisy_signal_out, gsl_rng_gfsr4, num_amt, llambda );
    //generate_poisson_noise( noisy_signal_out, gsl_rng_mt19937, num_amt, llambda );
    //generate_poisson_noise( noisy_signal_out, gsl_rng_ranlxd1, num_amt, llambda );
     


    return;
    
}

/************************************************************************************************/
