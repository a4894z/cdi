
/* 
 * to compile, I use: 
 * mex -O CFLAGS="\$CFLAGS -std=c99" random_uniform.c -lgsl -lgslcblas -lm  
 *
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "mex.h"

/* Input Arguments */

#define NUM_AMT_IN  prhs[0] 
#define LOW_LIMIT_IN  prhs[1] 
#define HIGH_LIMIT_IN prhs[2] 


/* Output Arguments */

#define	UNIF_OUT    plhs[0]

/************************************************************************************************/

static void random_unif(float rand_unif_out[], const gsl_rng_type * gsl_type, unsigned int num_amt, float low_limit, float high_limit ) {

    gsl_rng * rng;
      
    int i;
    unsigned int seed;
    //int * lam_poiss = malloc( num_amt  * sizeof(int));

    srand(time(NULL));                  // create random seed based on time
    seed  = rand(); 
  
    rng = gsl_rng_alloc( gsl_type );     // initialize RNG
    gsl_rng_set( rng, seed );   

    for (i = 0; i < num_amt; i++)
    {
      rand_unif_out[i] = (float) gsl_ran_flat( rng, (double) low_limit, (double) high_limit );
    }

    return;
} 

/************************************************************************************************/

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    float *rand_unif_out; 
    float low_limit, high_limit;  
    unsigned int num_amt; 
    
    ///*
    if (mxIsUint64(NUM_AMT_IN) != 1)
      mexErrMsgTxt(" Size of number of random unif numbers must be uint64!");
    //*/

    mwSize mLam, nLam;

    mLam = mxGetM(NUM_AMT_IN); // get rows 
    nLam = mxGetN(NUM_AMT_IN); // get cols

    if ( (nLam != 1) || (nLam != 1) )
	mexErrMsgTxt("Output array size of random unif dist must be scalar!");
  
    num_amt = mxGetScalar( NUM_AMT_IN ); // Assign pointer to input parameter 
  
   /* 
    if (num_amt < 0)
      mexErrMsgTxt("Output array size must be positive!");
    */

    mLam = mxGetM(LOW_LIMIT_IN); // get rows 
    nLam = mxGetN(LOW_LIMIT_IN); // get cols

    if ( (nLam != 1) || (nLam != 1) )
	mexErrMsgTxt("Input size of low limit must be scalar!");
  
    low_limit = mxGetScalar( LOW_LIMIT_IN );    // Assign pointer to input parameter 

    mLam = mxGetM(HIGH_LIMIT_IN); // get rows 
    nLam = mxGetN(HIGH_LIMIT_IN); // get cols

    if ( (nLam != 1) || (nLam != 1) ) 
	mexErrMsgTxt("Input size of high limit must be scalar!");
  
    high_limit = mxGetScalar( HIGH_LIMIT_IN );  // Assign pointer to input parameter  
            
    mwSignedIndex dims[2];
    dims[0] = num_amt;
    dims[1] = 1;
    UNIF_OUT = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    rand_unif_out  = mxGetData(UNIF_OUT);
  
    random_unif( rand_unif_out, gsl_rng_taus, num_amt, low_limit, high_limit );

    return;
    
}

/************************************************************************************************/
