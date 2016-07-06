/*******************************************************************
A MEX-adapted version of a C program "minfo.c" by Eric R. Weeks

Compiling:
mex -largeArrayDims -O -lut lag.cpp

Adaptations:

Jonathan Hadida, July 2016, University of Oxford, UK
Dmytro S. Lituiev, March 2013, University of Zurich, Switzerland

Handles a 2D arrays, with the signals running along the first dimension.
i.e. in an array [ T x N ] -- N time series, each with T time points

**************************************************
minfo.c -- Eric R. Weeks -- started 2/28/97

does the mutual information algorithm discussed by Fraser & Swinney
(Phys Rev A 33 (1986) p1134-1140)

v01:  2-28-97: taken from shell.c (7/1/96)
			quicksort routine taken from sane.c (original version)
v02:  2-28-97: revised sorting for s[] (different than q[])
			sped up math
v03:  2-28-97: add in lag loop
	 3-01-97: fix for variable number of input; add -b option
v04:  3-01-97: take out chi2 tests for substructure
v05:  3-01-97: realize that with chi2 tests taken out, number()
			function doesn't need to be called very often.  remove
			a[] and b[][] arrays!  Much faster now.

This program is public domain, although please leave my name and
email address attached.

email: weeks@physics.emory.edu
web: http://www.physics.emory.edu/~weeks/

explanation of how to use this program:
    http://www.physics.emory.edu/~weeks/software/minfo.html

*******************************************************************/

#include "mex.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>

// ------------------------------------------------------------------------

// Maximum signal length (longer inputs will be truncated with warning)
#define MAX_LENGTH 32768

// Maximum recursion depth for function ffunct (multiple of 3)
#define MAX_DEPTH 48

// Default maximum lag
#define DEFAULT_LAG_MAX 100

// Math constants
#define PI 3.14159265358979323846264338328
#define EE 2.71828182845904523536

// Assertion helpers
#define ASSERT(P,msg) if (!(P)) { mexErrMsgTxt(msg); exit(1); }
#define IS_DOUBLE_MATRIX(P) (!mxIsComplex(P) && mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_DOUBLE_SCALAR(P) (IS_DOUBLE_MATRIX(P) && mxGetNumberOfElements(P) == 1)

// The index type used throughout
typedef mwSize index_t;

// Detect keyboard interruptions
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
#else
    extern bool utIsInterruptPending();
#endif

// ------------------------------------------------------------------------

double ILOG2;
index_t n0, nbins;
index_t pow2[MAX_DEPTH];

std::vector<index_t> s, q, pindex;
std::vector<double> pop;

double  findisq();
double  ffunct( index_t kmarray[], index_t m );
index_t number( index_t kmarray[], index_t m );
void    saneqsort( index_t p, index_t r );
index_t qpartition_neurons( index_t p, index_t r );

void usage() { mexErrMsgTxt("Bad usage."); }

// ------------------------------------------------------------------------

template <class T, class I = index_t>
struct MatrixWrapper
{
    T *data;
    I nrows, ncols;

    MatrixWrapper()
        { clear(); }
    MatrixWrapper( T *data_, I nrows_, I ncols_ )
        { set(data_,nrows_,ncols_); }

    inline void clear()
        { data = NULL; nrows = ncols = 0; }

    inline void set( T *data_, I nrows_, I ncols_ )
        { data = data_; nrows = nrows_; ncols = ncols_; }

    inline T& operator() ( I r, I c ) const
        { return data[ r + nrows*c ]; }
};



//--------------------     ==========     --------------------//
//--------------------     **********     --------------------//



void mexFunction(	int nargout, mxArray *out[],
					int nargin, const mxArray *in[] )
{
    ILOG2 = 1.0/log(2.0);
    index_t i = 0;

    if ( nargin  <  1 ) { usage(); return; }
	if ( nargout != 1 ) { usage(); return; }

    // Read data
    ASSERT( IS_DOUBLE_MATRIX(in[0]), "Expects a double matrix as first input." );
    MatrixWrapper<const double> X( mxGetPr(in[0]), mxGetM(in[0]), mxGetN(in[0]) );

    // Maximum lag
    index_t lag_max = std::min( static_cast<index_t>(DEFAULT_LAG_MAX), X.nrows/2 );
    if ( nargin >= 2 )
    {
        ASSERT( IS_DOUBLE_SCALAR(in[1]), "Expected a scalar as second input." );
        lag_max = static_cast<index_t>(mxGetScalar( in[1] ));
    }
    ASSERT( lag_max < X.nrows, "Maximum lag cannot be larger than the number of timepoints." );

    // Number of bins
    nbins = -1; // maximum value for unsigned type
    if ( nargin >= 3 )
    {
        ASSERT( IS_DOUBLE_SCALAR(in[2]), "Expected a scalar as third input." );
        nbins = static_cast<index_t>(mxGetScalar( in[2] ));
    }

    // Create output matrix
    out[0] = mxCreateDoubleMatrix( lag_max+1, X.ncols, mxREAL );
    MatrixWrapper<double> B( mxGetPr(out[0]), lag_max+1, X.ncols );

    // ----------  =====  ----------

    // n0 = 2^floor(log2( X.nrows-lag_max ))
	for ( n0 = 1; n0 <= (X.nrows-lag_max); n0 *= 2 ) {} n0 /= 2;

    // If n0 is too large, reduce it
    if ( n0 > MAX_LENGTH )
    {
        mexWarnMsgTxt("Input signal is too long, using truncated version instead (see MAX_LENGTH flag).");
        n0 = MAX_LENGTH;
    }

    // Resize arrays
    s.resize(n0);
    q.resize(n0);
    pindex.resize(n0);
    pop.resize(n0);

    // Loop on columns of the input
    for ( index_t c = 0; c < X.ncols; c++ )
    {
        for ( i=0; i<MAX_DEPTH; i++ )
            pow2[i] = n0 >> i;

        for ( i=0; i<n0; i++ )
        {
            pop[i]    = X(i,c);
            pindex[i] = i;
        }
        saneqsort( 0, n0-1 );

        // WARNING!!!!!! Note that this definition is somewhat opposite
        // of what 'makes sense' but will make number() routine faster
        for ( i=0; i<n0; i++ )
            s[i] = pindex[i]; // s[pindex[i]] = i;

        for ( index_t lag=0; lag<=lag_max; lag++ )
        {
            /* now do lag offset for q[] index_t */
            for ( i=0; i<n0; i++ )
            {
                pop[i]    = X(i+lag,c);
                pindex[i] = i;
            }
            saneqsort( 0, n0-1 );

            // Inverse mapping
            for ( i=0; i<n0; i++ )
                q[pindex[i]] = i;

            /* assume at this time that s[], q[] contain integers from
            * 0 to 2^X.ncols-1.   q[] is based on the time lagged data from
            * X[].  s[] is just X[].
            *
            * example: s[] = 0 4 5 3 6 1 2 7
            *          q[] = 7 4 5 3 6 2 0 1
            *
            *     o.......      so the first partition is s: 0-3 | 4-7
            *     ......o.                                q: 0-3 | 4-7
            *     .....o..
            *     ....o...      with 3 in LL, 1 in UL, 3 in UR, 1 in LR
            *     ...o....      quadrants.
            *     .o......
            *     .......o      second partition is s: 01 | 23 | 45 | 67
            *     ..o.....                          q: 01 | 23 | 45 | 67
            *                   with distribution: 1001
            *				                    0020
            *								1100
            *								0101
            */

            // now find I(S,Q) according to formula (19)
            B(lag,c) = findisq();

            // check for interruptions
            ASSERT( !utIsInterruptPending(), "Keyboard interruption." );
        }
    }

}



//--------------------     ==========     --------------------//
//--------------------     **********     --------------------//



double findisq()
{
	index_t kmarray[MAX_DEPTH];
	double x, y, isq;

	kmarray[0] = 0;
	x = n0;
	y = ffunct(kmarray,0);
	isq = (1.0/x)*y - log(x)*ILOG2;

	return isq;
}

double ffunct( index_t kmarray[], index_t m )
{
	/* THIS FUNCTION CAN CALL ITSELF RECURSIVELY */
    index_t temparray[MAX_DEPTH];

    ASSERT( m < MAX_DEPTH, "Maximum recursion depth reached." );

    // Copy input array
	for ( index_t j=0; j<=m; j++ )
        temparray[j] = kmarray[j];

	index_t n = number(temparray,m);
	double  v = n;

	if (n<=1)  {
		v = 0.0;
	} else if (n==2)  {
		v = 4.0;
	} else if (m==nbins)  {
		/* no substructure */
		v *= log(v)*ILOG2;
	} else {
		/* assume substructure exists */
		v *= 2;
		for ( index_t j=0; j<=3; j++ )
        {
			temparray[m+1] = j;
			v += ffunct(temparray,m+1);
		}
	}

	return v;
}

index_t number( index_t kmarray[], index_t m )
{
	/* THIS FUNCTION IS NOT RECURSIVE */
    index_t i, j, los, his, loq, hiq, ivalue;

	if ( m > 0 )
    {
		los = loq = 0;
		his = hiq = n0;

		for ( i=1; i <= m; i++ )
        {
			if (kmarray[i]%2==0)  his -= pow2[i];
			else                  los += pow2[i];
			if (kmarray[i]<2)     hiq -= pow2[i];
			else                  loq += pow2[i];
		}
		ivalue = 0;
		for ( i=los; i<his; i++ )
        {
			j = q[s[i]];
			if ( (j>=loq) && (j<hiq) ) ++ivalue;
		}
	}
    else ivalue = n0;

	return ivalue;
}



// quicksort, taken from sane-ea.c by David E. Moriarty
// modified for this purpose
void saneqsort( index_t p, index_t r )
{
    index_t q;
    if ( p < r )
    {
        q = qpartition_neurons(p,r);
        saneqsort(p,q);
        saneqsort(q+1,r);
    }
}

// partition function for saneqsort
index_t qpartition_neurons( index_t p, index_t r )
{
    index_t i, j, tempi;
    double x, temp;

    x = pop[p];
    i = p - 1;
    j = r + 1;

    while (true)
    {
        do { --j; } while ( pop[j] < x );
        do { ++i; } while ( pop[i] > x );

        // here's where the in-place swap takes place
        // fortunately pindex[] stores the *original* location
        if ( i < j ) {
            temp = pop[i]; pop[i] = pop[j]; pop[j] = temp;
            tempi = pindex[i]; pindex[i] = pindex[j]; pindex[j] = tempi;
        }
        else return j;
    }
}
