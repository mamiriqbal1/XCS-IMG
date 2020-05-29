#include "xcsMacros.h"
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_01.hpp>



long seed = 143907; // a number between 1 and _M-1 */
const long _Q = _M/_A; // constant for the random number generator (=_M/_A).
const long _R = _M%_A; // constant for the random number generator (=_M mod _A).

void setSeed(long s)  // sets a random seed in order to randomize the pseudo random generator.
{
    seed=s;
}

long getSeed()
{
    return seed;
}

//boost::rand48 intgen(114665);
//boost::uniform_01<boost::rand48> gen(intgen);

//boost::mt19937 intgen(114665);
//boost::uniform_01<boost::mt19937> gen(intgen);


//double drand()
//{
//    return gen();
//}

double drand()  // returns a floating-point random number generated according to uniform distribution from [0,1)
{
    long hi   = seed / _Q;
    long lo   = seed % _Q;
    long test = _A*lo - _R*hi;
    if (test>0)
    {
        seed = test;
    }
    else
    {
        seed = test+_M;
    }
    return (double)(seed)/_M;
}

int irand(int n)  // returns a random number generated according to uniform distribution from [0,n-1]
{
    int num = (int)(drand()*(float)n);
    while(num == n)
    {
        num = (int)(drand()*(float)n);
    }
    return num;
}
