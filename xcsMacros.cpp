#include "xcsMacros.h"
#include <random>
#include <iostream>

std::mt19937_64 mt;
std::uniform_real_distribution<double> real_distribution(0.0, 1.0);
std::uniform_int_distribution<int> integer_distribution;


void initialize_random_number_generator(bool fixed_seed)
{
    unsigned int seed = 0;
    if(fixed_seed) {
        seed = 51973;
    }else{
        std::random_device rd;
        seed = rd();
    }
    mt.seed(seed);
}


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

double drand()
{
    return real_distribution(mt);
}

double drand_old()  // returns a floating-point random number generated according to uniform distribution from [0,1)
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


int irand(int n)
{
    std::uniform_int_distribution<int>::param_type parameters(0, n-1);
    integer_distribution.param(parameters);
    return integer_distribution(mt);
}

int irand_old(int n)  // returns a random number generated according to uniform distribution from [0,n-1]
{
    int num = (int)(drand()*(float)n);
    while(num == n)
    {
        num = (int)(drand()*(float)n);
    }
    return num;
}
