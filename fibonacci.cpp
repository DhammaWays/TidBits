// Fibonacci Numbers
// Lekhraj
 
#include <cmath>

#include "fibonacci.h"


long int fibonacci::fib(const int iPos)
{
    // recursively calculate fibonacci number
    if( iPos <= 0 ) return 0;
    else if ( iPos == 1 ) return 1;
    else
        return fib(iPos-2) + fib(iPos-1);
}

long int fibonacci::fib_iter(const int iPos)
{
    // Iterative version of fibonacci series
    // fib(n) = fib(n-2) + fib(n-1)
    
    long int iF0 = 0, iF1 = 1, iFn = 1;
    if( iPos <= 0 ) return 0;
    else {
        for( int i = 2; i <= iPos; i++ ) {
            iFn = iF0 + iF1;
            iF0 = iF1;
            iF1 = iFn;
        }
        return iFn;
    }
}

long int fibonacci::fib_direct(const int iPos)
{
    // Given golden ratios, one can calculate fibonacci number directly
    // Phi = (1+sqrt(5)/2) and Phi2 = (1-sqrt(5))/2
    // Fibonacci(n) = (Phi^n - Phi2^n)/sqrt(5)
    
    //long double fPhi = 1·61803398874989484820L;
    //long double fPhi2 = -0·61803398874989484820L;
    //long double fInvSqrt5 = 0.44721359549995794L;
    const double fPhi = 1.61803398875;
    const double fPhi2 = -0.61803398875;
    const double fInvSqrt5 = 0.447213595499958;
    
    
    return (long int)((pow(fPhi, iPos) -  pow(fPhi2, iPos)) * fInvSqrt5);   
}

bool isSqr(long int iNum)
{
    // Check if input is perfect square
    long int iSqrt = sqrt(iNum);
    
    return (iSqrt*iSqrt == iNum);
}

bool fibonacci::isFib(const long int iNum)
{
    // Check if given number is fibonacci number
    // n is fibonacci only if either 5n^2 + 4 or 5n^2 - 4 is perfect square
    // Watchout: Suqaring large number can overflow!
    long int iSqr = 5 * iNum * iNum;
    
    return (isSqr(iSqr+4) || isSqr(iSqr-4));    
}

int fibonacci::getFibIdx(const long int iNum)
{
    // Check if given input is fibonacci number, return its index else -1 if not
    // We try to calculate nearest position of given number in fibonacci series
    // iPos = (log(n) + log(sqrt(5))) / log((1+sqrt(5))/2)
    // Then we check if given number is either at "iPos" or at "iPos+1" in series
    // (to be on safe side due to rounding of errors)
    
    double fLog5by2 = 0.804718956217050;
    double fInvLogPhi = 2.07808692124;
    
    int iPos = (log(iNum) + fLog5by2) * fInvLogPhi;
    
    if( fib_direct(iPos) == iNum ) return iPos;
    else if( fib_direct(iPos+1) == iNum ) return iPos+1;
    else return -1;
}
 