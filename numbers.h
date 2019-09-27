#ifndef MY_NUMBERS_H
#define MY_NUMBERS_H

// Few experiments with numbers
// Lekhraj

#include <set>
#include <tuple>

//Simple macros to return MAX, MIN value
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define MIN(a,b) ((a) <= (b) ? (a) : (b))

//Returns max sum of max sum subsequence
// (istart, iEnd) is max sub subsequence
int max_sum_subarray(const int aiElem[], const int nSize, int& iStart, int& iEnd);

//Returns max backward difference in a sequence
// 
int maxDiff(const int aiElem[], const int nSize);

//Returns 3 numbers which when added in pair result in square number in passed in vector
int genSquareTriplets(const int nMax, std::set<std::tuple<int,int,int>>& setTriplet);

//Returns numbers who are mirror squares (e. 13 <-> 31)
int genMirrorSquares(const int nMax, std::set<std::pair<int,int>>& setMirrorSq);

//Returns sqrt
float bsqrt(const float n, const float EPS=0.001, const int NITER=25);

//Returns sqrt
float msqrt(const float n, const float EPS=0.001, const int NITER=25);

//Returns kth root of given number, for sqrt:k=2, cuberoot:k=3, etc
float kthRoot(const double n, const int k=2, const float EPS=0.001, const int NITER=100);

#endif