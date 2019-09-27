// Experiment with numbers
// Lekhraj

#include <limits>
#include <cmath>
#include <iostream>

#include "numbers.h"

//Returns max sum of max sum subsequence
// (start, end) returned is index for max sub subsequence
int max_sum_subarray(const int A[], const int nSize,  int& start, int& end) {
    // Inspired by Kadane algorithm.
    // Kadane's algorithm begins with a simple inductive question: if we know the maximum subarray
    // sum ending at position "i", what is the maximum subarray sum ending at position "i+1"?
    // The answer: Either the maximum subarray sum ending at position "i+1" includes
    // the maximum subarray sum ending at position "i" as a prefix, or it doesn't. 
    // => new_max_ending_here = MAX(current element, max_ending_here + current_element)

    int newstart = 0;
    int max_ending_here = 0, max_so_far = std::numeric_limits<int>::min(); // assigning to INT_MN to handle all negative case
    start = end = 0;
    for(int i=0, x; i < nSize; i++) {
        x = A[i];
        max_ending_here = MAX(x, max_ending_here + x);
        max_so_far = MAX(max_so_far, max_ending_here);
        
        // Rest of the logic is to track the start and end index of the max sum sub sequence
        if( max_ending_here == x ) // Found a possible new sequence
            newstart = i;
          
        if( max_so_far == max_ending_here ) // This sequence is current max, so end it here
            end = i;
          
        if( newstart > start ) // Found a new sequence further down the line
            start = newstart;  
        }
            
    if( start > end ) // Why?
      start = end;
   
    // Max sum subsequence is: A[start:end]
    return max_so_far;       
}

//Returns max backward difference in a sequence
// 
int maxDiff(const int A[], const int nSize) {
   // Trivial case: Treat empty list and single element list as no difference
    if( nSize <= 1 )
        return 0;
    
    // Since maximum element must follow minimum (j > i), we can just track min and diff
    // Idea is to iterate over the list maintaining min and maxdiff we have so far    
    int min_sofar = MIN(A[0], A[1]);
    int maxdiff_sofar = A[1] - A[0];
    for(int i=2; i < nSize; i++) {
        // maxdiff_sofar = MAX(A[i]-min_sofar, maxdiff_sofar)
        if( (A[i] - min_sofar) > maxdiff_sofar )
            maxdiff_sofar = A[i] - min_sofar;
        
        // min_sofar = MIN(A[i], min_sofar)
        if( A[i] < min_sofar ) // Keep a running min 
            min_sofar = A[i];
    }    
            
    return maxdiff_sofar;    
}

bool isSqr(const int iNum)
{
    // Check if input is perfect square
    int iSqrt = sqrt(iNum);
    
    return (iSqrt*iSqrt == iNum);
}


#define MAKE_SORTED_TUPLE(a, b, c, eTup) \
    std::get<0>(eTup) = MIN(MIN(a, b), c) ; \
    std::get<2>(eTup) = MAX(MAX(a, b), c) ; \
    std::get<1>(eTup) =  a + b + c - std::get<0>(eTup) - std::get<2>(eTup);
       

//Returns 3 numbers which when added in pair result in square number in passed in vector
int genSquareTriplets(const int nMax, std::set<std::tuple<int,int,int>>& setTriplet) {
    // We need to find three numbers n1, n2, n3 which when added up in pair make a perfect square
    // a = n1+n2, b = n1+n3, c = n2+n3, where a, b, c all should be perfect square
    
    int nTriplets = 0;
    std::tuple<int, int, int> eTup;
    for(int i=1; i <= nMax; i++) {
        for(int j=1; j <= nMax; j++) {
            for(int k=1; k <= nMax; k++) {
                // std::cout << "I:(" << i << ", " << j << ", " << k << ")\n";
                if( (i != j && i != k && j != k) && isSqr(i+j) && isSqr(i+k) && isSqr(j+k) ) {
                    //std::cout << "O:(" << i << ", " << j << ", " << k << ")\n";
                    MAKE_SORTED_TUPLE(i, j, k, eTup);
                    setTriplet.insert(eTup);
                }               
            }            
        }        
    }
    
    return nTriplets = setTriplet.size();
}

int revNum(const int n) {
    //return parseInt(n.toString().split('').reverse().join(''));
    int N=n, revN = 0, d=N;
    while( N > 0 ) {
        //revN = revN * 10 + N % 10; N /= 10;
        N /= 10;
        d = d - N * 10; // equivalent to N%10 but faster
        revN = revN * 10 + d;
        d = N;
    }
    
    return revN;
}

int ndigitsNum(const int n) {
  if( n < 10 )
    return 1;
  else if( n < 100 )
    return 2;
  else if( n < 1000 )
    return 3;
  else {
      int i = 0, N = n;
      while ( N > 0 ) {
          N /= 10;
          i++;
      }
      return i;      
  }      
}

#define MAKE_SORTED_PAIR(a, b, ePair) \
    if( a <= b ) { \
        ePair.first = a; ePair.second = b; \
    } \
    else { \
        ePair.first = b; ePair.second = a; \
    }


//Returns numbers who are mirror squares (e. 13 <-> 31)
int genMirrorSquares(const int nMax, std::set<std::pair<int,int>>& setMirrorSq) {
    // Mirror squares are those integers whose mirror number square is also a mirror!

int j, jSqr, iSqr;
std::pair<int, int> ePair;

for(int i = 12; i <= nMax; i++) { 
	j = revNum(i);
    if( i != j ) {
      jSqr = j *j; iSqr = i * i;
      if( (ndigitsNum(iSqr) == ndigitsNum(jSqr)) && jSqr == revNum(iSqr) ) {
            MAKE_SORTED_PAIR(i, j, ePair);
        	setMirrorSq.insert(ePair);
        	//printf("(", i, ":", i*i, ", ", j, ":", j*j, ")");
      }
    }
}

return setMirrorSq.size();        
}

#define MABS(x) \
	((x) >= 0 ? (x) : (-(x)))

//Returns sqrt, default parameters EPS=0.001, NITER=25
float bsqrt(const float n, const float EPS, const int NITER) {
    // Babylon method: iterate on x1 = 1/2 * (x0 + n/x0) to get desired accuracy
	//
	// Babylon method is just a special case of newton-raphson method to solve
	// roots of f(x) = 0, in case of sqrt, f(x) = x^2 - n 
	//
	// In Newton Rapshon method, for a function f(x), we calulate the equation
	// of line passing through initial point (x0, y0=f(x0)) and having a slope of
	// m = f'(x0), y = m(x-x0) + y0, we find the root of this line by putting 
	// y = 0 => x1 = x0  - f(x0)/f'(x0)
	// Since for sqrt f(x) = x^2 - n => f'(x) = 2x, so our above newton equation
	// becomes => x1 = x0 - (x0^2 - n)/2x0, which when further simplified leads to
	// babylon iterative equation => x1 = 1/2 * (x0 + n/x0)
	//
    
    float x = n * 0.5; // our initial guess
    float g = 0.0;
    int i = 1;
    while( MABS(x - g) > EPS && ++i <= NITER ) { // x != g, keep looping until guesses start to converge
    	g = x; // save previous guess
        x = (x + n/x) * 0.5; // babylon iteration for next guess
    }
    
    return x;    
}

bool isOdd(const int n) {
	return (n&01 != 0);
}

bool isSqrInt(int n) { 
    for (int i = 1; i * i <= n; i++) {   
        // Check (i * i == n) 
        if ((n % i == 0) && (n / i == i)) { 
            return true; 
        } 
    } 
    return false; 
}

// Bisection method:  default parameters EPS=0.001, NITER=25
float msqrt(const float n, const float EPS, const int NITER) { 
	// binary chop: bisection method
	float m = n;
    float mSqr = m * m, lo = 0, hi = n;

	bool sqrNum = false, oddNum = false;
	
	if( n == (float)((int)n) ) { // is given number full integer
	   sqrNum = isSqrInt((int)n), oddNum = isOdd((int)n);
    }

	// For numbers less then 1, sqrt will be between n..1 range
    if( n < 1 ) { lo = n; hi = 1; }
  
    for(int i=1; MABS(mSqr - n) > EPS && ++i <= NITER; mSqr = m*m) {
		// keep adjusting bisection range so that result falls in between lo..hi
    	if( mSqr > n ) hi = m; else lo = m;
		
		// next guess: midlle point of range lo..hi
        m = (lo+hi) * 0.5;
        
		// Round off to two decimal digits for smaller numbers	
	   	if( m < 1 ) {m = round(m * 100)/100.0;}
	   
	   	// Trick to get perfect sqrt for a perfect square
       	if( sqrNum ) {
        	if( oddNum ) m = round(m); else m = (int)m;
       }
    }
    
	// m ~= sqrt(n)
    return m;
}

//Returns kth root of given number, for sqrt:k=2, cuberoot:k=3, etc
//Default parmeters, k=2, EPS=0.001, NITER=100
float kthRoot(const double n, const int k, const float EPS, const int NITER) {
    // Generalized babylon method: Iterate on x1 = 1/k * ((k-1)*x0 + n/(x0^(k-1)) to get desired accuracy
	//
	// Generalized babylon method to fin kth root is just a special case of
	// newton-raphson method to solve roots of f(x) = 0, in case of kth root,
	// f(x) = x^k - n 
	//
	// In Newton Rapshon method, for a function f(x), we calulate the equation
	// of line passing through initial point (x0, y0=f(x0)) and having a slope of
	// m = f'(x0), y = m(x-x0) + y0, we find the root of this line by putting 
	// y = 0 => x1 = x0  - f(x0)/f'(x0)
	// Since for kth root f(x) = x^k - n => f'(x) = kx^(k-1), so our above newton equation
	// becomes => x1 = x0 - (x0^k - n)/kx0^(k-1), which when further simplified leads to
	// babylon iterative equation => x1 = 1/k * ((k-1)*x0 + n/x0^(k-1))
	//
    
	// Initial guess is quite important for newton method to converge especially for
	// higher kth roots, otherwise it may take large number of iterations.
	// Easy initial guess: n/k
	// Fast converging guess: 2^(log(n)/k)
    double x1 = 2 << (int)ceil(ceil(log2(n))/k); // our initial guess
    double x0 = 0.0;
    int i = 1;
    while( MABS(x1 - x0) > EPS && ++i <= NITER ) { // x != g, keep looping until guesses start to converge
    	x0 = x1; // save previous guess
        x1 = ((k-1) * x0 + n/pow(x0, k-1)) / k; // babylon iteration for next guess
    }
    
    return x1;    	
}




