// Numbers: Template implementation
// Lekhraj

#include <type_traits>
#include <cassert>


//Returns a root of continuous function in given interval
template <typename Func, typename Scalar>
Scalar rootFunction(Func&& F, Scalar a, Scalar b, const Scalar EPS, const int NITER) {
	// Use "bisection method" to find root of F(x) in given interval [a, b]
	//
	// We know that for a function f(x) which is continuous in interval [a,b] and
	// f(a) * f(b) < 0, then there exists a "c" in [a, b] where f(c) = 0
	// To find this "c" (root) where f(c)=0, we proceed iteratively by assigning
	// c to mid of interval [a, b], where a < b. We keep shortening interval based on sign of f(c).
	//
	// Number of iterations to converge within EPS:
	// Since in each iteration we halve the interval: bn - an = (b - a)/2^n
	// When we converge: bn - an <= EPS, therefore 2^n >= (b - a)/EPS
	// 	n >= log((b-a)/EPS)/log2
	
	
	// Pre-conditions
	
	if( F(a) * F(b) > 0 ) throw "f(a) and f(b) must be of different sign!";
	static_assert(std::is_floating_point<Scalar>::value, "Scalar must be floating point type!");
	
	// We want a < b
	if( a > b ) std::swap(a, b);
	int flipSign = F(a) < 0 ? 1 : -1;
	Scalar c = (a+b) * 0.5; // initial guess: mid point
	Scalar fV = F(c);
	for( int i=1; ABS(fV) > EPS && i <= NITER; i++ ) {
		assert( a < c && c < b );
		assert( SIGN(F(a)) != SIGN(F(b)) );
		
		// Keep shortening interval to ensure f(a) and f(b) lie on opposite side		
		if( flipSign * fV > 0 ) // SIGN(f(c)) == SIGN(f(b))
			b = c; 
		else // SIGN(f(c)) == SIGN(f(a))
			a = c;
			
		// Next guess: mid point	
		c = (a+b) * 0.5;
		fV = F(c);
	}
	
	return c;
}

//Returns integral of a continuous function in given interval
template <typename Func, typename Scalar>
Scalar integrateFunction(Func&& F, Scalar a, Scalar b, const Scalar EPS, const int NITER) {
	// Integral of F(X) in interval [a, b] = Area of F(X) in interval [a,b]
	// We use "quadrature method" to integrate F(x) in given interval [a, b]
	// -We divide given interval into "M" quads of base "h"
	// -Area of each quad of base "h" and height "F(mid point of base)" is h * F(mid_point)
	// -Now total area is just sum of these "M" quads
	//
	//  integral of F(x) in interval [a, b] = h * SUM(F(Xm)) for each m = 1, ..., b-a/h
	//  M = (b-a)/h
	//
	// Number of iterations "M" to converge within EPS:
	// EPS <= C * h^2, where C = ((b-a)/24) * max(F''(x)) in interval [a, b]
	
	
	// Pre-conditions
	
	static_assert(std::is_floating_point<Scalar>::value, "Scalar must be floating point type!");
	
	// We want a < b
	if( a > b ) std::swap(a, b);
	
	if( a == b ) // Zero Area
		return (Scalar)0;
		
	Scalar h = (b - a)/NITER;
	Scalar Xm = a + h * 0.5;
	Scalar area = 0;
	for(int i=0; i < NITER; i++) {
		area += F(Xm);
		Xm += h;
	}
	area *= h;
	
	return area;
}
	
	