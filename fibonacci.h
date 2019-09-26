#ifndef MY_FIBONACCI_H
#define MY_FIBONACCI_H


// Few fibonacii generation and query related experiments
// Lekhraj 

namespace fibonacci {

long int fib(const int iPos);       // Recursive version
long int fib_iter(const int iPos);  // Iterative version
long int fib_direct(const int iPos);// Direct calculation using golden ratio forumla
bool isFib(const long int iNum);    // Returns true if given number is in fibonacci series
int getFibIdx(const long int iNum); // Returns position (-1 if not found) of given number in fibonacci series
  
}

#endif