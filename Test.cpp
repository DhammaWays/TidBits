// Test Driver
// Lekhraj
 
#include <iostream>
#include <string>
#include <cstdlib>

#include "timer.h"
#include "fibonacci.h"
#include "sort.h"
#include "numbers.h"
#include "lstrings.h"

using namespace measure;
using namespace fibonacci;

using namespace std;


void TestFib() {
  // Test Fibonacci number functions  
  std::string strIn;
  int iPos;
  long int iNum;
  
  // Direct method starts to give inaacurat ereuslt beyond 56 due to real numbers precision 
  iPos = 41;
  
  timer tim;
  
  std::cout <<"\nTesting Fibonacci Numbers...\n";
  
  std::cout << "Fibonacci( " << iPos << " ):" << std::endl;
  
  //tim.start();
  //iNum = fib(iPos);
  //std::cout << "Recursive: " << iNum << " Time: " << tim.elasped() << std::endl;
 
  tim.start();
  iNum = fibonacci::fib_iter(iPos);
  std::cout << "Iterative: " << iNum << " Time: " << tim.elapsed() << std::endl;

  tim.start();
  iNum = fib_direct(iPos); 
  std::cout << "Direct: " << iNum << " Time: " << tim.elapsed() << std::endl;

  //std::cin >> iPos >> iNum;
  
  //std::cout << "Enter position for fibonacci series: ";
  //std::getline(std::cin, strIn); 
  //std::stoi(strIn, iPos);
  //std::cout << "Fibonacci: " << fib(iPos) << " " << fib_iter(iPos) << " " << fib_direct(iPos) << std::endl;
  
  
  //std::cout << "Enter number: ";
  //std::getline(std::cin, strIn); 
  //std::stoi(strIn, iNum);
  std::cout << "Is " << iNum << " a fibonacci number? " << isFib(iNum) << " its position in series: " << getFibIdx(iNum) << std::endl;
  std::cout << "Is " << iNum+1 << " a fibonacci number? " << isFib(iNum+1) << " its position in series: " << getFibIdx(iNum+1) << std::endl;
     
}

#define PRINT_ARRAY(arElem) \
    cout << "[ "; \
    for (const auto& value: (arElem)) \
		std::cout << value << ' '; \
	cout << "] ";
	
#define PRINT_ARRAY_RANGE(arElem, s, e) \
    cout << "[ "; \
    for (int i=s; i <= e; i++ ) \
		std::cout << (arElem)[i] << ' '; \
	cout << "] ";

#define FIND_VALUE(arElem, SIZE, Elem) \
   cout << "Value=" <<(Elem) << " is found at index=" << bsearch((arElem), SIZE, (Elem)) << endl;

#define ARRAY_SIZE(arElem) sizeof(arElem) / sizeof(arElem[0])

#define RANDOM_NUM(nBeg, nEnd) (rand() % ((nEnd)-(nBeg)+1) + (nBeg))


void TestSort() {
    timer tim;
    const int aiElem[] = {1, 7, 9, 12, 15, 21};
    const double afElem[] = {1.5, 7.0, 9.25, 12.1, 15.0, 21.7};
    const string asElem[] = {"A", "Abacus", "Acid", "Amul", "Ant", "Aeroplane"};
    const unsigned int iSIZE = ARRAY_SIZE(aiElem);
    
    std::cout <<"\nTesting Searching, Sorting...\n";
    PRINT_ARRAY(aiElem); FIND_VALUE(aiElem, iSIZE, aiElem[1]);
    PRINT_ARRAY(afElem); FIND_VALUE(afElem, iSIZE, afElem[iSIZE-1]);
    PRINT_ARRAY(asElem); FIND_VALUE(asElem, iSIZE, asElem[3]);
    
    int toSortAr[] = {21, 4, 1, 3, 9, 20, 25, 6, 21, 14, -5, 42, 30};
    const unsigned int iSORT_SIZE = ARRAY_SIZE(toSortAr);
    
    std::cout << "Sorting "; PRINT_ARRAY(toSortAr); std::cout << std::endl;
    qsort(toSortAr, iSORT_SIZE);
    PRINT_ARRAY(toSortAr); FIND_VALUE(toSortAr, iSORT_SIZE, toSortAr[7]);
    
    string toSortAr2[] = {"Shyam", "Radha", "Mohan", "Ram", "Sita", "Janaki", "Krishna"};
    const int iSize2 = ARRAY_SIZE(toSortAr2);
    
    std::cout << "Sorting "; PRINT_ARRAY(toSortAr2); std::cout << std::endl;
    isort(toSortAr2, iSize2);
    PRINT_ARRAY(toSortAr2); FIND_VALUE(toSortAr2, iSize2, toSortAr2[3]);
    
    const int nMAX = 100;
    int aiRandomElem[nMAX];
    int aiRandom2[nMAX];
    for(int i=0; i < nMAX; i++) {
        aiRandomElem[i] = RANDOM_NUM(1,nMAX);
        aiRandom2[i] = aiRandomElem[i];
    }
    
    cout << "Sorting "; PRINT_ARRAY(aiRandomElem); cout << endl;
    tim.start();
    isort(aiRandomElem, 15);
    cout << "Insertion Sort Time (small sequence): " << tim.elapsed(); cout << endl; 
    isort(aiRandomElem, nMAX);
    cout << "Insertion Sort Time (large sequence): " << tim.elapsed(); cout << endl; 
       
    PRINT_ARRAY(aiRandomElem); cout << endl;
    
    
    cout << "Sorting "; PRINT_ARRAY(aiRandom2); cout << endl;
    tim.start();
    qsort(aiRandom2, 15); 
    cout << "Quick Sort Time (small sequence): " << tim.elapsed(); cout << endl;   
    tim.start();
    qsort(aiRandom2, nMAX);
    cout << "Quick Sort Time (large sequence): " << tim.elapsed(); cout << endl;   
    PRINT_ARRAY(aiRandom2); cout << endl;
}

void TestNum() {
    int A[] = {-2, 1, -3, 4, -1, 2, 1, -5, 4}; // SubArray: [4, -1, 2, 1], Sum: 6
    int B[] = {10, -7, -4, 9, -2, -1, 5, -10, 7, 2};
    int C[] = {8, 3, 8, -5, 4, 3, -4, 3, 5};
    int D[] = {-4, -1, -3};    
    
    cout << "\nTesting numbers...\n";
    
    int sum, start, end;
    sum = max_sum_subarray(A, ARRAY_SIZE(A), start, end); 
    PRINT_ARRAY(A); cout << "Max Sum: " << sum << " Sub array: "; PRINT_ARRAY_RANGE(A, start, end); cout << endl;

    sum = max_sum_subarray(B, ARRAY_SIZE(B), start, end); 
    PRINT_ARRAY(B); cout << "Max Sum: " << sum << " Sub array: "; PRINT_ARRAY_RANGE(B, start, end); cout << endl;

    sum = max_sum_subarray(C, ARRAY_SIZE(C), start, end); 
    PRINT_ARRAY(C); cout << "Max Sum: " << sum << " Sub array: "; PRINT_ARRAY_RANGE(C, start, end); cout << endl;
    
    sum = max_sum_subarray(D, ARRAY_SIZE(D), start, end); 
    PRINT_ARRAY(D); cout << "Max Sum: " << sum << " Sub array: "; PRINT_ARRAY_RANGE(D, start, end); cout << endl;

    int E[] = {3, 15, 1, 3, 7, 14};
    int F[] = {1, 10, 2, 15, 1, -1, 10, 0, -2 , 3, 15};
    int G[] = {5, -2, -1, -1};
    
    int diff;
    diff = maxDiff(E, ARRAY_SIZE(E)); 
    PRINT_ARRAY(E); cout << "Max Diff: " << diff << endl;

    diff = maxDiff(F, ARRAY_SIZE(F)); 
    PRINT_ARRAY(F); cout << "Max Diff: " << diff << endl;
    
    diff = maxDiff(G, ARRAY_SIZE(G)); 
    PRINT_ARRAY(G); cout << "Max Diff: " << diff << endl;
    
    int nSet, n1, n2, n3;
    std::set<std::tuple<int,int,int>> vecSet;
    nSet = genSquareTriplets(100, vecSet);
    cout << "\nThree numbers when added in pairs sum to a perfect square are:\n" << nSet << endl;
    for(const auto& val: vecSet) {
        std::tie(n1, n2, n3) = val;
        cout << "(" << n1 << ", " << n2 << ", " << n3 << ")\n";
    }
    
    std::set<std::pair<int,int>> pairSet;
    nSet = genMirrorSquares(5000, pairSet);
    cout << "\nMagical square numbers are: " << nSet << endl;
    for(const auto& val: pairSet) {
        std::tie(n1, n2) = val;
        cout << n1 << ":" << n2 << " :: " << n1*n1  << ":" << n2*n2 << "\n";
    }
    
    float numAry[] = {0.81, 0.25, 0.1, 4, 9, 100, 169, 625, 1024, 2, 3, 5, 10, 1000, 999*999};

    cout << "\nSquare Roots:" << endl;
    for( int i = 0; i < ARRAY_SIZE(numAry); i++ ) {
    	cout << "sqrt(" << numAry[i] << ") = " << bsqrt(numAry[i]) << ":" << msqrt(numAry[i]) << endl;    	
    }


    
}

void TestStr() {
    string s1 = "Dhamma Path";
    string s2 = "Kalam";
    
    cout << "\nTesting strings...\n";
    cout << s1 << " (reverse)=> " << rev(s1) << endl;
    cout << rev(s1) << " (reverse)=> " << rev(rev(s1)) << endl;
    cout << s2 << " (reverse)=> " << rev(s2) << endl;
    cout << rev(s2) << " (reverse)=> " << rev(rev(s2)) << endl;
    
    string s3 = "xabbay";
    string s4 = "xyabbahelloZZollehzzz";
    string s5 = "xaba";
    string s6 = "abcbax";
    string s7 = "abbaabba"; 
    cout << s3 << " (longest palindrome)=> " << palindrome(s3) << endl;
    cout << s4 << " (longest palindrome)=> " << palindrome(s4) << endl;
    cout << s5 << " (longest palindrome)=> " << palindrome(s5) << endl;
    cout << s6 << " (longest palindrome)=> " << palindrome(s6) << endl;
    cout << s7 << " (longest palindrome)=> " << palindrome(s7) << endl;
}


// Test Driver
int main()
{ 
  TestFib();
  TestSort();
  TestNum();
  TestStr();
   
  return 0;
}