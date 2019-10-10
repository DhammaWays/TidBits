// Test Driver
// Lekhraj
 
#include <iostream>
#include <string>
#include <cstdlib>
#include <math.h>

#include "timer.h"
#include "fibonacci.h"
#include "sort.h"
#include "numbers.h"
#include "lstrings.h"
#include "geometry.h"

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

    double kAry[] = {64, 64, 81, 3125, 6*6*6*6*6*6, 7*7*7*7*7*7*7, 9.22744694e+15};
	
    cout << "\n kth Root:" << endl;
    for( int i = 0; i < ARRAY_SIZE(kAry); i++ ) {
    	cout << "kthRoot(" << kAry[i] << ", " << i+2 << ") = " << kthRoot(kAry[i], i+2) << endl;    	
    }
	
	vector<float> myPoly[] = {{1, 0, -625}, {1, 0, 0, -4913}, {1, -2, 1}, {9, 27, -162}, 
							  {2, 31, 137, 180}, { 1, 2, -5, 1}, {10, 0, 2, -1, 0, -3},
							  {5, 3, -2, -1, -10}, {7, 0, 5, 0, -3, 0, 1, 7777}};
	cout << "\n Roots of polynomial:" << endl;
	for( int i = 0; i < ARRAY_SIZE(myPoly); i++ ) {
		PRINT_ARRAY(myPoly[i]); cout << " : " << largestPolyRealRoot(myPoly[i]) << endl;
		}
		
	cout << "\n Roots of function:" << endl;
	auto F1 = [] (const double x) { return sin(x)*cos(x); };
	auto F2 = [] (const double x) { return 3*x + sin(x) - exp(x); };
	auto F3 = [] (const double x) { return exp(x) * cos(x) - x * sin(x); };
	auto F4 = [] (const double x) { return x * exp(-x*x); };
	auto F5 = [] (const double x) { return pow((x-0.5), 3); };
	auto F6 = [] (const double x) { return x*x*x - 125;};
	
	cout << "\n sin(x)*cos(x) has a zero at " << rootFunction(F1, 1., 2.) << endl;
	cout << "\n 3x + sin(x) - e^x has a zero at " << rootFunction(F2, 0.0, 0.5) << endl;
	cout << "\n e^x * cos(x) - x * sin(x) has a zero at " << rootFunction(F3, 0., 3.) << endl;
	cout << "\n x * e^(-x^2) has a zero at " << rootFunction(F4, -1., 1.) << endl;
	cout << "\n (x - 0.5)^3 has a zero at " << rootFunction(F5, -1., 1.) << endl;
	cout << "\n (x^3 - 125) has a zero at " << rootFunction(F6, 0., 10.) << endl;
	
	cout << "\n Next number in series:" << endl;
	std::vector<std::vector<int>> vecSeries = { {10, 20, 30, 40}, {2, 6, 12, 20},
		{1, 8, 27, 64}, {0, 3, 8, 15}, {2, 9, 28, 65}, {1, 4, 9, 16},
		{80, 63, 48, 35}, {1, 2, 6, 15}, {4, 18, 48, 100}, {30, 68, 130, 222},
		{17, 19, 25, 37}, {11, 31, 69, 131}, {4, 16, 36, 64},
		{720, 120, 24, 6}, {5, 9, 17, 33}, {3, 6, 12, 15}, {6, 10, 18, 34},
		{3, 5, 14, 48} };
	for(const auto& series: vecSeries) {
		PRINT_ARRAY(series); cout << ": " << nextNumInSeries(series) << endl;
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

void TestGeom() {
    std::array<double, 3> arVec = {0, 0, 1};
	double arV[] = {-1, 1};
	Vec v1, v2(0,0), v3(1,0,0), v4(0,1), v5(1,1,1), v6(arVec),
	    v7(reinterpret_cast<std::array<double, 2>&>(arV)),
		v8(Vec(0,1,0), 2);
	
	cout << "\nTesting geometry...\n";
	cout << v1 << "==" << v2 << " : " << (v1==v2) << endl;
	cout << v3 << " magnitude is " << (v3.mag()) << endl;
	cout << v3 << "." << v4 << " : " << (v3*v4) << endl;
	cout << v3 << "x" << v4 << " : " << (v3.cross(v4)) << endl;
	cout << v5 << " unit vector is " << (v5.unit()) << endl;
	cout << v3 << "+" << v4 << " : " << (v3+v4) << endl;
	cout << v3 << "-" << v4 << " : " << (v3-v4) << endl;
	cout << v3 << "+ 2 * " << v4 << " : " << (v3 + 2 * v4) << endl;
	cout << (v3+v4+v6) << "==" << v5 << " : " << ((v3+v4+v6) == v5) << endl;
	cout << "Is angle between " << v5 << " and " << v3 << " acute? : " << v5.isAcute(v3) << endl;
	cout << "Is angle between " << v7 << " and " << v3 << " acute? : " << v7.isAcute(v3) << endl;
	cout << "Projection of " << v8 << " on " << v3 << " : " << v8.project(v3) << endl;
	cout << v8 << " is in anticlockwise direction of " << v3 << " : " << v3.isAntiClock(v8) << endl;
	
	Line l1({0, 0, 0}, {1, 0, 0}), l2({0, 1, 0}, {1, 1, 0}), l3({0, -1, 0}, {1, 1, 0}),
	     l4({1,1,0}, {4,5,0}), l5({2,3,0}, {4,3,0}), l6({3,3,0}, {5,3,0}), l7({3,2,0}, {6,2,0}),
		 l8({0,0,0},{0,2,0}), l9({2.5,3,0}, {0,3,0});
	
	cout << l1 << " intersects " << l2 << " : " << l1.doesIntersect(l2) << endl;
	cout << l1 << " intersects " << l3 << " : " << l1.doesIntersect(l3) << endl;
	cout << l4 << " intersects " << l5 << " : " << l4.doesIntersect(l5) << endl;
	cout << l5 << " intersects " << l6 << " : " << l5.doesIntersect(l6) << endl;
	cout << l4 << " intersects " << l7 << " : " << l4.doesIntersect(l7) << endl;
	cout << l5 << " intersects " << l4 << " : " << l5.doesIntersect(l4) << endl;
	cout << l6 << " intersects " << l5 << " : " << l6.doesIntersect(l5) << endl;
	cout << l7 << " intersects " << l4 << " : " << l7.doesIntersect(l4) << endl;
	cout << l2 << " intersects " << l4 << " : " << l2.doesIntersect(l4) << endl;
	cout << l8 << " intersects " << l2 << " : " << l8.doesIntersect(l2) << endl;
	cout << l4 << " intersects " << l9 << " : " << l4.doesIntersect(l9) << endl;
	
	Point x1({0,0}), x2({0.5,0}), x3({0.5,0.5}), x4({1,1}), x5({3,3}), x6({2,3}), x7({2.5,3});
	cout << endl;
	cout << x1 << " is within " << l3 << " : " << l3.isInside(x1) << endl;
	cout << x2 << " is within " << l3 << " : " << l3.isInside(x2) << endl;
	cout << x3 << " is within " << l3 << " : " << l3.isInside(x3) << endl;
	cout << x4 << " is within " << l3 << " : " << l3.isInside(x4) << endl;
	cout << x5 << " is within " << l5 << " : " << l5.isInside(x5) << endl;
	cout << x6 << " is within " << l4 << " : " << l4.isInside(x6) << endl;
	cout << x7 << " is within " << l4 << " : " << l4.isInside(x7) << endl;
	
	Polygon p1({{0, 0, 0}, {4, 0, 0}, {3, 3, 0}}), p2({{5, 2, 0}, {10, 2, 0}, {12, 5, 0}, {6.1, 3.1, 0}}),
	        p3({{3, 2, 0}, {6, 3, 0}, {7, 7, 0}}), p4({{1,0,0}, {2,1,0}, {3,2,0}, {2,3,0}, {0,1,0}}),
			p5({{1,0,0}, {3,1,0}, {4,1,0}, {1,3,0}}), p6({{1,0,0}, {3,0,0}, {2,1,0}, {3,2,0}, {0,3,0}}),
			p7({{0,1,0}, {2,4,0}, {1,2,0}, {3,0,0}}), p8({{-1, -1,0}, {5,0,0}, {5, 4, 0}, {0, 4,0}}),
			p9({{1,-1,0}, {5,2,0}, {0,2,0}});
	Point p({-1, 1, 0}), q({2, 2, 0}), r({3, 2, 0}), s({-1,1,0}), t({1,1,0}), u({0,1,0}), v({2,1,0});
	
	cout << endl;
	cout << p << " is inside " << p1 << " : " << p1.isInside(p) << endl;
	cout << q << " is inside " << p1 << " : " << p1.isInside(q) << endl;
	cout << r << " is inside " << p1 << " : " << p1.isInside(r) << endl;
	cout << s << " is inside " << p4 << " : " << p4.isInside(s) << endl;
	cout << t << " is inside " << p4 << " : " << p4.isInside(t) << endl;
	cout << u << " is inside " << p5 << " : " << p5.isInside(u) << endl;
	cout << v << " is inside " << p5 << " : " << p5.isInside(v) << endl;

	cout << endl;
	cout << p1 << " intersects " << p2 << " : " << p1.doesIntersect(p2) << endl;
	cout << p1 << " intersects " << p3 << " : " << p1.doesIntersect(p3) << endl;
	cout << p3 << " intersects " << p2 << " : " << p3.doesIntersect(p2) << endl;
	cout << p1 << " overlaps " << p8 << " : " << p1.doesIntersect(p8) << endl;
	cout << p1 << " intersects " << p9 << " : " << p1.doesIntersect(p9) << endl;

	
	Line l10({-1,2,0}, {1,2,0}), l11({-1,2,0}, {2,2,0}), l12({-1,2,0}, {5,2,0});
	cout << endl;
	cout << l10 << " intersects " << p1 << " : " << l10.doesIntersect(p1) << endl;
	cout << l11 << " intersects " << p4 << " : " << l11.doesIntersect(p4) << endl;
	cout << l12 << " intersects " << p4 << " : " << l12.doesIntersect(p4) << endl;
		
	cout << endl;
	cout << p1 << " is convex : " << p1.isConvex() << endl;
	cout << p2 << " is convex : " << p2.isConvex() << endl;
	cout << p3 << " is convex : " << p3.isConvex() << endl;
	cout << p4 << " is convex : " << p4.isConvex() << endl;
	cout << p5 << " is convex : " << p5.isConvex() << endl;	
	cout << p6 << " is convex : " << p6.isConvex() << endl;	
	cout << p7 << " is convex : " << p7.isConvex() << endl;	

	cout << endl;
	cout << p1 << " area : " << p1.area() << endl;
	cout << p2 << " area : " << p2.area() << endl;
	cout << p3 << " area : " << p3.area() << endl;
	cout << p4 << " area : " << p4.area() << endl;
	cout << p5 << " area : " << p5.area() << endl;	
	cout << p6 << " area : " << p6.area() << endl;	
	cout << p7 << " area : " << p7.area() << endl;	
				
}


// Test Driver
int main()
{ 
  TestFib();
  TestSort();
  TestNum();
  TestStr();
  TestGeom(); 
  return 0;
}