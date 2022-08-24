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
#include "matrix.h"

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

#define PRINT_TRI(vecTri) \
	for( auto& tri : (vecTri) ) { \
		cout << "[" << std::get<0>(tri) << ", " << std::get<1>(tri) << ", " << std::get<2>(tri) << "], "; \
	} \
	cout << endl;
	
#define GEN_TRI(poly, vecTri, vecTriPnts) \
    vecTriPnts.clear(); \
	for(auto& e: vecTri) { \
		int i, j, k; \
		std::tie(i,j,k) = e; \
		vecTriPnts.push_back(std::make_tuple(poly.mVertices[i], poly.mVertices[j], poly.mVertices[k])); \
	} \
	std::swap(vecTriPnts.front(), vecTriPnts.back()); \
	std::swap(vecTriPnts.front(), vecTriPnts[vecTri.size()/2]); \
	std::swap(vecTriPnts.back(), vecTriPnts[vecTri.size()/3]);
	
#define PRINT_RECTS(rects, idx, flipped) \
	if(!(flipped)) \
		cout << "[" << (rects)[idx].w << ", " << (rects)[idx].h << "] "; \
	else \
		cout << "[" << (rects)[idx].h << ", " << (rects)[idx].w << "] "; 
	
	
double areaFromTri(const std::vector<std::tuple<Point,Point,Point>>& vecTriPnts) {
	double area = 0;
	for(auto& tri: vecTriPnts) {
		area += Polygon({std::get<0>(tri), std::get<1>(tri), std::get<2>(tri)}).Area();
	}
	
	return area;
}

void printPackedRects(const std::vector<Rect>& rects, const std::vector<std::tuple<int, Point, bool>>& packedRects, const double& W, const double& H, const double& density) {
	cout << endl << "W: " << W << " H: " << H << " Density: " << density << endl;
	int idx;
	Point loc;
	bool flipped;
	for(auto e: packedRects) {
		std::tie(idx, loc, flipped) = e;
		PRINT_RECTS(rects, idx, flipped);
		cout << "x: " << loc.x << " y: " << loc.y << endl;
	}
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
			p9({{1,-1,0}, {5,2,0}, {0,2,0}}), 
			p10({{-2,-1,0}, {1,1,0}, {2,-2,0}, {5,2,0}, {-1,2,0}, {3,4,0}, {-2,4,0}}),
			p11({{1,1,0}, {2,0,0}, {3,1,0}, {5,2,0}, {3,3,0}, {2,4,0}, {1,3,0}, {-1,2,0}}),
			p12({{1,0,0}, {2,0,0}, {2,2,0}, {3,2,0}, {1.5,4}, {0,2,0}, {1,2,0}}),
			p13({{0,0,0}, {4,0,0}, {4,2,0}, {3,1,0}, {3,3,0}, {4,2.5,0}, {4,4,0}, {0,4,0}, {2,3,0}, {2,1,0}});
	Point p({-1, 1, 0}), q({2, 2, 0}), r({3, 2, 0}), s({-1,1,0}), t({1,1,0}), u({0,1,0}), v({2,1,0});
	
	cout << endl;
	cout << "Are " << p4 << " and " << p6 << " same?: " << (p4 == p6) << endl;
	cout << "Are " << p4 << " and " << p6 << " congruent?: " << p4.isCongruent(p6) << endl;
	cout << "Are " << p4 << " and " << p4 << " congruent?: " << p4.isCongruent(p4) << endl;
	
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
	cout << p1 << " area : " << p1.Area() << endl;
	cout << p2 << " area : " << p2.Area() << endl;
	cout << p3 << " area : " << p3.Area() << endl;
	cout << p4 << " area : " << p4.Area() << endl;
	cout << p5 << " area : " << p5.Area() << endl;	
	cout << p6 << " area : " << p6.signedArea() << endl;	
	cout << p7 << " area : " << p7.signedArea() << endl;	
	
	std::vector<std::tuple<int,int,int>> vecTri;
	std::vector<std::tuple<Point,Point,Point>> vecTriPnts;
	//int i, j, k;
	cout << endl;
	p5.genTriangles(vecTri);
	cout << p5 << " triangulation is : "; PRINT_TRI(vecTri);	
	GEN_TRI(p5, vecTri, vecTriPnts);
	Polygon polyRecon(vecTriPnts);
	cout << "Reconstructed polygon from triangles is: " << polyRecon << endl;
	cout << "Does reconstructed polygon match original? : " << (p5 == polyRecon) << endl;
	cout << "Are they congruent?: " << p5.isCongruent(polyRecon) << endl;
	cout << "Polygon area = " << p5.Area() << ", Sum of area of all triangles = " << areaFromTri(vecTriPnts) << endl;
	
	p6.genTriangles(vecTri);
	cout << endl << p6 << " triangulation is : "; PRINT_TRI(vecTri);
	GEN_TRI(p6, vecTri, vecTriPnts);
	polyRecon = Polygon(vecTriPnts);
	cout << "Reconstructed polygon from triangles is: " << polyRecon << endl;
	cout << "Does reconstructed polygon match original? : " << (p6 == polyRecon) << endl;
	cout << "Are they congruent?: " << p6.isCongruent(polyRecon) << endl;
	cout << "Polygon area = " << p6.Area() << ", Sum of area of all triangles = " << areaFromTri(vecTriPnts) << endl;
	
	p10.genTriangles(vecTri);
	cout << endl << p10 << " triangulation is : "; PRINT_TRI(vecTri);
	GEN_TRI(p10, vecTri, vecTriPnts);
	polyRecon = Polygon(vecTriPnts);
	cout << "Reconstructed polygon from triangles is: " << polyRecon << endl;
	cout << "Does reconstructed polygon match original? : " << (p10 == polyRecon) << endl;
	cout << "Are they congruent?: " << p10.isCongruent(polyRecon) << endl;
	cout << "Polygon area = " << p10.Area() << ", Sum of area of all triangles = " << areaFromTri(vecTriPnts) << endl;
	
	p11.genTriangles(vecTri);
	cout << endl << p11 << " triangulation is : "; PRINT_TRI(vecTri);
	GEN_TRI(p11, vecTri, vecTriPnts);
	polyRecon = Polygon(vecTriPnts);
	cout << "Reconstructed polygon from triangles is: " << polyRecon << endl;
	cout << "Does reconstructed polygon match original? : " << (p11 == polyRecon) << endl;
	cout << "Are they congruent?: " << p11.isCongruent(polyRecon) << endl;
	cout << "Polygon area = " << p11.Area() << ", Sum of area of all triangles = " << areaFromTri(vecTriPnts) << endl;
	
	p12.genTriangles(vecTri);
	cout << endl << p12 << " triangulation is : "; PRINT_TRI(vecTri);
	GEN_TRI(p12, vecTri, vecTriPnts);
	polyRecon = Polygon(vecTriPnts);
	cout << "Reconstructed polygon from triangles is: " << polyRecon << endl;
	cout << "Does reconstructed polygon match original? : " << (p12 == polyRecon) << endl;
	cout << "Polygon area = " << p12.Area() << ", Sum of area of all triangles = " << areaFromTri(vecTriPnts) << endl;
	
	p13.genTriangles(vecTri);
	cout << endl << p13 << " triangulation is : "; PRINT_TRI(vecTri);
	GEN_TRI(p13, vecTri, vecTriPnts);
	polyRecon = Polygon(vecTriPnts);
	cout << "Reconstructed polygon from triangles is: " << polyRecon << endl;
	cout << "Does reconstructed polygon match original? : " << (p13 == polyRecon) << endl;
	cout << "Are they congruent?: " << p13.isCongruent(polyRecon) << endl;
	cout << "Polygon area = " << p13.Area() << ", Sum of area of all triangles = " << areaFromTri(vecTriPnts) << endl;
	
	Polygon p14({{0,0,0}, {2,0,0}, {2,4,0}, {0,4,0}}), p15({{0,0,0}, {4,0,0}, {4,2,0}, {0,2,0}}),
	        p14_1({{0,0,0}, {4,0,0}, {4,4,0}, {0,4,0}}),
	        p16({{0,0,0}, {4,0,0}, {0,3,0}}), p17({{0,0,0}, {0,3,0}, {-4,0,0}}),
			p18({{2,0,0}, {7,0,0}, {9,2,0}, {2,4,0}}), p19({{0,4,0}, {7,6,0}, {5,8,0}, {0,8,0}}),
			p20({{9.5,2,0}, {12,8,0}, {6,8,0}, {8,4,0}}),
			p21({{0,0,0}, {4,0,0}, {5,1,0}, {1,1,0}}), p22({{0,0,0}, {2,0,0}, {4,2,0}, {2,2,0}}),
			p23({{0,0,0}, {1,0,0}, {1,1,0}, {2,1,0}, {0,2,0}}), p24({{1,0,0}, {2,0,0},{2,2,0}, {0,1,0}, {1,1,0}}),
			p25({{0,0,0}, {2,0,0}, {2,1,0}, {1,1,0}, {1,3,0}, {0,3,0}}),
			p26({{-1,0,0}, {0,0,0}, {0,3,0}, {-2,3,0}, {-2,2,0}, {-1,2,0}}),
			p27({{-1,0,0}, {0,0,0}, {0,2,0}, {-2,2,0}, {-2,1,0}, {-1,1,0}});
			
			
	
	cout << endl << p14 << " : " << p15 << endl;
	cout << "Are the equal? : " << (p14 == p15) << endl;
	cout << "Are they congruent? : " << p14.isCongruent(p15) << endl;
	
	cout << endl << p14 << " : " << p14_1 << endl;
	cout << "Are the equal? : " << (p14 == p14_1) << endl;
	cout << "Are they congruent? : " << p14.isCongruent(p14_1) << endl;
	
	cout << endl << p21 << " : " << p22 << endl;
	cout << "Are the equal? : " << (p21 == p22) << endl;
	cout << "Are they congruent? : " << p21.isCongruent(p22) << endl;
	
	cout << endl << p16 << " : " << p17 << endl;
	cout << "Are the equal? : " << (p16 == p17) << endl;
	cout << "Are they congruent? : " << p16.isCongruent(p17) << endl;
				
	cout << endl << p18 << " : " << p19 << endl;
	cout << "Are the equal? : " << (p18 == p19) << endl;
	cout << "Are they congruent? : " << p18.isCongruent(p19) << endl;
	
	cout << endl << p18 << " : " << p20 << endl;
	cout << "Are the equal? : " << (p18 == p20) << endl;
	cout << "Are they congruent? : " << p18.isCongruent(p20) << endl;
	cout << "Are their areas same: " << p18.Area() << " == " << p20.Area() << endl;
	
	cout << endl << p23 << " : " << p24 << endl;
	cout << "Are the equal? : " << (p23 == p24) << endl;
	cout << "Are they congruent? : " << p23.isCongruent(p24) << endl;
	
	p24.genTriangles(vecTri);
	cout << endl << p24 << " triangulation is : "; PRINT_TRI(vecTri);
	GEN_TRI(p24, vecTri, vecTriPnts);
	polyRecon = Polygon(vecTriPnts);
	cout << "Reconstructed polygon from triangles is: " << polyRecon << endl;
	cout << "Does reconstructed polygon match original? : " << (p24 == polyRecon) << endl;
	cout << "Are they congruent?: " << p24.isCongruent(polyRecon) << endl;	
	
	cout << endl << p25 << " : " << p26 << endl;
	cout << "Are the equal? : " << (p25 == p26) << endl;
	cout << "Are they congruent? : " << p25.isCongruent(p26) << endl;
	
	cout << endl << p25 << " : " << p27 << endl;
	cout << "Are the equal? : " << (p25 == p27) << endl;
	cout << "Are they congruent? : " << p25.isCongruent(p27) << endl;
	
	Polygon p28({{2,1,0}, {4,1,0}, {4,2,0}, {3,3,0}, {4,5,0}, {2,3,0}, {0,5,0},  {0,3,0}, {-1,3,0}, {0,2,0}});
	Polygon p29({{4,5,0}, {-1,3,0}, {4,2,0}, {0,3,0}, {2,1,0}, {0,2,0}, {4,1,0}, {0,5,0}, {3,3,0}, {2,3,0}}, true);
	Polygon p30(p25.mVertices, true), p31(p23.mVertices,true), p32(p11.mVertices, true);
	Polygon p33({{0,-1,0}, {5,2,0}, {2,1,0}, {0,5,0}, {-2,1,0}, {-5,2,0}}),
	        p34({{-5,2,0}, {5,2,0}, {-2,1,0}, {0,-1,0}, {0,5,0}, {2,1,0}}, true);
	cout << endl << "Constructed from random order points: " << p29 << " == " << p28 << " : " << (p29 == p28) << endl;
	cout << endl << "Constructed from random order points: " << p30 << " == " << p25 << " : " << (p30 == p25) << endl;
	cout << endl << "Constructed from random order points: " << p31 << " == " << p23 << " : " << (p31 == p23) << endl;
	cout << endl << "Constructed from random order points: " << p32 << " == " << p11 << " : " << (p32 == p11) << endl;
	cout << endl << "Constructed from random order points: " << p34 << " == " << p33 << " : " << (p34 == p33) << endl;

	std::vector<Rect> rects({{7,5}, {10,15}, {3,5},{5,5}, {5,7}, {3,5}, {7,3}, {6,5}, {5,2}});
	std::vector<Rect> rects1({{5,7}, {10,15}, {3,5}, {5,5}, {5,7}, {3,5}, {3,7}, {5,6}, {5,2}});
	std::vector<Rect> rects2({{500, 200}, {250, 200}, 
	                          {50,50}, {50,50}, {50,50}, {50,50}, {50,50},
	                          {50,50}, {50,50}, {50,50}, {50,50}, {50,50},
                          	  {50,50}, {50,50}, {50,50}, {50,50}, {50,50},
	                          {50,50}, {50,50}, {50,50}, {50,50}, {50,50}});
	std::vector<Rect> rects3({ {400,50}, {400,50}, {300,50}, {300,50}, {300,50}, {300,50}, {300,50},
			{200,50}, {200,50}, {200,50}, {200,50}, {200,50}, {200,50}, {200,50}, {200,50}, {200,50}, {200,50},
			{100,50}, {100,50}, {100,50}, {100,50}, {100,50}, {100,50}, {100,50}, {100,50}, {100,50}, {100,50},
			{100,50}, {100,50}, {100,50}, {100,50}, {100,50}, {100,50}, {100,50}, {100,50}, {100,50}, {100,50},
			{50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50},
			{50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50},
			{50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50},
			{50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50} });
	std::vector<Rect> rects4({ {400,100}, {100,400}, {400,100}, {100,400}, {400,100}, {100,400} });
	std::vector<Rect> rects5({ 
			{50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50},
			{50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50}, {50,50},
			{47,31}, {47,31}, {47,31}, {47,31}, {47,31}, {47,31}, {47,31}, {47,31}, {47,31}, {47,31},
			{47,31}, {47,31}, {47,31}, {47,31}, {47,31}, {47,31}, {47,31}, {47,31}, {47,31}, {47,31},
			{23,17}, {23,17}, {23,17}, {23,17}, {23,17}, {23,17}, {23,17}, {23,17}, {23,17}, {23,17},
			{23,17}, {23,17}, {23,17}, {23,17}, {23,17}, {23,17}, {23,17}, {23,17}, {23,17}, {23,17},
			{109,42}, {109,42}, {109,42}, {109,42}, {109,42}, {109,42}, {109,42}, {109,42}, {109,42}, {109,42},
			{109,42}, {109,42}, {109,42}, {109,42}, {109,42}, {109,42}, {109,42}, {109,42}, {109,42}, {109,42},
			{42,109}, {42,109}, {42,109}, {42,109}, {42,109}, {42,109}, {42,109}, {42,109}, {42,109}, {42,109},
			{42,109}, {42,109}, {42,109}, {42,109}, {42,109}, {42,109}, {42,109}, {42,109}, {42,109}, {42,109},
			{17,33}, {17,33}, {17,33}, {17,33}, {17,33}, {17,33}, {17,33}, {17,33}, {17,33}, {17,33},
			{17,33}, {17,33}, {17,33}, {17,33}, {17,33}, {17,33}, {17,33}, {17,33}, {17,33}, {17,33},			 
			 });
	std::vector<Rect> rects6({ {400,100}, {200,150}, {50,250}, {200,50}, {50,200}, {200,50} });
	std::vector<Rect> rects7({ {100,400}, {150,200}, {50,250}, {50,200}, {50,200}, {50,200} });
			
	auto rectsInput = {rects, rects1, rects2, rects3, rects4, rects5, rects6, rects7};
			
	std::vector<std::tuple<int, Point, bool>> packedRects;
	double W = 0, H = 0, density;
	bool flipped;

	cout << endl << "Packing rectangles..." << endl;
	
	for(auto r: rectsInput) {
		W=0; H=0;
		density = packRect(packedRects, W, H, r);
		printPackedRects(r, packedRects, W, H, density);
		
		W=0; H=0;
		density = packRect(packedRects, W, H, r, true);
		printPackedRects(r, packedRects, W, H, density);

		density = packRectSqr(packedRects, W, H, r);
		printPackedRects(r, packedRects, W, H, density);
		
		W=0; H=0;
		density = packRectBest(packedRects, W, H, r, true, 4);
		printPackedRects(r, packedRects, W, H, density);	
	}
	
}

void TestMatrix() {	
	Matrix<double> mI;	
	Matrix<double> m1({{2, -1, 3}, {0, 5, 2}, {1, -1, -2}});
	//m1(0,0) = 2; m1(0,1) = -1;  m1(0,2) = 3;
	//m1(1,0) = 0; m1(1,1) = 5; m1(1,2) = 2;
	//m1(2,0) = 1; m1(2,1) = -1; m1(2,2) = -2;
	
	cout << "\nTesting matrix...\n";
	cout << m1 << endl;
	cout << "M * IdenityMatrix = \n" << m1 * mI << endl;
	
	cout << "Determinant = " << m1.determinant() << endl;
	cout << "Minor = \n" << m1.minor() << endl;
	cout << "Cofactor = \n" << m1.cofactor() << endl;
	cout << "Adjoint = \n" << m1.adjoint() << endl;
		
	//cout << m1.minor_sub(0,0) << m1.minor_sub(0,0).determinant() << endl;
	//cout << m1.minor_sub(0,1) << m1.minor_sub(0,1).determinant() << endl;
	//cout << m1.minor_sub(1,0) << m1.minor_sub(1,0).determinant() << endl;	
	//cout << m1.minor_sub(1,1) << m1.minor_sub(1,1).determinant() << endl;	
	
	cout << "Inverse = \n" << m1.inverse(true) << endl;
	
	cout << "M * Inv(M) = \n" << m1 * m1.inverse(true) << endl;
	
	
	cout << "Solving Equations: \n";
	
	Matrix<double> m2({{2, 5}, {1, 2}});
	Matrix<double> m2C({10, 2}, true);	
	cout << m2;
	cout << m2C;	
	cout << "Solution: \n";
	cout << m2.inverse() * m2C;

	Matrix<double> m3({{1, -3, 3}, {2, 3, -1}, {4, -3, -1}});
	Matrix<double> m3C({-4, 15, 19}, true);	
	cout << m3;
	cout << m3C;	
	cout << "Solution: \n";
	cout << m3.inverse() * m3C;
	
	Matrix<double> m4({{1, 1, 1, 1}, {2, 3, 0, -1}, {-3, 4, 1, 2}, {1, 2, -1, 1}});
	Matrix<double> m4C({13, -1, 10, 1}, true);
	cout << m4;
	cout << m4C;	
	cout << "Solution: \n";
	cout << m4.inverse() * m4C << endl;
	
	/*
	Problem to Solve:
	A jar contains 100 blue, green, red and yellow marbles. There are twice as many yellow marbles as blue;
	there are 10 more blue marbles than red; the sum of the red and yellow marbles is the same as the
	sum of the blue and green.
	How many marbles of each color are there: blue, green, red, yellow?
	*/
	m4 = Matrix<double>({{1, 1, 1, 1}, {2, 0, 0, -1}, {1, 0, -1, 0}, {-1, -1, 1, 1}});
	m4C = Matrix<double>({100, 0, 10, 0}, true);
	cout << m4;
	cout << m4C;	
	cout << "Solution: \n";
	cout << m4.inverse() * m4C << endl;
	
	// Using GEM rule
	cout << "\nSolution using gauss method: \n";

	std::vector<double> solVec = SolveLinearEquation(m4, m4C);
	PRINT_ARRAY(solVec);
	cout << "\n" << endl;
		
	
	Matrix<double> m7({{1, 1, 1, 1, 2, 1, 4},
	                   {2, 3, 0, -1, 0, 1, 2},
					   {-3, 4, 1, 2, 1, 0, -1},
					   {1, 2, -1, 1, 3, 3, 5},
					   {5, 1, 4, 2, -1, 1, 7},
					   {3, 5, -2, 5, 0, 3, 1},
					   {2, 1, 4, -1, 4, 1, -3}});
	Matrix<double> m7C({50, -2, 10, 65, 80, 35, 110}, true);
	cout << m7;
	cout << m7C;	
	cout << "\nSolution using inverse method: \n";
	solVec = SolveLinearEquation(m7, m7C, SolverMethod::Inverse);
	PRINT_ARRAY(solVec)
	cout << endl;
	
	cout << "\nSolution using cramer rule: \n";
	solVec = SolveLinearEquation(m7, m7C, SolverMethod::Cramer);
	PRINT_ARRAY(solVec)
	cout << endl;
	
	cout << "\nSolution using GEM method: \n";
	solVec = SolveLinearEquation(m7, m7C, SolverMethod::Gauss);
	PRINT_ARRAY(solVec)
	cout << endl;
	
	cout << "\nUpper Triangular reduction: \n";
	Matrix<double> m7U(m7);
	Matrix<double> m7CU(m7C);
	m7U.reduceTriangular();
	cout << m7U;
	cout << endl;
			
	cout << "\nLower Triangular reduction: \n";
	Matrix<double> m7L(m7);
	Matrix<double> m7CL(m7C);
	m7L.reduceTriangular(NULL, false);
	cout << m7L;
	cout << endl;
	
	cout << "\nDeterminant: " << m7.determinant() << ", " << m7.determinant(true) << endl;
	double detU = 1.0, detL = 1.0;
	for(int i=0; i < m7.nrows(); i++) {
		detU *= m7U(i,i);
		detL *= m7L(i,i);
	}
	cout << "U: " << detU << "  L: " << detL << endl;
	
	
	cout << "\nIdentity reduction: \n";
	Matrix<double> m7I(m7);
	Matrix<double> m7CI(m7.nrows());
	m7I.reduceIdentity(&m7CI);
	cout << m7I << endl;
	cout << m7CI << endl;
    cout << m7.inverse() << endl;
	cout << m7 * m7CI << endl;
	
    cout << "\nInvesre via identity reduction: \n";
	Matrix<double> m1U(m1);
	Matrix<double> m1I(m1.nrows());
	m1U.reduceIdentity(&m1I);
	cout << m1U << endl;
	cout << m1I << endl;
	cout << m1.inverse() << endl;
	cout << endl;
		
}


// Test Driver
int main()
{ 
  TestFib();
  TestSort();
  TestNum();
  TestStr();
  TestGeom();
  TestMatrix(); 
  return 0;
}