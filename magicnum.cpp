int MagicNum(const int n) {	
	
	int N=n, magN = 0, d=N, i=0;
	while( N > 0 ) {	
		N /= 10;
		d = d - N * 10; 
		magN = magN * 10 + d;
		d = N;
		i++;
	}
	
	N = n;
	while( i > 0 ) {
		 N *= 10;
		 i--;
	 }
	magN += N;
	
	return magN;
	}

#include <iostream>	
using namespace std;

main()
{
	#154451
	cout << MagicNum(154) << endl;
	cout << MagicNum(19) << endl;
	cout << MagicNum(9792) << endl;
	cout << MagicNum(5) << endl;
	cout << MagicNum(10) << endl;
}
