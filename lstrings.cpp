
// Experiments with strings
// Lekhraj

#include "lstrings.h"

template <typename T>
void swap(T a, T b) {
    T tmp = a;
    a = b;
    b = tmp;
}

// Reverse string
std::string rev(const std::string& s) {
    // Reverse in-place a copy of input string
    char tCh;
    std::string s1 = s;
    for(int i=0, e=s1.length()-1; i < e; i++, e--) { // keep swapping until pointers cross
        // Swap
        tCh = s1[i];
        s1[i] = s1[e];
        s1[e] = tCh;
    }
    
    return s1;    
}

// Helper methods: Find longest palindrome match in a string with given centers
std::pair<int, int> matchCenter(const std::string& s, const int n, const int i1, const int i2) {
    bool matchWhole = true;
    int i=i1, j=i2;
    for(; i > -1 && j < n;  i--, j++) {
        if( s[i] != s[j] ) {
            matchWhole = false;
            break;
        }
    }
    

    if( matchWhole ) {
        if( i < 0 ) { i = 0; j =j-1; }
        if( j >= n ) { i = i+1; j = j-1; }
        return std::make_pair(i, j);
    }
    else
      return std::make_pair(i+1, j-1);
  }


// Return length of string given its start and end index tuple
int lenStr( std::pair<int, int> curIdx ) {
    if(curIdx.second-curIdx.first+1 > 0)
        return (curIdx.second-curIdx.first+1);
    else
        return 0;
}

  
// Longest palindrome
std::string palindrome(const std::string& a) {
    /*
    Given a string "a", returns the longest palindromic substring contained in "a"
    ASSUMPTIONS: Single character string is also a palindrome
    ARGUMENTS: String "a" for which longest palindromic substring is to be extracted
    MODIFIES: None
    RETURN: Longest palindromic substring in string "a". In case of multiple longest
            palindromes, first one from start of string is returned.
    */
    
    // Trivial case
    int size = a.length();
    if( size == 0 )
        return "";
    
    // Iterate over whole string, keeping track of current maximum palindrome 
    std::pair<int, int> maxStr = std::make_pair(0,0);
    std::pair<int, int> curStr = std::make_pair(0,0);
    
    for(int i=0; i <= size-1; i++) {
        // Find all odd matches which have same center in a palindrome
        curStr = matchCenter(a, size, i, i);
        if( lenStr(curStr) > lenStr(maxStr) )
            maxStr = curStr;
      
        // Find all even matches which have two centers in a palindrome
        curStr = matchCenter(a, size, i, i+1);
        if( lenStr(curStr) > lenStr(maxStr) )
            maxStr = curStr;
        }

   
   // return a[maxStr[0]:maxStr[1]+1]
   return a.substr(maxStr.first, (maxStr.second-maxStr.first)+1);    
}
