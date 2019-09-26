// Implementation for "sort.h"
// Lekhraj

#define PRINT_AR(arElem, s_i, e_i) \
    std::cout << "[ "; \
    for(int i=s_i; i <= e_i; i++) \
		std::cout << (arElem)[i] << ' '; \
	std::cout << "] " << std::endl;
	
template <typename T>
int bsearch(T const arElem[], const int nSize, const T value) {
    // Binary Search: Search for value in a increasing order sorted array
    // Iterative version
    
    int s = 0;
    int e = nSize;
    int m;
    
    // Keep looping, narrowing the range for search until we find the given value.
    // Idea is to keep dividing serach range in two halfs, and keep choosing
    // the half, the given value is likely to lie in. Since we know sequence
    // is already sorted, we can just compare which side the value lies by comparing
    // with the middle value we have chosen.
    // Complexity O(long n)
    
    while( s < e ) {
        m = (s+e) / 2;
        if( value == arElem[m] ) // Found the value, we are done!
            return m;
        else if( value < arElem[m] )
            e = m;
        else
            s = m + 1;
        }        

    // OOPS: Did not find the given value
    return -1;
}

template <typename T>
void qsort_range(T array[], const int s, const int e) {
    // Helper Function for in-place sort
    if( e <= s ) {// just one element or empty, return as it is
        return;
    }
    else if( e-s == 1 ) { // two elements, just swap them if needed
        if( array[s] > array[e] ) {
            T temp = array[e];
            array[e] = array[s];
            array[s] = temp;
        }
    }
    else { // Quicksort: Partition sort
        // Idea is to partition the given sequence in two half around a pivot point. 
        // The first half should hav eall elements less than pivot element, and
        // the second half should have all elements higher than pivot element.
        
        // We keep partitioning and sorting until no more shuffle is needed, so we are done!
        bool shift = true;
        int i = s;
        int p = e;
        while( shift ) {
            while( (i < p) && (array[i] <= array[p]) ) // Find next element on left which is bigger then pivot 
                i += 1;
            
            if( i >= p ) { // Pivot is in right place, need to sort left and right side
                shift = false;
                qsort_range(array, s, p-1);
                qsort_range(array, p+1, e);
            }
            else { // Shift pivot to left
                // Move front element to back (to current pivot position)
                // Move last but one (pivot-1) element to front
                // Pivot moves to left: p--
                T temp = array[p];
                array[p] = array[i];
                array[i] = array[p-1];
                array[p-1] = temp;
                p -= 1;
                shift = true;
            }

        }
    }
                
}

template <typename T>
void qsort(T arElem[], const int nSize) {
    // Quick Sort
    // 
    qsort_range(arElem, 0, nSize-1);
}

template <typename T>
void isort_insert(T arElem[], const int iPos) {
    // Insert element at given position in right place in sequence [0..i-1]
    // Assumption is sequence [0..i-1] is already in right order.
    int i=iPos-1;
    T value = arElem[iPos];
    for(; i >= 0 && arElem[i] > value; i--) { // Find the right place and insert
        arElem[i+1] = arElem[i]; // Keep shifting values to right
    }
    
    // Insert at the right place
    arElem[i+1] = value;
}

template <typename T>
void isort(T arElem[], const int nSize) {
    // Insertion Sort
    // Like pack o fcars, we keep inserting elements from start to their right position.
    // At anypoint of time our sequence [0..i] is in order.
    // We are done when we have scanned the whole sequence and inserted all elements
    // to their correct position.
    // Complexity: O(n^2) but usually better if array if SMALL sequence is mostly in order, best is O(n)
    //
    for(int i=1; i < nSize; i++)
        isort_insert(arElem, i);
}



