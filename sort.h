#ifndef MY_SORT_H
#define MY_SORT_H

// Few sorting and searching related experiments
// Lekhraj 


template <typename T>
int bsearch(T const arElem[], const int nSize, const T value ); // binary serach

template <typename T>
void qsort(T arElem[], const int nSize); // quick sort

template <typename T>
void isort(T arElem[], const int nSize); // insertion sort

#include "sort.tpp"

#endif