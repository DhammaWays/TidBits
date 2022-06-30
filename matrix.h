#ifndef LS_MATRIX_H
#define LS_MATRIX_H

// Voyage down the matrix lane: Motivation was to solve few equations
// Lekhraj


#include <vector>

template <typename T=double> class Matrix {

 // Public Interface
 public:
 // Constructors
  Matrix(unsigned rows, unsigned cols, const T& initial=(T)0);
  Matrix(unsigned rowscols=3); // Creates Identity Matrix of given size
  Matrix(const Matrix<T>& rhs); 
  // Creates either column or row matrix given a vector: Matrix({1, 0, 1}, true)
  Matrix(const std::vector<T>& vect, bool colmatrix); 
  // Easy way to create/populate a new matrix: Matrix({{1,0}, {0,1}})
  Matrix(const std::vector<std::vector<T>>& rhs);
  
  // Destructors
  virtual ~Matrix();

  // Matrix assignment: A = B                                                                                                                                                          
  Matrix<T>& operator=(const Matrix<T>& rhs);
  
  // Special Matrices
  Matrix<T> identity(unsigned rowscols=3) const; // I

  // Matrix mathematical operations                                                                                                                                                                                               
  Matrix<T> operator+(const Matrix<T>& rhs); // C = A + B
  Matrix<T>& operator+=(const Matrix<T>& rhs); // A += B
  Matrix<T> operator-(const Matrix<T>& rhs); // C = A - B
  Matrix<T>& operator-=(const Matrix<T>& rhs); // A -= B
  Matrix<T> operator*(const Matrix<T>& rhs); // C = A * B
  Matrix<T>& operator*=(const Matrix<T>& rhs); // A *= B
  
  Matrix<T> transpose() const;
  // To get sub matrix removing given "row" and "col"
  Matrix<T> minor_sub(const unsigned& row, const unsigned& col) const;
  Matrix<T> minor() const;
  Matrix<T> cofactor() const;
  Matrix<T> adjoint() const;
  Matrix<T> inverse() const;
  T determinant() const;

  // Matrix/scalar operations                                                                                                                                                                                                     
  Matrix<T> operator+(const T& rhs); // C = A + k
  Matrix<T> operator-(const T& rhs); // C = A - k
  Matrix<T> operator*(const T& rhs); // C = k * A
  Matrix<T> operator/(const T& rhs); // C = A / k

  // Matrix/vector operations                                                                                                                                                                                                     
  std::vector<T> operator*(const std::vector<T>& rhs); // V1 = A * Vec
  std::vector<T> diag_vec() const;
  std::vector<T> get_row(const unsigned& row_idx) const; // Returns "i" row 
  std::vector<T> get_col(const unsigned& col_idx) const; // Returns "j" column
  Matrix<T> replace_row(const unsigned& row_idx, const std::vector<T>& vec) const; // Replaces "i" row 
  Matrix<T> replace_col(const unsigned& col_idx, const std::vector<T>& vec) const; // Replaces "j" column


  // Access the individual elements                                                                                                                                                                                               
  T& operator()(const unsigned& row, const unsigned& col);
  const T& operator()(const unsigned& row, const unsigned& col) const;

  // Access the row and column sizes                                                                                                                                                                                              
  unsigned nrows() const;
  unsigned ncols() const;
  
  // Private Implementation 
  private:
  std::vector<std::vector<T> > _mat; // Matrix is vector of vectors
  unsigned _rows;
  unsigned _cols;

};

// Linear Equation Solvers

// Find solution for "X" given: A * X = C, here A is coefficient matrix and C is constant matrix
// Example:
//  Solve: 2x + 5y = 10, x + 2y = 2
//  A = [ 2, 5 ], C = [ 10 ]
//      [ 1, 2 ]      [  2 ]
//
//  Solver should return: [ -10, 6 ], i.e. x = -10, y = 6 is the solution

template <typename T=double>
std::vector<T> SolveLinearEquation(const Matrix<T>& mA, const Matrix<T>& mC, bool cramerRule = true);

#include "matrix.cpp"

#endif