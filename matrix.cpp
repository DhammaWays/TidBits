// Matrix Implemention
// Lekhraj Sharma
// June 2022


// Handy print
template<typename T>
std::ostream& operator<<(std::ostream& o, const Matrix<T>& m) {
 for(unsigned i=0; i < m.nrows(); i++) {
  o << "[ ";
  for(unsigned j=0; j < m.ncols(); j++) {
  	o << m(i,j) << ", "; 
	}
  o << "]" << "\n";
  }
  o << std::endl;
  return (o);
}

// Constructors
 
 template<typename T>
 Matrix<T>::Matrix(unsigned rows, unsigned cols, const T& initial) : _rows(rows), _cols(cols) {
  _mat.resize(_rows);
  if( _cols == 0 ) _cols = _rows; // square matrix: rows x rows
  for (unsigned i=0; i < _mat.size(); i++) {
    _mat[i].resize(_cols, initial);
  }
 }
  
 //Identity Matrix constructor
 template<typename T>
 Matrix<T>::Matrix(unsigned rowscols) : _rows(rowscols), _cols(rowscols) {
  _mat.resize(_rows);
  for (unsigned i=0; i < _mat.size(); i++) {
    _mat[i].resize(_cols, (T)0);
	_mat[i][i] = (T)1;
  }
 }
 
 template<typename T>
 Matrix<T>::Matrix(const Matrix<T>& rhs) {
  _mat = rhs._mat;
  _rows = rhs.nrows();
  _cols = rhs.ncols();
 }
 
 template<typename T>
 Matrix<T>::Matrix(const std::vector<T>& vect, bool colmatrix) {
  if( colmatrix ) {
  	_rows = vect.size();
	_cols = 1;
  }
  else {
  	_rows = 1;
	_cols = vect.size();
  }
  
  _mat.resize(_rows);
  unsigned k=0;
  for (unsigned i=0; i < _rows; i++) {
    _mat[i].resize(_cols, (T)0);
	for(unsigned j = 0; j < _cols; j++) {
		_mat[i][j] = vect[k++];
	}
  }
	
 }
 
 template<typename T>
 Matrix<T>::Matrix(const std::vector<std::vector<T>>& rhs) {
  _rows = rhs.size();
  _cols = rhs[0].size();
  
  _mat.resize(_rows);
  for (unsigned i=0; i < _rows; i++) {
    _mat[i].resize(_cols, (T)0);
	for(unsigned j = 0; j < _cols; j++) {
		_mat[i][j] = rhs[i][j];
	}
  } 
 	
 }
 
 // Destructors
 template<typename T>
 Matrix<T>::~Matrix() {}

 // Assignment: A = B
 
 template<typename T>                                                                                                                                                          
 Matrix<T>& Matrix<T>::operator=(const Matrix<T>& rhs) {
  if (&rhs == this)
    return *this;

  _rows = rhs.nrows();
  _cols = rhs.ncols();

  _mat.resize(_rows);
  for (unsigned i=0; i < _mat.size(); i++) {
    _mat[i].resize(_cols);
	for(unsigned j=0; j<_cols; j++) {
      _mat[i][j] = rhs(i, j);
    }
  }

  return *this;
	
 }
  
 // Special Matrices
 
 template<typename T>
 Matrix<T> Matrix<T>::identity(unsigned rowscols) const {
 	Matrix<T> result(rowscols);
	return result;
 }

 // Matrix mathematical operations 
 
 template<typename T>                                                                                                                                                                                              
 bool Matrix<T>::operator==(const Matrix<T>& rhs) {
  if( _rows != rhs.nrows() || _cols != rhs.ncols() )
  	return false;

  for (unsigned i=0; i < _rows; i++) {
    for (unsigned j=0; j < _cols; j++) {
      if( this->_mat[i][j] != rhs(i,j) )
	  	return false;
     }
  }
	
  return true; 
 }
 
 template<typename T>                                                                                                                                                                                              
 bool Matrix<T>::operator!=(const Matrix<T>& rhs) {
 	return !(*this == rhs); 
 }
 
 
 template<typename T>                                                                                                                                                                                              
 Matrix<T> Matrix<T>::operator+(const Matrix<T>& rhs) {
  Matrix<T> result(_rows, _cols, (T)0);

  for (unsigned i=0; i < _rows; i++) {
    for (unsigned j=0; j < _cols; j++) {
      result(i,j) = this->_mat[i][j] + rhs(i,j);
     }
  }
	
  return result; 
 }
 
 template<typename T>  
 Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs) {
  unsigned _rows = rhs.nrows();
  unsigned _cols = rhs.ncols();

  for (unsigned i=0; i < _rows; i++) {
    for (unsigned j=0; j < _cols; j++) {
      this->_mat[i][j] += rhs(i,j);
    }
  }

  return *this;
 }
 
 template<typename T>                                                                                                                                                                                              
 Matrix<T> Matrix<T>::operator-(const Matrix<T>& rhs) {
  Matrix<T> result(_rows, _cols, (T)0);

  for (unsigned i=0; i < _rows; i++) {
    for (unsigned j=0; j < _cols; j++) {
      result(i,j) = this->_mat[i][j] - rhs(i,j);
    }
  }
	
  return result; 
 }
 
 template<typename T>  
 Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs) {
  unsigned _rows = rhs.nrows();
  unsigned _cols = rhs.ncols();

  for (unsigned i=0; i < _rows; i++) {
    for (unsigned j=0; j < _cols; j++) {
      this->_mat[i][j] -= rhs(i,j);
    }
  }

  return *this;
 }
 
 
 template<typename T>  
 Matrix<T> Matrix<T>::operator*(const Matrix<T>& rhs) {
  unsigned _rows = rhs.nrows();
  unsigned _cols = rhs.ncols();
  Matrix<T> result(_rows, _cols, (T)0);

  for (unsigned i=0; i < _rows; i++) {
    for (unsigned j=0; j < _cols; j++) { 
      for (unsigned k=0; k < _rows; k++) {
        result(i,j) += this->_mat[i][k] * rhs(k,j); // C = A'row x B'column
      }
    }
  }

  return result;	
 }
 
 template<typename T>  
 Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& rhs) {
  Matrix<T> result = (*this) * rhs; // just reuse our * operator
  (*this) = result;
  return *this;
 }
 

 // Matrix/scalar operations
 
 template<typename T>                                                                                                                                                                                                     
 Matrix<T> Matrix<T>::operator+(const T& rhs) {
  Matrix<T> result(_rows, _cols, (T)0);

  for (unsigned i=0; i < _rows; i++) {
    for (unsigned j=0; j < _cols; j++) {
      result(i,j) = this->_mat[i][j] + rhs;
    }
  }

  return result; 
  }
  
 template<typename T>                                                                                                                                                                                                     
 Matrix<T> Matrix<T>::operator-(const T& rhs) {
  Matrix<T> result(_rows, _cols, (T)0);

  for (unsigned i=0; i < _rows; i++) {
    for (unsigned j=0; j < _cols; j++) {
      result(i,j) = this->_mat[i][j] - rhs;
    }
  }

  return result; 
  }
  
 template<typename T>                                                                                                                                                                                                     
 Matrix<T> Matrix<T>::operator*(const T& rhs) {
  Matrix<T> result(_rows, _cols, (T)0);

  for (unsigned i=0; i < _rows; i++) {
    for (unsigned j=0; j < _cols; j++) {
      result(i,j) = this->_mat[i][j] * rhs;
    }
  }

  return result; 
  }
  
 template<typename T>                                                                                                                                                                                                     
 Matrix<T> Matrix<T>::operator/(const T& rhs) {
  Matrix<T> result(_rows, _cols, (T)0);

  for (unsigned i=0; i < _rows; i++) {
    for (unsigned j=0; j < _cols; j++) {
      result(i,j) = this->_mat[i][j] / rhs;
    }
  }

  return result; 
  }
 

 // Matrix/vector operations 
 
 // Misc vector utility 
 
 template<typename T>                                                                                                                                                                                                     
 std::vector<T> operator*(std::vector<T> vec1, std::vector<T> vec2) {
  std::vector<T> result(vec1);
  for(int i=0; i < vec1.size(); i++) {
  	result[i] *= vec2[i]; 
  }
  
  return result;
  }
  
 template<typename T>                                                                                                                                                                                                     
 std::vector<T> operator-(std::vector<T> vec1, std::vector<T> vec2) {
  std::vector<T> result(vec1);
  for(int i=0; i < vec1.size(); i++) {
  	result[i] -= vec2[i]; 
  }
  
  return result;
  }
  
 template<typename T>                                                                                                                                                                                                     
 std::vector<T> operator+(std::vector<T> vec1, std::vector<T> vec2) {
  std::vector<T> result(vec1);
  for(int i=0; i < vec1.size(); i++) {
  	result[i] += vec2[i]; 
  }
  
  return result;
  }
  
 template<typename T>                                                                                                                                                                                                     
 std::vector<T> operator/(std::vector<T> vec1, std::vector<T> vec2) {
  std::vector<T> result(vec1);
  for(int i=0; i < vec1.size(); i++) {
  	result[i] /= vec2[i]; 
  }
  
  return result;
  }
 
 template<typename T>                                                                                                                                                                                                     
 std::vector<T> operator*(std::vector<T> vec, const T& rhs) {
  std::vector<T> result(vec);
  for(int i=0; i < vec.size(); i++) {
  	result[i] *= rhs; 
  }
  
  return result;
  }
  
 template<typename T>                                                                                                                                                                                                     
 std::vector<T> operator-(std::vector<T> vec, const T& rhs) {
  std::vector<T> result(vec);
  for(int i=0; i < vec.size(); i++) {
  	result[i] -= rhs; 
  }
  
  return result;
  }
  
 
 template<typename T>                                                                                                                                                                                                     
 std::vector<T> operator+(std::vector<T> vec, const T& rhs) {
  std::vector<T> result(vec);
  for(int i=0; i < vec.size(); i++) {
  	result[i] += rhs; 
  }
  
  return result;
  }
  
  template<typename T>                                                                                                                                                                                                     
 std::vector<T> operator/(std::vector<T> vec, const T& rhs) {
  std::vector<T> result(vec);
  for(int i=0; i < vec.size(); i++) {
  	result[i] /= rhs; 
  }
  
  return result;
  }
 template<typename T>                                                                                                                                                                                                    
 std::vector<T> Matrix<T>::operator*(const std::vector<T>& rhs) {
  std::vector<T> result(_rows, (T)0);

  for (unsigned i=0; i < _rows; i++)
    for (unsigned j=0; j < _cols; j++)
      result[i] += this->_mat[i][j] * rhs[j];

  return result; 
 } 
  
 template<typename T>
 std::vector<T> Matrix<T>::diag_vec() const {
  std::vector<T> result(_rows, (T)0);

  for (unsigned i=0; i < _rows; i++) {
    result[i] = this->_mat[i][i];
  }

  return result;	
 }
 
 // Returns "i" row
 template<typename T>
 std::vector<T> Matrix<T>::get_row(const unsigned& row_idx) const {
  std::vector<T> result(_cols, (T)0);

  for (unsigned i=0; i < _cols; i++) 
    result[i] = this->_mat[row_idx][i];

  return result;		
 } 
  
 // Returns "j" column
 template<typename T>
 std::vector<T> Matrix<T>::get_col(const unsigned& col_idx) const {
  std::vector<T> result(_rows, (T)0);

  for (unsigned i=0; i < _rows; i++) 
    result[i] = this->_mat[i][col_idx];

  return result;	  	
  } 
  
 // Replaces "i" row
 template<typename T> 
 Matrix<T> Matrix<T>::replace_row(const unsigned& row_idx, const std::vector<T>& vec) const {
  Matrix<T> result(*this);
  
  for (unsigned i=0; i < _cols; i++) 
    result(row_idx,i) = vec[i];

  return result;		 	
 }  
 
 // Replaces "j" column
 template<typename T>
 Matrix<T> Matrix<T>::replace_col(const unsigned& col_idx, const std::vector<T>& vec) const {
  Matrix<T> result(*this);
  
  for (unsigned i=0; i < _rows; i++) 
    result(i, col_idx) = vec[i];

  return result;		 		
 }


  // Access the individual elements
  template<typename T>                                                                                                                                                                                               
  T& Matrix<T>::operator()(const unsigned& row, const unsigned& col) {
   return this->_mat[row][col];
  }
  
  template<typename T>
  const T& Matrix<T>::operator()(const unsigned& row, const unsigned& col) const {
  	return this->_mat[row][col]; 
  }

  // Access the row and column sizes 
  
  template<typename T>                                                                                                                                                                                             
  unsigned Matrix<T>::nrows() const {
  	return this->_rows;
  }
  
  template<typename T>
  unsigned Matrix<T>::ncols() const {
  	return this->_cols;
  }
  
 // Special matrix operations
 
 template<typename T>
 T Matrix<T>::determinant(bool bStandardWay) const {
 	//return (T)1;
	if( this->_rows == 1 && this->_cols == 1 ){
		return this->_mat[0][0];
	}
	else if( this->_rows == 2 && this->_cols == 2 )
	{
		return this->_mat[0][0] * this->_mat[1][1] - this->_mat[0][1] * this->_mat[1][0];
	}
	
	T det = (T)0;
	
	if( bStandardWay ){
		int sign = 1;
		// Summing up determinant of sub-minor matrix for first row
		for(unsigned j = 0; j < _cols; j++) {
			det += sign * this->_mat[0][j] * this->minor_sub(0, j).determinant();
			sign = -sign; // Keep flipping sign to simulate (-1) ^ (i+j) 
		}
	} else {
		Matrix<T> result(*this);
		
		// Reduce the matrix to upper triangular matrix
		result.reduceTriangular();
		
		// Determinant of a triangular matrix is product of its diagonal elements
	    det = (T)1;
		for(int i=0; i < result.nrows(); i++)
			det *= result(i,i);
	}
	
	return det;
 }
 
 template<typename T>
 Matrix<T> Matrix<T>::transpose() const {
  Matrix<T> result(_rows, _cols, (T)0);

  for (unsigned i=0; i < _rows; i++) {
    for (unsigned j=0; j < _cols; j++) {
      result(i,j) = this->_mat[j][i];
    }
  }

  return result;	
 }
 
 template<typename T>
 Matrix<T> Matrix<T>::minor_sub(const unsigned& row, const unsigned& col) const {
  Matrix<T> result(_rows-1, _cols-1, (T)0);
  
  // Need to find value for each element
  unsigned k=0,l=0;
  for(unsigned i=0; i < _rows; i++) {
  	for(unsigned j=0; j < _cols; j++) {
		if( i != row && j != col ) {
			result(k, l) = this->_mat[i][j];
			l++;
			if( l == _cols - 1 ) { // this row is done; reset for next row
				l = 0;
				k++;
			}
		}
	}
  }
  
  return result;  	
 }
 
  template<typename T>
  Matrix<T> Matrix<T>::minor() const {
  Matrix<T> result(_rows, _cols, (T)0);
  Matrix<T> result_sub(_rows-1, _cols-1, (T)0);

  // Need to find value for each element
  unsigned k=0,l=0;
  for(unsigned i=0; i < _rows; i++) {
  	for(unsigned j=0; j < _cols; j++) {
		result_sub = this->minor_sub(i, j);
		result(i, j) = result_sub.determinant();
	}
  }
  
  return result;  	
 }
 
  
 template<typename T>
 Matrix<T> Matrix<T>::cofactor() const {
  Matrix<T> result(_rows, _cols, (T)0);
  
  result = this->minor();
  // Multiply by appropriate signs
  int sign = 1;
  for (unsigned i=0; i < _rows; i++) {
    for (unsigned j=0; j < _cols; j++) {
	  sign = ( ((i+j)&0x1) == 0 ? 1 : -1 ); // Simulate (-1)^(i+j)
      result(i,j) = sign * result(i, j);
    }
  }
  
  return result;  		
 }
 
 template<typename T>
 Matrix<T> Matrix<T>::adjoint() const {
  Matrix<T> result(_rows, _cols, (T)0);
  
  result = this->cofactor().transpose();
  
  return result;  		
 }
 
 template<typename T>
 Matrix<T> Matrix<T>::inverse(bool bStandardWay) const {
  
  if( bStandardWay ) {
	  Matrix<T> result(_rows, _cols, (T)0);
	  
	  T det_val = this->determinant();
	  assert(det_val != (T)0);
	  result = this->adjoint() * (1/det_val);
	  
	  return result;
  }
  else {
    Matrix<T> result(*this);
  	Matrix<T> mI(_rows);  
  	result.reduceIdentity(&mI);
	
	return mI;
  }
  
  // Should not get here!
  return Matrix<T>(_rows);  			
 }
 
// Reduce matrix to upper/lower triangular matrix (lower/upper half becomes zeros)
// Reduces it in-place, so orginal matrix along with given adjacent matrix is directly reduced
template <typename T>
Matrix<T>& Matrix<T>::reduceTriangular(Matrix<T>* pmAdj, bool bUpperTriangular) {
	// Reduction: Reduce coeff matrix into an upper/lower diagonal matrix 
	//            Also perform same operations on given (if any) adjacent matrix
	// Complexity: O(n^3)
	//

	// Only valid for square matrix	
	assert(this->nrows() == this->ncols() );
	if( pmAdj != NULL ) assert(pmAdj->nrows() == this->nrows());
	
	const unsigned N = this->ncols();
	T redFactor; // reduction factor to for each row/element
		
	// Reduction
	if( bUpperTriangular ) { // Recuce to upper triangular matrix, i.e. lower half is zeros
		for(int k = 0; k < N - 1; k++) {
			for(int i = k+1; i < N; i++) {
				assert(this->_mat[k][k] != (T)0 );
				redFactor = this->_mat[i][k] / this->_mat[k][k];
				// Row(i) = Row(i) - reductionFactor * Row(k)
				// A = A.replace_row(i, A.get_row(i) - A.get_row(k) * redFactor);
				
				// Loop is from k+1..N as A[0..k-1] should already be zero by now!
				this->_mat[i][k] = 0; // Our reduction factor gaurantees that "k" elemnet of "i" row would be zero after subtracting!
				for(int j = k+1; j < N; j++) {
					this->_mat[i][j] -= redFactor * this->_mat[k][j];
				}
				//mAdj(i, 0) -= redFactor * mAdj(k, 0); 
				if( pmAdj != NULL )
					pmAdj->_mat[i] = pmAdj->_mat[i] - pmAdj->_mat[k] * redFactor;
			}
		}
	}
	else { // Reduce to Lower Triangular Matrix; i.e. upper half is zeros
		for(int k = N-1; k > 0; k--) {
			for(int i = k-1; i >= 0; i--) {
				assert(this->_mat[k][k] != (T)0 );
				redFactor = this->_mat[i][k] / this->_mat[k][k];
					
				// Loop is from k..1 as A[k+1...N] should already be zero by now!
				this->_mat[i][k] = 0; // Our reduction factor gaurantees that "k" elemnet of "i" row would be zero after subtracting!
				for(int j = k-1; j >= 0; j--) {
					this->_mat[i][j] -= redFactor * this->_mat[k][j];
				}
				//mAdj(i, 0) -= redFactor * mAdj(k, 0); 
				if( pmAdj != NULL )
					pmAdj->_mat[i] = pmAdj->_mat[i] - pmAdj->_mat[k] * redFactor;
			}
		}
		
	}
		
	return *this;	
}

// Reduce matrix to identity matrix
// Reduces it in-place, so orginal matrix along with given adjacent matrix is directly reduced
template <typename T>
Matrix<T>& Matrix<T>::reduceIdentity(Matrix<T>* pmAdj) {
	// Reduction: Reduce coeff matrix into an identity matrix 
	//            Also perform same operations on given (if any) adjacent matrix
	// Complexity: O(n^3)
	//	
	
	//Reduce first to upper and then lower triangular to get to a diagonal matrix
	this->reduceTriangular(pmAdj);
	this->reduceTriangular(pmAdj, false);

	// Now we convert diagonal matrix into identoy by just divding each row by diagonal value
	for(int i=0; i < this->nrows(); i++) {
		if( pmAdj != NULL )
			pmAdj->_mat[i] = pmAdj->_mat[i] / this->_mat[i][i];
		
		this->_mat[i][i] = (T)1; // By design it will become 1 as we are dividing by same value!
	}
	
	// Our originalmatrix should become an identity matrix by now!
	assert( *this == this->identity(this->nrows()) );
	return *this;	
}

 
 
// Solve linear equations using GEM (Gauss Elimination Method) 
template <typename T>
std::vector<T> SolveLinEqGEM(const Matrix<T>& mA, const Matrix<T>& mC) {
	// GEM Algorithm:
	// 1. Reduction: Reduce coeff matrix into an upper diagonal matrix ( also perform same operations on constant matrix)
	// 2. Back subsitution: Working backwards from (last row/variable) calculate value of
	//    each variable, i.e. Xn, Xn-1, Xn-2, ..., X1
	//
	// Assumption: No diagonal element of coefficient matrix is zero!
	// Complexity: O(n^3), more like 2/3 * n^3 - 3/2 * n^2 - 7/6 * n
	//
	
	const unsigned N = mA.ncols();
	std::vector<T> solVec;
	solVec.resize(N);
	Matrix<T> Ag(mA); // make copy of coef matrix; we will reduce it
	Matrix<T> Cg(mC); // make copy of constant matrix as well
	T redFactor; // reduction factor to for each row/element
	
	// std::cout << "GEM: N: " << N << " Given coeff matrix: " << std::endl << Ag;
	
	// Reduction
	/*
	for(unsigned k = 0; k < N - 1; k++) {
		for(unsigned i = k+1; i < N; i++) {
			assert(Ag(k, k) != (T)0 );
			redFactor = Ag(i, k) / Ag(k, k);
			// Row(i) = Row(i) - reductionFactor * Row(k)
			// Ag = Ag.replace_row(i, Ag.get_row(i) - Ag.get_row(k) * redFactor);
			
			// Loop is from k+1..N as Ag[0..k-1] should already be zero by now!
			Ag(i, k) = 0; // Our reduction factor gaurantees that "k" elemnet of "i" row would be zero after subtracting!
			for(unsigned j = k+1; j < N; j++) {
				Ag(i, j) -= redFactor * Ag(k, j);
			}
			Cg(i, 0) -= redFactor * Cg(k, 0); 
		}
	}
	*/
	Ag.reduceTriangular(&Cg);
	// Assert: Reduced Coefficient matrix should be upper diagonal matrix now!
	//         Therefore ALL lower diagonal elements should be zero!
	// std::cout << "Reduced coeff matrix: " << std::endl << Ag;

	
	// Back substitution
	// Calculate solution working backwards
	T sum;
	solVec[N-1] = Cg(N-1,0) / Ag(N-1, N-1); // Last variable
	// Work backwards for remaining variables: Xn-2, Xn-3, ..., X0
	for(int i = N - 2; i >= 0; i--) { // cannot have "unsigned" since "i" would below zero in check!
	    // Substitute all caluclated variables in current row: Xn-1, Xn-2, ..., Xi+1
		// For given row, find the forward sum of remaining part of the equation
		// for already calculated variables.
		sum = 0;
		for(int j = i + 1; j <= N-1; j++) {
			sum += Ag(i, j) * solVec[j];
		}
		// Calculate value of next variable
		assert(Ag(i, i) != (T)0);
		solVec[i] = (Cg(i, 0) - sum) / Ag(i, i);
	}
	
	// std::cout << "GEM Solution: " << Matrix<double>(solVec, false) << std::endl;	
	return solVec;	
}

// Solve linear equations using Cramer Rule 
template <typename T>
std::vector<T> SolveLinEqCramer(const Matrix<T>& mA, const Matrix<T>& mC) {
	// Cramer rule:
	// 1. Ai = For "i" variable solution, substitue "i" column in coeff matrix with constants
	// 2. Xi = det|Ai| / det|A|
	//
	// Complexity: O(n!), closer to 3 * factorial(n+1), so it is super expensive!
	//
	// Why does cramer rule work?
	// Let us see in action in 2 variable case: a1 * x + b1 * y = c1, a2 * x + b2 * y = c2
	//  [ a1 b1 ]   [ x ]   [ c1 ]
	//  [ a2 b2 ] * [ y ] = [ c2 ]
	// can be rewritten as:
	//  [ a1 b1 ]   [ x  0]   [ c1 b1]
	//  [ a2 b2 ] * [ y  1] = [ c2 b2]
	//
	//  So here, leftmost matrix is our original coeff matrix "A", and rightmost matrix "Ai" is our
	//  substituted coeff matrix where first column has been replaced by constant values.
	//
	//  if we now take determinant on both sides:
	//    det|A| * x = det|Ai|
	//    x = det|Ai| / det|A|
	//
	//  Similiarly, we can rewrite original equation as:
	//  [ a1 b1 ]   [ 1  x]   [ a1 c1]
	//  [ a2 b2 ] * [ 0  y] = [ a2 c2]
	//  where we can see now our substituted coeff matrix "Ai" now has its second column replaced
	//  with constants.
	//   det|A| * y = det|Ai|
	//   y = det|Ai| / det|A|
	//
	
	std::vector<T> solVec;
	solVec.resize(mA.ncols());
	Matrix<T> Ai(mA);
	std::vector<T> constVec = mC.get_col(0);
	T detA = mA.determinant();
	for(unsigned i = 0; i < mA.ncols(); i++) {
		Ai = mA.replace_col(i, constVec);
		solVec[i] = Ai.determinant()/ detA;
	}
	
	return solVec;		
}

 
 // Linear Equation Solvers
template <typename T>
std::vector<T> SolveLinearEquation(const Matrix<T>& mA, const Matrix<T>& mC, const SolverMethod solveMethod) {
	if( solveMethod == SolverMethod::Gauss ) {
		return SolveLinEqGEM(mA, mC);
	}
	else if( solveMethod == SolverMethod::Cramer ) {
		return SolveLinEqCramer(mA, mC);
	}
	else { // Use standard inverse method
		//  A * X = C => X = Inv(A) * C		
		return (mA.inverse() * mC).get_col(0);		
	}
}

 
 
