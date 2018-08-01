/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_MATRIX_H_
#define TANGRAM_MATRIX_H_

#include <cassert>
#include <vector>
#include <array>
#include <iostream>
#include <type_traits>

#include "Vector.h"   // Tangram::Vector

/*!
  @file Matrix.h
  @brief Matrix class for Tangram
*/

namespace Tangram {


class Matrix {
 public:
  Matrix() : Rows_(0), Columns_(0) {}
  
  Matrix(int const Rows, int const Columns) :
      Rows_(Rows), Columns_(Columns) {
    A_.resize(Rows_*Columns_);    // uninitialized
  }
  
  // Initialize to some constant value
  explicit Matrix(int const Rows, int const Columns,
                  double initval) :
      Rows_(Rows), Columns_(Columns) {
    A_.resize(Rows_*Columns_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] = initval;
  }
  
  // Initialize from vector of vectors - assume that each row vector
  // is of the same size - if its not, then some values will be
  // uninitialized
  
  explicit Matrix(std::vector<std::vector<double>> const& invals) {
    Rows_ = invals.size();
    Columns_ = invals[0].size();
    A_.resize(Rows_*Columns_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] = invals[i][j];
  }

  Matrix(Matrix const& M) : Rows_(M.rows()), Columns_(M.columns()) {
    A_.resize(Rows_*Columns_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] = M[i][j];
  }

  Matrix& operator=(Matrix const& M) {
    Rows_ = M.rows();
    Columns_ = M.columns();
    A_.resize(Rows_*Columns_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] = M[i][j];
    return *this;
  }

  /// Add the Matrix @c rhs to this Matrix.
  Matrix& operator+=(const Matrix& rhs) {
    assert((Rows_ == rhs.rows()) && (Columns_ == rhs.columns()));

    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] += rhs[i][j];
  
    return *this;
  }

  /// Subtract the Matrix @c rhs from this Matrix.
  Matrix& operator-=(const Matrix& rhs) {
    assert((Rows_ == rhs.rows()) && (Columns_ == rhs.columns()));

    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] -= rhs[i][j];
  
    return *this;
  }  

  /// Multiply this Matrix by a scalar
  Matrix& operator*=(const double rhs) {
    for (int i = 0; i < A_.size(); ++i)
      A_[i] *= rhs;
  
    return *this;
  }

  // Destructor
  ~Matrix() {}

  //! Number of rows
  int rows() const {return Rows_;}

  //! Number of columns
  int columns() const {return Columns_;}

  //! Return a row of values
  //
  // The main utility of this operator is to facilitate the use
  // of the [][] notation to refer to elements of the matrix
  
  double * operator[](int const RowIndex) {
    return &(A_[RowIndex*Columns_]);
  }
  
  //! Return a row of values that cannot be modified
  //
  // The main utility of this operator is to facilitate the use
  // of the [][] notation to refer to elements of a const matrix
  
  double const * operator[](int const RowIndex) const {
    return &(A_[RowIndex*Columns_]);
  }
  
  //! Get a transpose
  Matrix transpose() const {
    Matrix AT(Columns_, Rows_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        AT[j][i] = A_[i*Columns_+j];
    return AT;
  }

  //! Get Inverse - only if its a square matrix

  Matrix inverse() const {
    if (Rows_ != Columns_) {
      std::cerr << "Matrix is rectangular" << std::endl;
      throw std::exception();
    }
    
    Matrix Ainv(Rows_, Rows_, 0.0);

    // Create a temporary matrix with twice as many columns
    Matrix D(Rows_, 2*Rows_);
    
    // Initialize the reduction matrix
    int Rows2 = 2*Rows_;
    for (int i = 0; i < Rows_; i++) {
      for (int j = 0; j < Rows_; j++) {
        D[i][j] = A_[i*Columns_+j];
        D[i][Rows_+j] = 0.0;
      }
      D[i][Rows_+i] = 1.0;
    }
    
    // Do the reduction
    for (int i = 0; i < Rows_; i++) {
      double alpha = D[i][i];
      if (alpha == 0.0) {
        std::cerr << "invert: Singular Matrix" << std::endl;
        return Ainv;
      }
      
      for (int j = 0; j < Rows2; j++)
        D[i][j] = D[i][j]/alpha;
      
      for (int k = 0; k < Rows_; k++) {
        if ((k - i) == 0)
          continue;
        
        double beta = D[k][i];
        for (int j = 0; j < Rows2; j++)
          D[k][j] = D[k][j] - beta*D[i][j];
      }
    }
    
    // Copy result into output matrix
    for (int i = 0; i < Rows_; i++)
      for (int j = 0; j < Rows_; j++)
        Ainv[i][j] = D[i][j + Rows_];

    return Ainv;
  }
  
  /*!
    @brief  Matrix-Vector multiply with std::vector
    @param[in] X The vector to post-multiply with
    
  */
  
  std::vector<double> operator*(std::vector<double> const& X) {
    assert(Columns_ == X.size());

    std::vector<double> AX(Rows_);
    for (int i = 0; i < Rows_; ++i) {
      AX[i] = 0.0;
      for (int j = 0; j < Columns_; ++j)
        AX[i] += A_[i*Columns_+j]*X[j];
    }
    return AX;
  }
  
  /*!
    @brief  Matrix-Vector multiply with Portage::Vector
    @param[in] X The vector to post-multiply with
    
  */
  
  template<long D>
  Vector<D> operator*(Vector<D> const& X) {
    assert(Rows_ == D && Columns_ == D);

    Vector<D> AX;
    for (int i = 0; i < Rows_; ++i) {
      AX[i] = 0.0;
      for (int j = 0; j < Columns_; ++j)
        AX[i] += A_[i*Columns_+j]*X[j];
    }
    return AX;
  }
  
  /*!
    @brief  Matrix-Matrix multiply
    @param[in] B   matrix to post-multiply with
  */
  
  Matrix operator*(Matrix const& B) {
    assert(Columns_ == B.rows());
    
    Matrix AB(Rows_, B.columns(), 0.0);
    
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < B.columns(); ++j)
        for (int k = 0; k < Columns_; ++k)
          AB[i][j] += A_[i*Columns_+k]*B[k][j];
    
    return AB;
  }


 private:
  int Rows_, Columns_;
  std::vector<double> A_;

  friend class MatrixRow;
};

// Add two matrices.
inline
const Matrix operator+(const Matrix& A, const Matrix& B) {
  assert((A.rows() == B.rows()) && (A.columns() == B.columns()));

  Matrix Sum(A);
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.columns(); ++j)
      Sum[i][j] += B[i][j];
  
  return Sum;
}

/// Subtract two materices.
inline
const Matrix operator-(const Matrix& A, const Matrix& B) {
  assert((A.rows() == B.rows()) && (A.columns() == B.columns()));

  Matrix Diff(A);
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.columns(); ++j)
      Diff[i][j] -= B[i][j];
  
  return Diff;
}

/*!
  @brief  Multiply a Matrix by a scalar
  @param[in] x The scaling factor 
*/
inline  
const Matrix operator*(const Matrix& A, const double& s) {
  Matrix As(A);
  
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.columns(); ++j)
      As[i][j] *= s;
  
  return As;
}

/*!
  @brief  Multiply a Matrix by a scalar
  @param[in] x The scaling factor
*/
inline  
const Matrix operator*(const double& s, const Matrix& A) {
  Matrix sA(A);
  
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.columns(); ++j)
      sA[i][j] *= s;
  
  return sA;
}

// Multiply the first vector by the transpose of the second vector
template<long D>
inline
Matrix operator*(const Vector<D>& a, const Vector<D>& b) {
  Matrix prod(D, D);
  for (int i = 0; i < D; i++) 
    for (int j = 0; j < D; j++)
      prod[i][j] = a[i]*b[j];
  return prod;
}

}  // namespace Tangram

#endif  // TANGRAM_MATRIX_H_
