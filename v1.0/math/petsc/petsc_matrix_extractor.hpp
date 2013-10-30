// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_extractor.hpp
 * @author William A. Perkins
 * @date   2013-10-09 13:24:22 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_matrix_extractor_hpp_
#define _petsc_matrix_extractor_hpp_

#include <boost/assert.hpp>
#include "gridpack/utilities/uncopyable.hpp"
#include "implementation_visitor.hpp"

namespace gridpack {
namespace math {

class Matrix;
class PETScMatrixImplementation;

// -------------------------------------------------------------
//  class PETScMatrixExtractor
// -------------------------------------------------------------

class PETScMatrixExtractor 
  : public ImplementationVisitor
{
public:

  /// Default constructor.
  PETScMatrixExtractor(void)
    : matrix_(NULL)
  {}

  /// Destructor
  ~PETScMatrixExtractor(void)
  {}

  /// 
  void visit(PETScMatrixImplementation& petsc_impl)
  {
    matrix_ = petsc_impl.get_matrix();
  }

  Mat *matrix(void) const
  {
    return matrix_;
  }

protected:
  /// Where the matrix goes if it's found
  Mat *matrix_;

};


// -------------------------------------------------------------
//  class PETScConstMatrixExtractor
// -------------------------------------------------------------
class PETScConstMatrixExtractor 
  : public ConstImplementationVisitor
{
public:

  /// Default constructor.
  PETScConstMatrixExtractor(void)
    : matrix_(NULL)
  {}

  /// Destructor
  ~PETScConstMatrixExtractor(void)
  {}

  /// 
  void visit(const PETScMatrixImplementation& petsc_impl)
  {
    matrix_ = petsc_impl.get_matrix();
  }

  const Mat *matrix(void) const
  {
    return matrix_;
  }

protected:
  /// Where the matrix goes if it's found
  const Mat *matrix_;

};

// -------------------------------------------------------------
// PETScMatrix
// -------------------------------------------------------------
/// Get a PETSc matrix from a Matrix
inline Mat *
PETScMatrix(Matrix& A)
{
  Mat *result(NULL);
  PETScMatrixExtractor extract;
  A.accept(extract);
  result = extract.matrix();

  // a null pointer means the Matrix was not implemented in PETSc -- a
  // programming error
  BOOST_ASSERT(result);

  return result;
}

/// Get a (const) PETSc matrix from a Matrix
inline const Mat *
PETScMatrix(const Matrix& A)
{
  const Mat *result(NULL);
  PETScConstMatrixExtractor extract;
  A.accept(extract);
  result = extract.matrix();

  // a null pointer means the Matrix was not implemented in PETSc -- a
  // programming error
  BOOST_ASSERT(result);

  return result;
}
  

} // namespace math
} // namespace gridpack

#endif
