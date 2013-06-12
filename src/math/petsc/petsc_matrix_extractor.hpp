// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_extractor.hpp
 * @author William A. Perkins
 * @date   2013-06-11 14:17:34 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April 17, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

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
