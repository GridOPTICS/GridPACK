// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_extractor.hpp
 * @author William A. Perkins
 * @date   2013-05-16 08:31:57 d3g096
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

#include "gridpack/utilities/uncopyable.hpp"
#include "gridpack/math/implementation_visitor.hpp"

namespace gridpack {
namespace math {

class PETScMatrixImplementation;

// -------------------------------------------------------------
//  class PETScMatrixExtractor
// -------------------------------------------------------------

class PETScMatrixExtractor 
  : public ImplementationVisitor,
    private utility::Uncopyable
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
  : public ImplementationVisitor,
    private utility::Uncopyable
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

} // namespace math
} // namespace gridpack

#endif
