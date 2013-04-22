// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   petsc_matrix_extractor.hpp
 * @author William A. Perkins
 * @date   Wed Apr 17 13:57:20 2013
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

#include "gridpack/utility/uncopyable.hpp"
#include "implementation_visitor.hpp"
#include "petsc_matrix_implementation.hpp"
#include "petscmat.h"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScMatrixExtractor
// -------------------------------------------------------------

class PETScMatrixExtractor 
  : public ImplementationVisitor,
    private utility::UnCopyable
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
  void visit(PETScMatrixImplentation& petsc_impl)
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

#endif
