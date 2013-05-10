// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   petsc_vector_extractor.hpp
 * @author William A. Perkins
 * @date   2013-05-10 09:03:48 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  1, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _petsc_vector_extractor_hpp_
#define _petsc_vector_extractor_hpp_

#include "gridpack/utilities/uncopyable.hpp"
#include "gridpack/math/implementation_visitor.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScVectorExtractor
// -------------------------------------------------------------
class PETScVectorExtractor 
  : public ImplementationVisitor
{
public:

  /// Default constructor.
  PETScVectorExtractor(void)
    : ImplementationVisitor(), p_vector(NULL)
  { }

  /// Destructor
  ~PETScVectorExtractor(void) {}

  /// Get the vector
  void visit(PETScVectorImplementation& petsc_impl)
  {
    p_vector = petsc_impl.get_vector();
  }

  Vec *vector(void)
  {
    return p_vector;
  }

protected:

  /// Where the (non-const) vector goes if it's found
  Vec *p_vector;

};


// -------------------------------------------------------------
//  class PETScConstVectorExtractor
// -------------------------------------------------------------
class PETScConstVectorExtractor 
  : public ConstImplementationVisitor
{
public:

  /// Default constructor.
  PETScConstVectorExtractor(void)
    : ConstImplementationVisitor(), p_vector(NULL)
  { }

  /// Destructor
  ~PETScConstVectorExtractor(void) {}

  /// Get the vector
  void visit(const PETScVectorImplementation& petsc_impl) 
  {
    p_vector = petsc_impl.get_vector();
  }

  const Vec *vector(void)
  {
    return p_vector;
  }

protected:

  /// Where the (const) vector goes if it's found
  const Vec *p_vector;

};

} // namespace math
} // namespace gridpack



#endif
