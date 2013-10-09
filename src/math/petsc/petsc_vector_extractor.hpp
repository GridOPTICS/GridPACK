// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_vector_extractor.hpp
 * @author William A. Perkins
 * @date   2013-10-09 13:25:32 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_vector_extractor_hpp_
#define _petsc_vector_extractor_hpp_

#include "gridpack/utilities/uncopyable.hpp"
#include "vector.hpp"
#include "petsc/petsc_vector_implementation.hpp"
#include "implementation_visitor.hpp"

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
    p_vector = petsc_impl.getVector();
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
    p_vector = petsc_impl.getVector();
  }

  const Vec *vector(void)
  {
    return p_vector;
  }

protected:

  /// Where the (const) vector goes if it's found
  const Vec *p_vector;

};

// -------------------------------------------------------------
// PETScVector
// -------------------------------------------------------------
/// Get a PETSc vector from a Vector
inline Vec *
PETScVector(Vector& A)
{
  Vec *result(NULL);
  PETScVectorExtractor extract;
  A.accept(extract);
  result = extract.vector();

  // a null pointer means the Vector was not implemented in PETSc -- a
  // programming error
  BOOST_ASSERT(result);

  return result;
}

/// Get a (const) PETSc vector from a Vector
inline const Vec *
PETScVector(const Vector& A)
{
  const Vec *result(NULL);
  PETScConstVectorExtractor extract;
  A.accept(extract);
  result = extract.vector();

  // a null pointer means the Vector was not implemented in PETSc -- a
  // programming error
  BOOST_ASSERT(result);

  return result;
}


} // namespace math
} // namespace gridpack



#endif
