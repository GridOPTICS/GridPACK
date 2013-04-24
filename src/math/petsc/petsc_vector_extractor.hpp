// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   petsc_vector_extractor.hpp
 * @author William A. Perkins
 * @date   Mon Apr  1 11:31:00 2013
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

#include "gridpack/utility/uncopyable.hpp"
#include "implementation_visitor.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class PETScVectorExtractor
// -------------------------------------------------------------
class PETScVectorExtractor 
  : public ImplementationVisitor,
    private utility::UnCopyable
{
public:

  /// Default constructor.
  PETScVectorExtractor(void);

  /// Destructor
  ~PETScVectorExtractor(void);

  /// 
  void visit(PETScVectorImplentation& petsc_impl)
  {
    vector_ = petsc_impl.get_vector();
  }

  Vec *vector(void) const
  {
    return vector_;
  }

protected:
  /// Where the vector goes if it's found
  Vec *vector_;

};

} // namespace math
} // namespace gridpack



#endif
