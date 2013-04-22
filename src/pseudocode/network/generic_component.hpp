/**
 * @file   generic_component.hpp
 * @author Kevin A. Glass
 * @date   Mon Apr  19 13:51 2013
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  1, 2013 by Kevin A. Glass
// Last Change:
// -------------------------------------------------------------

#ifndef _generic_component_hpp_
#define _generic_component_hpp_

#include "gridpack/utility/noncopyable.hpp"
#include "gridpack/utility/configurable.hpp"
#include "gridpack/parallel/distributed.hpp"
#include "gridpack/network/component_visitor.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class GenericComponent
//  OutputInterface holds the information required by the visitor
// -------------------------------------------------------------
class GenericComponent
  : public parallel::Distributed,
    public utility::Configurable,
    private utility::UnCopyable,
    public network::OutputInterface
{
public:
  
  /// Default constructor.
  GenericComponent(const parallel::Distribution& dist,
                             const Matrix& A);
  
  /// Destructor
  virtual ~GenericComponent(void);
  
  /// Solve w/ the specified RHS, put result in specified vector
  void solve(const Vector& b, Vector& x) const
  {
    this->solve_(b, x);
  }

  /// Allow visits by implemetation visitor
  void size(ComponentVisitor& visitor)
  {
    this->size_(visitor);
  }
  /// Allow visits by implemetation visitor
  void mapData(ComponentVisitor& visitor)
  {
    this->mapData_(visitor);
  }


protected:

};

} // namespace math
} // namespace gridpack



#endif
