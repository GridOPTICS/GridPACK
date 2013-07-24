/**
 * @file   linear_solver_implementation.hpp
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

#ifndef _bus_component_hpp_
#define _bus_component_hpp_

#include "gridpack/network/matrix_tnterface.hpp"
#include "gridpack/math/matrix.hpp"
#include "gridpack/math/vector.hpp"
#include "gridpack/utility/noncopyable.hpp"
#include "gridpack/utility/configurable.hpp"
#include "gridpack/parallel/distributed.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class BusComponent
//  OutputInterface holds the information required by the visitor
// -------------------------------------------------------------
class BusComponent : public MatrixInterface
{
public:
  
  /// Default constructor.
  BusComponent(const parallel::Distribution& dist) : MatrixInterface(1, 1, 1, 1);
  
  /// Destructor
  virtual ~BusComponent(void);
  

protected:

  virtual void accept_(BusCountVisitor & visitor){visitor->count();};
  virtual void accept_(BusSizeVisitor & visitor){setSize(&visitor.nx);};
  virtual void accept_(BusDataVisitor & visitor){setData(visitor);};
};

} // namespace math
} // namespace gridpack



#endif
