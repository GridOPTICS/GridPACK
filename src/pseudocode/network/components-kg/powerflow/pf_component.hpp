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

#ifndef _pf_component_hpp_
#define _pf_component_hpp_

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
// The BusComponent accepts Bus Count and Data visitors
// -------------------------------------------------------------

class PFComponent : public Index
{
public:
	PFComponent(int n, int m) : n_(n),m_(m){};
protected:
    void increment(MapSize * map){map->size(n_, m_);};

private:
    int n_;
    int m_;
};

} // namespace math
} // namespace gridpack

#endif
