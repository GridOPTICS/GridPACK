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

#ifndef _map_constructor_hpp_
#define _map_constructor_hpp_

#include "gridpack/network/matrix_tnterface.hpp"
#include "gridpack/math/matrix.hpp"
#include "gridpack/math/vector.hpp"
#include "gridpack/utility/noncopyable.hpp"
#include "gridpack/utility/configurable.hpp"
#include "gridpack/parallel/distributed.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class BranchComponent
// -------------------------------------------------------------

class MapConstructor
{
public:
	MapConstructor() : n_(0), m_(0){};
    virtual ~MapConstructor(void){};

    void setSize(int n, int m)
    {
        n_ = n;
        m_ = m;
        matrix_ = new math::Matrix(n_, m_);
    }
    void setData(int i, int j, complex_type * data);
    math::Matrix * getMatrix(){return matrix_;};
protected:
    int n_;
    int m_;
    math::Matrix   * matrix_;
}
} // namespace math
} // namespace gridpack

#endif
