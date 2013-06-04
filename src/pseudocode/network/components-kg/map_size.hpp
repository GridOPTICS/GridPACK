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

#ifndef _map_size_hpp_
#define _map_size_hpp_

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

class MapSize : public MapConstructor
{
public:
	MapSize() : n_(0), m_(0){};
    virtual ~MapSize(void){};

    virtual void setSize(int n, int m){n_ += n; m_ += m;};
    virtual void getSize(int * size){size[0] = n_, size[1] = m_;};
protected:
    int n_;
    int m_;
`} // namespace math
} // namespace gridpack

#endif
