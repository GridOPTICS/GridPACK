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

#ifndef _map_data_hpp_
#define _map_data_hpp_

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

class MapData : public MapConstructor
{
public:
    MapData(math::Matrix * matrix) : matrix_(matrix), n_(0), m_(0){};
    virtual ~MapData(void){};

    virtual void mapData(int region_n, int region_n, complex_double * region);
protected:
    math::Matrix * matrix_;
    int n_;
    int m_;
`} // namespace math
} // namespace gridpack

#endif
