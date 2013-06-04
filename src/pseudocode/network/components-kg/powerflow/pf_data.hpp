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

#ifndef _pf_data_hpp_
#define _pf_data_hpp_

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

class PFData : public MapData
{
public:
    PFData(){};
    virtual ~MapData(void){};

    // gets information from component regarding
    // I, Y and V. These values are created and sent
    // to the problem
    virtual void mapData();
    math::Matrix * getY(){return Y;};
    math::Vector * getI(){return I;};
    math::Vector * getV(){return V;};
protected:
    math::Vector *I;
    math::Vector *V;
    math::Matrix *Y;
    int n_;
    int m_;
`} // namespace math
} // namespace gridpack

#endif
