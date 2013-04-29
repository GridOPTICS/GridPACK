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
// The BusComponent accepts Bus Count and Data visitors
// -------------------------------------------------------------

class PFComponent : public Index
{
public:
    /// Default constructor.
    PFComponent(const parallel::Distribution& dist, int yn, int ym, int cin, int cim) :
        y_(1, 1, 1, 1), ci_(1,1,1,1){};

    /// Destructor
    virtual ~PFComponent(void);
    virtual void accept(ComponentVisitor & visitor) = 0;

    void setYSourceIndices(int ul, int lr){
       y_->setSourceIndices(ul,lr);
    }

    void setCISourceIndices(int i, int j){
       ci_->setSourceIndices(i,j);
    }

    void setVisitorYValues(ComponentYVisitor & visitor);
    void setVisitorCIValues(ComponentCIVisitor & visitor);
protected:

private:
    MatrixInterface           y_;
    MatrixInterface           ci_;
};

} // namespace math
} // namespace gridpack

#endif
