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

#ifndef _pf_bus_component_hpp_
#define _pf_bus_component_hpp_

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

class PFBranchComponent : public PFComponent
{
public:
    /// Default constructor, reader maps data to Y and CI
    PFBranchComponent(const parallel::Distribution& dist, const reader::ParameterReader);

    /// Destructor
    virtual ~PFBranchComponent(void);

    virtual void connectInputBus(PFBusComponent * bus) {input = bus;};
    virtual void connectOutputBus(PFBusComponent * bus){output = bus;};
    virtual void accept(ComponentVisitor & visitor){accept_(visitor);};

protected:
    virtual void accept_(ComponentCountVisitor & visitor){};
    virtual void accept_(BranchCountVisitor & visitor){visitor.increment();};

    virtual void accept_(ComponentYVisitor & visitor){};
    virtual void accept_(BranchYVisitor & visitor){setVisitorYValues(visitor);};

    virtual void accept_(ComponentCIVisitor & visitor){};
    virtual void accept_(BranchCIVisitor & visitor){setVisitorCIValues(visitor);};
};

    PFBusComponent    * input;
    PFBusComponent    * output;
} // namespace math
} // namespace gridpack

#endif
