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
//  class BusComponent
// -------------------------------------------------------------

class PFBusComponent : public PFComponent
{
public:
    /// Default constructor, reader maps data to Y and CI
    PFBusComponent(const parallel::Distribution& dist);

    /// Destructor
    virtual ~PFBusComponent(void);

    virtual void addInputBranch(PFBranchComponent * branch){inputBranches_.push_back(branch);};
    virtual void addOutputBranch(PFBranchComponent * branch){outputBranches_.push_back(branch);};

    virtual void accept(ComponentVisitor & visitor){accept_(visitor);};

protected:
    virtual void accept_(ComponentVisitor & visitor){};

    virtual void accept_(ComponentCountVisitor & visitor){};
    virtual void accept_(BusCountVisitor & visitor){visitor.increment();};

    virtual void accept_(ComponentYVisitor & visitor){};
    virtual void accept_(BusYVisitor & visitor){setVisitorYValues(visitor);};

    virtual void accept_(ComponentCIVisitor & visitor){};
    virtual void accept_(BusCIVisitor & visitor){setVisitorCIValues(visitor);};

private:
    std::vector<PFBranchComponent *>   inputBranches_;
    std::vector<PFBranchComponent *>   outputBranches_;

};

} // namespace math
} // namespace gridpack

#endif
