/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   optimization_ifc.hpp
 * @author Yilin Fang and Bruce Palmer
 * @date   2015-06-12
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _optimization_ifc_h_
#define _optimization_ifc_h_

#include <vector>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "boost/smart_ptr/weak_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/expression/variable.hpp"
#include "gridpack/expression/expression.hpp"

#include <boost/serialization/export.hpp>

namespace gridpack{
namespace component{

// -------------------------------------------------------------
//  class OptimizationInterface:
//    A general interface for defining variables, constraints,
//    bounds and contributions to the objective function from
//    network components
// -------------------------------------------------------------
class OptimizationInterface {
  public:

    /**
     * Constructor
     */
    OptimizationInterface(void);

    /**
     * Destructor
     */
    virtual ~OptimizationInterface(void);

    /**
     * Return a vector of optimization variables associated witht this interface
     * @return list of variables
     */
    virtual std::vector<gridpack::optimization::Variable*> getVariables();

    /**
     * Return contribution from bus to a global constraint
     * @param tag string that can be parsed by bus to determine which constraint
     * contribution is being requested
     * @param flag bool that returns false if there is no contribution to constaint
     * @return contribution to global constraint
     */
    virtual gridpack::optimization::Expression* getGlobalConstraint(const char*
        tag, bool *flag);
  private:

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar;
  }


};
}    // component
}    // gridpack

#endif
