/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   optimization_ifc.hpp
 * @author Yilin Fang and Bruce Palmer
 * @date   2015-10-12 13:13:54 d3g096
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
     * Return a vector of optimization variables associated with this interface
     * @return list of variables
     */
    virtual std::vector<boost::shared_ptr<gridpack::optimization::Variable> >
      getVariables();

    /**
     * Return a vector of auxiliary variables associated with this interface.
     * These are variables that are used in expressions but may not be
     * defined by this network
     * @return list of variables
     */
    virtual std::vector<boost::shared_ptr<gridpack::optimization::Variable> >
      getAuxVariables();

    /**
     * Return contribution from bus to a global constraint
     * @param tag string that can be parsed by bus to determine which constraint
     * contribution is being requested
     * @return contribution to global constraint. If no contribution, return
     * null pointer
     */
    virtual boost::shared_ptr<gridpack::optimization::Expression> 
      getGlobalConstraint(const char* tag);

    /**
     * Return a list of local constraints from component
     * @return list of constraints
     */
    virtual std::vector<boost::shared_ptr<gridpack::optimization::Constraint> >
      getLocalConstraints();

    /**
     * Return contribution to objective function
     * @return expression representing contribution to objective function. If no
     * contribution, return null pointer
     */
    virtual boost::shared_ptr<gridpack::optimization::Expression>
      getObjectiveFunction();
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
