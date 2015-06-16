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

    enum OptType{OPT_INT, OPT_DBL, OPT_BOOL};

    /**
     * Constructor
     */
    OptimizationInterface(void);

    /**
     * Destructor
     */
    virtual ~OptimizationInterface(void);

    /**
     * Get the number of optimization variables contributed from component
     * @return number of optimization variables
     */
    virtual int numOptVariables(void);

    /**
     * Return the data type of the optimization variable
     * @param idx index of variable inside network component
     * @return data type of network variable
     */
    virtual int optVariableType(int idx);

    /**
     * Return bounds on variable. If the lower bound is a null pointer the lower
     * bound is -inf, if the upper bound is a null pointer, the upper bound is
     * inf
     * @param idx index of variable inside network component
     * @param lo lower bound
     * @param hi upper bound
     */
    virtual void optVariableBounds(int idx, double *lo, double *hi);
    virtual void optVariableBounds(int idx, int *lo, int *hi);

    /**
     * Return the objective function contribution from the component
     * @return contribution to objective function coming from component
     */
    virtual double objectiveFunction(void);

    /**
     * Return the solution from the bus 
     */
    virtual bool solution(void);
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
