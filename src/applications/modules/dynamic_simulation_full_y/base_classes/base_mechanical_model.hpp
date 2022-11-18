/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_mechanical_model.hpp
 * @author Shuangshuang Jin
 * @Created on:    Nov 15, 2022
 * 
 * @brief  
 * 
 * 
 */

#ifndef _base_mechanical_model_h_
#define _base_mechanical_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"

namespace gridpack {
namespace dynamic_simulation {
class BaseMechanicalModel
{
  public:
    /**
     * Basic constructor
     */
    BaseMechanicalModel();

    /**
     * Basic destructor
     */
    virtual ~BaseMechanicalModel();

    /**
     * Load parameters from DataCollection object into mechanical model
     * @param data collection of mechanical parameters from input files
     * @param index of mechanical on bus
     * TODO: might want to move this functionality to BaseMechanicalModel
     */
    virtual void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Initialize mechanical model before calculation
     * @param mag voltage magnitude
     * @param ang voltage angle
     * @param ts time step
     */
    virtual void init(double mag, double ang, double ts);

    /**
     * Predict new state variables for time step
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    virtual void predictor(double t_inc, bool flag);

    /**
     * Correct state variables for time step
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    virtual void corrector(double t_inc, bool flag);



  private:

};
}  // dynamic_simulation
}  // gridpack
#endif
