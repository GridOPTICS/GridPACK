/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_plant_model.hpp
 * @author Bruce Palmer
 * @Last modified:   May 18, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _base_plant_model_h_
#define _base_plant_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"

namespace gridpack {
namespace dynamic_simulation {
class BasePlantControllerModel
{
  public:
    /**
     * Basic constructor
     */
    BasePlantControllerModel();

    /**
     * Basic destructor
     */
    virtual ~BasePlantControllerModel();

    /**
     * Load parameters from DataCollection object into exciter model
     * @param data collection of exciter parameters from input files
     * @param index of exciter on bus
     * TODO: might want to move this functionality to BaseExciterModel
     */
    virtual void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Initialize exciter model before calculation
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

    virtual void setGenPQV(double P, double Q, double Vt);

	virtual void setBusFreq(double freq);
	
	virtual void setPrefQext(double Pref, double Qext);
	
	virtual void setExtBusNum(int ExtBusNum);
	
	virtual void setExtGenId(std::string ExtGenId);
	
	virtual double getPref( );
	
	virtual double getQext( );

  private:
    
    //double Vterminal, w;

};
}  // dynamic_simulation
}  // gridpack
#endif
