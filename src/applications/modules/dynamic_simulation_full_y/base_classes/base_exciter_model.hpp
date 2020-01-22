/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_exciter_model.hpp
 * @author Bruce Palmer
 * @Last modified:   May 18, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _base_exciter_model_h_
#define _base_exciter_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"

namespace gridpack {
namespace dynamic_simulation {
class BaseExciterModel
{
  public:
    /**
     * Basic constructor
     */
    BaseExciterModel();

    /**
     * Basic destructor
     */
    virtual ~BaseExciterModel();

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

    /**
     * Set the field voltage parameter inside the exciter
     * @param fldv value of the field voltage
     */
    virtual void setFieldVoltage(double fldv);

    /**
     * Set the field current parameter inside the exciter
     * @param fldc value of the field current
     */
    virtual void setFieldCurrent(double fldc);

    /** 
     * Get the value of the field voltage parameter
     * @return value of field voltage
     */
    virtual double getFieldVoltage();

    /** 
     * Get the value of the field current parameter
     * @return value of field current
     */
    virtual double getFieldCurrent();

    /**
     * Set the value of the Vterminal
     * @return value of field current
     */
    virtual void setVterminal(double mag);

    /** 
     * Set the value of the Omega 
     * @return value of field current
     */
    virtual void setOmega(double omega);

    /**
     * Set the value of the Vcomp
     * @return value of teh Vcomp
     */
    virtual void setVcomp(double vtmp);
	
	virtual void setVstab(double Vstab);
	
	virtual void setWideAreaFreqforPSS(double freq);	

  private:
    
    //double Vterminal, w;

};
}  // dynamic_simulation
}  // gridpack
#endif
