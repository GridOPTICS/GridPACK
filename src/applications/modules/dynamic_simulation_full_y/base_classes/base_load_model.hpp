/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_load_model.hpp
 * @author Bruce Palmer
 * @Last modified:   September 23, 2016
 * 
 * @brief  
 * 
 * 
 */

#ifndef _base_load_model_h_
#define _base_load_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"

namespace gridpack {
namespace dynamic_simulation {
class BaseLoadModel
{
  public:
    /**
     * Basic constructor
     */
    BaseLoadModel();

    /**
     * Basic destructor
     */
    virtual ~BaseLoadModel();

    /**
     * Load parameters from DataCollection object into load model
     * @param data collection of load parameters from input files
     * @param index of load on bus
     */
    virtual void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx, double loadP, double loadQ, int ibCMPL);

    /**
     * Initialize load model before calculation
     * @param mag voltage magnitude
     * @param ang voltage angle
     * @param ts time step
     */
    virtual void init(double mag, double ang, double ts);

    /**
     * Return contribution to Norton current
     * @return contribution to Norton vector
     */
    virtual gridpack::ComplexType INorton();

    /**
     * Return Norton impedence
     * @return value of Norton impedence
     */
    virtual gridpack::ComplexType NortonImpedence();

    /**
     * Predict new state variables for time step
     * @param flag initial step if true
     */
    virtual void predictor_currentInjection(bool flag);

    /**
     * Correct state variables for time step
     * @param flag initial step if true
     */
    virtual void corrector_currentInjection(bool flag);

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
	* post process for time step
	* @param t_inc time step increment
	* @param flag initial step if true
	*/
	virtual void dynamicload_post_process(double t_inc, bool flag);
	
	 /**
     * Set voltage on each generator
     */
    virtual void setVoltage(gridpack::ComplexType voltage);
	
	/**
     * get intialized reactive power of the dynamic load model
     */
    virtual double getInitReactivePower(void);
	
	/**
 * Set terminal voltage frequency on each load
 */
	virtual void setFreq(double dFreq);

    /**
     * Write output from loads to a string.
     * @param string (output) string with information to be printed out
     * @param bufsize size of string buffer in bytes
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if bus is contributing string to output, false otherwise
     */
    virtual bool serialWrite(char *string, const int bufsize,
        const char *signal);

    /**
     * Write out load state
     * @param signal character string used to determine behavior
     * @param string buffer that contains output
     */
    virtual void write(const char* signal, char* string);
	
	void setDynLoadP(double pl);
	void setDynLoadQ(double ql);
	void setDynLoadID(std::string load_id);
	
	double getDynLoadP();
	double getDynLoadQ();
	std::string getDynLoadID();
	

    /**
     * Set watch parameter
     * @param flag set watch to true or false
     */
    void setWatch(bool flag);

    /**
     * Get current value of watch parameter
     * @return current value of watch parameter
     */
    bool getWatch();

  private:
	
	double dyn_p;   // initial value of the dynamic load model real power P
	double dyn_q;	 // initial value of the dynamic load model real power Q
	std::string dyn_load_id; // load ID

    bool p_watch;
};
}  // dynamic_simulation
}  // gridpack
#endif
