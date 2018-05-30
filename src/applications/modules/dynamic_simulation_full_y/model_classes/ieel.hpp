/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ieel.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   Oct 25, 2016
 * 
 * @brief  
 * 
 * 
 */

#ifndef _ieel_h_
#define _ieel_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_load_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class IeelLoad : public BaseLoadModel
{
  public:
    /**
     * Basic constructor
     */
    IeelLoad();

    /**
     * Basic destructor
     */
    virtual ~IeelLoad();

    /**
     * Load parameters from DataCollection object into load model
     * @param data collection of load parameters from input files
     * @param index of load on bus
     * TODO: might want to move this functionality to BaseLoadModel
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx, double dloadP, double dloadQ, int ibCMPL);

    /**
     * Initialize load model before calculation
     * @param mag voltage magnitude
     * @param ang voltage angle
     * @param ts time step 
     */
    void init(double mag, double ang, double ts);

    /**
     * Return contribution to Norton current
     * @return contribution to Norton vector
     */
    gridpack::ComplexType INorton();

    /**
     * Return Norton impedence
     * @return value of Norton impedence
     */
    gridpack::ComplexType NortonImpedence();

    /**
     * Predict part calculate current injections
     * @param flag initial step if true
     */
    void predictor_currentInjection(bool flag);

    /**
     * Corrector part calculate current injections
     * @param flag initial step if true
     */
    void corrector_currentInjection(bool flag);

    /**
     * Predict new state variables for time step
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    void predictor(double t_inc, bool flag);

    /**
     * Correct state variables for time step
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    void corrector(double t_inc, bool flag);
	
	/**
	* Set voltage on each load
	*/
	void setVoltage(gridpack::ComplexType voltage);
	
/**
 * Set terminal voltage frequency on each load
 */
	void setFreq(double dFreq);

    /**
     * Write out load state
     * @param string buffer that contains output
     * @param bufsize size of string buffer
     * @param signal character string used to determine behavior
     */
    bool serialWrite(char* string, const int bufsize, const char* signal);

  private:

    //double p_sbase;
    double p_pl, p_ql;
    //int p_status;
    
    // model definition parameters
    double a1, a2, a3, a4, a5, a6, a7, a8;
    double n1, n2, n3, n4, n5, n6; 
 
    // Internal temporary parameters
    double P0; // initial load P value
    double Q0; // initial load Q value
    double P; // actual load P at each time step
    double Q; // actual load Q at each time step

    gridpack::ComplexType p_INorton, nortonY;
    double vt_init;
    double Vsmin_pu; // constant P will be scaled down if voltage is less than Vsmin_pu
    double Vimin_pu; // constant I will be scaled down if voltage is less than Vimin_pu
 
    std::string p_loadid;
    int p_bus_id;
	
	double presentMag, presentAng, presentFreq;  //presentFreq is p.u (1.0)
	gridpack::ComplexType vt_complex;  // load terminal voltage (complex format, pu)

    friend class boost::serialization::access;

    template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar &
          & p_pl & p_ql
          & p_INorton
          & p_bus_id;
      }

};
}  // dynamic_simulation
}  // gridpack
#endif
