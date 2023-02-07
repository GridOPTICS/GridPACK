/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   epria1.hpp
 * @author Shrirang Abhyankar
 * @Added: February 7 2023
 * 
 *
 **/


#ifndef _epria1_h_
#define _epria1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_generator_model.hpp"

extern "C" {
#include "GFLIBR.h"
}

namespace gridpack {
namespace dynamic_simulation {
class Epria1Generator : public BaseGeneratorModel
{
  public:
    /**
     * Basic constructor
     */
    Epria1Generator();

    /**
     * Basic destructor
     */
    virtual ~Epria1Generator();

    /**
     * Load parameters from DataCollection object into generator model
     * @param data collection of generator parameters from input files
     * @param index of generator on bus
     * TODO: might want to move this functionality to BaseGeneratorModel
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);
   
    /**
     * Initialize generator model before calculation
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
     * Set voltage on each generator
     */
    void setVoltage(gridpack::ComplexType voltage);
	
	/**
     * Set frequency on each generator, frequency is perunit
     */
    void setFreq(double dFreq);

    /**
	* return true if trip generator successfully
	* 
	*/
    bool tripGenerator();
	
    /*
     * Write out generator state
     * @param string buffer that contains output
     * @param bufsize size of string buffer
     * @param signal character string used to determine behavior
     */
    bool serialWrite(char* string, const int bufsize, const char* signal);

    /**
     * return a vector containing any generator values that are being
     * watched
     * @param vals vector of watched values
     */
    void getWatchValues(std::vector<double> &vals);
		
  private:

    double p_sbase;
    double p_mbase;
    double p_pg, p_qg;
    int p_status;
    bool p_tripped;

    gridpack::ComplexType p_INorton;
    gridpack::ComplexType p_Norton_Ya, Y_a;
  
    double Vt, theta, VR, VI, busfreq;

    std::string p_gen_id;
    int p_bus_num;

    IEEE_Cigre_DLLInterface_Instance model;

    // Parameters
    MyModelParameters modelparams;
    // Inputs
    MyModelInputs     modelinputs;
    // Outputs
    MyModelOutputs    modeloutputs;
    // states
    double            modelstates[15];
  
    friend class boost::serialization::access;

    template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar &
          & p_sbase
	  & p_mbase
          & p_pg & p_qg
          & p_status
          & p_INorton
          & p_bus_num;
      }

};
}  // dynamic_simulation
}  // gridpack
#endif
