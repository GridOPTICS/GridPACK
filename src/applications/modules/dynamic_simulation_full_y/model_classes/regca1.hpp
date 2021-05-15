/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   regca1.hpp
 * @author Renke Huang
 * @Last modified:   May. 16, 2021
 * 
 * @brief  
 * 
 * 
 */


#ifndef _regca1_h_
#define _regca1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_generator_model.hpp"
#include "base_plant_model.hpp"
#include "DelayBlock.hpp"
#include "DelayBlockwithLimit.hpp"
#include "PIBlockwithLimit.hpp"
#include "LeadLagBlock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Regca1Generator : public BaseGeneratorModel
{
  public:
    /**
     * Basic constructor
     */
    Regca1Generator();

    /**
     * Basic destructor
     */
    virtual ~Regca1Generator();

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
	
	/**
	* return true if modify the generator parameters successfully
	* input controlTyp: 0: GFI mp adjust; 1: GFI mq adjust; 2: GFI Pset adjust; 3: GFI Qset adjust; others: invalid 
    * input newParValScaletoOrg:  GFI new parameter scale factor to the very initial parameter value at the begining of dynamic simulation
	* 
	*/
    bool applyGeneratorParAdjustment(int controlType, double newParValScaletoOrg);

    /**
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
		
	//three state blocks
	gridpack::dynamic_simulation::DelayBlock delayblock_x3vmeas;
	gridpack::dynamic_simulation::DelayBlockwithLimit delayblocklimit_x1ip;
	gridpack::dynamic_simulation::DelayBlockwithLimit delayblocklimit_x2iq;


  private:

    double p_sbase;
    double p_pg, p_qg;
    int p_status;
	bool p_tripped;
	bool bmodel_debug;
    
    double MVABase;
	int lvplsw;
    double tg, rrpwr, brkpt, zerox, lvpl1, volim, lvpnt1, lvpnt0, lolim, tfltr, khv, iqrmax, iqrmin, accel;

    //double Vterm, Theta, Ir, Ii;
    
    double x1ip_0, x2iq_0, x3vmeas_0;
    double x1ip_1, x2iq_1, x3vmeas_1;
    double dx1ip_0, dx2iq_0, dx3vmeas_0;
    double dx1ip_1, dx2iq_1, dx3vmeas_1;

    double ipcmd, iqcmd, pref, qext, ip, iq, IrNorton, IiNorton, genP, genQ, busfreq;  // busfreq is perunit, 1.0
	
    boost::shared_ptr<BaseExciterModel> p_exciter;
	boost::shared_ptr<BasePlantControllerModel> p_plant;
    
    gridpack::ComplexType p_INorton;
	gridpack::ComplexType p_Norton_Ya;

    double presentMag, presentAng;
	double Vterm, Theta;
	//double B, G;

    std::string p_ckt;
    int p_bus_id;

    friend class boost::serialization::access;

    template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar &
          & p_sbase
          & p_pg & p_qg
          & p_status
          & MVABase 
          & p_INorton
          & p_bus_id;
      }

};
}  // dynamic_simulation
}  // gridpack
#endif
