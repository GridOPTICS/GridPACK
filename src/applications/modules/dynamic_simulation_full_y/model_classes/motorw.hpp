/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   motorw.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   Dec 12, 2016
 * 
 * @brief  
 * 
 * 
 */

#ifndef _motorw_h_
#define _motorw_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_load_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class MotorwLoad : public BaseLoadModel
{
  public:
    /**
     * Basic constructor
     */
    MotorwLoad();

    /**
     * Basic destructor
     */
    virtual ~MotorwLoad();

    /**
     * Load parameters from DataCollection object into load model
     * @param data collection of load parameters from input files
     * @param index of load on bus
     * TODO: might want to move this functionality to BaseLoadModel
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx, double dloadP, double dloadQ, int ibCMPL);

    /**
     * Update parameters in DataCollection object with current values from
     * load
     * @param data collection object for bus that hosts load
     * @param index of generator on bus
     */
    void updateData(boost::shared_ptr<gridpack::component::DataCollection> data,
        int idx);

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

    double getMotorQ();

    /**
     * get intialized reactive power of the dynamic load model
     */
    double getInitReactivePower(void);
	
	/**
     * get the variable Fonline for each load
     */
    double getFonline(void);
	
	/**
	 * return true if load change is enabled
	 * @param percentageFactor: the fraction (percentage) of load that is changed. Negative: load reduction, Positive: load increase
	 */
	bool changeLoad(double percentageFactor);
	
	/**
     * Set same bus static load p and q for load shedding action usage
     */
    void setSameBusStaticLoadPQ(double static_pl, double static_ql, double mag);

  private:

    double p_sbase;
    double p_pl, p_ql;
    int p_status;
	bool bdebugprint;
	
    // parameters
    double loadFactor, Pul, rs, lls, lm, lmr, rr1, llr1, rr2, llr2;
    double H, MVABase, A, B, C0, D, E;
    double Ls, lmp, Lp, lmpp, Lpp, tpo, tppo;

    // boundary variables
    double volt, freq, Id, Iq;
    
    // state varialbes for predictor
    double epq0, epd0, eppq0, eppd0, slip0;

    // state variables for corrector
    double epq, epd, eppq, eppd, slip;

    // derivatives of state vriables for predictor
    double depq_dt0, depd_dt0, deppq_dt0, deppd_dt0, dslip_dt0;

    // derivatives of state variables for corrector
    double depq_dt, depd_dt, deppq_dt, deppd_dt, dslip_dt;

    // oter varialbes
    double w0, TL, Tm0, p, q, Pmotor, Qmotor, Qmotor_init, sysMVABase;
	gridpack::ComplexType nortonImpedance_sysMVA;
	double Fonline;
	double motorwload_perc; //motorw load percentage to the total load;
	double samebus_static_loadP;  // system base mva pu value
    double samebus_static_loadQ;  // system base mva pu value
	gridpack::ComplexType samebus_static_equivY_sysMVA;

    gridpack::ComplexType p_INorton;

    std::string p_loadid;
    int p_bus_id;
	
	double presentMag, presentAng, presentFreq;  //presentFreq is p.u (1.0)
	gridpack::ComplexType vt_complex;  // load terminal voltage (complex format, pu)

    friend class boost::serialization::access;

    template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar &
          & p_sbase
          & p_pl & p_ql
          & p_status
          & MVABase 
          & p_INorton
          & p_bus_id;
      }

};
}  // dynamic_simulation
}  // gridpack
#endif
