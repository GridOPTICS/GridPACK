/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   acmotor.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   Oct 21, 2016
 * 
 * @brief  
 * 
 * 
 */

#ifndef _acmotor_h_
#define _acmotor_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_load_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class AcmotorLoad : public BaseLoadModel
{
  public:
    /**
     * Basic constructor
     */
    AcmotorLoad();

    /**
     * Basic destructor
     */
    virtual ~AcmotorLoad();

    /**
     * Load parameters from DataCollection object into load model
     * @param data collection of load parameters from input files
     * @param index of load on bus
     * TODO: might want to move this functionality to BaseLoadModel
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx, double loadP, double loadQ, int ibCMPL);

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
	* post process for time step
	* @param t_inc time step increment
	* @param flag initial step if true
	*/
	void dynamicload_post_process(double t_inc, bool flag);
	
	/**
	* Set voltage on each load
	*/
	void setVoltage(gridpack::ComplexType voltage);
	
	/**
 * Set terminal voltage frequency on each load
 */
	void setFreq(double dFreq);
	
	/**
     * get intialized reactive power of the dynamic load model
     */
    double getInitReactivePower(void);

    /**
     * Write out load state
     * @param string buffer that contains output
     * @param bufsize size of string buffer
     * @param signal character string used to determine behavior
     */
    bool serialWrite(char* string, const int bufsize, const char* signal);

  private:

    double p_sbase;
    double p_pl, p_ql;
    int p_status;
    
    // AC Motor parameters used in definition
    double Tstall, Trst, Tv, Tf, CompLF, CompPF, Vstall, Rstall, Xstall, LFadj;
    double Kp1, Np1, Kq1, Nq1, Kp2, Np2, Kq2, Nq2;
    double Vbrk, Frst, Vrst, CmpKpf, CmpKqf, Vc1off, Vc2off;
    double Vc1on, Vc2on, Tth, Th1t, Th2t, Fuvr, Uvtr1, Ttr1, Uvtr2, Ttr2;
 
    // Network parameters used in the model
    double volt;
    double freq;

    // State variables and derivative of state variables
    double volt_measured0, freq_measured0, temperatureA0, temperatureB0;
    double volt_measured, freq_measured, temperatureA, temperatureB;
    double dv_dt0, dfreq_dt0, dThA_dt0, dThB_dt0;
    double dv_dt, dfreq_dt, dThA_dt, dThB_dt;

    // Internal temporary parameters
    double MVABase;
    double Pinit, Qinit, Pinit_pu, Qinit_pu, P0, Q0, PA, QA, PB, QB, Pmotor, Qmotor;
    double Vstallbrk;
    double Gstall, Bstall;
    gridpack::ComplexType equivY, equivY_sysMVA;
    gridpack::ComplexType INorton_sysMVA;
    double Imotor_init;
    int statusA, statusB;
    double stallTimer;
    double restartTimer;
    double FthA, FthB;
    double TthA, TthB;
    double thEqnA, thEqnB;
    double Kuv, Kcon, fcon_trip;
    int isContractorActioned;
    double UVTimer1, UVTimer2;
    gridpack::ComplexType equivYpq_motorBase;
    double I_conv_factor_M2S;
   
 
    gridpack::ComplexType p_INorton;

    double presentMag, presentAng, presentFreq;  //presentFreq is p.u (1.0)
	gridpack::ComplexType vt_complex;  // load terminal voltage (complex format, pu)

    //double Efdinit, Pmechinit;
            
    std::string p_loadid;
    int p_bus_id;

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
