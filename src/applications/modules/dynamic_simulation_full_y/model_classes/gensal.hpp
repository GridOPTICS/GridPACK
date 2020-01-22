/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   gensal.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   June 1, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _gensal_h_
#define _gensal_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_generator_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class GensalGenerator : public BaseGeneratorModel
{
  public:
    /**
     * Basic constructor
     */
    GensalGenerator();

    /**
     * Basic destructor
     */
    virtual ~GensalGenerator();

    /**
     * Load parameters from DataCollection object into generator model
     * @param data collection of generator parameters from input files
     * @param index of generator on bus
     * TODO: might want to move this functionality to BaseGeneratorModel
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Saturation function
     * @ param x
     */
    double Sat(double x);
    
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
     * Get the value of the field voltage parameter
     * @return value of field voltage
     */
    double getFieldVoltage();
	
	void setWideAreaFreqforPSS(double freq);

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

  private:

    double p_sbase;
    double p_pg, p_qg;
    int p_status;
    
    double MVABase;
    double H, D, Ra, Xd, Xq, Xdp, Xdpp, Xl;
    double Tdop, Tdopp, Tqopp, S10, S12;
    
    double Vterm, Theta, Ir, Ii;
    
    double x1d_0, x2w_0, x3Eqp_0, x4Psidp_0, x5Psiqpp_0;
    double x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqpp_1;
    double dx1d_0, dx2w_0, dx3Eqp_0, dx4Psidp_0, dx5Psiqpp_0;
    double dx1d_1, dx2w_1, dx3Eqp_1, dx4Psidp_1, dx5Psiqpp_1;
    double Id, Iq;
    
    double Efd, LadIfd, Pmech, Vstab;
    
    double B, G;

    double IrNorton, IiNorton;
    
    gridpack::ComplexType p_INorton;

    boost::shared_ptr<BaseGovernorModel> p_governor;
    boost::shared_ptr<BaseExciterModel> p_exciter;
	boost::shared_ptr<BasePssModel> p_pss;
    double presentMag, presentAng;

    double Efdinit, Pmechinit;
            
    ////////////////////////////////////////////////////

    /*gridpack::ComplexType p_pelect, p_volt;
    gridpack::ComplexType p_mac_ang_s0, p_mac_spd_s0;
    gridpack::ComplexType p_mac_ang_s1, p_mac_spd_s1;
    gridpack::ComplexType p_dmac_ang_s0, p_dmac_spd_s0;
    gridpack::ComplexType p_dmac_ang_s1, p_dmac_spd_s1;
    gridpack::ComplexType p_eqprime, p_mech;
    gridpack::ComplexType p_eprime_s0, p_eprime_s1;
    gridpack::ComplexType p_INorton;*/
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
          & MVABase & Ra
          & p_INorton
          & p_bus_id;
      }

};
}  // dynamic_simulation
}  // gridpack
#endif
