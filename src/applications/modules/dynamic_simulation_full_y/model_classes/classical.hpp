/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   classical.hpp
 * @author Bruce Palmer
 * @Last modified:   May 18, 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef _classical_gen_h_
#define _classical_gen_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_generator_model.hpp"

namespace gridpack {
namespace dynamic_simulation {
class ClassicalGenerator : public BaseGeneratorModel
{
  public:
    /**
     * Basic constructor
     */
    ClassicalGenerator();

    /**
     * Basic destructor
     */
    virtual ~ClassicalGenerator();

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
     * Set voltage on each generator
     */
    void setVoltage(gridpack::ComplexType voltage);

    /**
     * Write out generator state
     * @param signal character string used to determine behavior
     * @param string buffer that contains output
     */
    //void write(const char* signal, char* string);

    /**
     * Write output from generators to a string.
     * @param string (output) string with information to be printed out
     * @param bufsize size of string buffer in bytes
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if bus is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const int bufsize, const char *signal);

    /**
     * Get the roter angle
     */
    double getAngle();

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
    double p_mva, p_r, p_dstr, p_dtr;
    double p_d0, p_h;
    double p_PI;

    gridpack::ComplexType p_pelect, p_volt;
    gridpack::ComplexType p_mac_ang_s0, p_mac_spd_s0;
    gridpack::ComplexType p_mac_ang_s1, p_mac_spd_s1;
    gridpack::ComplexType p_dmac_ang_s0, p_dmac_spd_s0;
    gridpack::ComplexType p_dmac_ang_s1, p_dmac_spd_s1;
    gridpack::ComplexType p_eqprime, p_pmech;
    gridpack::ComplexType p_eprime_s0, p_eprime_s1;
    gridpack::ComplexType p_INorton;
    std::string p_ckt;
    int p_bus_id;

    double IrNorton, IiNorton;

    friend class boost::serialization::access;

    template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar &
          & p_sbase
          & p_pg & p_qg
          & p_status
          & p_mva & p_r
          & p_dstr & p_dtr
          & p_d0 & p_h
          & p_pelect & p_volt
          & p_mac_ang_s0 & p_mac_spd_s0
          & p_mac_ang_s1 & p_mac_spd_s1
          & p_dmac_ang_s0 & p_dmac_spd_s0
          & p_dmac_ang_s1 & p_dmac_spd_s1
          & p_eqprime & p_pmech
          & p_eprime_s0 & p_eprime_s1
          & p_INorton
          & p_ckt
          & p_bus_id;
      }

};
}  // dynamic_simulation
}  // gridpack
#endif
