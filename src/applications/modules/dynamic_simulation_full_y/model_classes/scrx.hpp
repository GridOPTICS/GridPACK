/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   scrx.hpp
 * @author Yousu Chen
 * @date   2026-03-27
 *
 * @brief  SCRX — Simple Controllable Rectifier Exciter
 *
 *  Topology (identical to SEXS, plus CSWITCH bus-fed option):
 *
 *  CSWITCH = 0 : bus-fed  — limits scale with terminal voltage Vt
 *  CSWITCH = 1 : solid-fed — fixed limits (same as SEXS)
 *
 *  DYR format:
 *    I, 'SCRX', ID, TATB, TB, K, TE, EMIN, EMAX, CSWITCH
 */

#ifndef _scrx_h_
#define _scrx_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"
#include <string>

#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class ScrxModel : public BaseExciterModel
{
  public:
    /**
     * Basic constructor
     */
    ScrxModel();

    /**
     * Basic destructor
     */
    virtual ~ScrxModel();

    /**
     * Load parameters from DataCollection object into exciter model
     * @param data collection of exciter parameters from input files
     * @param index of exciter on bus
     */
    void load(boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Initialize exciter model before calculation
     * @param mag voltage magnitude
     * @param ang voltage angle
     * @param ts time step
     */
    void init(double mag, double ang, double ts);

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
     * Set the initial field voltage value
     */
    void setFieldVoltage(double fldv);

    /**
     * Get the value of the field voltage parameter
     */
    double getFieldVoltage();

    /**
     * Set the value of the terminal voltage
     */
    void setVterminal(double mag);

    /**
     * Set the voltage stabilizer signal
     */
    void setVstab(double vstab);

    /**
     * Set the exciter bus number
     */
    void setExtBusNum(int ExtBusNum);

    /**
     * Set the exciter generator id
     */
    void setExtGenId(std::string ExtGenId);

    /**
     * Set internal state parameter in exciter
     */
    bool setState(std::string name, double value);

    /**
     * Get internal state parameter in exciter
     */
    bool getState(std::string name, double *value);

  private:

    bool zero_TE; // If TE == 0 the filter block is replaced by a gain block

    // SCRX parameters
    double TA_OVER_TB, TA, TB, K, TE, EMIN, EMAX;
    int    CSWITCH; // 0 = bus-fed (limits × Vt), 1 = solid-fed (fixed limits)

    // Linear control blocks
    LeadLag    leadlagblock;
    Filter     filterblock;
    GainLimiter gainblock;

    // Model inputs
    double Ec;   // Terminal voltage
    double Vref; // Reference voltage
    double Vs;   // Stabilizing voltage signal

    // Model output
    double Efd;  // Field voltage

    std::string p_ckt;
    int p_bus_id;
};
}  // dynamic_simulation
}  // gridpack
#endif
