/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   sexs.hpp
 * 
 * @brief  


                                                        EMAX
                                                       ------
                                                      /
                                                     /
                            -----------        ------------ 
 *                          |          |       |           |
 *                          | 1 + TAs  |       |     K     |
 *   Vref - Ec + Vs ------> |--------- |------>| --------- |-----> Efd
                            | 1 + TBs  |       |  1 + TEs  |
 *                          |          |       |           |
 *                          ------------       ------------
                                                     / 
                                                    /
                                               -----
                                                EMIN
*/                 

#ifndef _sexs_h_
#define _sexs_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_exciter_model.hpp"
#include <string>

namespace gridpack {
namespace dynamic_simulation {
class SexsModel : public BaseExciterModel
{
  public:
    /**
     * Basic constructor
     */
    SexsModel();

    /**
     * Basic destructor
     */
    virtual ~SexsModel();

    /**
     * Load parameters from DataCollection object into exciter model
     * @param data collection of exciter parameters from input files
     * @param index of exciter on bus
     * TODO: might want to move this functionality to BaseExciterModel
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
     * @param initial value of the field voltage
     */
    void setFieldVoltage(double fldv);

    /** 
     * Get the value of the field voltage parameter
     * @return value of field voltage
     */
    double getFieldVoltage();

    /** 
     * Set the value of the Vterminal
     * 
     */
    void setVterminal(double mag);

    /**
     * Set the initial voltage stabilizer signal
     * @param stabilizing signal input
    **/
    void setVstab(double vstab);
	
    /** 
     * Set the exciter bus number
     * @return value of exciter bus number
     */
    void setExtBusNum(int ExtBusNum);
	
    /** 
     * Set the exciter generator id
     * @return value of generator id
     */
    void setExtGenId(std::string ExtGenId);

  private:

    // Internal variables
    bool zero_TE; // If TE == 0 then the filter block is replaced by gain block

    // Exciter SEXS parameters from dyr
    double TA_OVER_TB, TA, TB, K, TE, EMIN, EMAX;

    // Linear control blocks
    LeadLag leadlagblock;
    Filter  filterblock;
    Gain    gainblock;

    // Model inputs
    double Ec; // Terminal voltage
    double Vref; // Reference voltage
    double Vs;  // Stabilizing voltage signal
  
    // Model outputs 
    double Efd;     // Field Voltage

    std::string p_ckt; // id of the generator where the exciter is installed on
    int p_bus_id;  // bus number of the generator 

};
}  // dynamic_simulation
}  // gridpack
#endif
