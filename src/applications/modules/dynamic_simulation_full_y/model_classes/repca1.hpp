/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   repca1.hpp
 * @author Shrirang Abhyankar
 * @Updated November 20, 2022
 * 
 * @brief  
 * 
 * 
 */

#ifndef _repca1_h_
#define _repca1_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_plant_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Repca1Model : public BasePlantControllerModel
{
  public:
    /**
     * Basic constructor
     */
    Repca1Model();

    /**
     * Basic destructor
     */
    virtual ~Repca1Model();

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

    void setGenPQV(double P, double Q, double Vt);

    void setBusFreq(double freq);
	
    void setPrefQext(double Pref, double Qext);
	
    void setExtBusNum(int ExtBusNum);
  
    void setExtGenId(std::string ExtGenId);
	
    double getPref( );
	
    double getQext( );
	
  private:

  // Model parameters
  double ireg;   // Bus number of regulated bus
  double fbus, tbus; // To and from bus of the monitored branch
  std::string  br_ckt; // Branch circuit id
  int VCompFLAG; // droop flag
  int RefFLAG;  // Q control or V control
  int FreqFLAG;   // Freq control
  
  double Rc, Xc; // line compensation parameters
  double Kc;     // Gain for branch reactive power measurement
  double Tfltr;  // Filter time constant
  double dbd1,dbd2; // Deadband limits
  double Emax,Emin;  // Limits for voltage
  double Kp, Ki;     // PI controller gains for reactive power
  double Qmax,Qmin;  // Limits on Q PI controller
  double Tft, Tfv;    // Parameters for Qext lead lag block
  double Tp;        // Pbranch filter time constant
  double fdbd1,fdbd2; // Deadband limits for frequency error
  double Ddn, Dup;   // Droop constant
  double Kpg, Kig;   // Pext PI contoller gains
  double femax,femin; // Limits on Frequency error
  double Pmax, Pmin;  // Non wind-up limits for Pext PI controller
  double Tlag;      // Pext filter time constant
  double Tg;
  double Vfrz;      // State freeze voltage

  std::string p_gen_id;
  int p_bus_num;
  
  // Model inputs
  double Pbranch, Qbranch; // Line flow on designated branch (on system MVAbase)
  double Vreg;             // Voltage of regulated bus
  double Ibranch;          // Line current on designated branch (on system Ibase)
  double Freq;             // Frequency at POI

  // Model outputs
  double Pext,Qext;        // Output reference signals from the plant.
  double Pg, Qg;           // Generator active and reactive power (on machine Mbase

  // Internal variables
  double Vt;
  double p_sbase,p_mbase; // System base and machine Mbase
};
}  // dynamic_simulation
}  // gridpack
#endif
