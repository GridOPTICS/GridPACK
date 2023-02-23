/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   gdform.hpp
 * @author Renke Huang
 * 
 * @brief  
 * Grid-Forming inverter model 
 * 
 * @updated
 * @author Shrirang Abhyankar
 * Model re-implementation with cblocks
 * Jan. 6, 2023
 *
 * Updated version of model including Qlimit dynamics and Vflag
 */


#ifndef _gdform_h_
#define _gdform_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_generator_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class GridFormingGenerator : public BaseGeneratorModel
{
  public:
    /**
     * Basic constructor
     */
    GridFormingGenerator();

    /**
     * Basic destructor
     */
    virtual ~GridFormingGenerator();

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
     * Set frequency
     */
    void setFreq(double dFreq);

    // Called by both predictor and corrector to compute the model
  void computeModel(double t_inc, IntegrationStage int_flag,bool flag);

  // Update the limits on E for current limiting

  gridpack::ComplexType CurrentLimitLogic(gridpack::ComplexType I);
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

  private:

  double p_sbase;
  double p_mbase;
  double p_pg, p_qg;
  int p_status;
  bool p_tripped;

  // parameters
  double Ra, Xl, mq, kpv, kiv, mp, kppmax, kipmax, Pmax, Pmin;
  double Emax, Emin, Tpf, Imax, Qmax, Qmin;
  double kpqmax,kiqmax,Tqf,Tvf;
  int    Vflag;

  // Inputs set at steady-state (t=0)
  double Vset,Pset;
  gridpack::ComplexType Zsource;

  // Constants
  double omega0; // base angular velocity

  // Blocks
  Filter P_filter_blk;
  double Pinv; // Output of P filter block

  Filter Q_filter_blk;
  double Qinv; // Output of Q filter block

  Filter V_filter_blk;
  double Vmeas; // Output of V filter block

  PIControl Edroop_PI_blk; // PI control for E

  GainLimiter Edroop_limiter_blk; // Limiter for E when Vflag = 0

  double Edroop; // Output of Edroop PI block OR Edroop limiter block
  
  PIControl Pmax_PI_blk; // PI controller for Pmax limit
  double Pmax_PI_blk_out; // Output of PI controller for Pmax limit

  PIControl Pmin_PI_blk; // PI controller for Pmin limit
  double Pmin_PI_blk_out; // Output of PI controller for Pmin limit

  PIControl Qmax_PI_blk; // PI controller for Qmax limit
  double Qmax_PI_blk_out; // Output of PI controller for Qmax limit

  PIControl Qmin_PI_blk; // PI controller for Qmin limit
  double Qmin_PI_blk_out; // Output of PI controller for Qmin limit

  Integrator Delta_blk; // Integrator block for PLL
  double delta;         // Output of integrator block

  // Internal variables
  double busfreq;
  std::string p_gen_id;
  int p_bus_num;
  double Vt, theta, VR, VI, Im;
  gridpack::ComplexType E,V,I,S;
  double omega;
  gridpack::ComplexType p_Norton_Ya;
  double B, G;
  double Edroop_max,Edroop_min; // Only set when current exceeds Imax
  bool   zero_Tpf,zero_Tqf,zero_Tvf;
  double Vthresh; // Threshold below which states are frozen

  // Output
  gridpack::ComplexType p_INorton;
  
  friend class boost::serialization::access;

    template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar &
          & p_sbase
          & p_pg & p_qg
          & p_status
          & p_mbase 
          & p_INorton
          & p_bus_num;
      }

};
}  // dynamic_simulation
}  // gridpack
#endif
