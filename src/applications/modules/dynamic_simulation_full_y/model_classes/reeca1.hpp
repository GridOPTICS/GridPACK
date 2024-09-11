/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   reeca1.hpp
 * @author Shrirang Abhyankar
 * @Created November 17, 2022
 *
 * @brief
 * Renewable energy electrical controller model REECA1
 *
 */

#ifndef _reeca1_h_
#define _reeca1_h_

#include "base_generator_model.hpp"
#include "boost/smart_ptr/shared_ptr.hpp"
#include "cblock.hpp"
#include "dblock.hpp"

namespace gridpack {
namespace dynamic_simulation {
class Reeca1Model : public BaseExciterModel {
public:
  /**
   * Basic constructor
   */
  Reeca1Model();

  /**
   * Basic destructor
   */
  virtual ~Reeca1Model();

  /**
   * Load parameters from DataCollection object into exciter model
   * @param data collection of exciter parameters from input files
   * @param index of exciter on bus
   * TODO: might want to move this functionality to BaseExciterModel
   */
  void load(boost::shared_ptr<gridpack::component::DataCollection> data,
            int idx);

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
   * Set the value of the Vterminal
   */
  void setVterminal(double mag);

  /** 
   * Set the value of the omega_g
   */
  void setOmega(double omega);

  /* 
   * This function is only used during initialization
   */
  void setIpcmdIqcmd(double ipcmd, double iqcmd);

  void setPrefQext(double pref, double qref);

  double getIpcmd();
  double getIqcmd();

  /*
   * Return active power output Pord - Used by pitch controller model
   */
  double getPord();

  /**
   * Set internal state parameter in exciter
   * @param name character string corresponding to state variable
   * @param value new value for state parameter
   * @return false if no variable corresponding to name is found
   */
  bool setState(std::string name, double value);

  /**
   * Get internal state parameter in exciter
   * @param name character string corresponding to state variable
   * @param value current value for state parameter
   * @return false if no variable corresponding to name is found
   */
  bool getState(std::string name, double *value);
	
  
  // The next two methods should not be in this class. Since Pref and Qext is
  // not an output of this model. This is simply used because the plant
  // controller initialization gets done in the generator controller model. That
  // initialization should be done here instead.
  double getPref();
  double getQext();

  // Yuan added below 2020-6-23
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

  /**
   * Set generator active and reactive power (machine MVA base)
   * 
   */
  void setGeneratorPower(double Pg, double Qg);

  /**
     Current limit logic
  **/
  void CurrentLimitLogic(int PQFLAG,double Vt_filter, double Ipcmd, double Iqcmd, double *Ipmin, double *Ipmax, double *Iqmin, double *Iqmax);

  /**
   *  ComputeModel - Used for both predictor and corrector methods. Single method that computes the model.
   **/ 
  void computeModel(bool Voltage_dip, int Iqinj_sw,double t_inc, IntegrationStage int_flag);
  
  /**
   * Write output from exciter to a string.
   * @param string (output) string with information to be printed out
   * @param bufsize size of string buffer in bytes
   * @param signal an optional character string to signal to this
   * routine what about kind of information to write
   * @return true if governor is contributing string to output, false otherwise
   */
  bool serialWrite(char *string, const int bufsize, const char *signal);

private:
  double Vt;               // Terminal voltage magnitude
  double p_sbase, p_mbase; // System and machine MVA base

  // Model inputs
  double Pref, Qref;
                       // initialization
  double Pgen, Qgen; // Generator real and reactive power on machine MVAbase

  // Model outputs
  double Ipcmd, Iqcmd;

  // Model Parameters
  double Vdip, Vup; // Lower and upper thresholds for reactive current injection
                    // activation
  double Trv;       // time constant for voltage filter
  double Vp1, Ip1, Vp2, Ip2, Vp3, Ip3, Vp4,
      Ip4; // Piecewise linear V-I points for real power control
  double Vq1, Iq1, Vq2, Iq2, Vq3, Iq3, Vq4,
      Iq4; // Piecewise linear V-I points for real power control
  double dbd1, dbd2, Kqv, Imax, Tiq, dPmax, dPmin, Pmax, Pmin, Tpord;
  double Iqh1, Iql1, Vref0, Iqfrz;
  double Thld, Thld2, Tp;
  double Qmax, Qmin, Vmax, Vmin;
  double Kqp, Kqi, Kvp, Kvi, Vbias;

  // Model control flags
  int PFFLAG; // Power factor control (1 power factor control, 0 = Q control)
  int QFLAG;  // Reactive power control (1 voltage or Q control, constant pf
                 // or Q control)
  int VFLAG;  // Voltage control (1 if Q control, 0 voltage control)
  int PFLAG;  // active current has speed dependency
  int PQFLAG; // 0 - Q priority, 1 - P priority

  // Limits on active and reactive power current
  double Ipmin, Ipmax, Iqmin, Iqmax;

  // Blocks
  Filter Vt_filter_blk; // Vt filter block
  double Vt_filter;     // Output of Vt filter block

  Deadband V_err_deadband;
  double V_err; // Output of voltage error block

  GainLimiter Iqv_limit_blk; // Limiter for Iqv
  double Iqv;           // Input for limiter block
  double Iqinj;         // Output of limiter block

  GainLimiter Iqcmd_limit_blk; // Limiter for Iqcmd

  Filter Pe_filter_blk; // Electrical power filter
  double Pe_filter;     // Output of electrical power filter

  GainLimiter Qlim_blk;   // Q limiter block with limits Qmin and Qmax
  double Qlim_out; // output of Q limiter block

  PIControl Q_PI_blk; // PI control for reactive power
  double V_PI;        // Output of Q_PI_blk

  GainLimiter Vlim_blk;   // Voltage limiter block with limits Vmin and Vmax
  double Vlim_out; // Output of voltage limiter block

  PIControl Verr_PI_blk; // PI Control for voltage error
  double Verr_PI_out;    // Output of Verr PI control block

  Filter Iq_lag_blk; // Lag for Iq current block
  double Iq_lag_out; // Output of Iq lag block

  GainLimiter Vt_filter_lowcap_blk; // Lower cap for Vt_filter used in division (limit
                             // 0.01)
  double Vt_filter_lowcap_out; // Output of Vt filter lower cap

  PiecewiseSlope VDL1;
  double VDL1_out;

  PiecewiseSlope VDL2;
  double VDL2_out;

  GainLimiter Pref_limit_blk;   // Pref limiter
  double Pref_limit_out; // Output of Pref limiter block

  Filter Pord_blk; // Pord filter block
  double Pord;     // Output of Pord filter block

  GainLimiter Ipcmd_limit_blk; // Limiter for Iqcmd

  bool p_has_drivetrain; // Is the model connected to drive train?
  double omega_g; // Input from drive train (1.0 if not available)

  int Iqinj_sw; // Iqinj switch

  bool Voltage_dip; // True if voltage dip
  bool Voltage_dip_prev; // Variable to hold Voltage_dip value at last time-step. Used in Iqinj_sw logic
  
  double thld_timer; // Timer for thld
  
  std::string p_gen_id;
  int p_bus_num;

  bool getVoltageDip(double);
};
} // namespace dynamic_simulation
} // namespace gridpack
#endif
