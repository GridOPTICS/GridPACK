/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   reeca1.hpp
 * 
 * @brief  Reeca1 header template file
 * 
 * 
 */

#ifndef _reeca1_h_
#define _reeca1_h_

#include <base_exc_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "cblock.hpp"
#include "dblock.hpp"

class Reeca1 : public BaseEMTExcModel
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  Reeca1();
  
  /**
   * Basic destructor
   */
   ~Reeca1();

  /**
     Prestep function
  */
  void preStep(double time ,double timestep);

  /**
     Poststep function
  */
  void postStep(double time);

  /**
    Number of variables
  */ 
  void getnvar(int *nvar);

  /**
     Set the reference power inputs
  **/
   void setPrefQext(double Pref, double Qext) { }

  /**
     Get the power order - used by pitch controller model
  */
   double getPord() { return 0.0; }

  /**
     Set omega_g - from drive train model
  */
   void setOmega(double omega) { }

  /**
     Note: This is a custom version of the load method from the BaseComponent Class. It takes in an extra argument idx to specify which component is being read. Ideally, this method should be moved to the MatVecInterface

   * Load data from DataCollection object into corresponding
   * component. This needs to be implemented by every component
   * @param data data collection associated with component
   */
   void load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx);

  /**
   * Initialize exciter model before calculation
   * @param [output] values - array where initialized exciter variables should be set
   */
   void init(gridpack::RealType *values);
  
  /**
   * Write output from exciters to a string.
   * @param string (output) string with information to be printed out
   * @param bufsize size of string buffer in bytes
   * @param signal an optional character string to signal to this
   * routine what about kind of information to write
   * @return true if bus is contributing string to output, false otherwise
   */
   bool serialWrite(char *string, const int bufsize,
			   const char *signal);
  
  /**
   * Write out exciter state
   * @param signal character string used to determine behavior
   * @param string buffer that contains output
   */
   void write(const char* signal, char* string);

  /**
     set exciter status
  **/
  void setStatus(int estatus) {status = estatus;}
  
  /**
   * return the bolean indicating whether the exciter is ON or OFF
   */
  bool getStatus() {return status;}
  
  /**
   * Set bus voltage
   */
  void setVoltage(double busVD, double busVQ) {VD = busVD; VQ = busVQ; }

  /**
   * Copy over initial bus voltage from the bus (power flow solution)
   */
  void setInitialVoltage(double inVm,double inVa) {p_Vm0 = inVm; p_Va0 = inVa;}

  void setEvent(gridpack::math::RealDAESolver::EventManagerPtr);

  /**
   * Update the event function values
   */
  void eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues);

  /**
   * Event handler function 
   */
  void eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state);
  
  /**
   * Get number of matrix values contributed by generator
   * @return number of matrix values
   */
   int matrixNumValues();

  /**
   * Return values from Jacobian matrix
   * @param nvals: number of values to be inserted
   * @param values: pointer to matrix block values
   * @param rows: pointer to matrix block rows
   * @param cols: pointer to matrix block cols
   */
   void matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols);


  /**
   * Return vector values from the generator model 
   * @param values - array of returned values
   *
   * Note: This function is used to return the entries in vector,
   * for e.g., the entries in the residual vector from the generator
   * object
   */
   void vectorGetValues(gridpack::RealType *values);

  /**
   * Pass solution vector values to the generator object
   * @param values - array of returned values
   *
   * Note: This function is used to pass the entries in vector
   * to the generator object,
   * for e.g., the state vector values for this generator
   */
   void setValues(gridpack::RealType *values);

  /**
     Get the current command references
  **/
  void getIpcmdIqcmd(double *Ipcmdout, double *Iqcmdout);

  
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

  void CurrentLimitLogic(int PQFLAG,double Vt_filter, double Ipcmd, double Iqcmd,double *Ipmin_out, double *Ipmax_out, double *Iqmin_out, double *Iqmax_out);


};

// Class for defining events for Reeca1 model
class Reeca1Event
  :public gridpack::math::RealDAESolver::Event
{
public:

  // Default constructor
  Reeca1Event(Reeca1 *reeca1):gridpack::math::RealDAESolver::Event(0),p_exc(reeca1)
  {
    std:fill(p_term.begin(),p_term.end(),false);

    std::fill(p_dir.begin(),p_dir.end(),gridpack::math::CrossZeroNegative);

  }

  // Destructor
  ~Reeca1Event(void) {}
protected:
  Reeca1 *p_exc;

  void p_update(const double& t, gridpack::RealType *state);

  void p_handle(const bool *triggered, const double& t, gridpack::RealType *state);
};

#endif
