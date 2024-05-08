/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   gdform.hpp
 *
 * WECC generic grid forming inverter model
 **/


#ifndef _gdform_h_
#define _gdform_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "base_gen_model.hpp"
//#include "base_plant_model.hpp"
#include "cblock.hpp"
#include "dblock.hpp"

class Gdform : public BaseEMTGenModel
{
  public:
    /**
     * Basic constructor
     */
    Gdform();

    /**
     * Basic destructor
     */
    virtual ~Gdform();

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
     * @param [output] values - array where initialized generator variables should be set
     */
  void init(gridpack::RealType *values);

    /**
     * Write output from generators to a string.
     * @param string (output) string with information to be printed out
     * @param bufsize size of string buffer in bytes
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if bus is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const int bufsize,
        const char *signal);

    /**
     * Write out generator state
     * @param signal character string used to determine behavior
     * @param string buffer that contains output
     */
    void write(const char* signal, char* string);

  /**
   * Return the generator current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
  void getCurrent(double *ia, double *ib, double *ic);

    /**
   * Return the generator real and reactive power
   * @param [input] time - current time
   * @param [output] Pg - generator real power
   * @param [output] Qg - generator reactive power
   *
   * Note: Power is on system MVA base
   */
  virtual void getPower(double time, double *Pg, double *Qg);

  /**
   * Return the generator frequency (pu)
   * @param [output] freq - machine frequency
   *
   * Note: Frequency is per unit. Steady-state frequency is 1.0
  */
  double getFreq();


  /**
   * Return the generator initial real and reactive power
   * @param [output] Pg(t0) - generator real power
   * @param [output] Qg(t0) - generator reactive power
   *
   * Note: Power is pu on system MVA base
   */
  virtual void getInitialPower(double *Pg, double *Qg);

  /**
   * Return the global location for the generator current injection 
   * @param [output] i_gloc - global location for the first current variable
   */
  void getCurrentGlobalLocation(int *i_gloc);

  /**
   * Return the rotor speed deviation
   * @param 
   */
  double getSpeedDeviation() { return 0.0; }

  /**
   * Return the speed deviation and its global location 
   * @param[output] rotor speed deviation
   * @param[output] rotor speed deviation global location
   */
  double getSpeedDeviation(int *dw_gloc)
  {
    return 0.0;
  }

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
   * Returns the initial field voltage (Efd(t0))
   * @param [out] Efd0 - Initial field voltage
   */
  double getInitialFieldVoltage() { return 0.0; }

  /**
   * Returns the initial mechanical power (Pmech(t0))
   * @param [out] Pmech0 - Initial mechanical power
   */
  double getInitialMechanicalPower() { return 0.0; }

  /**
   * Return the machine angle
   * @param [output] delta - machine angle
   */
  double getAngle()
  {
    return theta;
  }

  /**
   * Return the machine angle and its global location
   * @param [output] delta - machine angle
   */
  double getAngle(int *delta_gloc)
  {
    return theta;
    *delta_gloc = -1;
  }
  private:
    // parameters
    int lvplsw;
    double tg, rrpwr, brkpt, zerox, lvpl1, volim, lvpnt1, lvpnt0, lolim, tfltr, khv, iqrmax, iqrmin, accel;

    // Blocks
    Filter Ip_blk;
    double Ip; // Output of Ip_blk

    Filter Iq_blk;
    double Iq; // Output of Iq_blk

    Filter Vt_filter_blk;
    double Vt_filter; // Output of Vt_filter blk

    Slope  Lvpnt_blk;
    double Lvpnt_out; // Output of Lvpnt_blk

    Slope Lvpl_blk;
    double Lvpl_out; // Output of Lvpl_blk

    GainLimiter    Iqlowlim_blk;

    PIControl Pll_block;
    double    dw; // Output of PLL block

    Integrator angle_block;
    double    theta; // Output of angle_block
  
    double  Ipout,Iqout; // inverter current output in inverter reference frame
  double  Irout, Iiout; // inverter current output in network reference frame
  
  double Ipcmd, Iqcmd, busfreq;  // busfreq is perunit, 1.0
  // parameters
  double Ra, Xl, mq, kpv, kiv, mp, kppmax, kipmax, Pmax, Pmin;
  double Emax, Emin, Tpf, Imax, Qmax, Qmin;
  double kpqmax,kiqmax,Tqf,Tvf;
  int    Vflag;

  // Inputs set at steady-state (t=0)
  double Vset,Pset;

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

  gridpack::ComplexType Zsource;
  double L;

  // Internal variables
  double Vt, VR, VI, Im;
  gridpack::ComplexType E;
  double omega;
  double Edroop_max,Edroop_min; // Only set when current exceeds Imax
  bool   zero_Tpf,zero_Tqf,zero_Tvf;
  double Vthresh; // Threshold below which states are frozen

  double iabc[3], diabc[3];

  // phase voltages
  double vabc[3], eabc[3], vdq0[3], idq0[3], iout[3];

};

#endif
