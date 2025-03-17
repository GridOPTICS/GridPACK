/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file  wttqa1.hpp
 * 
 * @brief  Header file for WTTQA1 torque controller
 * 
 * 
 */

#ifndef _wttqa1_h_
#define _wttqa1_h_

#include <base_mechanical_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "cblock.hpp"
#include "dblock.hpp"

class Wttqa1 : public BaseEMTRMechModel
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  Wttqa1();
  
  /**
   * Basic destructor
   */
   ~Wttqa1();

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
  
   void resetEventFlags(void) {}

    /**
   *  Set Pref0 - Sets reference power
   *  @param Pref0 - reference power
   *  Set by plant controller model
   **/
   void setPref0(double Pref0);

  /**
   *  setPelec - Electrical power Pelec
   *  @param Pelec - Electrical power input
   **/
   void setPelec(double Pelec);

  /**
   *  setGeneratorSpeedDeviation - Set the speed deviation
   *  @param - domega : speed deviation
   *  From drive train model
   **/
   void setGeneratorSpeedDeviation(double domega);

  /**
   * setVdip - Voltage dip flag
   * @param vdip - flag to indicate voltage dip
   * From elecrical controller
   **/
   void setVdip(bool Vdip);
  
  /**
   * getPref - Output of torque controller
   * @param  - Pref : reference power
   **/
   double getPref();

  /**
   *  getOmegaref - Output of torque controller
   * @param - omega_ref : reference speed
   **/
   double getOmegaref();

private:
    // Parameters
  int    Tflag;
  double Kpp, Kip, Tp, Twref, Temax, Temin;
  double p1, spd1, p2, spd2, p3, spd3, p4, spd4;
  double Trate;
  
  // Inputs:
  double Pelec; // Generator real power ouput
  double Pref0; // Plant controller reference
  double domega_g; // speed deviation
  bool   Vdip;   // Voltage dip flag
  
  // Outputs
  double Pref; // Reference power output
  double omega_ref; // Reference speed (pu)

  // Internal variables
  double MBase; // Machine base

  // Control blocks
  Filter Pelec_filter_blk; // Pelec filter block
  double Pelec_filter_blk_out; // output of Pelec filter block
  
  Filter wref_filter_blk; // wref filter block

  PIControl Tref_pi_blk; // PI controller for Tref
  double    Tref_pi_blk_out; // Output of Tref PI control block

  PiecewiseSlope Pomega_blk; // Power speed piecewise linear block
  double         Pomega_blk_out; // Output of Pomega_blk
};

// Class for defining events for Reeca1 model
class Wttqa1Event
  :public gridpack::math::RealDAESolver::Event
{
public:

  // Default constructor
  Wttqa1Event(Wttqa1 *wttqa1):gridpack::math::RealDAESolver::Event(0),p_torquecon(wttqa1)
  {
    std:fill(p_term.begin(),p_term.end(),false);

    std::fill(p_dir.begin(),p_dir.end(),gridpack::math::CrossZeroNegative);

  }

  // Destructor
  ~Wttqa1Event(void) {}
protected:
  Wttqa1 *p_torquecon;

  void p_update(const double& t, gridpack::RealType *state);

  void p_handle(const bool *triggered, const double& t, gridpack::RealType *state);
};


#endif
