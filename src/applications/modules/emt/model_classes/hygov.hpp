/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   hygov.hpp
 *
 * @brief HYGOV governor model header file  
 * 
 * 
 */

#ifndef _hygov_h_
#define _hygov_h_

#include <base_gov_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "cblock.hpp"
#include "dblock.hpp"

class Hygov: public BaseEMTGovModel
{
public:
  /**
   * Basic constructor
   */
  Hygov();
  
  /**
   * Basic destructor
     */
  ~Hygov();

  /**
    Number of variables
  */ 
  void getnvar(int *nvar);

  /**
     Prestep function
  */
  void preStep(double time ,double timestep);

  /**
     Poststep function
  */
  void postStep(double time);
  
  /**
   * Load parameters from DataCollection object into governor model
   * @param data collection of governor parameters from input files
   * @param index of governor on bus
   * TODO: might want to move this functionality to BaseGoviterModel
   */
  void load(const boost::shared_ptr<gridpack::component::DataCollection>
	    data, int idx);

  /**
   * Initialize governor model before calculation
   * @param [output] values - array where initialized governor variables should be set
   */
  void init(gridpack::RealType *values);
  
  /**
   * Write output from governors to a string.
   * @param string (output) string with information to be printed out
   * @param bufsize size of string buffer in bytes
   * @param signal an optional character string to signal to this
   * routine what about kind of information to write
   * @return true if bus is contributing string to output, false otherwise
   */
  bool serialWrite(char *string, const int bufsize,
		   const char *signal);
  
  /**
   * Write out governor state
   * @param signal character string used to determine behavior
   * @param string buffer that contains output
   */
  void write(const char* signal, char* string);
  
  /**
   * Set the internal values of the voltage magnitude and phase angle. Need this
   * function to push values from vectors back onto governors
   * @param values array containing governor state variables
   */
  void setValues(gridpack::RealType*);
  
  /**
   * Return the values of the governor vector block
   * @param values: pointer to vector values
   */
  void vectorGetValues(gridpack::RealType *values);

  /**
   * Get number of matrix values contributed by generator
   * @return number of matrix values
   */
  int matrixNumValues();

  /**
 * Return values from a matrix block
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
  void matrixGetValues(int *nvals,gridpack::RealType *values,
      int *rows, int *cols);

  /**
   * Set the mechanical power during initialization inside the governor
   * @param pmech value of the mechanical power 
   */
  void setInitialMechanicalPower(double pmech);
  
  /** 
   * Get the value of the mechanical power parameter
   * @return value of the mechanical power 
   */
  double getMechanicalPower();

  /** 
   * Get the value of the mechanical power and its global location
   * @return value of the mechanical power
   *
   * Note: Used in Jacobian calculation
   */
  double getMechanicalPower(int *Pmech_gloc);

  
  /**
   * Partial derivatives of Mechanical Power Pmech w.r.t. governor variables
   * @param xgov_loc locations of governor variables
   * @param dPmech_dxgov partial derivatives of mechanical power Pmech w.r.t governor variables
   */
  bool getMechanicalPowerPartialDerivatives(int *xgov_loc,double *dPmech_dxgov);
  
  /**
   * Set the value of the Vcomp
   * @return value of teh Vcomp
   */
  void setVcomp(double vtmp);
  
    /**
   * Set Event 
   */
  void setEvent(gridpack::math::RealDAESolver::EventManagerPtr);

  /**
   * Update the event function values
   */
  void eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues);

  /**
   * Event handler function 
   */
  void eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state);

private:
  
    // Governor Hygov Parameters read from dyr
    double R, r, TR, TF, TG, VELM, GMAX, GMIN, TW, AT, Dt, qNL;

    // Hygov control blocks
    Filter    filter_block; // output e
    PIControl gate_block; // output c
    Filter    opening_block; // output g
    Integrator turbine_flow_block; // output q

    double opening_block_out, gate_block_out, turbine_flow_block_out, filter_block_in;
    // Outputs: Mechnical Power
    double Pmech;
	double xout;

    // Inputs
    double nref;     // reference
    double delta_w;  // speed deviation

  //// Hygov state variables
  //double x1;  // Turbine power 
  //double x2;  // Valve position

  //// Hygov output
  //double xout;
  //
  //// Hygov state-variable derivatives
  //double dx1, dx2;

  //// Flags for limiter
  //bool x1_at_min,x1_at_max;

};

// Class for defining events for hygov model
class HygovEvent
  :public gridpack::math::RealDAESolver::Event
{
public:

  // Default constructor
  HygovEvent(Hygov *gov):gridpack::math::RealDAESolver::Event(2),p_gov(gov)
  {
    std:fill(p_term.begin(),p_term.end(),false);

    std::fill(p_dir.begin(),p_dir.end(),gridpack::math::CrossZeroNegative);

  }

  // Destructor
  ~HygovEvent(void) {}
protected:
  Hygov *p_gov;

  void p_update(const double& t, gridpack::RealType *state);

  void p_handle(const bool *triggered, const double& t, gridpack::RealType *state);
};

#endif
