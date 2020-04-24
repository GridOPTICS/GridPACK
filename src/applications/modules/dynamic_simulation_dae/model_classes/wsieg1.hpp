/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wsieg1.hpp
 * @author Shuangshuang Jin 
 * @author Shrirang Abhyankar
 * @Last modified:   04/22/20
 *
 * @brief WSIEG1 governor model header file  
 * 
 * 
 */

#ifndef _wsieg1_h_
#define _wsieg1_h_

#include <base_gov_model.hpp>
#include <gridpack/include/gridpack.hpp>

class Wsieg1Gov: public BaseGovModel
{
public:
  /**
   * Basic constructor
   */
  Wsieg1Gov();
  
  /**
   * Basic destructor
     */
  ~Wsieg1Gov();
  
  /**
   * Load parameters from DataCollection object into governor model
   * @param data collection of governor parameters from input files
   * @param index of governor on bus
   * TODO: might want to move this functionality to BaseGoviterModel
   */
  void load(const boost::shared_ptr<gridpack::component::DataCollection>
	    data, int idx);

  /**
   * Set Jacobian block
   * @param values a 2-d array of Jacobian block for the bus
   */
  bool setJacobian(gridpack::ComplexType **values);
  
  /**
   * Initialize governor model before calculation
   * @param [output] values - array where initialized governor variables should be set
   */
  void init(gridpack::ComplexType *values);
  
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
   *  Set the number of variables for this governor model
   *  @param [output] number of variables for this model
   */
  bool vectorSize(int *nvar) const;
  
  /**
   * Set the internal values of the voltage magnitude and phase angle. Need this
   * function to push values from vectors back onto governors
   * @param values array containing governor state variables
   */
  void setValues(gridpack::ComplexType*);
  
  /**
   * Return the values of the governor vector block
   * @param values: pointer to vector values
   * @return: false if governor does not contribute
   *        vector element
   */
  bool vectorValues(gridpack::ComplexType *values);
  
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
  
private:
  
  // Governor WSIEG1 Parameters read from dyr
  double K, T1, T2, T3, Uo, Uc, Pmax, Pmin;
  double T4, K1, K2, T5, K3, K4, T6, K5, K6, T7, K7, K8;
  double Db1, Err, Db2;
  double Gv1, PGv1, Gv2, PGv2, Gv3, PGv3, Gv4, PGv4, Gv5, PGv5;
  double Iblock;
  
  // WSIEG1 state variables
  double xLL; // Lead-lag block state
  double xGV; // Governor output
  double xT1; // First turbine integrator output
  double xT2; // Second turbine integrator output
  double xT3; // Third turbine integrator output
  double xT4; // Fourth turbine integrator output

  // WSIEG1 state-variable derivatives
  double dxLL, dxGV, dxT1, dxT2, dxT3, dxT4;

  // WSIEG1 previous step solution
  double xLLprev, xGVprev, xT1prev, xT2prev, xT3prev, xT4prev;

  // Outputs: Mechnical Power Gen1 and Gen 2
  double Pmech1, Pmech2;

  bool SecondGenExists;

  // Inputs
  double Pref;

  int iseq_diff[6];   
};

#endif
