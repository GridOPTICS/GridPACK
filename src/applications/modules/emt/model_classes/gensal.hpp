/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   gensal.hpp
 * 
 * @brief: Header file for GENSAL model  
 * 
 * 
 */

#ifndef _gensal_h_
#define _gensal_h_

#include <base_gen_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "boost/smart_ptr/shared_ptr.hpp"

class GensalGen: public BaseGenModel
{
public:
  /**
   * Basic constructor
   */
  GensalGen();
  
  /**
   * Basic destructor
   */
  ~GensalGen();
  
  /**
   * Load parameters from DataCollection object into generator model
   * @param data collection of generator parameters from input files
   * @param index of generator on bus
   * TODO: might want to move this functionality to BaseGeneratorModel
   */
  void load(const boost::shared_ptr<gridpack::component::DataCollection>
	    data, int idx);
  
  /**
   * Saturation function
   * @ param x
   */
  double Sat(double x); 
  
  /**
   * Initialize generator model before calculation
   * @param [output] values - array where initialized generator variables should be set
   */
  void init(gridpack::ComplexType *values);
  
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
  
  double getAngle();
  
  /**
   * Write out generator state
   * @param signal character string used to determine behavior
   * @param string buffer that contains output
   */
  void write(const char* signal, char* string);
  
  /**
   *  Set the number of variables for this generator model
   *  @param [output] number of variables for this model
   */
  bool vectorSize(int *nvar) const;
  
  /**
   * Set the internal values of the voltage magnitude and phase angle. Need this
   * function to push values from vectors back onto generators
   * @param values array containing generator state variables
   */
  void setValues(gridpack::ComplexType*);
  
  /**
   * Return the values of the generator vector block
   * @param values: pointer to vector values
   * @return: false if generator does not contribute
   *        vector element
   */
  bool vectorValues(gridpack::ComplexType *values);
  
  /**
   *  Set Jacobian values
   *  @param values a 2-d array of Jacobian block for the bus
   */
  bool setJacobian(gridpack::ComplexType **values);
  
  /**
   * Return the generator current injection (in rectangular form) 
   * @param [output] IGD - real part of the generator current
   * @param [output] IGQ - imaginary part of the generator current
   */
  void getCurrent(double *IGD, double *IGQ);
  
  /* Return the field current */ 
  double getFieldCurrent();
  
  /* Return rotor speed deviation */ 
  double getRotorSpeedDeviation();
  
  /* return rotor speed deviation variable location.
     The variable location is w.r.t. variables in the bus array
  */
  int getRotorSpeedDeviationLocation();
  
private:
  // Machine parameters
  double H, D, Ra, Xd, Xq, Xdp, Xdpp, Xl;
  double Tdop, Tdopp, Tqopp, S10, S12, Tqop;
  double Xqp, Xqpp;
  
  double LadIfd; // Field current

  // Inputs
  double Efd; // Field voltage
  double Pmech; // Mechanical power

  // constants
  double B, G;
  
  // Variables and their derivatives
  double delta, dw, Eqp, Psidp, Psiqpp; 
  double ddelta, ddw, dEqp, dPsidp, dPsiqpp;
  // previous step values of the variables
  double deltaprev, dwprev, Eqpprev, Psidpprev, Psiqppprev; 
  
  double sat_A,sat_B; // Saturation function coefficients
};

#endif
