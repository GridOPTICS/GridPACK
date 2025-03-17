/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file  wtara1.hpp
 * 
 * @brief  Header file for Aerodynamic controller
 * 
 * 
 */

#ifndef _wtara1_h_
#define _wtara1_h_

#include <base_mechanical_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "cblock.hpp"
#include "dblock.hpp"

class Wtara1 : public BaseEMTRMechModel
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  Wtara1();
  
  /**
   * Basic destructor
   */
   ~Wtara1();

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
   * setTheta - sets pitch angle
   * @param Theta - pitch angle
   **/
   void setTheta(double Theta);

  /**
   * setTurbineSpeedDeviation - sets the turbine speed deviation
   * @param domega_turb : turbine speed deviation
   * From drive train model
   **/
   void setTurbineSpeedDeviation(double domega_turb);

  /**
   *  setPmech - Mechanical power Pmech
   *  @param Pmech - Mechanical power input
   **/
  void setPmech(double Pmech);

  /**
   * getTaero - returns the aero-dynamic torque
   * @output Taero - aero dynamic torque, output of aerodynamic model
   **/
   double getTaero();

  /**
   * getTheta - Get initial value of pitch controller
   * @output theta0 - initial value of pitch controller
   **/
  double getTheta();

private:

  // Parameters
  double Ka, Theta0;
    
  // Model inputs
  double Theta; // External pitch angle
  double domega_t; // Turbine speed deviation
  double Pmech0;

  // Model outputs
  double Taero; // Aerodynamic torque
};

// Class for defining events for Reeca1 model
class Wtara1Event
  :public gridpack::math::RealDAESolver::Event
{
public:

  // Default constructor
  Wtara1Event(Wtara1 *wtara1):gridpack::math::RealDAESolver::Event(0),p_aerocon(wtara1)
  {
    std:fill(p_term.begin(),p_term.end(),false);

    std::fill(p_dir.begin(),p_dir.end(),gridpack::math::CrossZeroNegative);

  }

  // Destructor
  ~Wtara1Event(void) {}
protected:
  Wtara1 *p_aerocon;

  void p_update(const double& t, gridpack::RealType *state);

  void p_handle(const bool *triggered, const double& t, gridpack::RealType *state);
};


#endif
