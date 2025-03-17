/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file  wtdta1.hpp
 * 
 * @brief  Header file for WTDTA1 drive train controller
 * 
 * 
 */

#ifndef _wtdta1_h_
#define _wtdta1_h_

#include <base_mechanical_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "cblock.hpp"
#include "dblock.hpp"

class Wtdta1 : public BaseEMTRMechModel
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  Wtdta1();
  
  /**
   * Basic destructor
   */
   ~Wtdta1();

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
   * setTmech - sets mechanical torque
   * @param Tmech - mechanical torque
   * From aerodynamic model
   **/
   void setTmech(double Tmech);

  /**
   * getTmech - gets initial mechanical torque
   * @returnh Tmech- initial mechanical torque
   **/
   double getTmech();

  /**
   * setTelec - sets electrical torque
   * @param Telec - electrical torque
   * From generator?
   **/
   void setTelec(double Tmech);

  /**
   * getTurbineSpeedDeviation - gets the turbine speed deviation
   * @ouput turbine speed deviation
   **/
   double getTurbineSpeedDeviation();

  /**
   * getGenSpeedDeviation - gets the generator speed deviation
   * @ouput generator speed deviation
   **/
   double getGeneratorSpeedDeviation();

  /**
   * getGenRotorAngleDeviation - gets the rotor angle deviation
   * @ouput rotor angle deviation
   **/
   double getRotorAngleDeviation();

  /**
   *  setOmegaref - Output of torque controller
   * @param - omega_ref : reference speed
   * Only set during initialization
   **/
  void setOmegaref(double omega_ref);

private:
  // Parameters
  double H, D, Hfrac;
  double Freq1, Dshaft;
  
  // Model inputs
  double Tm; // Mechanical torque
  double Te; // Electrical torque
  double s0; // Initial slip

  // Model outputs
  double domega_g; // speed deviation
  double dtheta_g; // rotor angle speed deviation
  double domega_t; // Turbine speed deviation

  // Blocks
  Integrator domegag_blk;
  Integrator domegat_blk;
  Integrator dthetag_blk;
  Integrator Tshaft_blk;

  // Internal variables
  double Hg; // generator inertia
  double Ht; // turbine inertia
  double Kshaft;
  bool   double_mass; // True -> also includes turbine piece

  int p_bus_num;
};

// Class for defining events for Reeca1 model
class Wtdta1Event
  :public gridpack::math::RealDAESolver::Event
{
public:

  // Default constructor
  Wtdta1Event(Wtdta1 *wtdta1):gridpack::math::RealDAESolver::Event(0),p_dtrain(wtdta1)
  {
    std:fill(p_term.begin(),p_term.end(),false);

    std::fill(p_dir.begin(),p_dir.end(),gridpack::math::CrossZeroNegative);

  }

  // Destructor
  ~Wtdta1Event(void) {}
protected:
  Wtdta1 *p_dtrain;

  void p_update(const double& t, gridpack::RealType *state);

  void p_handle(const bool *triggered, const double& t, gridpack::RealType *state);
};


#endif
