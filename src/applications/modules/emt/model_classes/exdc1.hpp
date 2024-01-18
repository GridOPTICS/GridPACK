/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   exdc1.hpp
 * 
 * @brief EXDC1 exciter model header file 
 * 
 * 
 */

#ifndef _exdc1_h_
#define _exdc1_h_

#include <base_exc_model.hpp>
#include <gridpack/include/gridpack.hpp>

// This model is also termed as 'ESDC1A' in PSSE

class Exdc1: public BaseEMTExcModel
{
public:
  /**
   * Basic constructor
   */
  Exdc1();
  
  /**
   * Basic destructor
   */
  ~Exdc1();
  
  /**
   * Load parameters from DataCollection object into exciter model
   * @param data collection of exciter parameters from input files
   * @param index of exciter on bus
   * TODO: might want to move this functionality to BaseExciterModel
   */
    void load(const boost::shared_ptr<gridpack::component::DataCollection>
	      data, int idx);
  
  /**
   * Set Jacobian block
   * @param values a 2-d array of Jacobian block for the bus
   */
  bool setJacobian(gridpack::ComplexType **values);
  
  /**
   * Initialize exciter model before calculation
   * @param [output] values - array where initialized exciter variables should be set
   */
  void init(gridpack::ComplexType *values);
  
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
  void matrixGetValues(int *nvals, gridpack::ComplexType *values, int *rows, int *cols);

  /**
   * Return vector values from the generator model 
   * @param values - array of returned values
   *
   * Note: This function is used to return the entries in vector,
   * for e.g., the entries in the residual vector from the generator
   * object
   */
  void vectorGetValues(gridpack::ComplexType *values);

  /**
   * Pass solution vector values to the generator object
   * @param values - array of returned values
   *
   * Note: This function is used to pass the entries in vector
   * to the generator object,
   * for e.g., the state vector values for this generator
   */
  void setValues(gridpack::ComplexType *values);

  /** 
   * Get the value of the field voltage parameter
   * @return value of field voltage
   */
  double getFieldVoltage();

  /** 
   * Get the value of the field voltage parameter
   * and its global location
   * @return value of field voltage
   */
  double getFieldVoltage(int *Efd_gloc);

  /**
   * Partial derivatives of field voltage Efd w.r.t. exciter variables
   * @param xexc_loc locations of exciter variables
   * @param dEfd_dxexc partial derivatives of field voltage w.r.t exciter variables
   * @param dEfd_dxgen partial derivatives of field voltage Efd w.r.t generator variables
   */
  bool getFieldVoltagePartialDerivatives(int *xexc_loc,double *dEfd_dxexc,double *dEfd_dxgen);

  /**
   * Set Event 
   */
  void setEvent(gridpack::math::DAESolver::EventManagerPtr);

  /**
   * Update the event function values
   */
  void eventFunction(const double&t,gridpack::ComplexType *state,std::vector<std::complex<double> >& evalues);

  /**
   * Event handler function 
   */
  void eventHandlerFunction(const bool *triggered, const double& t, gridpack::ComplexType *state);

  /**
   * Updated limiter flags after event has occured. Only called when the network is resolved
   */
  void resetEventFlags(void);
  
private:
  
  // Exciter exdc1 parameters from dyr
  double TR, VRmax, VRmin, TC, TB;
  double KA, TA, KE, TE;
  double KF, TF;
  int    SWITCH;
  
  // EXDC1 state variables
  double Vmeas; // Measured voltage by transducer
  double xLL;   // Lead-lag block state variable
  double VR;    //   Voltage regulator output
  double Efd;   // Exciter output
  double xf;   // Feedback block state variable
  
  // EXDC1 derivatives
  double dVmeas, dxLL, dVR, dEfd, dxf;    
  
  // Saturation function data points
  double E1, SE1, E2, SE2;

  // Efd threshold for satuarion (Efd < Efdthresh => SE = 0)
  double Efdthresh;

  // Saturation function parameters (SE = B*(Efd - A)^2/Efd)
  double satA, satB;

  // EXDC1 inputs
  double Ec; // Terminal voltage
  double Vothsg; // Voltage signal from stabilizer
  double Vuel; // Under excitation limiter voltage
  double Voel; // Over excitation limiter voltage
  double Vcomp, Vterm, Vstab;
  
  // Voltage regulator reference
  double Vref;
  
  bool VR_at_min,VR_at_max;
  bool sat_zero; // Saturation function SE is zero
};

// Class for defining events for Exdc1 model
class Exdc1Event
  :public gridpack::math::DAESolver::Event
{
public:

  // Default constructor
  Exdc1Event(Exdc1 *exc):gridpack::math::DAESolver::Event(2),p_exc(exc)
  {
    std:fill(p_term.begin(),p_term.end(),false);

    std::fill(p_dir.begin(),p_dir.end(),gridpack::math::CrossZeroNegative);

  }

  // Destructor
  ~Exdc1Event(void) {}
protected:
  Exdc1 *p_exc;

  void p_update(const double& t, gridpack::ComplexType *state);

  void p_handle(const bool *triggered, const double& t, gridpack::ComplexType *state);
};


#endif
