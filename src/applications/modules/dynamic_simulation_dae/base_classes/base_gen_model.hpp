/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_gen_model.hpp
 * @author Shrirang Abhyankar
 * @Last modified:   02/17/19
 * 
 * @brief  
 * 
 * 
 */

#ifndef _base_gen_model_h_
#define _base_gen_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"
#include <constants.hpp>
#include <base_exc_model.hpp>
#include <base_gov_model.hpp>

class BaseExcModel; // Forward declaration for BaseExcModel
class BaseGovModel; // Forward declaration for BaseGovModel

class BaseGenModel : public gridpack::component::BaseComponent
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  BaseGenModel();
  
  /**
   * Basic destructor
   */
  virtual ~BaseGenModel();
  
  /**
   * Initialize generator model before calculation
   * @param [output] values - array where initialized generator variables should be set
   */
  virtual void init(gridpack::ComplexType *values);
  
  /**
   * Write output from generators to a string.
   * @param string (output) string with information to be printed out
   * @param bufsize size of string buffer in bytes
   * @param signal an optional character string to signal to this
   * routine what about kind of information to write
   * @return true if bus is contributing string to output, false otherwise
   */
  virtual bool serialWrite(char *string, const int bufsize,
			   const char *signal);
  
  virtual double getAngle();
  
  /**
   * Write out generator state
   * @param signal character string used to determine behavior
   * @param string buffer that contains output
   */
  virtual void write(const char* signal, char* string);
  
  /**
   * return the bolean indicating whether the gen is ON or OFF
   */
  bool getGenStatus() {return status;}
  
  /**
   * Set bus voltage
   */
  void setVoltage(double busVD, double busVQ) {VD = busVD; VQ = busVQ; }
  
  /**
   * Set TSshift: This parameter is passed by PETSc and is to be used in the Jacobian calculation only.
   */
  void setTSshift(double inshift) {shift = inshift;}
  
  /**
   * Return the generator current injection (in rectangular form) 
   * @param [output] IGD - real part of the generator current
   * @param [output] IGQ - imaginary part of the generator current
   */
  virtual void getCurrent(double *IGD, double *IGQ);
  
  /**
     Note: This is a custom version of the load method from the BaseComponent Class. It takes in an extra argument idx to specify which component is being read. Ideally, this method should be moved to the MatVecInterface

   * Load data from DataCollection object into corresponding
   * component. This needs to be implemented by every component
   * @param data data collection associated with component
   */
  virtual void load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx);


/**
 * Set Jacobian block
 * @param values a 2-d array of Jacobian block for the bus
 */
  virtual bool setJacobian(gridpack::ComplexType **values);

  /****************************************************
 The following methods are inherited from the BaseComponent class and are 
to be overwritten by the implementation */
  
  
  /**
   * Set an internal variable that can be used to control the behavior of the
   * component. This function doesn't need to be implemented, but if needed,
   * it can be used to change the behavior of the component in different phases
   * of the calculation. For example, if a different matrix needs to be
   * generated at different times, the mode of the calculation can changed to
   * get different values from the MatVecInterface functions
   * @param mode integer indicating which mode should be used
   */
  void setMode(int mode) { p_mode = mode;}

  /**
   * Return size of vector block contributed by component
   * @param isize number of vector elements
   * @return false if network component does not contribute
   *        vector element
   */
  bool vectorSize(int *isize) const;

  /**
   * Return the values of the vector block
   * @param values pointer to vector values
   * @return false if network component does not contribute
   *        vector element
   */
  bool vectorValues(gridpack::ComplexType *values);

  /**
   * Set values in the component based on values in a vector or
   * matrix
   * @param values values in vector or matrix
   */
  void setValues(gridpack::ComplexType *values);

  /**
   * Set the field current parameter inside the exciter
   * @param fldc value of the field current
   */
  virtual double getFieldCurrent(void);

  /**
   * Return the rotor speed deviation
   * @param rotor speed deviation 
   */
  virtual double getRotorSpeedDeviation();

  /**
   * Set the offset for first variable for the generator in the array for all bus variables 
   * @param offset offset
   */
  void setBusOffset(int offset) {offsetb = offset;}

  /***************************************/


  void setExciter(boost::shared_ptr<BaseExcModel> &p_exciter);

  boost::shared_ptr<BaseExcModel> getExciter();
  
  bool getphasExciter();

  void setGovernor(boost::shared_ptr<BaseGovModel> &p_governor);

  boost::shared_ptr<BaseGovModel> getGovernor();
  
  bool getphasGovernor();

 protected:
  double        pg; /**< Generator active power output */
  double        qg; /**< Generator reactive power output */
  double        mbase; /**< MVA base of the machine */
  int           status; /**< Machine status */
  double        sbase;  /** The system MVA base */
  double        shift; // shift (multiplier) used in the Jacobian calculation.
  double        VD, VQ;
  bool          p_hasExciter;
  bool          p_hasGovernor;
  boost::shared_ptr<BaseExcModel> p_exciter; //SJin: new pointer to exciter
  boost::shared_ptr<BaseGovModel> p_governor; //SJin: new pointer to governor
  int           offsetb; /**< offset for the first variable for the generator in the array for all bus variables */
};

#endif
