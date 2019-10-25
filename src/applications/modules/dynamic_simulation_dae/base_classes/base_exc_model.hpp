/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_exc_model.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   10/15/19
 * 
 * @brief  
 * 
 * 
 */

#ifndef _base_exc_model_h_
#define _base_exc_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"
#include <constants.hpp>

class BaseExcModel : public gridpack::component::BaseComponent
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  BaseExcModel();
  
  /**
   * Basic destructor
   */
  virtual ~BaseExcModel();
  
  /**
   * Initialize exciter model before calculation
   * @param [output] values - array where initialized exciter variables should be set
   */
  virtual void init(gridpack::ComplexType *values);
  
  /**
   * Write output from exciters to a string.
   * @param string (output) string with information to be printed out
   * @param bufsize size of string buffer in bytes
   * @param signal an optional character string to signal to this
   * routine what about kind of information to write
   * @return true if bus is contributing string to output, false otherwise
   */
  virtual bool serialWrite(char *string, const int bufsize,
			   const char *signal);
  
  /**
   * Write out exciter state
   * @param signal character string used to determine behavior
   * @param string buffer that contains output
   */
  virtual void write(const char* signal, char* string);
  
  /**
   * return the bolean indicating whether the exciter is ON or OFF
   */
  bool getExcStatus() {return status;}
  
  /**
   * Set bus voltage
   */
  void setVoltage(double busVD, double busVQ) {VD = busVD; VQ = busVQ; }
  
  /**
   * Set TSshift: This parameter is passed by PETSc and is to be used in the Jacobian calculation only.
   */
  void setTSshift(double inshift) {shift = inshift;}
  
  /**
   * Return the matrix entries
   * @param [output] nval - number of values set
   * @param [output] row - row indices for matrix entries
   * @param [output] col - col indices for matrix entries
   * @param [output] values - matrix entries
   * return true when matrix entries set

   Note: It becomes cumbersome to use matrixDiagValues for large matrix blocks. matrixDiagEntries makes this easier by also allowing to set the row and col indices for the entries. Eventually, we would want to move this function to the MatVecInterface and then have the loadXXXData methods in the mapper class call it for setting the matrix values 
   */

  virtual bool matrixDiagEntries(int *nval,int *row, int *col, gridpack::ComplexType *values);

  /**
     Note: This is a custom version of the load method from the BaseComponent Class. It takes in an extra argument idx to specify which component is being read. Ideally, this method should be moved to the MatVecInterface

   * Load data from DataCollection object into corresponding
   * component. This needs to be implemented by every component
   * @param data data collection associated with component
   */
  virtual void load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx);
  
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

  /***************************************/

  /**
   * Set the initial field voltage (at t = tstart) for the exciter
   * @param fldv value of the field voltage
   */
  virtual void setInitialFieldVoltage(double fldv);

  /**
   * Set the field current parameter inside the exciter
   * @param fldc value of the field current
   */
  virtual void setFieldCurrent(double fldc);

  /** 
   * Get the value of the field voltage parameter
   * @return value of field voltage
   */
  virtual double getFieldVoltage();

  /** 
   * Get the value of the field current parameter
   * @return value of field current
   */
  virtual double getFieldCurrent();
  
  /**
   * Set initial time step
   */
  void setInitialTimeStep(double timestep);
 
protected:
  double        VD, VQ;
  int           status; /**< Exciter status */
  double        shift; // shift (multiplier) used in the Jacobian calculation.i
  double        dt0;
};

#endif
