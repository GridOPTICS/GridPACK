/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_gov_model.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   10/02/19
 * 
 * @brief  
 * 
 * 
 */

#ifndef _base_gov_model_h_
#define _base_gov_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"
#include <constants.hpp>

class BaseGovModel : public gridpack::component::BaseComponent
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  BaseGovModel();
  
  /**
   * Basic destructor
   */
  virtual ~BaseGovModel();
  
  /**
   * Initialize governor model before calculation
   * @param [output] values - array where initialized governor variables should be set
   */
  virtual void init(gridpack::ComplexType *values);
  
  /**
   * Write output from governors to a string.
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
   * Write out governor state
   * @param signal character string used to determine behavior
   * @param string buffer that contains output
   */
  virtual void write(const char* signal, char* string);
  
  /**
   * return the bolean indicating whether the governor is ON or OFF
   */
  bool getGovStatus() {return status;}
  
  /**
   * Set bus voltage
   */
  void setVoltage(double busVD, double busVQ) {VD = busVD; VQ = busVQ; }
  
  /**
   * Set TSshift: This parameter is passed by PETSc and is to be used in the Jacobian calculation only.
   */
  void setTSshift(double inshift) {shift = inshift;}
  
  /**
   * Return the governor current injection (in rectangular form) 
   * @param [output] IGD - real part of the governor current
   * @param [output] IGQ - imaginary part of the governor current
   */
  virtual void getCurrent(double *IGD, double *IGQ);
  
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
     * Set the mechanical power parameter inside the governor
     * @param pmech value of the mechanical power 
     */
    virtual void setMechanicalPower(double pmech);

    /**
     * Set the rotor speed deviation parameter inside the governor
     * @param delta_o value of the rotor speed deviation 
     */
    virtual void setRotorSpeedDeviation(double delta_o);

    /** 
     * Get the value of the mechanical power parameter
     * @return value of the mechanical power 
     */
    virtual double getMechanicalPower();

    /** 
     * Get the value of the rotor speed deviation parameter
     * @return value of the rotor speed deviation 
     */
    virtual double getRotorSpeedDeviation();

    /**
     * Set the value of the Vcomp
     * @return value of the Vcomp
     */
    virtual void setVcomp(double vtmp);

    /**
     * Set the value of the time step
     * @return value of the time step
     */
    virtual void setTimestep(double timestep);
 
    /**
     * Set the value of the time increment 
     * @return value of the time increment
     */
    virtual void setTimeincrement(double timeincrement);
 
protected:
  //double        pg; /**< Generator active power output */
  //double        qg; /**< Generator reactive power output */
  //double        mbase; /**< MVA base of the machine */
  //int           status; /**< Machine status */
  //double        sbase;  /** The system MVA base */
  //double        shift; // shift (multiplier) used in the Jacobian calculation.
  //bool          hasGovernor;
  //double        VD, VQ;
  double        VD, VQ;
  int           status; /**< Machine status */
  double        shift; // shift (multiplier) used in the Jacobian calculation.
  double        ts, t_inc;
};

#endif
