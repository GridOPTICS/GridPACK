/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_gov_model.hpp
 * 
 * @brief  Base governor class header file
 * 
 * 
 */

#ifndef _base_gov_model_h_
#define _base_gov_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"
#include <gridpack/applications/modules/emt/constants.hpp>
#include <gridpack/applications/modules/emt/base_classes/base_gen_model.hpp>
#include <gridpack/math/dae_solver.hpp>

class BaseEMTGenModel; // Forward declaration for BaseGenModel

class BaseEMTGovModel : public gridpack::component::BaseComponent
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  BaseEMTGovModel();
  
  /**
   * Basic destructor
   */
  virtual ~BaseEMTGovModel();
  
  /**
     Note: This is a custom version of the load method from the BaseComponent Class. It takes in an extra argument idx to specify which component is being read. Ideally, this method should be moved to the MatVecInterface

   * Load data from DataCollection object into corresponding
   * component. This needs to be implemented by every component
   * @param data data collection associated with component
   */
  virtual void load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx);

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

  /**
     set governor status
  **/
  void setStatus(int gstatus) {status = gstatus;}
  
  /**
   * return the bolean indicating whether the exciter is ON or OFF
   */
  bool getStatus() {return status;}

  /**
   * Write out governor state
   * @param signal character string used to determine behavior
   * @param string buffer that contains output
   */
  virtual void write(const char* signal, char* string);
  
  /**
   * Set bus voltage
   */
  void setVoltage(double busVD, double busVQ) {VD = busVD; VQ = busVQ; }

  /**
   * Copy over voltage from the bus
   */
  void setVoltage(double inva, double invb,double invc) {p_va = inva; p_vb = invb; p_vc = invc;}

  
  /**
   * Set TSshift: This parameter is passed by PETSc and is to be used in the Jacobian calculation only.
   */
  void setTSshift(double inshift) {shift = inshift;}

  /**
   * set current time
   */
  void setTime(double time) {p_time = time; }

  /**
   * return number of variables/states
   */
  void getnvar(int *nvar) {*nvar = nxgov; }

  virtual void setEvent(gridpack::math::DAESolver::EventManagerPtr);
    
  /**
   * Set Jacobian block
   * @param values a 2-d array of Jacobian block for the bus
   */
  virtual bool setJacobian(gridpack::ComplexType **values);

    /**
   * Return vector values from the governor model 
   * @param values - array of returned values
   *
   * Note: This function is used to return the entries in vector,
   * for e.g., the entries in the residual vector from the generator
   * object
   */
  virtual void vectorGetValues(gridpack::ComplexType *values);

  /**
   * Set values in the component based on values in a vector or
   * matrix
   * @param values values in vector or matrix
   */
  virtual void setValues(gridpack::ComplexType *values);

  /**
   * Get number of matrix values contributed by generator
   * @return number of matrix values
   */
  virtual int matrixNumValues();

  /**
 * Return values from a matrix block
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
  virtual void matrixGetValues(int *nvals,gridpack::ComplexType *values,
      int *rows, int *cols);

  /**
   * Set the mechanical power parameter inside the governor
   * @param pmech value of the mechanical power 
   */
  virtual void setInitialMechanicalPower(double pmech);

  /** 
   * Get the value of the mechanical power parameter
   * @return value of the mechanical power 
   */
  virtual double getMechanicalPower();

  /** 
   * Get the value of the mechanical power and its global location
   * @return value of the mechanical power and its global location
   */
  virtual double getMechanicalPower(int *Pmech_gloc);

  
  /**
   * Partial derivatives of Mechanical Power Pmech w.r.t. governor variables
   * @param xgov_loc locations of governor variables
   * @param dPmech_dxgov partial derivatives of mechanical power Pmech w.r.t governor variables
   */
  virtual bool getMechanicalPowerPartialDerivatives(int *xgov_loc,double *dPmech_dxgov);

  /**
   * Set the value of the Vcomp
   * @return value of the Vcomp
   */
  virtual void setVcomp(double vtmp);
  
  void setGenerator(BaseEMTGenModel* generator) {p_gen = generator; };

  BaseEMTGenModel* getGenerator(void) { return p_gen; }

  /**
   * Set the offset for first variable for the governor in the array for all bus variables 
   * @param offset offset
   */
  void setBusOffset(int offset) {offsetb = offset;}

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

  void setBusLocalOffset(int offset) {p_busoffset = offset;}

  /*
    set the location for the first variable in the solution vector
  */
  void setGlobalLocation(int gloc) {p_gloc = gloc;}

  /**
   * return offset in the local vector 
   */
  int getLocalOffset()
  {
    return p_busoffset + offsetb;
  }

  virtual void resetEventFlags(void) {}

protected:
  double        VD, VQ;
  double        p_va, p_vb, p_vc; // Instantaneous bus voltages
  int           status; /**< Machine status */
  double        shift; // shift (multiplier) used in the Jacobian calculation.
  double        p_time; // Current time

  BaseEMTGenModel *p_gen;
  int           offsetb; /**< offset for the first variable for the generator in the array for all bus variables */
  int           p_gloc; /* Location of the first variable in the solution vector */
  int           nxgov; // Number of variables for the model (set by the derived class)
  int           p_busoffset; // Starting location of the bus variables in the local vector */

  double        mbase,sbase; // Machine and system MVA bases
  int           p_nrows;  // number of rows (equations) contributed by this governor
  int           p_ncols;  // number of columns (variables) contributed by this governor
  std::vector<int>   p_rowidx; // global index for rows
  std::vector<int>   p_colidx; // global index for columns
};

#endif
