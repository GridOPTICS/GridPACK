/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_exc_model.hpp
 * 
 * @brief  Base exciter class header file
 * 
 * 
 */

#ifndef _base_exc_model_h_
#define _base_exc_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"
#include <gridpack/applications/modules/emt/constants.hpp>
#include <gridpack/applications/modules/emt/emtutilfunctions.hpp>
#include <gridpack/math/matrix.hpp>
#include <gridpack/applications/modules/emt/base_classes/base_gen_model.hpp>
#include <gridpack/math/dae_solver.hpp>

class BaseEMTGenModel; // Forward declaration for BaseGenModel

class BaseEMTExcModel : public gridpack::component::BaseComponent
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  BaseEMTExcModel();
  
  /**
   * Basic destructor
   */
  virtual ~BaseEMTExcModel();

  /**
     Prestep function
  */
  virtual void preStep(double time ,double timestep) { }
  
  /**
     Poststep function
  */
  virtual void postStep(double time) { }

  /**
    Number of variables
  */ 
  virtual void getnvar(int *nvar) { *nvar = nxexc; }

  /**
     Set the reference power inputs
  **/
  virtual void setPrefQext(double Pref, double Qext) { }

  /**
     Get the current command references
  **/
  virtual void getIpcmdIqcmd(double *Ipcmd, double *Iqcmd) { }

  /**
     Get the power order - used by pitch controller model
  */
  virtual double getPord() { return 0.0; }

  /**
     Set omega_g - from drive train model
  */
  virtual void setOmega(double omega) { }

  /**
     Note: This is a custom version of the load method from the BaseComponent Class. It takes in an extra argument idx to specify which component is being read. Ideally, this method should be moved to the MatVecInterface

   * Load data from DataCollection object into corresponding
   * component. This needs to be implemented by every component
   * @param data data collection associated with component
   */
  virtual void load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx);

  /**
   * Initialize exciter model before calculation
   * @param [output] values - array where initialized exciter variables should be set
   */
  virtual void init(gridpack::RealType *values);
  
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
     set exciter status
  **/
  void setStatus(int estatus) {status = estatus;}
  
  /**
   * return the bolean indicating whether the exciter is ON or OFF
   */
  bool getStatus() {return status;}
  
  /**
   * set current time
   */
  void setTime(double time) {p_time = time; }

  /**
   * Set TSshift: This parameter is passed by PETSc and is to be used in the Jacobian calculation only.
   */
  void setTSshift(double inshift) {shift = inshift;}

  /**
   * Set bus voltage
   */
  void setVoltage(double busVD, double busVQ) {VD = busVD; VQ = busVQ; }

  /**
   * Copy over initial bus voltage from the bus (power flow solution)
   */
  void setInitialVoltage(double inVm,double inVa) {p_Vm0 = inVm; p_Va0 = inVa;}

  /**
   * Get the initial field voltage (at t = tstart) for the exciter
   * @param fldv value of the field voltage
   */
  double getInitialFieldVoltage();
  
  virtual void setEvent(gridpack::math::RealDAESolver::EventManagerPtr);

  /**
   * Set the offset for first variable for the exciter in the array for all bus variables 
   * @param offset offset
   */
  void setBusOffset(int offset) {offsetb = offset;}


  /**
   * Get number of matrix values contributed by generator
   * @return number of matrix values
   */
  virtual int matrixNumValues();

  /**
   * Return values from Jacobian matrix
   * @param nvals: number of values to be inserted
   * @param values: pointer to matrix block values
   * @param rows: pointer to matrix block rows
   * @param cols: pointer to matrix block cols
   */
  virtual void matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols);


  /**
   * Return vector values from the generator model 
   * @param values - array of returned values
   *
   * Note: This function is used to return the entries in vector,
   * for e.g., the entries in the residual vector from the generator
   * object
   */
  virtual void vectorGetValues(gridpack::RealType *values);

  /**
   * Pass solution vector values to the generator object
   * @param values - array of returned values
   *
   * Note: This function is used to pass the entries in vector
   * to the generator object,
   * for e.g., the state vector values for this generator
   */
  virtual void setValues(gridpack::RealType *values);
  
  /** 
   * Get the value of the field voltage parameter
   * @return value of field voltage
   */
  virtual double getFieldVoltage();

  /**
   * Get the current command references during initialization
   */
  void getInitialIpcmdIqcmd(double *Ipcmd0, double *Iqcmd0);

  
  /** 
   * Get the value of the field voltage parameter
   * and return the global location of field voltage variable
   * @return value of field voltage
   */
  virtual double getFieldVoltage(int *Efd_gloc);

  /**
   * Partial derivatives of field voltage Efd w.r.t. exciter variables
   * @param xexc_loc locations of exciter variables
   * @param dEfd_dxexc partial derivatives of field voltage w.r.t exciter variables
   * @param dEfd_dxgen partial derivatives of field voltage w.r.t. generator variables
   */
  virtual bool getFieldVoltagePartialDerivatives(int *xexc_loc,double *dEfd_dxexc,double *dEfd_dxgen);

  /**
   * Set the generator associated with this exciter
   **/
  void setGenerator(BaseEMTGenModel* generator) {p_gen = generator; }

  /**
   * Get the generator associated with this exciter
   */
  BaseEMTGenModel* getGenerator() { return p_gen; }

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

  /*
    set the global location for first voltage variable for this bus
  */
  void setVoltageGlobalLocation(int glocvoltage) {p_glocvoltage = glocvoltage;}


  /**
   * return offset in the local vector 
   */
  int getLocalOffset()
  {
    return p_busoffset + offsetb;
  }

  virtual void resetEventFlags(void) {}

  /**
   * Copy over voltage from the bus
   */
  void setVoltage(double inva, double invb,double invc) {p_va = inva; p_vb = invb; p_vc = invc;}

  /**
   * Set Machine angle
   */
  void setMachineAngle(double deltain) {p_delta = deltain; }

  /**
   * Set type of integration algorithm
   */
  void setIntegrationType(EMTMachineIntegrationType type) {integrationtype = type; }

protected:
  double        VD, VQ;
  int           status; /**< Exciter status */
  double        p_time;   /** Current time */
  double        p_Vm0, p_Va0; /** Initial voltage magnitude and angle **/
  double        p_va,p_vb,p_vc; /** Bus voltage **/
  double        p_delta;   /** Machine angle */
  double        shift; // shift (multiplier) used in the Jacobian calculation
  
  EMTMachineIntegrationType integrationtype; // Integration type 

  BaseEMTGenModel* p_gen; // Generator model
  int           offsetb; /**< offset for the first variable for the generator in the array for all bus variables */
  int           p_gloc; // Global location of the first variable for the generator
  int           p_glocvoltage; // Global location for the first voltage variable for the bus

  int           nxexc;    /** Number of variables for the exciter model */
  int           p_busoffset; /** Starting location for bus variables in the local state vector */
  int           p_nrows;  // number of rows (equations) contributed by this excitor
  int           p_ncols;  // number of columns (variables) contributed by this excitor
  std::vector<int>   p_rowidx; // global index for rows
  std::vector<int>   p_colidx; // global index for columns

};

#endif
