/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_gen_model.hpp
 * 
 * @brief  Base generator class header file
 * 
 */

#ifndef _base_emt_gen_model_h_
#define _base_emt_gen_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"
#include <constants.hpp>
#include <emtutilfunctions.hpp>
#include <gridpack/math/matrix.hpp>
#include <base_exc_model.hpp>
#include <base_gov_model.hpp>

class BaseExcModel; // Forward declaration for BaseExcModel
class BaseGovModel; // Forward declaration for BaseGovModel

class BaseEMTGenModel : public gridpack::component::BaseComponent
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  BaseEMTGenModel();
  
  /**
   * Basic destructor
   */
  virtual ~BaseEMTGenModel();

  /**
     Note: This is a custom version of the load method from the BaseComponent Class. It takes in an extra argument idx to specify which component is being read.
     
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
  
  /**
   * Write out generator state
   * @param signal character string used to determine behavior
   * @param string buffer that contains output
   */
  virtual void write(const char* signal, char* string);
  
  /**
   * return the bolean indicating whether the gen is ON or OFF
   */
  int getStatus() {return status;}

  /**
   * set the bolean indicating whether the gen is ON or OFF
   */
  void setStatus(int gstatus) {status = gstatus;}

  /**
   * set current time
   */
  void setTime(double time) {p_time = time; }
  
  /**
   * Copy over voltage from the bus
   */
  void setVoltage(double inva, double invb,double invc) {p_va = inva; p_vb = invb; p_vc = invc;}

  /*
    Set the bus voltage global location
  */
  void setVoltageGlobalLocation(int v_gloc) { p_glocvoltage = v_gloc; }

  /**
   * Copy over initial bus voltage from the bus (power flow solution)
   */
  void setInitialVoltage(double inVm,double inVa) {p_Vm0 = inVm; p_Va0 = inVa;}

  
  /**
   * Set TSshift: This parameter is passed by PETSc and is to be used in the Jacobian calculation only.
   */
  void setTSshift(double inshift) {shift = inshift;}
  
  /**
   * Return the generator current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
  virtual void getCurrent(double *ia, double *ib, double *ic);

  /**
   * Return the global location for the generator current injection 
   * @param [output] i_gloc - global location for the first current variable
   */
  virtual void getCurrentGlobalLocation(int *i_gloc);

  
  /**
   * Return the number of variables
   * @param [output] nvar - number of variables
   */
  void getnvar(int *nvar) {*nvar = nxgen;}

  /**
   * Get number of matrix values contributed by generator
   * @return number of matrix values
   */
  virtual int matrixNumValues();

  /**
   * Return values from a matrix block
   * @param matrix - the Jacobian matrix
   */
  virtual void matrixGetValues(gridpack::math::Matrix &matrix);

  /**
   * Return vector values from the generator model 
   * @param values - array of returned values
   *
   * Note: This function is used to return the entries in vector,
   * for e.g., the entries in the residual vector from the generator
   * object
   */
  virtual void vectorGetValues(gridpack::ComplexType *values);

  /**
   * Pass solution vector values to the generator object
   * @param values - array of returned values
   *
   * Note: This function is used to pass the entries in vector
   * to the generator object,
   * for e.g., the state vector values for this generator
   */
  virtual void setValues(gridpack::ComplexType *values);
  
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
   * Return the location of speed rotor speed deviation variable in the bus array
   * @param rotor speed deviation location
   */
  virtual int getRotorSpeedDeviationLocation();

  /**
   * Set the offset for first variable for the generator in the array for all bus variables 
   * @param offset offset
   */
  void setBusOffset(int offset) {offsetb = offset;}

  void setExciter(boost::shared_ptr<BaseExcModel> &p_exciter);

  boost::shared_ptr<BaseExcModel> getExciter();
  
  bool hasExciter();

  void setGovernor(boost::shared_ptr<BaseGovModel> &p_governor);

  boost::shared_ptr<BaseGovModel> getGovernor();
  
  bool hasGovernor();

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
  double        pg; /**< Generator active power output */
  double        qg; /**< Generator reactive power output */
  double        mbase; /**< MVA base of the machine */
  int           busnum; /**< Bus number */
  int           status; /**< Machine status */
  double        sbase;  /** The system MVA base */
  double        p_time; /** Current time */
  double        shift; // shift (multiplier) used in the Jacobian calculation.
  double        p_Vm0,p_Va0; // Initial bus voltage and angle
  double        p_va, p_vb, p_vc; // phase voltages
  bool          p_hasExciter; // Flag indicating whether this generator has exciter
  bool          p_hasGovernor; // Flag indicating whether this generator has governor
  boost::shared_ptr<BaseExcModel> p_exciter; // Exciter
  boost::shared_ptr<BaseGovModel> p_governor; // Governor
  int           offsetb; /**< offset for the first variable for the generator in the array for all bus variables */
  int           p_gloc; // Global location of the first variable for the generator
  
  int           nxgen; /* Number of variables for the generator model */
  int           p_busoffset; /** Offset for the bus variables in the local vector. Used only for events */
  int           p_glocvoltage; /* Global location of the first bus voltage variable. This is set by the bus */

  // Arrays used in coupling blocks between generator and exciter. These should be allocated and destroyed by the derived class
  int           *xexc_loc;   // locations for exciter variables in the bus variables array
  double        *dEfd_dxexc; // Partial derivatives of field voltage Efd w.r.t. exciter variables (size = nxexc)
  double        *dEfd_dxgen; // Partial derivatives of field voltage Efd w.r.t. generator variables (size = nxgen)

  // Arrays used in coupling blocks between generator and governor. These should be allocated and destroyed by the derived class
  int           *xgov_loc;   // locations for governor variables in the bus variables array
  double        *dPmech_dxgov; // Partial derivatives of mechanical power Pmech w.r.t. governor variables (size = nxgov)

  std::vector<int>   p_rowidx; // global index for rows
  std::vector<int>   p_colidx; // global index for columns
};

#endif
