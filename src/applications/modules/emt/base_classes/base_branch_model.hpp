/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_branch_model.hpp
 * 
 * @brief  Base branch class header file
 * 
 */

#ifndef _base_branch_model_h_
#define _base_branch_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"
#include <gridpack/applications/modules/emt/constants.hpp>
#include <gridpack/applications/modules/emt/emtutilfunctions.hpp>
#include <gridpack/applications/modules/emt/emtnetwork.hpp>
#include <gridpack/math/matrix.hpp>

class EmtBus; // Forward declaration for EmtBus

class BaseEMTBranchModel : public gridpack::component::BaseComponent
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  BaseEMTBranchModel();
  
  /**
   * Basic destructor
   */
  virtual ~BaseEMTBranchModel();

  /**
     Note: This is a custom version of the load method from the BaseComponent Class. It takes in an extra argument idx to specify which component is being read.
     
     * Load data from DataCollection object into corresponding
     * component. This needs to be implemented by every component
     * @param data data collection associated with component
     */
  virtual void load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx);

  /**
   * Set up branch model before calculation
   */
  virtual void setup() {}
  
  /**
   * Initialize branch model before calculation
   * @param [output] values - array where initialized branch variables should be set
   */
  virtual void init(gridpack::RealType *values);
  
  /**
   * Write output from branchs to a string.
   * @param string (output) string with information to be printed out
   * @param bufsize size of string buffer in bytes
   * @param signal an optional character string to signal to this
   * routine what about kind of information to write
   * @return true if bus is contributing string to output, false otherwise
   */
  virtual bool serialWrite(char *string, const int bufsize,
			   const char *signal);
  
  /**
   * Write out branch state
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
     Prestep function
  */
  virtual void preStep(double time ,double timestep);

  /**
     Poststep function
  */
  virtual void postStep(double time);

  /*
    Set the bus voltage global location
  */
  void setVoltageGlobalLocation(int v_gloc) { p_glocvoltage = v_gloc; }

  /**
   * Set TSshift: This parameter is passed by PETSc and is to be used in the Jacobian calculation only.
   */
  void setTSshift(double inshift) {shift = inshift;}
  
  /**
   * Return the branch current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
  virtual void getCurrent(double *ia, double *ib, double *ic);


  /**
   * Return the global location for the branch current injection 
   * @param [output] i_gloc - global location for the first current variable
   */
  virtual void getCurrentGlobalLocation(int *i_gloc);

  /**
   * Return the location for the current in the local branch array
   * @param [output] i_gloc - location for the first current variable in the local branch array
   */
  virtual void getCurrentLocalLocation(int *i_loc) { *i_loc = 0; }


  /**
   * Return the number of variables
   * @param [output] nvar - number of variables
   */
  virtual void getnvar(int *nvar) {*nvar = nxbranch;}

  /**
   * Get number of matrix values contributed by branch
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
   * Return vector values from the branch model 
   * @param values - array of returned values
   *
   * Note: This function is used to return the entries in vector,
   * for e.g., the entries in the residual vector from the branch
   * object
   */
  virtual void vectorGetValues(gridpack::RealType *values);

  /**
   * Pass solution vector values to the branch object
   * @param values - array of returned values
   *
   * Note: This function is used to pass the entries in vector
   * to the branch object,
   * for e.g., the state vector values for this branch
   */
  virtual void setValues(gridpack::RealType *values);
  
  /**
   * Set the offset for first variable for the branch in the array for all bus variables 
   * @param offset offset
   */
  void setBranchOffset(int offset) {offsetb = offset;}

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

  void setBranchLocalOffset(int offset) {p_branchoffset = offset;}

  /*
    set the location for the first variable in the solution vector
  */
  void setGlobalLocation(int gloc) {p_gloc = gloc;}

  /**
   * return offset in the local vector 
   */
  int getLocalOffset()
  {
    return p_branchoffset + offsetb;
  }

  /**
   * Set the type of machine integration algorithm
   * @param [in] type - the integration type
   */
  void setIntegrationType(EMTMachineIntegrationType type) {integrationtype = type;}

  /**
   * Set the from bus associated with this branch
   * @param [in] frombus - from bus
   */
  void setFromBus(EmtBus *frombus) {fbus = frombus;}

    /**
   * Set the to bus associated with this branch
   * @param [in] tobus - to bus
   */
  void setToBus(EmtBus *tobus) {tbus = tobus;}

 protected:
  int           status; /**< Branch status */
  std::string   cktid; // circuit id
  double        sbase;  /** The system MVA base */
  double        p_time; /** Current time */
  double        shift; // shift (multiplier) used in the Jacobian calculation.

  EMTMachineIntegrationType integrationtype;

  int           offsetb; /**< offset for the first variable for the branch in the array for all branch variables */
  int           p_gloc; // Global location of the first variable for the branch
  
  int           nxbranch; /* Number of variables for the branch model */
  int           p_branchoffset; /** Offset for the bus variables in the local vector. Used only for events */
  int           p_glocvoltage; /* Global location of the first bus voltage variable. This is set by the bus */

  EmtBus *fbus; // From bus
  EmtBus *tbus; // To bus

  std::vector<int>   p_rowidx; // global index for rows
  std::vector<int>   p_colidx; // global index for columns
};

#endif
