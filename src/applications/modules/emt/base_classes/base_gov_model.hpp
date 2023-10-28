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
#include <constants.hpp>
#include <base_gen_model.hpp>
#include <gridpack/math/dae_solver.hpp>

class BaseEMTGenModel; // Forward declaration for BaseGenModel

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
  virtual void init(gridpack::RealType *values);
  
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
   * Set TSshift: This parameter is passed by PETSc and is to be used in the Jacobian calculation only.
   */
  void setTSshift(double inshift) {shift = inshift;}

  virtual void setEvent(gridpack::math::RealDAESolver::EventManagerPtr);
  
  /**
     Note: This is a custom version of the load method from the BaseComponent Class. It takes in an extra argument idx to specify which component is being read. Ideally, this method should be moved to the MatVecInterface

   * Load data from DataCollection object into corresponding
   * component. This needs to be implemented by every component
   * @param data data collection associated with component
   */
  virtual void load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx);
  
  /**
   * Set the number of rows contributed by this generator
   * @param nrows number of rows
   */
  virtual void matrixSetNumRows(int nrows);

  /**
   * Set the number of columns contributed by this generator
   * @param ncols number of columns
   */
  virtual void matrixSetNumCols(int ncols);

  /**
   * Number of rows (equations) contributed to by generator
   * @return number of rows
   */
  int matrixNumRows();

  /**
   * Number of rows (equations) contributed to by generator
   * @return number of rows
   */
  int matrixNumCols();

  /** 
   * Set global row index
   * @param irow local row index
   * @param global row index
   */
  void matrixSetRowIndex(int irow, int idx);

  /** 
   * Set global column index
   * @param icol local column index
   * @param global column index
   */
  void matrixSetColIndex(int icol, int idx);

  /**
   * Return global row index given local row index
   * @param irow local row index
   * @return global row index
   */
  int matrixGetRowIndex(int irow);

  /**
   * Return global column index given local column index
   * @param icol local column index
   * @return global column index
   */
  virtual int matrixGetColIndex(int icol);

  /**
   * Get number of matrix values contributed by generator
   * @return number of matrix values
   */
  int matrixNumValues();

  /**
 * Return values from a matrix block
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
  void matrixGetValues(int *nvals,gridpack::RealType *values,
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
  
  void setGenerator(BaseEMTGenModel* generator);

  BaseEMTGenModel* getGenerator(void);

  /**
   * Set the offset for first variable for the governor in the array for all bus variables 
   * @param offset offset
   */
  void setBusOffset(int offset) {offsetb = offset;}

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
  bool vectorValues(gridpack::RealType *values);

  /**
   * Set values in the component based on values in a vector or
   * matrix
   * @param values values in vector or matrix
   */
  void setValues(gridpack::RealType *values);

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
  int           status; /**< Machine status */
  double        shift; // shift (multiplier) used in the Jacobian calculation.

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
