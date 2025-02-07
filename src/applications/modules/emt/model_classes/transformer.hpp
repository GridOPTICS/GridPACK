/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   transformer.hpp
 * 
 * @brief  Transformer header file
 * 
 */

#ifndef _transformer_h_
#define _transformer_h_

#include <base_branch_model.hpp>
#include <gridpack/include/gridpack.hpp>

class Transformer : public BaseEMTBranchModel
{
public:
  /**
   * Basic constructor
   */
  Transformer();
  
  /**
   * Basic destructor
   */
  ~Transformer();

  /**
     Note: This is a custom version of the load method from the BaseComponent Class. It takes in an extra argument idx to specify which component is being read.
     
     * Load data from DataCollection object into corresponding
     * component. This needs to be implemented by every component
     * @param data data collection associated with component
     */
  void load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx);

  /**
   * Set up branch model before calculation
   */
  void setup();

  /**
   * Initialize branch model before calculation
   * @param [output] values - array where initialized branch variables should be set
   */
  void init(gridpack::RealType *values);
  
  /**
   * Write output from branchs to a string.
   * @param string (output) string with information to be printed out
   * @param bufsize size of string buffer in bytes
   * @param signal an optional character string to signal to this
   * routine what about kind of information to write
   * @return true if bus is contributing string to output, false otherwise
   */
  bool serialWrite(char *string, const int bufsize,
			   const char *signal);
  
  /**
   * Write out branch state
   * @param signal character string used to determine behavior
   * @param string buffer that contains output
   */
  void write(const char* signal, char* string);
  
  /**
     Prestep function
  */
  void preStep(double time ,double timestep);

  /**
     Poststep function
  */
  void postStep(double time);

  /**
   * Return the branch from bus current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
  void getFromBusCurrent(double *ia, double *ib, double *ic);

  /**
   * Return the branch to bus current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
  void getToBusCurrent(double *ia, double *ib, double *ic);

  /**
   * Return the global location for the branch current injection 
   * @param [output] i_gloc - global location for the first current variable
   */
  void getCurrentGlobalLocation(int *i_gloc);

  /**
   * Return the location for the current in the local branch array
   * @param [output] i_gloc - location for the first current variable in the local branch array
   */
  void getCurrentLocalLocation(int *i_loc);

  /**
   * Return the number of variables
   * @param [output] nvar - number of variables
   */
  void getnvar(int *nvar);

  /**
   * Get number of matrix values contributed by branch
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
   * Return vector values from the branch model 
   * @param values - array of returned values
   *
   * Note: This function is used to return the entries in vector,
   * for e.g., the entries in the residual vector from the branch
   * object
   */
  void vectorGetValues(gridpack::RealType *values);

  /**
   * Pass solution vector values to the branch object
   * @param values - array of returned values
   *
   * Note: This function is used to pass the entries in vector
   * to the branch object,
   * for e.g., the state vector values for this branch
   */
  void setValues(gridpack::RealType *values);
  
 protected:
  // Branch parameters
  double R, X;
  double tap;

  double R1,L1;
  double R0,L0;

  double p_R[3][3], p_L[3][3];

  
  bool p_hasResistance;
  bool p_hasInductance;

  // Some temporary arrays
  // From bus and to bus current and their derivatives
  double ibr_from[3],ibr_to[3];
  double dibrf_dt[3],dibrt_dt[3];
};

#endif
