/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   constantimpedance.hpp
 * 
 * @brief Constant Impedance Load Model 
 * 
 * 
 */

#ifndef _constantimpedance_model_h_
#define _constantimpedance_model_h_

#include <base_load_model.hpp>
#include <gridpack/include/gridpack.hpp>


class Constantimpedance: public BaseEMTLoadModel
{
   public:
  /**
     * Basic constructor
     */
    Constantimpedance();

    /**
     * Basic destructor
     */
    ~Constantimpedance();

    /**
     * Load parameters from DataCollection object into load model
     * @param data collection of load parameters from input files
     * @param index of load on bus
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);


    /**
     *  Set Jacobian values
     *  @param values a 2-d array of Jacobian block for the bus
     */
    bool setJacobian(gridpack::ComplexType **values);

    /**
     * Initialize load model before calculation
     * @param [output] values - array where initialized load variables should be set
     */
  void init(gridpack::ComplexType *values);

    /**
     * Write output from loads to a string.
     * @param string (output) string with information to be printed out
     * @param bufsize size of string buffer in bytes
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if bus is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const int bufsize,
        const char *signal);

    /**
     * Write out load state
     * @param signal character string used to determine behavior
     * @param string buffer that contains output
     */
    void write(const char* signal, char* string);

  /**
   * Return the load current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
  void getCurrent(double *ia, double *ib, double *ic);

  /**
   * Return the global location for the load current 
   * @param [output] i_gloc - global location for the first current variable
   */
  void getCurrentGlobalLocation(int *i_gloc);

  /**
   * Get number of matrix values contributed by load
   * @return number of matrix values
   */
  int matrixNumValues();

  /**
   * Return values from a matrix block
   * @param matrix - the Jacobian matrix
   */
  void matrixGetValues(gridpack::math::Matrix &matrix);

  /**
   * Return vector values from the load model 
   * @param values - array of returned values
   *
   * Note: This function is used to return the entries in vector,
   * for e.g., the entries in the residual vector from the load
   * object
   */
  void vectorGetValues(gridpack::ComplexType *values);

  /**
   * Pass solution vector values to the load object
   * @param values - array of returned values
   *
   * Note: This function is used to pass the entries in vector
   * to the load object,
   * for e.g., the state vector values for this load
   */
  void setValues(gridpack::ComplexType *values);

  private:
    // Load parameters
    double p_R[3]; // Resistance part of constant impedance 
    double p_L[3];  // Inductance part of constant impedance

    // Load variables
  double p_i[3]; // Inductor current
  double p_idot[3]; // Current derivative
};

#endif
