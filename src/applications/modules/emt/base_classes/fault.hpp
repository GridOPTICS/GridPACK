/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   fault.hpp
 * 
 * @brief Three-phase fault Model 
 * 
 * 
 */

#ifndef _fault_h_
#define _fault_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include <gridpack/include/gridpack.hpp>
#include "gridpack/component/base_component.hpp"
#include <gridpack/applications/modules/emt/constants.hpp>
#include <gridpack/applications/modules/emt/emtutilfunctions.hpp>
#include <gridpack/applications/modules/emt/emtnetwork.hpp>
#include <gridpack/math/matrix.hpp>

class EmtBus; // Forward declaration to resolve circular dependency

class Fault: public gridpack::component::BaseComponent
{
   public:
  /**
     * Basic constructor
     */
    Fault();

    /**
     * Basic destructor
     */
    ~Fault();

    /**
     * Load parameters from DataCollection object into load model
     * @param data collection of load parameters from input files
     * @param index of load on bus
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);


  /**
     Set fault parameters

    @param [in] ton - fault on time
    @param [in] toff - fault off time
    @param [in] type - fault type (SLG, LLG, or 3Phase)
    @param [in] phases - faulted phases
    @param [in] Ron - fault resistance
    @param [in] Rgnd - ground resistance
  */
  void setparams(EmtBus *emtbus,double ton, double toff, std::string type, std::string phases, double Ron, double Rgnd);

    /**
     *  Set Jacobian values
     *  @param values a 2-d array of Jacobian block for the bus
     */
    bool setJacobian(gridpack::ComplexType **values);


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
     * Write out fault current state
     * @param signal character string used to determine behavior
     * @param string buffer that contains output
     */
    void write(const char* signal, char* string);

  /**
   * Return the fault current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
  void getCurrent(double *ia, double *ib, double *ic);

  /**
   * Return the global location for the fault current 
   * @param [output] i_gloc - global location for the first current variable
   */
  void getCurrentGlobalLocation(int *i_gloc);

  /**
   * Get number of matrix values contributed by fault
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
  void matrixGetValues(int *nvals, gridpack::ComplexType *values, int *rows, int *cols);

  /**
   * Return vector values from the fault model 
   * @param values - array of returned values
   *
   * Note: This function is used to return the entries in vector,
   * for e.g., the entries in the residual vector from the fault
   * object
   */
  void vectorGetValues(gridpack::ComplexType *values);

  /**
   * Pass solution vector values to the fault object
   * @param values - array of returned values
   *
   * Note: This function is used to pass the entries in vector
   * to the fault object,
   * for e.g., the state vector values for this fault
   */
  void setValues(gridpack::ComplexType *values);

  private:
    EmtBus *bus;
  
    // Fault parameters
    double Ron; // Fault resistance
    double Rgnd;  // Ground resistance
    int    faulttype; // 1 for SLG, 2 for LL, 3 for 3 phase
    int    faultedphases[3]; // Which phases are faulted (1 for faulted, zero for unfaulted
    double ton, toff; // fault on and fault off times
    bool   faulton;
    int    p_gloc; // Global location
    int    nxfault;
    int    offsetb; 

    // Fault variables
    double ifault[3]; // Fault current
};

#endif
