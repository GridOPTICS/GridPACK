/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   gencls.hpp
 * 
 * @brief Classical generator model 
 * 
 * 
 */

#ifndef _gencls_model_h_
#define _gencls_model_h_

#include <base_gen_model.hpp>
#include <gridpack/include/gridpack.hpp>


class Gencls: public BaseGenModel
{
   public:
  /**
     * Basic constructor
     */
    Gencls();

    /**
     * Basic destructor
     */
    ~Gencls();

    /**
     * Load parameters from DataCollection object into generator model
     * @param data collection of generator parameters from input files
     * @param index of generator on bus
     * TODO: might want to move this functionality to BaseGeneratorModel
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);


    /**
     *  Set Jacobian values
     *  @param values a 2-d array of Jacobian block for the bus
     */
    bool setJacobian(gridpack::ComplexType **values);

    /**
     * Initialize generator model before calculation
     * @param [output] values - array where initialized generator variables should be set
     */
  void init(gridpack::ComplexType *values);

    /**
     * Write output from generators to a string.
     * @param string (output) string with information to be printed out
     * @param bufsize size of string buffer in bytes
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if bus is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const int bufsize,
        const char *signal);

    /**
     * Write out generator state
     * @param signal character string used to determine behavior
     * @param string buffer that contains output
     */
    void write(const char* signal, char* string);

  /**
   * Return the generator current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
  void getCurrent(double *ia, double *ib, double *ic);

  /**
   * Return the number of variables
   * @param [output] nvar - number of variables
   */
  void getnvar(double *nvar);
  
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
  void matrixGetValues(int *nvals,gridpack::ComplexType *values,
      int *rows, int *cols);

  
  private:
    // Machine parameters
    double p_Rs; // Machine stator resistance
    double p_L;  // Machine transient inductance
    double p_H;      // Machine Inertia constant
    double p_D;      // Machine damping coefficient

    // Internal constants
    double p_Pm;  // Mechanical power input
    double p_Ep;  // Internal emf

    // Generator variables and their derivatives
    double p_delta,p_dw;
    double p_deltadot,p_dwdot;

    // Previous step values of the variables
    double p_deltaprev, p_dwprev;

    int bid;
};

#endif
