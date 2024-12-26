/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   epria1.hpp
 * 
 * @brief EPRIA1 model wrapper
 * 
 * 
 */

#ifndef _epria1_model_h_
#define _epria1_model_h_

#include <base_gen_model.hpp>
#include <gridpack/include/gridpack.hpp>

extern "C" {
#include "GFLIBR.h"
}


class Epria1: public BaseEMTGenModel
{
   public:
  /**
     * Basic constructor
     */
    Epria1();

    /**
     * Basic destructor
     */
    ~Epria1();

    /**
     * Load parameters from DataCollection object into generator model
     * @param data collection of generator parameters from input files
     * @param index of generator on bus
     * TODO: might want to move this functionality to BaseGeneratorModel
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);


  /**
   * return number of variables
   */
  void getnvar(int *nvar);

  /**
     Prestep function
  */
  void preStep(double time ,double timestep);

  /**
     Poststep function
  */
  void postStep(double time);


    /**
     * Initialize generator model before calculation
     * @param [output] values - array where initialized generator variables should be set
     */
  void init(gridpack::RealType *values);

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
   * Return the global location for the generator current injection 
   * @param [output] i_gloc - global location for the first current variable
   */
  void getCurrentGlobalLocation(int *i_gloc);

  /**
   * Get number of matrix values contributed by generator
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
   * Return vector values from the generator model 
   * @param values - array of returned values
   *
   * Note: This function is used to return the entries in vector,
   * for e.g., the entries in the residual vector from the generator
   * object
   */
  void vectorGetValues(gridpack::RealType *values);

  /**
   * Pass solution vector values to the generator object
   * @param values - array of returned values
   *
   * Note: This function is used to pass the entries in vector
   * to the generator object,
   * for e.g., the state vector values for this generator
   */
  void setValues(gridpack::RealType *values);

  private:
    // Machine parameters
    double p_Rs; // Machine stator resistance
    double p_Xdp;  // Machine transient reactance
    double p_L;  // Machine transient inductance
    double p_H;      // Machine Inertia constant
    double p_D;      // Machine damping coefficient

    double VD, VQ; // Dq axis voltages

    // Internal constants
    double Ed,Eq;  // Internal emf
    double phi_PLL;

    // Some temporary arrays
  double p_vdq0[3]; // voltage in dq0 reference frame
  double p_idq0[3]; // current in dq0 reference frame
  double p_vabc[3]; // voltage in abc reference frame
  double p_eabc[3]; // Internal voltage output by model
    // Generator variables and their derivatives
  double p_iabc[3],p_iLabc[3];
  double p_iLdot[3];

  IEEE_Cigre_DLLInterface_Instance model;
  
  // Parameters
  MyModelParameters modelparams;
  // Inputs
  MyModelInputs     modelinputs;
  // Outputs
  MyModelOutputs    modeloutputs;
  // states
  double            modelstates[15];

  gridpack::ComplexType Zsource;

};

#endif
