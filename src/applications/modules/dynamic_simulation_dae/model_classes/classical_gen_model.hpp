/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   classical_generator_model.hpp
 * @author Shrirang Abhyankar
 * @Last modified:   02/17/19
 * 
 * @brief  
 * 
 * 
 */

#ifndef _classical_gen_model_h_
#define _classical_gen_model_h_

#include <base_gen_model.hpp>
#include <gridpack/include/gridpack.hpp>


class ClassicalGen: public BaseGenModel
{
   public:
  /**
     * Basic constructor
     */
    ClassicalGen();

    /**
     * Basic destructor
     */
    ~ClassicalGen();

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

    double getAngle();

    /**
     * Write out generator state
     * @param signal character string used to determine behavior
     * @param string buffer that contains output
     */
    void write(const char* signal, char* string);

    /**
     *  Set the number of variables for this generator model
     *  @param [output] number of variables for this model
     */
    bool vectorSize(int *nvar) const;

    /**
     * Set the internal values of the voltage magnitude and phase angle. Need this
     * function to push values from vectors back onto generators
     * @param values array containing generator state variables
     */
     void setValues(gridpack::ComplexType*);

    /**
     * Return the values of the generator vector block
     * @param values: pointer to vector values
     * @return: false if generator does not contribute
     *        vector element
     */
    bool vectorValues(gridpack::ComplexType *values);

    /**
     * Return the generator current injection (in rectangular form) 
     * @param [output] IGD - real part of the generator current
     * @param [output] IGQ - imaginary part of the generator current
     */
    void getCurrent(double *IGD, double *IGQ);

  private:
    // Machine parameters
    double p_Rs; // Machine stator resistance
    double p_Xdp;  // Machine transient reactance
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
