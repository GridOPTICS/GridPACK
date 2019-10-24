/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   genrou.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   07/25/19
 * 
 * @brief  
 * 
 * 
 */

#ifndef _genrou_h_
#define _genrou_h_

#include <base_gen_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "boost/smart_ptr/shared_ptr.hpp"

class GenrouGen: public BaseGenModel
{
   public:
  /**
     * Basic constructor
     */
    GenrouGen();

    /**
     * Basic destructor
     */
    ~GenrouGen();

    /**
     * Load parameters from DataCollection object into generator model
     * @param data collection of generator parameters from input files
     * @param index of generator on bus
     * TODO: might want to move this functionality to BaseGeneratorModel
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Saturation function
     * @ param x
     */
    double Sat(double x); 

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

    /**
     * Return the matrix entries
     * @param [output] nval - number of values set
     * @param [output] row - row indices for matrix entries
     * @param [output] col - col indices for matrix entries
     * @param [output] values - matrix entries
     * return true when matrix entries set
     */
    bool matrixDiagEntries(int *nval,int *row, int *col, gridpack::ComplexType *values);

    //void setExciter(boost::shared_ptr<BaseExcModel> &p_exciter);

    //boost::shared_ptr<BaseExcModel> getExciter();
    
    //void setphasExciter(bool info);

  private:
    /*// Machine parameters
    double p_Rs; // Machine stator resistance
    double p_Xdp;  // Machine transient reactance
    double p_H;      // Machine Inertia constant
    double p_D;      // Machine damping coefficient

    // Internal constants
    double p_Pm;  // Mechanical power input
    double p_Ep;  // Internal emf

    // Generator variables and their derivatives
    double p_delta,p_dw;
    double p_deltadot,p_dwdot;*/

    // Machine parameters
    double H, D, Ra, Xd, Xq, Xdp, Xdpp, Xl, Xqp, Xqpp;
    double Tdop, Tdopp, Tqopp, S10, S12, Tqop;

    // Internal constants
    double Efd; // Field voltage
    double LadIfd; // Field current
    double Pmech; // Mechanical power
  
    // Generator variables and their derivatives
    double Vterm, Theta; // Terminal voltage magnitude and angle
    double Ir, Ii; // Terminal current
    double x1d, x2w, x3Eqp, x4Psidp, x5Psiqp, x6Edp; 
    double dx1d, dx2w, dx3Eqp, dx4Psidp, dx5Psiqp, dx6Edp;
    double Id, Iq; // Generator current on d and q axis

    double IrNorton, IiNorton;
    gridpack::ComplexType INorton;

    bool p_hasExciter;
    boost::shared_ptr<BaseExcModel> p_exciter;
    bool p_hasGovernor;
    boost::shared_ptr<BaseGovModel> p_governor;

    double B, G, Vrterm, Viterm;
    double Vd, Vq;
    
    double presentMag, presentAng;
};

#endif