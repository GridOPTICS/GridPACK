/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   wsieg1.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   10/02/19
 * 
 * @brief  
 * 
 * 
 */

#ifndef _wsieg1_h_
#define _wsieg1_h_

#include <base_gov_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "GainBlockClass.hpp"
#include "BackLashClass.hpp"
#include "DBIntClass.hpp"

class Wsieg1Gov: public BaseGovModel
{
   public:
   /**
     * Basic constructor
     */
    Wsieg1Gov();

    /**
     * Basic destructor
     */
    ~Wsieg1Gov();

    /**
     * Load parameters from DataCollection object into governor model
     * @param data collection of governor parameters from input files
     * @param index of governor on bus
     * TODO: might want to move this functionality to BaseGoviterModel
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Initialize governor model before calculation
     * @param [output] values - array where initialized governor variables should be set
     */
  void init(gridpack::ComplexType *values);

    /**
     * Write output from governors to a string.
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
     * Write out governor state
     * @param signal character string used to determine behavior
     * @param string buffer that contains output
     */
    void write(const char* signal, char* string);

    /**
     *  Set the number of variables for this governor model
     *  @param [output] number of variables for this model
     */
    bool vectorSize(int *nvar) const;

    /**
     * Set the internal values of the voltage magnitude and phase angle. Need this
     * function to push values from vectors back onto governors
     * @param values array containing governor state variables
     */
     void setValues(gridpack::ComplexType*);

    /**
     * Return the values of the governor vector block
     * @param values: pointer to vector values
     * @return: false if governor does not contribute
     *        vector element
     */
    bool vectorValues(gridpack::ComplexType *values);

    /**
     * Return the governor current injection (in rectangular form) 
     * @param [output] IGD - real part of the governor current
     * @param [output] IGQ - imaginary part of the governor current
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

    /**
     * Set the mechanical power parameter inside the governor
     * @param pmech value of the mechanical power 
     */
    virtual void setMechanicalPower(double pmech);

    /**
     * Set the rotor speed deviation parameter inside the governor
     * @param delta_o value of the rotor speed deviation 
     */
    virtual void setRotorSpeedDeviation(double delta_o);

    /** 
     * Get the value of the mechanical power parameter
     * @return value of the mechanical power 
     */
    virtual double getMechanicalPower();

    /** 
     * Get the value of the rotor speed deviation parameter
     * @return value of the rotor speed deviation 
     */
    virtual double getRotorSpeedDeviation();

    /**
     * Set the value of the Vcomp
     * @return value of teh Vcomp
     */
    virtual void setVcomp(double vtmp);

    /**
     * Set the value of the time step
     * @return value of the time step
     */
    //virtual void setTimestep(double timestep);
 
    /**
     * Set the value of the time increment 
     * @return value of the time increment
     */
    //virtual void setTimeincrement(double timeincrement);


  private:
    // Governor WSIEG1 Parameters read from dyr
    double K, T1, T2, T3, Uo, Uc, Pmax, Pmin;
    double T4, K1, K2, T5, K3, K4, T6, K5, K6, T7, K7, K8;
    double Db1, Err, Db2;
    //double Gv1, PGv1, Gv2, PGv2, Gv3, PGv3, Gv4, PGv4, Gv5, PGv5;
    double Iblock;

    // WSIEG1 state variables
    double x1LL, x2GovOut, x3Turb1, x4Turb2, x5Turb3, x6Turb4;
    double x1LL_1, x2GovOut_1, x3Turb1_1, x4Turb2_1, x5Turb3_1, x6Turb4_1;
    double dx1LL, dx2GovOut, dx3Turb1, dx4Turb2, dx5Turb3, dx6Turb4;
    double dx1LL_1, dx2GovOut_1, dx3Turb1_1, dx4Turb2_1, dx5Turb3_1, dx6Turb4_1;

    // Outputs: Mechnical Power Gen1 and Gen 2
    double Pmech1, Pmech2;

    bool SecondGenExists, OptionToModifyLimitsForInitialStateLimitViolation;

    GainBlockClass GainBlock;
    BackLashClass BackLash;
    DBIntClass DBInt;

    double Pref;
    double w;

    //bool flag2, flag3, flag4, flag5; //flags for residual function conditions
    //double A, B; // Sat function variables
    
};

#endif
