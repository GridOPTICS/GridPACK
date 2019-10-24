/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   esst1a.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   09/20/19
 * 
 * @brief  
 * 
 * 
 */

#ifndef _esst1a_h_
#define _esst1a_h_

#include <base_exc_model.hpp>
#include <gridpack/include/gridpack.hpp>


class Esst1aExc: public BaseExcModel
{
   public:
   /**
     * Basic constructor
     */
    Esst1aExc();

    /**
     * Basic destructor
     */
    ~Esst1aExc();

    /**
     * Load parameters from DataCollection object into exciter model
     * @param data collection of exciter parameters from input files
     * @param index of exciter on bus
     * TODO: might want to move this functionality to BaseExciterModel
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int idx);

    /**
     * Saturation function
     * @ param x
     */
    double Sat(double x); 

    /**
     * Initialize exciter model before calculation
     * @param [output] values - array where initialized exciter variables should be set
     */
  void init(gridpack::ComplexType *values);

    /**
     * Write output from exciters to a string.
     * @param string (output) string with information to be printed out
     * @param bufsize size of string buffer in bytes
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if bus is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const int bufsize,
        const char *signal);

    /**
     * Write out exciter state
     * @param signal character string used to determine behavior
     * @param string buffer that contains output
     */
    void write(const char* signal, char* string);

    /**
     *  Set the number of variables for this exciter model
     *  @param [output] number of variables for this model
     */
    bool vectorSize(int *nvar) const;

    /**
     * Set the internal values of the voltage magnitude and phase angle. Need this
     * function to push values from vectors back onto exciters
     * @param values array containing exciter state variables
     */
     void setValues(gridpack::ComplexType*);

    /**
     * Return the values of the exciter vector block
     * @param values: pointer to vector values
     * @return: false if exciter does not contribute
     *        vector element
     */
    bool vectorValues(gridpack::ComplexType *values);

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
     * Set the initial field voltage (at t = tstart) for the exciter
     * @param fldv value of the field voltage
     */
    virtual void setInitialFieldVoltage(double fldv);

    /**
     * Set the field current parameter inside the exciter
     * @param fldc value of the field current
     */
    virtual void setFieldCurrent(double fldc);

    /** 
     * Get the value of the field voltage parameter
     * @return value of field voltage
     */
    virtual double getFieldVoltage();

    /** 
     * Get the value of the field current parameter
     * @return value of field current
     */
    virtual double getFieldCurrent();

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

    // Exciter esst1a parameters from dyr
    double UEL, VOS, Tr, Vimax, Vimin, Tc, Tb;
    double Tc1, Tb1, Ka, Ta, Vamax, Vamin;
    double Vrmax, Vrmin, Kc, Kf, Tf, Klr, Ilr;
    double Ta1;
    
    // ESST1A state variables
    double x1Va, x2Vcomp, x3LL1, x4LL2, x5Deriv;       
    double dx1Va, dx2Vcomp, dx3LL1, dx4LL2, dx5Deriv;    
   
    // ESST1A inputs
    double Vcomp, Vterm, Vstab;

    // Field Voltage Output
    double Efd;

    // Field Current Output
    double LadIfd;

    double Vref;

    double presentMag, presentAng;
    
    bool OptionToModifyLimitsForInitialStateLimitViolation;

    //bool flag2, flag3, flag4, flag5; //flags for residual function conditions
    
};

#endif
