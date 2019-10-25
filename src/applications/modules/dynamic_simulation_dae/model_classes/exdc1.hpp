/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   exdc1.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   08/06/19
 * 
 * @brief  
 * 
 * 
 */

#ifndef _exdc1_h_
#define _exdc1_h_

#include <base_exc_model.hpp>
#include <gridpack/include/gridpack.hpp>


class Exdc1Exc: public BaseExcModel
{
   public:
   /**
     * Basic constructor
     */
    Exdc1Exc();

    /**
     * Basic destructor
     */
    ~Exdc1Exc();

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

  private:
    // Exciter EXDC1 parameters from dyr
    double TR, KA, TA, TB, TC, Vrmax, Vrmin;
    double KE, TE, KF, TF, SWITCH; // TF?
    double E1, SE1, E2, SE2;

    // EXDC1 state variables
    double x1, x2, x3, x4, x5;
    double dx1, dx2, dx3, dx4, dx5;

    // Field Voltage Output
    double Efd;

    // Field Current Output
    double LadIfd;

    double Vref;

    double Vterminal;

    bool flag2, flag3, flag4, flag5; //flags for residual function conditions
    double A, B; // Sat function variables

};

#endif
