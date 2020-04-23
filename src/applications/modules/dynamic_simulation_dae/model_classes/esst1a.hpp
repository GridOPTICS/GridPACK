/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   esst1a.hpp
 * @author Shuangshuang Jin
 * @author Shrirang Abhyankar
 * @Last modified:   04/22/20
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
     * Set Jacobian block
     * @param values a 2-d array of Jacobian block for the bus
     */
    bool setJacobian(gridpack::ComplexType **values);

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
     * Set the initial field voltage (at t = tstart) for the exciter
     * @param fldv value of the field voltage
     */
    void setInitialFieldVoltage(double fldv);

    /** 
     * Get the value of the field voltage parameter
     * @return value of field voltage
     */
    double getFieldVoltage();

    /**
     * Partial derivatives of field voltage Efd w.r.t. exciter variables
     * @param xexc_loc locations of exciter variables
     * @param dEfd_dxexc partial derivatives of field voltage w.r.t exciter variables
     * @param dEfd_dxgen partial derivatives of field voltage Efd w.r.t generator variables
     */
  bool getFieldVoltagePartialDerivatives(int *xexc_loc,double *dEfd_dxexc,double *dEfd_dxgen);

  private:

    // Exciter esst1a parameters from dyr
    int    UEL,VOS;
    double Tr, Vimax, Vimin, Tc, Tb;
    double Tc1, Tb1, Ka, Ta, Vamax, Vamin;
    double Vrmax, Vrmin, Kc, Kf, Tf, Klr, Ilr;
    
    // ESST1A state variables
    double Vmeas; // Measured voltage by transducer
    double xLL1; // State variable for first lead-lag block
    double xLL2; // State variable for second lead-lag block
    double Va;   // Voltage regulator output
    double xf;   // Feedback block state variable
    double x1Va, x2Vcomp, x3LL1, x4LL2, x5Deriv;

    // ESST1A derivatives
    double dVmeas, dxLL1, dxLL2, dVa, dxf;    

    // ESST1A previous step solution
    double Vmeasprev,xLL1prev,xLL2prev,Vaprev,xfprev;

    // ESST1A inputs
    double Ec; // Terminal voltage
    double Vothsg; // Voltage signal from stabilizer
    double Vuel; // Under excitation limiter voltage
    double Voel; // Over excitation limiter voltage
    double Vcomp, Vterm, Vstab;

    // Field Voltage Output
    double Efd;

    // Field Current Input
    double LadIfd;

    // Voltage regulator reference
    double Vref;

    //bool flag2, flag3, flag4, flag5; //flags for residual function conditions

    // Flag to denote whether each equation is algebraic or differential.
    // iseq_diff[i] = 1 if equation is differential, 0 otherwise.
    int iseq_diff[5];

    int bid;
};

#endif
