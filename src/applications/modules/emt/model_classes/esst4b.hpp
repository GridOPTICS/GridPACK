/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   esst4b.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   Aug 31, 2024
 * 
 * @brief  
 * 
 * 
 */

#ifndef _esst4b_h_
#define _esst4b_h_

#include <base_exc_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "cblock.hpp"
#include "dblock.hpp"

class Esst4bExc : public BaseEMTExcModel
{
  public:
    /**
     * Basic constructor
     */
    Esst4bExc();

    /**
     * Basic destructor
     */
    ~Esst4bExc();

    /**
     * Load parameters from DataCollection object into exciter model
     * @param data collection of exciter parameters from input files
     * @param index of exciter on bus
     * TODO: might want to move this functionality to BaseExciterModel
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection>
	      data, int idx);

    /**
     Number of variables
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
     * Saturation function
     * @ param x
     */
    double Sat(double x);

    double sqr(double x);

    /**
     * FEX function
     * @ param IN
     */
    double FEX(double IN);

    /**
     * CalculateVb function
     * @ param Vterm, Theta, Ir, Ii, LadIfd
     */
    double CalculateVb(double Vterm, double Theta, double Ir, double Ii, double LadIfd);

    /**
     * Initialize exciter model before calculation
     * @param [output] values - array where initialized exciter variables should be set
     */
    void init(gridpack::RealType *values);

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

    /** 
     * Get the value of the field voltage parameter
     * @return value of field voltage
     */
    double getFieldVoltage();
    /** 
     * Get the value of the field voltage parameter
     * and its global location
     * @return value of field voltage
     */
    double getFieldVoltage(int *Efd_gloc);

    /**
     * Partial derivatives of field voltage Efd w.r.t. exciter variables
     * @param xexc_loc locations of exciter variables
     * @param dEfd_dxexc partial derivatives of field voltage w.r.t exciter variables
     * @param dEfd_dxgen partial derivatives of field voltage Efd w.r.t generator variables
     */
    bool getFieldVoltagePartialDerivatives(int *xexc_loc,double *dEfd_dxexc,double *dEfd_dxgen);

    /**
     * Set Event 
     */
    void setEvent(gridpack::math::RealDAESolver::EventManagerPtr);

    /**
     * Update the event function values
     */
    void eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues);

    /**
     * Event handler function 
     */
    void eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state);

    /**
     * Set the initial field voltage value
     * @param initial value of the field voltage
     */
    void setFieldVoltage(double fldv);

    /**
     * Set the field current parameter inside the exciter
     * @param fldc value of the field current
     */
    void setFieldCurrent(double fldc);

    /** 
     * Get the value of the field current parameter
     * @return value of field current
     */
    double getFieldCurrent();

    /** 
     * Set the value of the Vterminal
     * @return value of field current
     */
    void setVterminal(double mag);

    /** 
     * Set the value of the Omega 
     * @return value of field current
     */
    void setOmega(double omega);

    void setVuel(double vtmp);

    void setVs(double vtmp);

    /** 
     * Set the value of the terminal current
     */
    void setIri(double vIr, double vIi);

    void setVoel(double vtmp);
    
    void setVcomp(double vtmp);

  private:

    //double S10, S12; 

    // Exciter ESST4B parameters from dyr
    double Tr, Kpr, Kir, Vrmax, Vrmin, Ta, Kpm;
    double Kim, Vmmax, Vmmin, Kg, Kp, KI;
    double Vbmax, Kc, Xl, Kpang, Vgmax, Thetap;
    
    // ESST4B state variables
    double x1Vm, x2Vcomp, x3Va, x4Vr;    
    double x1Vm_1, x2Vcomp_1, x3Va_1, x4Vr_1;
    double dx1Vm, dx2Vcomp, dx3Va, dx4Vr;    
    double dx1Vm_1, dx2Vcomp_1, dx3Va_1, dx4Vr_1;

    bool zero_TR; // Time constant TR for measurement block is zero, no transfer function
    bool zero_TA; // Time constant TA for measurement block is zero, no transfer function  
    
    bool zero_KIM; // Integrator gain KIM for PI controller is zero, no transfer function
    bool zero_KIR; // Integrator gain KIR for PI controller is zero, no transfer function

    Filter Filter_blkR;
    PIControl PIControl_blkR;
    Filter Filter_blkA;
    PIControl PIControl_blmM;
    LVGate LVGate_blk;

    double Vuel, Vs;
    double Vref; // Reference voltage
    double Vmeas; // Ouptut of voltage measurement block
    double Voel;
   
    // ESST4B inputs
    double Ec; // Terminal voltage
    double Vcomp, Vterm, Theta, Ir, Ii, LadIfd, Vstab;
 
    // Field Voltage Output
    double Efd;

    // Keep around 
    //double Vref;
    double Kpvr, Kpvi, Kpir, Kpii;

    double presentMag, presentAng;
  
    bool OptionToModifyLimitsForInitialStateLimitViolation;

};

// Class for defining events for ESST4b model
class Esst4bExcEvent
  :public gridpack::math::RealDAESolver::Event
{
public:

  // Default constructor
  Esst4bExcEvent(Esst4bExc *exc):gridpack::math::RealDAESolver::Event(2),p_exc(exc)
  {
    std:fill(p_term.begin(),p_term.end(),false);

    std::fill(p_dir.begin(),p_dir.end(),gridpack::math::CrossZeroNegative);

  }

  // Destructor
  ~Esst4bExcEvent(void) {}

protected:
  Esst4bExc *p_exc;

  void p_update(const double& t, gridpack::ComplexType *state);

  void p_handle(const bool *triggered, const double& t, gridpack::ComplexType *state);
};

#endif