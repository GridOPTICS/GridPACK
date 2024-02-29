/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   sexs.hpp
 * 
 * @brief  


                                                        EMAX
                                                       ------
                                                      /
                                                     /
                            -----------        ------------ 
 *                          |          |       |           |
 *                          | 1 + TAs  |       |     K     |
 *   Vref - Ec + Vs ------> |--------- |------>| --------- |-----> Efd
                            | 1 + TBs  |       |  1 + TEs  |
 *                          |          |       |           |
 *                          ------------       ------------
                                                     / 
                                                    /
                                               -----
                                                EMIN
*/                 

#ifndef _sexs_h_
#define _sexs_h_

#include <base_exc_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "cblock.hpp"
#include "dblock.hpp"

class Sexs : public BaseEMTExcModel
{
  public:
    /**
     * Basic constructor
     */
    Sexs();

    /**
     * Basic destructor
     */
    ~Sexs();

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
     * Set the value of the Vterminal
     * 
     */
    //void setVterminal(double mag);

    /**
     * Set the initial voltage stabilizer signal
     * @param stabilizing signal input
    **/
    //void setVstab(double vstab);
	
  private:

    // Internal variables
    bool zero_TE; // If TE == 0 then the filter block is replaced by gain block

    // Exciter SEXS parameters from dyr
    double TA_OVER_TB, TA, TB, K, TE, EMIN, EMAX;

    // SEXS state variables
    double Vmeas; // First integrator
    double xLL; // Second integrator

    // SEXS derivatives
    double dVmeas;
    double dxLL;

    // Control blocks (when using explicit integration)
    // Model transfer blocks
    LeadLag leadlagblock;
    Filter  filterblock;
    GainLimiter    gainblock;

    // Model inputs
    double Ec; // Terminal voltage
    double Vref; // Reference voltage
    double Vs;  // Stabilizing voltage signal
  
    // Model outputs 
    double Efd;     // Field Voltage

    //std::string p_ckt; // id of the generator where the exciter is installed on
    //int p_bus_id;  // bus number of the generator 

};

// Class for defining events for Ieeet1 model
class SexsEvent
  :public gridpack::math::RealDAESolver::Event
{
public:

  // Default constructor
  SexsEvent(Sexs *exc):gridpack::math::RealDAESolver::Event(2),p_exc(exc)
  {
    std:fill(p_term.begin(),p_term.end(),false);

    std::fill(p_dir.begin(),p_dir.end(),gridpack::math::CrossZeroNegative);

  }

  // Destructor
  ~SexsEvent(void) {}
protected:
  Sexs *p_exc;

  void p_update(const double& t, gridpack::RealType *state);

  void p_handle(const bool *triggered, const double& t, gridpack::RealType *state);
};


#endif
