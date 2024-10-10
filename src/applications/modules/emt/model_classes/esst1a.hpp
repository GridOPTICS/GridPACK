/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   esst1a.hpp
 * 
 * @brief ESST1 exciter model header file 
 * @last updated by Shuangshuang Jin on Aug 23, 2024
 * 
 * 
 */

#ifndef _esst1a_h_
#define _esst1a_h_

#include <base_exc_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "cblock.hpp"
#include "dblock.hpp"

class Esst1aExc: public BaseEMTExcModel
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
   * Set the initial field voltage (at t = tstart) for the exciter
   * @param fldv value of the field voltage
   */
  //void setInitialFieldVoltage(double fldv);
  
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
   * Updated limiter flags after event has occured. Only called when the network is resolved
   */
  //void resetEventFlags(void);

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
  double Vterm, Vstab;
  
  // Initial Field Voltage (at t= 0)
  double Efd0;
  
  // Field Current Input
  double LadIfd;
  
  // Voltage regulator reference
  double Vref;
  
  //bool flag2, flag3, flag4, flag5; //flags for residual function conditions
  
  // Flag to denote whether each equation is algebraic or differential.
  // iseq_diff[i] = 1 if equation is differential, 0 otherwise.
  int iseq_diff[5];

  bool Efd_at_min,Efd_at_max;
  bool Vi_at_min,Vi_at_max;
  bool Va_at_min,Va_at_max;

  Filter Filter_blkR;
  HVGate HVGate_blk1; 
  LeadLag Leadlag_blkBC;
  LeadLag Leadlag_blkBC1;
  Filter Regulator_blk;
  GainLimiter Regulator_gain_blk;
  HVGate HVGate_blk2; 
  LVGate LVGate_blk; 
  Cblock Feedback_blkF;

  double Vf; // Output of Feedback block
  bool   zero_TA;      // Time constant TA for regulator block zero, no transfer function
  bool   zero_TR;      // Time constant TR for measurement block is zero, no transfer function

  bool zero_TF;   // Time constant TF for feedback block, if too small, no transfer function
  bool zero_TB;   // Time constant TB for first lead lag block, if too small, no transfer function
  bool zero_TB1;   // Time constant TB1 for second lead lag block, if too small, no transfer function
  bool OptionToModifyLimitsForInitialStateLimitViolation;

  double VA; // Output of Regulator blk
  double VLL1; // Output of LeadLag blk BC1
  double VLL; // Output of LeadLag blk BC
  //double Vref; // Reference voltage
  //double Vmeas; // Output of voltage measurement block

  double Efd;
};


// Class for defining events for ESST1a model
class Esst1aExcEvent
  :public gridpack::math::RealDAESolver::Event
{
public:

  // Default constructor
  Esst1aExcEvent(Esst1aExc *exc):gridpack::math::RealDAESolver::Event(2),p_exc(exc)
  {
    std:fill(p_term.begin(),p_term.end(),false);

    std::fill(p_dir.begin(),p_dir.end(),gridpack::math::CrossZeroNegative);

  }

  // Destructor
  ~Esst1aExcEvent(void) {}

protected:
  Esst1aExc *p_exc;

  void p_update(const double& t, gridpack::ComplexType *state);

  void p_handle(const bool *triggered, const double& t, gridpack::ComplexType *state);
};

#endif
