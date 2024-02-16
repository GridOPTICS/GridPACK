/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file  repca1.hpp
 * 
 * @brief  Header file for REPCA1 plant controller
 * 
 * 
 */

#ifndef _repca1_h_
#define _repca1_h_

#include <base_plant_model.hpp>
#include <gridpack/include/gridpack.hpp>
#include "cblock.hpp"
#include "dblock.hpp"

class Repca1 : public BaseEMTPlantControllerModel
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  Repca1();
  
  /**
   * Basic destructor
   */
   ~Repca1();

  /**
     Prestep function
  */
  void preStep(double time ,double timestep);
  
  /**
     Poststep function
  */
  void postStep(double time);

  /**
    Number of variables
  */ 
  void getnvar(int *nvar);

  /**
     Get the plant controller output
  */
  void getPrefQext(double *Pref, double *Qext);

  /**
     Note: This is a custom version of the load method from the BaseComponent Class. It takes in an extra argument idx to specify which component is being read. Ideally, this method should be moved to the MatVecInterface

   * Load data from DataCollection object into corresponding
   * component. This needs to be implemented by every component
   * @param data data collection associated with component
   */
   void load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx);

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
  
   void resetEventFlags(void) {}

private:
  // Model parameters
  double ireg;   // Bus number of regulated bus
  double fbus, tbus; // To and from bus of the monitored branch
  std::string  br_ckt; // Branch circuit id
  int VCompFLAG; // droop flag
  int RefFLAG;  // Q control or V control
  int FreqFLAG;   // Freq control
  
  double Rc, Xc; // line compensation parameters
  double Kc;     // Gain for branch reactive power measurement
  double Tfltr;  // Filter time constant
  double dbd1,dbd2; // Deadband limits
  double Emax,Emin;  // Limits for voltage
  double Kp, Ki;     // PI controller gains for reactive power
  double Qmax,Qmin;  // Limits on Q PI controller
  double Tft, Tfv;    // Parameters for Qext lead lag block
  double Tp;        // Pbranch filter time constant
  double fdbd1,fdbd2; // Deadband limits for frequency error
  double Ddn, Dup;   // Droop constant
  double Kpg, Kig;   // Pext PI contoller gains
  double femax,femin; // Limits on Frequency error
  double Pmax, Pmin;  // Non wind-up limits for Pext PI controller
  double Tlag;      // Pext filter time constant
  double Tg;
  double Vfrz;      // State freeze voltage

  std::string p_gen_id;
  int p_bus_num;

  // Model blocks
  Filter V_filter_blk;       // V filter block
  double V_filter_blk_out;   // Output of V_filter block

  Filter Qbranch_filter_blk; // Q branch filter block
  double Qbranch_filter_blk_out; // Output of Q branch filter block

  Deadband VQerr_deadband; // Deadband for V or Q error
  double VQerr_deadband_out; // Output of VQerr_deadband block

  GainLimiter VQerr_limiter; // VQ err limiter block
  double VQerr_limiter_out; // Output of VQerr limiter block

  PIControl Qref_PI_blk; // PI control for Qref
  double Qref_PI_blk_out; // Output of PI control for Qref

  LeadLag Qref_leadlag_blk; // Lead lag block for Qref

  Deadband Freqerr_deadband; // Frequency error deadband

  Filter   Pbranch_filter_blk; // Pbranch filter block
  double   Pbranch_filter_blk_out; // Output of Pbranch filter block

  GainLimiter     Freqerr_limiter;   // Frequency error limiter
  double   Freqerr_limiter_out; // Output of frequency error limiter
  PIControl Pref_PI_blk;
  double    Pref_PI_blk_out; // Output of Pref PI block

  Filter    Pref_filter_blk;  // Pref filter block
  
  // Model inputs
  double Pbranch, Qbranch; // Line flow on designated branch (on system MVAbase)
  double Vreg;             // Voltage of regulated bus
  double Ibranch;          // Line current on designated branch (on system Ibase)
  double Freq;             // Frequency at POI

  double Freq_ref;         // Reference frequency
  // Model outputs
  double Pref,Qref;        // Output signal
  double Pg, Qg;           // Generator active and reactive power (on machine Mbase
  double Plant_ref;       // Reference plant active power output

  // Internal variables
  double Vt;
  double p_sbase,p_mbase; // System base and machine Mbase
  double Vref; // Reference voltage used in QV control
  bool   Vfreeze; // Flag to check if Vfreeze state is active (Vt < Vfrz)


};

// Class for defining events for Reeca1 model
class Repca1Event
  :public gridpack::math::RealDAESolver::Event
{
public:

  // Default constructor
  Repca1Event(Repca1 *repca1):gridpack::math::RealDAESolver::Event(0),p_plant(repca1)
  {
    std:fill(p_term.begin(),p_term.end(),false);

    std::fill(p_dir.begin(),p_dir.end(),gridpack::math::CrossZeroNegative);

  }

  // Destructor
  ~Repca1Event(void) {}
protected:
  Repca1 *p_plant;

  void p_update(const double& t, gridpack::RealType *state);

  void p_handle(const bool *triggered, const double& t, gridpack::RealType *state);
};


#endif
