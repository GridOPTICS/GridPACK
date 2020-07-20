/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_gen_model.hpp
 * @author Shrirang Abhyankar
 * @Last modified:   04/24/20
 * 
 * @brief  Base generator class header file
 * 
 * 
 */

#ifndef _base_gen_model_h_
#define _base_gen_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"
#include <constants.hpp>
#include <base_exc_model.hpp>
#include <base_gov_model.hpp>

class BaseExcModel; // Forward declaration for BaseExcModel
class BaseGovModel; // Forward declaration for BaseGovModel

class BaseGenModel : public gridpack::component::BaseComponent
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  BaseGenModel();
  
  /**
   * Basic destructor
   */
  virtual ~BaseGenModel();
  
  /**
   * Initialize generator model before calculation
   * @param [output] values - array where initialized generator variables should be set
   */
  virtual void init(gridpack::ComplexType *values);
  
  /**
   * Write output from generators to a string.
   * @param string (output) string with information to be printed out
   * @param bufsize size of string buffer in bytes
   * @param signal an optional character string to signal to this
   * routine what about kind of information to write
   * @return true if bus is contributing string to output, false otherwise
   */
  virtual bool serialWrite(char *string, const int bufsize,
			   const char *signal);
  
  virtual double getAngle();
  
  /**
   * Write out generator state
   * @param signal character string used to determine behavior
   * @param string buffer that contains output
   */
  virtual void write(const char* signal, char* string);
  
  /**
   * return the bolean indicating whether the gen is ON or OFF
   */
  bool getGenStatus() {return status;}
  
  /**
   * Set bus voltage
   */
  void setVoltage(double busVD, double busVQ) {VD = busVD; VQ = busVQ; }
  
  /**
   * Set TSshift: This parameter is passed by PETSc and is to be used in the Jacobian calculation only.
   */
  void setTSshift(double inshift) {shift = inshift;}
  
  /**
   * Return the generator current injection (in rectangular form) 
   * @param [output] IGD - real part of the generator current
   * @param [output] IGQ - imaginary part of the generator current
   */
  virtual void getCurrent(double *IGD, double *IGQ);
  
  /**
     Note: This is a custom version of the load method from the BaseComponent Class. It takes in an extra argument idx to specify which component is being read. Ideally, this method should be moved to the MatVecInterface

   * Load data from DataCollection object into corresponding
   * component. This needs to be implemented by every component
   * @param data data collection associated with component
   */
  virtual void load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx);


  /**
   * Set Jacobian block
   * @param values a 2-d array of Jacobian block for the bus
   */
  virtual bool setJacobian(gridpack::ComplexType **values);

  /**
   * Set the field current parameter inside the exciter
   * @param fldc value of the field current
   */
  virtual double getFieldCurrent(void);

  /**
   * Return the rotor speed deviation
   * @param rotor speed deviation 
   */
  virtual double getRotorSpeedDeviation();

  /**
   * Return the location of speed rotor speed deviation variable in the bus array
   * @param rotor speed deviation location
   */
  virtual int getRotorSpeedDeviationLocation();

  /**
   * Set the offset for first variable for the generator in the array for all bus variables 
   * @param offset offset
   */
  void setBusOffset(int offset) {offsetb = offset;}

  void setExciter(boost::shared_ptr<BaseExcModel> &p_exciter);

  boost::shared_ptr<BaseExcModel> getExciter();
  
  bool hasExciter();

  void setGovernor(boost::shared_ptr<BaseGovModel> &p_governor);

  boost::shared_ptr<BaseGovModel> getGovernor();
  
  bool hasGovernor();

  /****************************************************
 The following methods are inherited from the BaseComponent class and are 
to be overwritten by the implementation */
  
  
  /**
   * Set an internal variable that can be used to control the behavior of the
   * component. This function doesn't need to be implemented, but if needed,
   * it can be used to change the behavior of the component in different phases
   * of the calculation. For example, if a different matrix needs to be
   * generated at different times, the mode of the calculation can changed to
   * get different values from the MatVecInterface functions
   * @param mode integer indicating which mode should be used
   */
  void setMode(int mode) { p_mode = mode;}

  /**
   * Return size of vector block contributed by component
   * @param isize number of vector elements
   * @return false if network component does not contribute
   *        vector element
   */
  bool vectorSize(int *isize) const;

  /**
   * Return the values of the vector block
   * @param values pointer to vector values
   * @return false if network component does not contribute
   *        vector element
   */
  bool vectorValues(gridpack::ComplexType *values);

  /**
   * Set values in the component based on values in a vector or
   * matrix
   * @param values values in vector or matrix
   */
  void setValues(gridpack::ComplexType *values);

  void setBusLocalOffset(int offset) {p_busoffset = offset;}

  /**
   * return offset in the local vector 
   */
  int getLocalOffset()
  {
    return p_busoffset + offsetb;
  }

 protected:
  double        pg; /**< Generator active power output */
  double        qg; /**< Generator reactive power output */
  double        mbase; /**< MVA base of the machine */
  int           status; /**< Machine status */
  double        sbase;  /** The system MVA base */
  double        shift; // shift (multiplier) used in the Jacobian calculation.
  double        VD, VQ;
  bool          p_hasExciter; // Flag indicating whether this generator has exciter
  bool          p_hasGovernor; // Flag indicating whether this generator has governor
  boost::shared_ptr<BaseExcModel> p_exciter; // Exciter
  boost::shared_ptr<BaseGovModel> p_governor; // Governor
  int           offsetb; /**< offset for the first variable for the generator in the array for all bus variables */
  int           nxgen; /* Number of variables for the generator model */
  int           p_busoffset; /** Offset for the bus variables in the local vector. Used only for events */

  // Arrays used in coupling blocks between generator and exciter. These should be allocated and destroyed by the derived class
  int           *xexc_loc;   // locations for exciter variables in the bus variables array
  double        *dEfd_dxexc; // Partial derivatives of field voltage Efd w.r.t. exciter variables (size = nxexc)
  double        *dEfd_dxgen; // Partial derivatives of field voltage Efd w.r.t. generator variables (size = nxgen)

  // Arrays used in coupling blocks between generator and governor. These should be allocated and destroyed by the derived class
  int           *xgov_loc;   // locations for governor variables in the bus variables array
  double        *dPmech_dxgov; // Partial derivatives of mechanical power Pmech w.r.t. governor variables (size = nxgov)
};

#endif
