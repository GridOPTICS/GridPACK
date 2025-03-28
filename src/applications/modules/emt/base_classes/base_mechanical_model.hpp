/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   base_emt_mechanical_model.hpp
 * 
 * @brief  Base model for renewable energy mechanical part
 *         This base model serves as a common base model
 *         for drive train, aerodynamics, and torque controller
 */

#ifndef _base_emt_mechanical_model_h_
#define _base_emt_mechanical_model_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"
#include <gridpack/applications/modules/emt/constants.hpp>
#include <gridpack/applications/modules/emt/emtutilfunctions.hpp>
#include <gridpack/math/matrix.hpp>
#include <gridpack/applications/modules/emt/base_classes/base_gen_model.hpp>
#include <gridpack/math/dae_solver.hpp>

class BaseEMTGenModel; // Forward declaration for BaseGenModel
class BaseEMTExcModel; // Forward declaration for BaseExcModel
class BaseEMTPlantControllerModel; // Forward declaration for Base Plant Controller Model

class BaseEMTRMechModel : public gridpack::component::BaseComponent
{
  /* Inheriting the BaseComponent class allows use of functions
     for loading data and accessing/setting values in the vector/matrix
  */
public:
  /**
   * Basic constructor
   */
  BaseEMTRMechModel();
  
  /**
   * Basic destructor
   */
  virtual ~BaseEMTRMechModel();

  /**
     Prestep function
  */
  virtual void preStep(double time ,double timestep) { }
  
  /**
     Poststep function
  */
  virtual void postStep(double time) { }

  /**
    Number of variables
  */ 
  virtual void getnvar(int *nvar) { *nvar = nxrmech; }

  /**
     Get the plant controller output
  */
  virtual void getPrefQext(double *Pref, double *Qext) { }

  /**
     Note: This is a custom version of the load method from the BaseComponent Class. It takes in an extra argument idx to specify which component is being read. Ideally, this method should be moved to the MatVecInterface

   * Load data from DataCollection object into corresponding
   * component. This needs to be implemented by every component
   * @param data data collection associated with component
   */
  virtual void load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx);

  /**
   * Initialize exciter model before calculation
   * @param [output] values - array where initialized exciter variables should be set
   */
  virtual void init(gridpack::RealType *values);
  
  /**
   * Write output from exciters to a string.
   * @param string (output) string with information to be printed out
   * @param bufsize size of string buffer in bytes
   * @param signal an optional character string to signal to this
   * routine what about kind of information to write
   * @return true if bus is contributing string to output, false otherwise
   */
  virtual bool serialWrite(char *string, const int bufsize,
			   const char *signal);
  
  /**
   * Write out exciter state
   * @param signal character string used to determine behavior
   * @param string buffer that contains output
   */
  virtual void write(const char* signal, char* string);

  /**
   * Set bus voltage
   */
  void setVoltage(double busVD, double busVQ) {VD = busVD; VQ = busVQ; }
  
  /**
     set exciter status
  **/
  void setStatus(int estatus) {status = estatus;}
  
  /**
   * return the bolean indicating whether the exciter is ON or OFF
   */
  bool getStatus() {return status;}
  
  /**
   * set current time
   */
  void setTime(double time) {p_time = time; }

  /**
   * Set TSshift: This parameter is passed by PETSc and is to be used in the Jacobian calculation only.
   */
  void setTSshift(double inshift) {shift = inshift;}

  virtual void setEvent(gridpack::math::RealDAESolver::EventManagerPtr);

  /**
   * Set the offset for first variable for the exciter in the array for all bus variables 
   * @param offset offset
   */
  void setBusOffset(int offset) {offsetb = offset;}


  /**
   * Get number of matrix values contributed by generator
   * @return number of matrix values
   */
  virtual int matrixNumValues();

  /**
   * Return values from Jacobian matrix
   * @param nvals: number of values to be inserted
   * @param values: pointer to matrix block values
   * @param rows: pointer to matrix block rows
   * @param cols: pointer to matrix block cols
   */
  virtual void matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols);


  /**
   * Return vector values from the generator model 
   * @param values - array of returned values
   *
   * Note: This function is used to return the entries in vector,
   * for e.g., the entries in the residual vector from the generator
   * object
   */
  virtual void vectorGetValues(gridpack::RealType *values);

  /**
   * Pass solution vector values to the generator object
   * @param values - array of returned values
   *
   * Note: This function is used to pass the entries in vector
   * to the generator object,
   * for e.g., the state vector values for this generator
   */
  virtual void setValues(gridpack::RealType *values);
  
  /**
   * Set the generator associated with this exciter
   **/
  void setGenerator(BaseEMTGenModel* generator) {p_gen = generator; }

  /**
   * Set the electrical controller associated
   */
  void setElectricalController(BaseEMTExcModel* econtroller) {p_econ = econtroller; }

  /**
   * Set the plant controller associated
   */
  void setPlantController(BaseEMTPlantControllerModel* pcontroller) {p_pcon = pcontroller; }


  /**
   * Get the generator associated with this plant controller
   */
  BaseEMTGenModel* getGenerator() { return p_gen; }

  /**
   * Get the electrical controller associated
   */
  BaseEMTExcModel* getElectricalController() { return p_econ; }

  /**
   * Get the plant controller associated
   */
  BaseEMTPlantControllerModel* getPlantController() { return p_pcon; }

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

  void setBusLocalOffset(int offset) {p_busoffset = offset;}

  /*
    set the location for the first variable in the solution vector
  */
  void setGlobalLocation(int gloc) {p_gloc = gloc;}

  /**
   * return offset in the local vector 
   */
  int getLocalOffset()
  {
    return p_busoffset + offsetb;
  }

  virtual void resetEventFlags(void) {}

  /**
   * Set type of integration algorithm
   */
  void setIntegrationType(EMTMachineIntegrationType type) {integrationtype = type; }

    /**
   *  Set Pref0 - Sets reference power
   *  @param Pref0 - reference power
   *  Set by plant controller model
   **/
  virtual void setPref0(double Pref0) {}

  /**
   *  setPelec - Electrical power Pelec
   *  @param Pelec - Electrical power input
   **/
  virtual void setPelec(double Pelec) {}

  /**
   *  setPmech - Mechanical power Pmech
   *  @param Pmech - Mechanical power input
   **/
  virtual void setPmech(double Pmech) {}


  /**
   *  setGeneratorSpeedDeviation - Set the speed deviation
   *  @param - domega : speed deviation
   *  From drive train model
   **/
  virtual void setGeneratorSpeedDeviation(double domega) {}

  /**
   * setVdip - Voltage dip flag
   * @param vdip - flag to indicate voltage dip
   * From elecrical controller
   **/
  virtual void setVdip(bool Vdip) {}

  /**
   * getPref - Output of torque controller
   * @param  - Pref : reference power
   **/
  virtual double getPref() {}

  /**
   *  getOmegaref - Output of torque controller
   * @return - omega_ref : reference speed
   **/
  virtual double getOmegaref() {return 1.0;}

  /**
   *  setOmegaref - Output of torque controller
   * @param - omega_ref : reference speed
   **/
  virtual void setOmegaref(double omega_ref) {}


  /**   Pitch controller **/

  /**
   * setTurbineSpeedDeviation - sets the turbine speed deviation
   * @param domega_turb : turbine speed deviation
   * From drive train model
   **/
  virtual void setTurbineSpeedDeviation(double domega_turb) {}

  /**
   * setPord - sets Pord
   * @param Pord - electric Pord
   * From electrical controller model
   **/
  virtual void setPord(double Pord) {}

  /**
   *  setPord0 - Sets initial power order
   *  @param Pord0 - Initial value of reference power
   *
   **/
  virtual void setPord0(double Pord0) {}

  /**
   * getTheta - Get output of pitch controller
   * @output theta - output of pitch controller
   **/
  virtual double getTheta() {return 0.0;}

  /** Aerodynamic model **/

  /**
   * setTheta - sets pitch angle
   * @param Theta - pitch angle
   **/
  virtual void setTheta(double Theta) {}

  /**
   * getTaero - returns the aero-dynamic torque
   * @output Taero - aero dynamic torque, output of aerodynamic model
   **/
  virtual double getTaero() {return 0.0;}

  /** Drive train model **/

  /**
   * setTmech - sets mechanical torque
   * @param Tmech - mechanical torque
   * From aerodynamic model
   **/
  virtual void setTmech(double Tmech) {}

  /**
   * getTmech - sets mechanical torque
   * @return Tmech - mechanical torque
   **/
  virtual double getTmech() {return 0.0; }

  /**
   * setTelec - sets electrical torque
   * @param Telec - electrical torque
   * From generator?
   **/
  virtual void setTelec(double Tmech) {}

  /**
   * getTurbineSpeedDeviation - gets the turbine speed deviation
   * @ouput turbine speed deviation
   **/
  virtual double getTurbineSpeedDeviation() {return 0.0; }

  /**
   * getGenSpeedDeviation - gets the generator speed deviation
   * @ouput generator speed deviation
   **/
  virtual double getGeneratorSpeedDeviation() {return 0.0; }

  /**
   * getGenRotorAngleDeviation - gets the rotor angle deviation
   * @ouput rotor angle deviation
   **/
  virtual double getRotorAngleDeviation() {return 0.0; }

  void setPitchController(boost::shared_ptr<BaseEMTRMechModel> &pcontroller) {p_pitchcon = pcontroller; }

  boost::shared_ptr<BaseEMTRMechModel> getPitchController() { return p_pitchcon; }

  void setTorqueController(boost::shared_ptr<BaseEMTRMechModel> &tcontroller) {p_torquecon = tcontroller; }

  boost::shared_ptr<BaseEMTRMechModel> getTorqueController() { return p_torquecon; }

  void setAeroDynamicController(boost::shared_ptr<BaseEMTRMechModel> &acontroller) {p_aerocon = acontroller; }

  boost::shared_ptr<BaseEMTRMechModel> getAeroDynamicController() { return p_aerocon; }

  void setDriveTrainController(boost::shared_ptr<BaseEMTRMechModel> &dtcontroller) {p_drivetraincon = dtcontroller; }

  boost::shared_ptr<BaseEMTRMechModel> getDriveTrainController() { return p_drivetraincon; }

protected:
  int           status; /**< Plant status */
  double        mbase,sbase; /** Machine base and system MVA base */
  int           busnum; /** Bus number */
  std::string   id;
  double        p_time = 0.0;   /** Current time */
  double        shift; // shift (multiplier) used in the Jacobian calculation
  double        VD,VQ;
  
  EMTMachineIntegrationType integrationtype; // Integration type 

  BaseEMTGenModel* p_gen; // Generator model
  BaseEMTExcModel* p_econ; // Electrical Controller model
  BaseEMTPlantControllerModel* p_pcon; // Plant Controller Model
  boost::shared_ptr<BaseEMTRMechModel> p_pitchcon; // Pitch controller
  boost::shared_ptr<BaseEMTRMechModel> p_torquecon; // Torque controller
  boost::shared_ptr<BaseEMTRMechModel> p_aerocon; // Aerodynamic controller
  boost::shared_ptr<BaseEMTRMechModel> p_drivetraincon; // Drive Train controller

  int           offsetb; /**< offset for the first variable for the generator in the array for all bus variables */
  int           p_gloc; // Global location of the first variable for the generator

  int           nxrmech;    /** Number of variables for the model */
  int           p_busoffset; /** Starting location for bus variables in the local state vector */
  int           p_nrows;  // number of rows (equations) contributed by this excitor
  int           p_ncols;  // number of columns (variables) contributed by this excitor
  std::vector<int>   p_rowidx; // global index for rows
  std::vector<int>   p_colidx; // global index for columns

};

#endif
