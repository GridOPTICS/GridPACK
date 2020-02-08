/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_components.hpp
 * @author Bruce Palmer
 * @date   2016-07-14 13:27:00 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _pf_components_h_
#define _pf_components_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/applications/components/y_matrix/ymatrix_components.hpp"

namespace gridpack {
namespace powerflow {

enum PFMode{YBus, Jacobian, RHS, S_Cal, State};

class PFBus
  : public gridpack::ymatrix::YMBus
{
  public:
    /**
     *  Simple constructor
     */
    PFBus(void);

    /**
     *  Simple destructor
     */
    ~PFBus(void);

    /**
     * Return size of matrix block on the diagonal contributed by
     * component
     * @param isize, jsize: number of rows and columns of matrix block
     * @return: false if network component does not contribute matrix
     *        element
     */
    bool matrixDiagSize(int *isize, int *jsize) const;

    /**
     * Return the values of for a diagonal matrix block. The values are
     * returned in row-major order
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute
     *        matrix element
     */
    bool matrixDiagValues(ComplexType *values);
    bool matrixDiagValues(RealType *values);

    /**
     * Return size of vector block contributed by component
     * @param isize: number of vector elements
     * @return: false if network component does not contribute
     *        vector element
     */
    bool vectorSize(int *isize) const;

    /**
     * Return the values of the vector block
     * @param values: pointer to vector values
     * @return: false if network component does not contribute
     *        vector element
     */
    bool vectorValues(ComplexType *values);
    bool vectorValues(RealType *values);

    /**
     * Set the internal values of the voltage magnitude and phase angle. Need this
     * function to push values from vectors back onto buses 
     * @param values array containing voltage magnitude and angle
     */
    void setValues(gridpack::ComplexType *values);
    void setValues(gridpack::RealType *values);

    /**
     * Return the size of the buffer used in data exchanges on the network.
     * For this problem, the voltage magnitude and phase angle need to be exchanged
     * @return size of buffer
     */
    int getXCBufSize(void);

    /**
     * Assign pointers for voltage magnitude and phase angle
     */
    void setXCBuf(void *buf);

    /**
     * Set values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    void setYBus(void);

    /**
     * Get values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    gridpack::ComplexType getYBus(void);

    /**
     * Load values stored in DataCollection object into PFBus object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       bus that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);

    /**
     * Set the mode to control what matrices and vectors are built when using
     * the mapper
     * @param mode: enumerated constant for different modes
     */
    void setMode(int mode);

    /**
     * Reset voltage and phase angle to initial values
     */
    void resetVoltage(void);

    /**
     * Set voltage limits on bus
     * @param vmin lower value of voltage
     * @param vmax upper value of voltage
     */
    void setVoltageLimits(double vmin, double vmax);

    /**
     * Check voltage for violation
     * @return false if there is a voltage violation
     */
    bool checkVoltageViolation(void);

    /**
     * Return the value of the voltage magnitude on this bus
     * @return voltage magnitude
     */
    double getVoltage(void);

    /**
     * Return the complex voltage on this bus
     * @return the complex voltage
     */
    ComplexType getComplexVoltage(void);

    /**
     * Return the value of the phase angle on this bus
     * @return: phase angle
     */
    double getPhase(void);

    /**
     * Get generator status
     * @param gen_id generator ID
     * @return current status of generator
     */
    bool getGenStatus(std::string gen_id);

    /**
     * Get list of generator IDs
     * @return vector of generator IDs
     */
    std::vector<std::string> getGenerators();

    /**
     * Get list of load IDs
     * @return vector of load IDs
     */
    std::vector<std::string> getLoads();

    /**
     * Return whether or not the bus is a PV bus (V held fixed in powerflow
     * equations)
     * @return true if bus is PV bus
     */
    bool isPV(void);

    /**
     * Set voltage value
     */
    void setVoltage(void);

    /**
     * Set phase angle value
     */
    void setPhase(void);

    /**
     * Set generator status
     * @param gen_id generator ID
     * @param status generator status
     */
    void setGenStatus(std::string gen_id, bool status);

    /**
     * Set isPV status
     * @param status isPV status
     */
    void setIsPV(int status);

    /**
     * Reset isPV status
     */
    void resetIsPV();

    /**
     * setSBus
    BUS = (CG*(GEN(ON,PG) + J*GEN(ON,QG)-(PD+J*QD))/BASEMVA
    */
    void setSBus(void);

    /**
     * Update pg of specified bus element based on their genID
     * @param busID
     * @param genID
     * @param value
     **/
//    void updatePg(int busID, std::string genID, double value);

    /**
     * Write output from buses to standard out
     * @param string (output) string with information to be printed out
     * @param signal an optional character string to signal to this
     * @param bufsize size of string buffer in bytes
     * routine what about kind of information to write
     * @return true if bus is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const int bufsize, const char *signal = NULL);

    /**
     * chkQlim
     check QLIM violations
    */
    bool chkQlim(void);

    /**
     * Clear changes that were made for Q limit violations and reset
     * bus to its original state
     */
    void clearQlim();

    /**
     * Save state variables inside the component to a DataCollection object.
     * This can be used as a way of moving data in a way that is useful for
     * creating output or for copying state data from one network to another.
     * @param data data collection object into which new values are inserted
     */
    void saveData(boost::shared_ptr<gridpack::component::DataCollection>
          data);

    /**
     * Modify parameters inside the bus module. This is designed to be
     * extensible
     * @param name character string describing parameter to be modified
     * @param value new value of parameter
     * @param idx index (if necessary) of variable to be modified
     */
    void setParam(std::string name, int busID, std::string genID, double value);
    void setParam(int busID, std::string genID, double value);

    /**
     * Access parameters inside the bus module. This is designed to be
     * extensible
     * @param name character string describing parameter to be accessed
     * @param value value of parameter
     * @param idx index (if necessary) of variable to be accessed
     */
    void getParam(std::string &name, double *value, int idx);
    void getParam(std::string &name, int *value, int idx);

    /**
     * Get index of internal bus element based on character string identifier
     * @param name character string describing element
     * @param tag character string specifying bus element
     * @return index of element
     */
    int getElementIndex(std::string &name, std::string &tag);

    /**
     * Set parameter to ignore voltage violations
     * @param flag value of ignore parameter
     */
    void setIgnore(bool flag);

    /**
     * Get parameter to ignore voltage violations
     * @return value of ignore parameter
     */
    bool getIgnore();

    /**
     * Get area parameter for bus
     * @return bus area index
     */
    int getArea();

    /**
     * Get zone parameter for bus
     * @return bus zone index
     */
    int getZone();

    /**
     * Evaluate diagonal block of Jacobian for power flow calculation and return
     * result as an array of real values
     * @param rvals values of Jacobian block
     * @return number of values returned
     */
    int diagonalJacobianValues(double *rvals);

    /**
     * Evaluate RHS values for powerflow equation and return result as an array
     * of real values
     * @param rvals values of Jacobian block
     * @return number of values returned
     */
    int rhsValues(double *rvals);

    /**
     * Push p_isPV values from exchange buffer to p_isPV variable
     */
    void pushIsPV();

    /**
     * Get vector containing generator participation
     * @return vector of generator participation factors
     */
    std::vector<double> getGeneratorParticipation();

    /**
     * Set value of real power on individual generators
     * @param tag generator ID
     * @param value new value of generator real power
     * @param data data collection object associated with bus
     */
    void setGeneratorRealPower(std::string tag, double value,
        gridpack::component::DataCollection *data);

    /**
     * Scale value of real power on all generators
     * @param character ID for generator
     * @param value scale factor for real power
     */
    void scaleGeneratorRealPower(std::string tag, double value);

    /**
     * Set value of real power on individual loads
     * @param tag load ID
     * @param value new value of load real power
     * @param data data collection object associated with bus
     */
    void setLoadRealPower(std::string tag, double value,
        gridpack::component::DataCollection *data);

    /**
     * Scale value of real and reactive power on loads
     * @param character ID for load
     * @param value scale factor for real power
     */
    void scaleLoadPower(std::string tag, double value);

    /**
     * Reset power for generators and loads back to original values
     */
    void resetPower();

    /**
     * Get available margin for generator
     * @param tag character ID for generator
     * @param current initial generation
     * @param pmin minimum allowable generation
     * @param pmax maximum allowable generation
     * @param status current status of generator
     */
    void getGeneratorMargins(std::vector<std::string> &tag,
        std::vector<double> &current, std::vector<double> &pmin,
        std::vector<double> &pmax, std::vector<int> &status);

    /**
     * Get current value of loads
     * @param tag character ID for load
     * @param pl initial value of load real power
     * @param ql initial value of load reactive power
     * @param status current status of load
     */
    void getLoadPower(std::vector<std::string> &tag, std::vector<double> &pl,
        std::vector<double> &ql, std::vector<int> &status);

    /**
     * Label bus as a source for real time path rating
     * @param flag identify bus as source
     */
    void setSource(bool flag);

    /**
     * Label bus as a sink for real time path rating
     * @param flag identify bus as sink
     */
    void setSink(bool flag);

    /**
     * Store scale factor
     * @param scale factor for scaling generation or loads
     */
    void setScale(double scale);

  private:
    double p_shunt_gs;
    double p_shunt_bs;
    bool p_shunt;
    bool p_load;
    int p_mode;
    bool p_ignore;

    // p_v and p_a are initialized to p_voltage and p_angle respectively,
    // but may be subject to change during the NR iterations
    double p_v, p_a;
    double p_theta; //phase angle difference
    double p_ybusr, p_ybusi;
    double p_P0, p_Q0; //double p_sbusr, p_sbusi;
    double p_angle;   // initial bus angle read from parser
    double p_voltage; // initial bus voltage read from parser
    // newly added priavate variables:
    std::vector<double> p_pg, p_qg, p_pFac;
    std::vector<double> p_savePg;
    std::vector<int> p_gstatus;
    std::vector<int> p_gstatus_save;
    std::vector<double> p_qmax,p_qmin;
    std::vector<double> p_qmax_orig, p_qmin_orig, p_pFac_orig;
    std::vector<double> p_vs;
    std::vector<std::string> p_gid;
    std::vector<double> p_pt;
    std::vector<double> p_pb;
    std::vector<double> p_pl, p_ql,p_ip,p_iq,p_yp,p_yq;
    std::vector<double> p_savePl;
    std::vector<double> p_saveQl;
    std::vector<int> p_lstatus;
    std::vector<std::string> p_lid;
    double p_sbase;
    double p_Pinj, p_Qinj;
    double p_vmin, p_vmax;
    bool p_isPV, p_saveisPV, p_save2isPV;
    bool *p_PV_ptr;
    int p_ngen;
    int p_nload;
    int p_type;
    int p_area;
    int p_zone;
    bool p_source;
    bool p_sink;
    double p_rtpr_scale;
    bool p_original_isolated;

    /**
     * Variables that are exchanged between buses
     */
    double* p_vMag_ptr;
    double* p_vAng_ptr;
    
    /**
     * Cache a pointer to DataCollection object
     */
    gridpack::component::DataCollection *p_data;

private:


  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar  & boost::serialization::base_object<gridpack::ymatrix::YMBus>(*this)
      & p_shunt_gs
      & p_shunt_bs
      & p_shunt
      & p_load
      & p_mode
      & p_ignore
      & p_v & p_a & p_theta
      & p_ybusr & p_ybusi
      & p_P0 & p_Q0
      & p_angle & p_voltage
      & p_pg & p_qg & p_pFac & p_qmin & p_qmax
      & p_qmin_orig & p_qmax_orig & p_pFac_orig
      & p_gstatus
      & p_vs & p_gid
      & p_pt & p_pb
      & p_pl & p_ql & p_ip & p_iq & p_yp & p_yq
      & p_savePl & p_saveQl
      & p_lstatus & p_lid
      & p_sbase
      & p_Pinj & p_Qinj
      & p_vmin & p_vmax
      & p_isPV
      & p_saveisPV
      & p_ngen & p_type & p_nload
      & p_area & p_zone
      & p_source & p_sink
      & p_rtpr_scale;
  }  

};

class PFBranch
  : public gridpack::ymatrix::YMBranch {
  public:
    /**
     *  Simple constructor
     */
    PFBranch(void);

    /**
     *  Simple destructor
     */
    ~PFBranch(void);

    /**
     * Return size of off-diagonal matrix block contributed by the component
     * for the forward/reverse directions
     * @param isize, jsize: number of rows and columns of matrix block
     * @return: false if network component does not contribute matrix element
     */
    bool matrixForwardSize(int *isize, int *jsize) const;
    bool matrixReverseSize(int *isize, int *jsize) const;

    /**
     * Return the values of the forward/reverse matrix block. The values are
     * returned in row-major order
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute matrix element
     */
    bool matrixForwardValues(ComplexType *values);
    bool matrixReverseValues(ComplexType *values);
    bool matrixForwardValues(RealType *values);
    bool matrixReverseValues(RealType *values);

    /**
     * Set values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    void setYBus(void);

    /**
     * Get values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    gridpack::ComplexType getYBus(void);

    /**
     * Load values stored in DataCollection object into PFBranch object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       branch that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);

    /**
     * Return the contribution to the Jacobian for the powerflow equations from
     * a branch
     * @param bus: pointer to the bus making the call
     * @param values: an array of 4 doubles that holds return metrix elements
     * @return: contribution to Jacobian matrix from branch
     */
    void getJacobian(PFBus *bus, double *values);

    /**
     * Return contribution to constraints
     * @param p: real part of constraint
     * @param q: imaginary part of constraint
     */
    void getPQ(PFBus *bus, double *p, double *q);

    /**
     * Set the mode to control what matrices and vectors are built when using
     * the mapper
     * @param mode: enumerated constant for different modes
     */
    void setMode(int mode);

    /**
     * Return complex power for line element
     * @param tag describing line element on branch
     * @return complex power
     */
    gridpack::ComplexType getComplexPower(std::string tag);

    /**
     * Write output from branches to standard out
     * @param string (output) string with information to be printed out
     * @param bufsize size of string buffer in bytes
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if branch is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const int bufsize, const char *signal = NULL);

    /**
     * Get the status of the branch element
     * @param tag character string identifying branch element
     * @return status of branch element
     */
    bool getBranchStatus(std::string tag);

    /**
     * Set the status of the branch element
     * @param tag character string identifying branch element
     * @param status status of branch element
     */
    void setBranchStatus(std::string tag, bool status);

    /**
     * get branch rating A value
     * @param tag transmission element ID
     * @return branch rating value
     */
    double getBranchRatingA(std::string tag);

    /**
     * get branch rating B value
     * @param tag transmission element ID
     * @return branch rating value
     */
    double getBranchRatingB(std::string tag);

    /**
     * get branch rating C value
     * @param tag transmission element ID
     * @return branch rating value
     */
    double getBranchRatingC(std::string tag);

    /**
     * Get list of line IDs
     * @return list of line identifiers
     */
    std::vector<std::string> getLineIDs();

    /**
     * Set parameter to ignore voltage violations
     * @param tag identifier of line element
     * @param flag value of ignore parameter
     */
    void setIgnore(std::string tag, bool flag);

    /**
     * Get parameter to ignore voltage violations
     * @param tag identifier of line element
     * @return value of ignore parameter
     */
    bool getIgnore(std::string tag);

    /**
     * Evaluate off-diagonal block of Jacobian for power flow calculation
     * and return result as an array of real values
     * @param rvals values of Jacobian block
     * @return number of values returned
     */
    int forwardJacobianValues(double *rvals);
    int reverseJacobianValues(double *rvals);

  private:
    std::vector<bool> p_ignore;
    std::vector<double> p_reactance;
    std::vector<double> p_resistance;
    std::vector<double> p_tap_ratio;
    std::vector<double> p_phase_shift;
    std::vector<double> p_charging;
    std::vector<double> p_shunt_admt_g1;
    std::vector<double> p_shunt_admt_b1;
    std::vector<double> p_shunt_admt_g2;
    std::vector<double> p_shunt_admt_b2;
    std::vector<bool> p_xform, p_shunt;
    std::vector<double> p_rateA;
    std::vector<double> p_rateB;
    std::vector<double> p_rateC;
    std::vector<bool> p_branch_status;
    std::vector<std::string> p_ckt;
    int p_mode;
    double p_ybusr_frwd, p_ybusi_frwd;
    double p_ybusr_rvrs, p_ybusi_rvrs;
    double p_theta;
    double p_sbase;
    int p_elems;
    bool p_active;

private:


  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar  & boost::serialization::base_object<gridpack::ymatrix::YMBranch>(*this)
      & p_ignore
      & p_reactance
      & p_resistance
      & p_tap_ratio
      & p_phase_shift
      & p_charging
      & p_shunt_admt_g1
      & p_shunt_admt_b1
      & p_shunt_admt_g2
      & p_shunt_admt_b2
      & p_xform & p_shunt 
      & p_rateA
      & p_branch_status
      & p_ckt
      & p_mode
      & p_ybusr_frwd & p_ybusi_frwd
      & p_ybusr_rvrs & p_ybusi_rvrs
      & p_theta
      & p_sbase
      & p_elems
      & p_active;
  }  

};


}     // powerflow
}     // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::powerflow::PFBus)
BOOST_CLASS_EXPORT_KEY(gridpack::powerflow::PFBranch)


#endif
