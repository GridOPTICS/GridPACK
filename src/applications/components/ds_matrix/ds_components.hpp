/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_components.hpp
 * @author Shuangshuang Jin 
 * @date   2016-07-14 13:42:39 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _ds_components_h_
#define _ds_components_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/applications/components/y_matrix/ymatrix_components.hpp"

namespace gridpack {
namespace dynamic_simulation_r {

enum DSMode{YBUS, YA, YB, YC, updateYbus,
            DAE_init, onFY, posFY, Eprime0, Eprime1, Current};

class DSBus
  : public gridpack::ymatrix::YMBus
{
  public:
    /**
     *  Simple constructor
     */
    DSBus(void);

    /**
     *  Simple destructor
     */
    ~DSBus(void);

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

    void setValues(ComplexType *values);

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
     * Load values stored in DataCollection object into DSBus object. The
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
     * Return the value of the voltage magnitude on this bus
     * @return: voltage magnitude
     */
    double getVoltage(void);

    /**
     * Return the value of the phase angle on this bus
     * @return: phase angle
     */
    double getPhase(void);

    /**
     * Return the number of generators on this bus
     * @return number of generators on bus
     */
    int getNumGen(void);

    /**
     * Return whether or not a bus is isolated
     * @return true if bus is isolated
     */
    bool isIsolated(void) const;

    /**
     * Set values of the IFunction on this bus (gen)
     */
    void setIFunc(void);

    /**
     * Set values of the IJaco on this bus (gen)
     */
    void setIJaco(void);

    /**
     * Check to see if a fault event applies to this bus and set an internal
     * flag marking the bus as the "from" or "to" bus for the event
     * @param from_idx index of "from" bus for fault event
     * @param to_idx index of "to" bus for fault event
     * @param branch_ptr pointer to branch on which fault occurs
     */
    void setEvent(int from_idx, int to_idx,
        gridpack::component::BaseBranchComponent* branch_ptr);

    /**
     * Clear fault event from bus
     */
    void clearEvent();

    /**
     * Write output from buses to standard out
     * @param string (output) string with information to be printed out
     * @param bufsize size of string buffer in bytes
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if bus is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const int bufsize, const char *signal);

    /**
     * Set an internal parameter that specifies that the rotor speed and angle
     * for the generator corresponding to the string tag are to be printed to
     * output
     * @param tag 2-character identifier of generator
     * @param flag set to true to monitor generator
     */
    void setWatch(std::string tag, bool flag);

    /**
     * Initialize dynamic simulation data structures
     */
    void setDSParams();

    /**
     * Evaluate first part of a dynamic simulation step
     * @param flag false if step is not initial step
     */
    void initDSStep(bool flag);

    /**
     * Evaluate predictor part of dynamic simulation step
     * @param t_inc time increment
     */
    void predDSStep(double t_inc);

    /**
     * Evaluate corrector part of dynamic simulation step
     * @param t_inc time increment
     */
    void corrDSStep(double t_inc);

  private:
    double p_shunt_gs;
    double p_shunt_bs;
    bool p_shunt;
    int p_mode;
    double p_theta; // phase angle difference
    double p_ybusr, p_ybusi;
    double p_angle, p_voltage;
    bool p_load;
    double p_pl, p_ql;
    double p_sbase;
    bool p_isGen;
    std::vector<double> p_pg, p_qg;
    std::vector<int> p_gstatus;
    std::vector<double> p_mva, p_r, p_dtr;
    int p_ngen;
    int p_type;
    gridpack::ComplexType p_permYmod;
    bool p_from_flag, p_to_flag;
    std::vector<std::string> p_genid;
    std::vector<bool> p_watch;
    double p_time;

    // DAE related variables
    //double user_eqprime, user_pmech, user_gen_d0, user_gen_h; // User app context variables
    //int user_ngen; // User app context variables
    std::vector<double> p_h, p_d0;
    //std::vector<double> x, xdot; // DAE variables

    std::vector<gridpack::ComplexType> p_pelect, p_eprime_s0, p_eprime_s1;
    std::vector<gridpack::ComplexType> p_mac_ang_s0, p_mac_spd_s0;
    std::vector<gridpack::ComplexType> p_mac_ang_s1, p_mac_spd_s1;
    std::vector<gridpack::ComplexType> p_dmac_ang_s0, p_dmac_spd_s0;
    std::vector<gridpack::ComplexType> p_dmac_ang_s1, p_dmac_spd_s1;
    std::vector<gridpack::ComplexType> p_eqprime, p_mech, p_curr;

    gridpack::component::BaseBranchComponent* p_branch;

    friend class boost::serialization::access;

    template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar &
          boost::serialization::base_object<gridpack::ymatrix::YMBus>(*this)
          & p_shunt_gs
          & p_shunt_bs
          & p_shunt
          & p_mode
          & p_theta
          & p_ybusr & p_ybusi
          & p_angle & p_voltage
          & p_load
          & p_pl & p_ql
          & p_sbase
          & p_isGen
          & p_pg & p_qg
          & p_gstatus
          & p_mva & p_r & p_dtr
          & p_ngen & p_type & p_permYmod
          & p_from_flag & p_to_flag
          & p_h & p_d0
          & p_pelect
          & p_pelect & p_eprime_s0 & p_eprime_s1
          & p_mac_ang_s0 & p_mac_spd_s0
          & p_mac_ang_s1 & p_mac_spd_s1
          & p_dmac_ang_s1 & p_dmac_spd_s1
          & p_eqprime & p_mech & p_curr
          & p_genid
          & p_watch
          & p_time;
      }

};

class DSBranch
  : public gridpack::ymatrix::YMBranch {
  public:
    // Small utility structure to encapsulate information about fault events
    struct Event{
      std::string name;
      double start,end;
      int from_idx, to_idx;
      double step;
    };

    /**
     *  Simple constructor
     */
    DSBranch(void);

    /**
     *  Simple destructor
     */
    ~DSBranch(void);

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

    /**
     * Set values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    void setYBus(void);

    /**
     * Load values stored in DataCollection object into DSBranch object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       branch that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);

    /**
     * Return the complex admittance of the branch
     * @return: complex addmittance of branch
     */
    gridpack::ComplexType getAdmittance(void);

    /**
     * Return transformer contribution from the branch to the calling
     * bus
     * @param bus: pointer to the bus making the call
     * @return: contribution from transformers to Y matrix
     */
    gridpack::ComplexType getTransformer(DSBus *bus);

    /**
     * Return the contribution to a bus from shunts
     * @param bus: pointer to the bus making the call
     * @return: contribution to Y matrix from shunts associated with branches
     */
    gridpack::ComplexType getShunt(DSBus *bus);

    /**
     * Set the mode to control what matrices and vectors are built when using
     * the mapper
     * @param mode: enumerated constant for different modes
     */
    void setMode(int mode);

    /**
     * Return the updating factor that will be applied to the ybus matrix at
     * the clear fault phase
     * @return: value of update factor
     */
    gridpack::ComplexType getPosfy11YbusUpdateFactor(int sw2_2, int sw3_2);
    gridpack::ComplexType getUpdateFactor();

    /**
     * Check to see if an event applies to this branch and set appropriate
     * internal parameters
     * @param event a struct containing parameters that describe a fault
     * event in a dyanamic simulation
     */
    void setEvent(const Event &event);

  private:
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
    int p_mode;
    double p_ybusr_frwd, p_ybusi_frwd;
    double p_ybusr_rvrs, p_ybusi_rvrs;
    double p_theta;
    double p_sbase;
    std::vector<int> p_branch_status;
    int p_elems;
    bool p_active;
    bool p_event;

    friend class boost::serialization::access;

    template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar &
          boost::serialization::base_object<gridpack::ymatrix::YMBranch>(*this)
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
          & p_mode
          & p_ybusr_frwd & p_ybusi_frwd
          & p_ybusr_rvrs & p_ybusi_rvrs
          & p_theta & p_sbase
          & p_branch_status
          & p_elems & p_active & p_event;
      }
};

/// The type of network used in the dynamic_simulation application
typedef network::BaseNetwork<DSBus, DSBranch > DSNetwork;

}     // dynamic_simulation
}     // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::dynamic_simulation_r::DSBus)
BOOST_CLASS_EXPORT_KEY(gridpack::dynamic_simulation_r::DSBranch)

#endif
