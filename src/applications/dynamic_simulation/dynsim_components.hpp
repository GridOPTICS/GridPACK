// -------------------------------------------------------------
/**
 * @file   dynsim_components.hpp
 * @author Shuangshuang Jin 
 * @date   September 19, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _dynsim_components_h_
#define _dynsim_components_h_

/**
 * Some preprocessor string declarations. These will need to be put in an
 * include file someplace else. Just declare them here for the time being.
 */

/*
#define BRANCH_REACTANCE   "branch_reactance"
#define BRANCH_RESISTANCE  "branch_resistance"
#define BRANCH_TAP_RATIO   "branch_tap_ratio"
#define BRANCH_PHASE_SHIFT "branch_phase_shift"
#define BRANCH_CHARGING    "branch_charging"
#define BUS_SHUNT_GS    "branch_shunt_gs"
#define BUS_SHUNT_BS    "branch_shunt_bs"
*/

/* These are defined in dictionary.hpp. */

/* #define BRANCH_SHUNT_ADMTTNC_G1 "branch_shunt_admttnc_g1" */
/* #define BRANCH_SHUNT_ADMTTNC_B1 "branch_shunt_admttnc_b1" */
/* #define BRANCH_SHUNT_ADMTTNC_G2 "branch_shunt_admttnc_g2" */
/* #define BRANCH_SHUNT_ADMTTNC_B2 "branch_shunt_admttnc_b2" */

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/network/base_network.hpp"

namespace gridpack {
namespace dynsim {

//enum DynSimMode{YBUS, JACOBIAN, GENERATOR};
enum DynSimMode{YL, PERM, YA, YB};

class DynSimBus
  : public gridpack::component::BaseBusComponent {
  public:
    /**
     *  Simple constructor
     */
    DynSimBus(void);

    /**
     *  Simple destructor
     */
    ~DynSimBus(void);

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
    bool matrixDiagValues(void *values);

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
    bool vectorValues(void *values);

    /**
     * Set values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    void setYBus(void);

    /**
     * Load values stored in DataCollection object into DynSimBus object. The
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
     * Return whether or not the bus is a generator
     * @return: true if bus is geneartor
     */
    bool isGen();

  private:
    //double p_shunt_gs;
    //double p_shunt_bs;
    //bool p_shunt;
    int p_mode;
    //double p_v, p_theta;
    double p_ybusr, p_ybusi;
    //double p_P0, p_Q0;
    double p_voltage;
    bool p_load;
    double p_pl, p_ql;
    double p_sbase;
    bool p_isGen;
    std::vector<double> p_pg, p_qg;
    std::vector<int> p_gstatus;
    std::vector<double> p_mva, p_r, p_dstr, p_dtr;
    int p_ngen;
    int ngen;
};

class DynSimBranch
  : public gridpack::component::BaseBranchComponent {
  public:
    /**
     *  Simple constructor
     */
    DynSimBranch(void);

    /**
     *  Simple destructor
     */
    ~DynSimBranch(void);

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
    bool matrixForwardValues(void *values);
    bool matrixReverseValues(void *values);

    /**
     * Set values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    void setYBus(void);

    /**
     * Load values stored in DataCollection object into DynSimBranch object. The
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
    //gridpack::ComplexType getAdmittance(void);

    /**
     * Return transformer contribution from the branch to the calling
     * bus
     * @param bus: pointer to the bus making the call
     * @return: contribution from transformers to Y matrix
     */
    //gridpack::ComplexType getTransformer(DynSimBus *bus);

    /**
     * Return the contribution to a bus from shunts
     * @param bus: pointer to the bus making the call
     * @return: contribution to Y matrix from shunts associated with branches
     */
    //gridpack::ComplexType getShunt(DynSimBus *bus);

    /**
     * Return the contribution to the Jacobian for the powerflow equations from
     * a branch
     * @param bus: pointer to the bus making the call
     * @param values: an array of 4 doubles that holds return metrix elements
     * @return: contribution to Jacobian matrix from branch
     */
    //void getJacobian(DynSimBus *bus, double *values);

    /**
     * Return contribution to constraints
     * @param p: real part of constraint
     * @param q: imaginary part of constraint
     */
    //void getPQ(DynSimBus *bus, double *p, double *q);

    /**
     * Set the mode to control what matrices and vectors are built when using
     * the mapper
     * @param mode: enumerated constant for different modes
     */
    //void setMode(int mode);

  private:
    double p_reactance;
    double p_resistance;
    double p_tap_ratio;
    double p_phase_shift;
    double p_charging;
    double p_shunt_admt_g1;
    double p_shunt_admt_b1;
    double p_shunt_admt_g2;
    double p_shunt_admt_b2;
    bool p_xform, p_shunt;
    int p_mode;
    double p_ybusr, p_ybusi;
    double p_theta;
};

}     // dynsim
}     // gridpack
#endif
