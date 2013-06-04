// -------------------------------------------------------------
/**
 * @file   pf_components.hpp
 * @author Bruce Palmer
 * @date   June 4, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

/**
 * Some preprocessor string declarations. These will need to be put in an
 * include file someplace else. Just declare them here for the time being.
 */

#define BRANCH_REACTANCE   "branch_reactance"
#define BRANCH_RESISTANCE  "branch_resistance"
#define BRANCH_TAP_RATIO   "branch_tap_ratio"
#define BRANCH_PHASE_SHIFT "branch_phase_shift"
#define BRANCH_CHARGING    "branch_charging"
#define BRANCH_SHUNT_ADMTTNC_G1 "branch_shunt_admttnc_g1"
#define BRANCH_SHUNT_ADMTTNC_B1 "branch_shunt_admttnc_b1"
#define BRANCH_SHUNT_ADMTTNC_G2 "branch_shunt_admttnc_g2"
#define BRANCH_SHUNT_ADMTTNC_B2 "branch_shunt_admttnc_b2"
#define BUS_SHUNT_GS    "branch_shunt_gs"
#define BUS_SHUNT_BS    "branch_shunt_bs"

#ifdef _pf_components_h_
#define _pf_components_h_

#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"

namespace gridpack {
namespace powerflow {

class PFBus
  : public BaseBusComponent {
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
     *  Return size of matrix block contributed by the component
     *  @param isize, jsize: number of rows and columns of matrix block
     *  @return: false if network component does not contribute matrix element
     */
    bool matrixSize(int *isize, int *jsize) const;

    /**
     * Return the values of the matrix block. The values are
     * returned in row-major order.
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute matrix element
     */
    bool matrixValues(void *values);

    /**
     * Load values stored in DataCollection object into PFBus object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       bus that were read in when network was initialized
     */
    void load(shared_ptr<gridpack::component::DataCollection> data);

  private:
    double p_shunt_gs;
    double p_shunt_bs;
    bool p_shunt;

};

class PFBranch
  : public BaseBranchComponent {
    /**
     *  Simple constructor
     */
    PFBranch(void);

    /**
     *  Simple destructor
     */
    ~PFBranch(void);

    /**
     *  Return size of matrix block contributed by the component
     *  @param isize, jsize: number of rows and columns of matrix block
     *  @return: false if network component does not contribute matrix element
     */
    bool matrixSize(int *isize, int *jsize) const;

    /**
     * Return the values of the matrix block. The values are
     * returned in row-major order.
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute matrix element
     */
    bool matrixValues(void *values);

    /**
     * Load values stored in DataCollection object into PFBranch object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       branch that were read in when network was initialized
     */
    void load(shared_ptr<gridpack::component::DataCollection> data);

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
    gridpack::ComplexType getTransformer(boost::shared_ptr<PFBus> bus);

    /**
     * Return the contribution to a bus from shunts
     * @param bus: pointer to the bus making the call
     * @return: contribution to Y matrix from shunts associated with branches
     */
    gridpack::ComplexType getShunt(boost::shared_ptr<PFBus> bus);

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
};

}
}     // gridpack
#endif
