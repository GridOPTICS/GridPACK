/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#ifndef _uc_components_h_
#define _uc_components_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"

namespace gridpack {
namespace ucCommitment {

class UCBus
  : public gridpack::component::BaseBusComponent {
  public: 
    /**
     *  Simple constructor
     */
    UCBus(void); 
    /**
     *  Simple destructor
     */
    ~UCBus(void);

    double objectiveFunction(void);

    /**
     * Return the size of the buffer used in data exchanges on the network.
     * For this problem, the number of plant units need to be exchanged
     * @return size of buffer
     */
//  void getXCBufSize(void);
    /**
     * Assign pointers for plant units
     */
//  void setXCBuf(void *buf);

    /**
     * Load values stored in DataCollection object into PFBus object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       bus that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);
    /**
     * Set values of constant cost parameters. These can then be used in subsequent
     * calculations of cost function
     */
  void setCostConst(void);
    /**
     * Get values of constant cost parameters. These can then be used in subsequent
     * calculations
     */
  void getCostConst(void);

    /**
     * Set values of linear cost parameters. These can then be used in subsequent
     * calculations of cost function
     */
  void setLinearConst(void);
    /**
     * Get values of linear cost parameters. These can then be used in subsequent
     * calculations
     */
  void getLinearConst(void);
    /**
     * Set values of quadratic cost parameters. These can then be used in subsequent
     * calculations of cost function
     */
  void setQuadConst(void);
    /**
     * Get values of quadratic cost parameters. These can then be used in subsequent
     * calculations
     */
  void getQuadConst(void);

  private:
// unit power and status
    int p_numUnits;
    int p_numHorizons;
    std::vector<double> p_up;
    std::vector<int> p_ustatus;
    std::vector<double> p_demand;
    std::vector<double> p_minPower;
    std::vector<double> p_maxPower;
    std::vector<double> p_minUptime;
    std::vector<double> p_minDntime;
    std::vector<double> p_constConst;
    std::vector<double> p_constLinear;
    std::vector<double> p_constQuad;
};

class UCBranch
  : public gridpack::component::BaseBranchComponent {
  public: 
    /** 
     *  Simple constructor
     */ 
    UCBranch(void); 
    /** 
     *  Simple destructor
     */
    UCBranch(void);
  private:

};
}     // ucCommitment
}     // gridpack
    
BOOST_CLASS_EXPORT_KEY(gridpack::ucCommitment::UCBus);
BOOST_CLASS_EXPORT_KEY(gridpack::ucCommitment::UCBranch);


#endif

