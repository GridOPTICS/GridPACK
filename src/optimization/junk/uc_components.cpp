/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#include <vector>
#include <iostream>
#include <cstring>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "uc_components.hpp"
#include "gridpack/parser/dictionary.hpp"

/**
 *  Simple constructor
 */
gridpack::ucCommitment::UCBus::UCBus(void)
{
//
  p_numUnits = 0;
  p_numHorizons = 0;
}

/**
 *  Simple destructor
 */
gridpack::ucCommitment::UCBus::~UCBus(void)
{
//
}


double objectiveFunction(void) 
{
    double obj;
    obj = 2.0;
    return obj;
}

/**
     * Return the size of the buffer used in data exchanges on the network.
     * For this problem, the number of plant units need to be exchanged
     * @return size of buffer
*/
int gridpack::ucCommitment::UCBus::getXCBufSize(void)
{
  return numUnits*sizeof(double);
}
/**
     * Assign pointers for plant units
*/
void gridpack::powerflow::PFBus::setXCBuf(void *buf)
{
}

/**
 * Load values stored in DataCollection object into YMBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::ucCommitment::UCBus::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  data->getValue(GENERATOR_NUMBER, &p_numUnits);

}

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

};

class UCBranch
{
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
    

