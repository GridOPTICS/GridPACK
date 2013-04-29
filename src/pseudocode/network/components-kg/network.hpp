// -------------------------------------------------------------
/**
 * @file   network.hpp
 * @author Kevin A. Glass
 * @date   Fri Apr  19 13:36:28 2013
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  19, 2013 by Kevin A. Glass
// Last Change:
// -------------------------------------------------------------

#ifndef _network_hpp_
#define _network_hpp_

namespace gridpack {
namespace network {

class MatrixInterface;
#include <iostream>
#include <vector>

// -------------------------------------------------------------
//  class PFNetwork
// -------------------------------------------------------------
/**
 * 
 * To be safe, these should be used simultaneously on all processes.  
 */

class PFNetwork {
public:

  /// Default constructor.
    PFNetwork() : nBuses(0), nBranches(0), nMeasurements(0){};
    {
    }

  /// Destructor
  virtual ~PFNetwork(void){
      // delete components in buses, branches and measurements
  };

  // factory calls this method, it assumes the bus-branch connections ar already made
  virtual void addBranch(PFBranchComponent  * branch)
  {
      ++nBranches;
      branches_.push_back(branch);
  }

  virtual void addBuses(PFBusComponent  * bus) {
      ++nBuses;
      buses_.push_back(bus);
  }

  virtual void addMeasurement(PFMeasurementComponent  * measurement) {
      ++nMeasurements;
      measurements_.push_back(measurement);
  }

  // loop through components to retrieve matrix value and position information
  virtual MatrixInterface * getYInterface() const;

  // loop through components to retrieve vector value and position information
  virtual MatrixInterface * getCIInterface() const;

protected:

private:
  std::vector<PFComponent *>           buses_;
  std::vector<PFComponent *>           branches_;
  std::vector<PFComponent *>           measurements_;
  int                    nBuses;
  int                    nBranches;
  int                    nMeasurements;
};

} // namespace math
} // namespace gridpack

#endif
