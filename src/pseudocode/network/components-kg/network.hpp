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
    PFNetwork();
  /// Destructor
  virtual ~PFNetwork(void){
      // delete components in buses, branches and measurements
  };


  virtual void getYMatrix(math::Matrix * matrix);
  virtual void getHMatrix(math::Matrix * matrix);

  }
protected:
  virtual void getMatrixSize(network::ComponentNetwork & network) {
      BusCountVisitor        visitor;
      while()
  }
private:
  std::vector<PFComponent *>           buses;
  std::vector<PFComponent *>           branches;
  std::vector<PFComponent *>           measurements;
  int                                nRows;
  int                                nCols;
};

} // namespace math
} // namespace gridpack

#endif
