// -------------------------------------------------------------
/**
 * @file   test_network.hpp
 * @author Bruce Palmer
 * @date   May 31, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _test_network_h_
#define _test_network_h_

// -------------------------------------------------------------
//  class TestNetwork:
//  Trivial class that instantiates network so we can check for compile time
//  errors in BaseNetwork class.
// -------------------------------------------------------------
namespace gridpack {
class TestNetwork {
public:

  /**
   * Default constructor.
   */
  TestNetwork(void);

  /**
   * Default destructor.
   */
  ~TestNetwork(void);

};
}  //namespace gridpack

#endif
