// -------------------------------------------------------------
/**
 * @file   base_factory.cpp
 * @author Bruce Palmer
 * @date   2013-08-08 10:06:51 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include "gridpack/component/base_component.hpp"
#include "gridpack/factory/base_factory.hpp"


namespace gridpack {

// just to force instantiation of a BaseFactory<>
static 
void a_function(void) {

  boost::mpi::communicator world;
  typedef network::BaseNetwork<component::BaseBusComponent,
                               component::BaseBranchComponent> ANetwork;
  typedef factory::BaseFactory<ANetwork> AFactory;
  AFactory::NetworkPtr a_network(new ANetwork(world));
  AFactory a_factory(a_network);
}

} // namespace gridpack
