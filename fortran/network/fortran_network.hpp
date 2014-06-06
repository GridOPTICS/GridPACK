#include <boost/smart_ptr/shared_ptr.hpp>
#include "gridpack/network/base_network.hpp"
#include "../component/fortran_component.hpp"

typedef gridpack::network::BaseNetwork<gridpack::fortran_component::FortranBusComponent,
        gridpack::fortran_component::FortranBranchComponent > FortranNetwork;

#define MAX_NETWORKS 5

struct p_network {
  bool active;
  boost::shared_ptr<FortranNetwork> network;
};

std::vector<p_network> p_networks;
