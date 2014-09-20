#include <boost/smart_ptr/shared_ptr.hpp>
#include "gridpack/factory/base_factory.hpp"
#include "../component/fortran_component.hpp"

typedef gridpack::network::BaseNetwork<
  gridpack::fortran_component::FortranBusComponent,
  gridpack::fortran_component::FortranBranchComponent>
  FortranNetwork;

typedef gridpack::factory::BaseFactory<FortranNetwork> FortranFactory;

struct networkWrapper {
  boost::shared_ptr<FortranNetwork> network;
};

/**
 * Create a new factory
 * @param factory pointer to Fortran factory object
 * @param network pointer to Fortran network object
 */
extern "C" void factory_create(FortranFactory **factory,
    networkWrapper *wnetwork)
{

  *factory = new FortranFactory(wnetwork->network);
}

/**
 * Destroy old factory
 * @param factory pointer to Fortran factory object
 */
extern "C" void factory_destroy(FortranFactory **factory)
{
  delete (*factory);
  *factory = NULL;
}

/**
 * Set pointers in each bus and branch component so that it points to
 * connected buses and branches. This routine operates on the generic
 * BaseBusComponent and BaseBranchComponent interfaces. It also sets some
 * indices in MatVecInterface for each component.
 * @param factory pointer to Fortran factory object
 */
extern "C" void factory_set_components(FortranFactory *factory)
{
  factory->setComponents();
}

/**
 * Generic method that invokes the "load" method on all branches and buses
 * to move data from the DataCollection objects on the network into the
 * corresponding buses and branches
 * @param factory pointer to Fortran factory object
 */
extern "C" void factory_load(FortranFactory *factory)
{
  factory->load();
}

/**
 * Set up the exchange buffers so that they work correctly. This should only
 * be called after the network topology has been specified
 * @param factory pointer to Fortran factory object
 */
extern "C" void factory_set_exchange(FortranFactory *factory)
{
  factory->setExchange();
}

/**
 * Set the mode for all BaseComponent objects in the network.
 * @param factory pointer to Fortran factory object
 * @param mode integer representing desired mode
 */
extern "C" void factory_set_mode(FortranFactory *factory, int mode)
{
  factory->setMode(mode);
}

/**
 * A convenience function that checks to see if something is true on all
 * processors
 * @param factory pointer to Fortran factory object
 * @param flag boolean flag on each processor
 * @return true if flag is true on all processors, false otherwise
 */
extern "C" bool factory_check_true(FortranFactory *factory, bool flag)
{
  return factory->checkTrue(flag);
}
