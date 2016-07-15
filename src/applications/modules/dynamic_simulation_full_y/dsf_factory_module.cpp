/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsf_factory_module.cpp
 * @author Shuangshuang Jin 
 * @date   2016-07-14 13:50:07 d3g096
 * @last modified date   May 13, 2015
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "dsf_factory_module.hpp"

namespace gridpack {
namespace dynamic_simulation {

// Dynamic simulation factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
DSFullFactory::DSFullFactory(DSFullFactory::NetworkPtr network)
  : gridpack::factory::BaseFactory<DSFullNetwork>(network)
{
  p_network = network;

  int i;
  p_numBus = p_network->numBuses();
  p_numBranch = p_network->numBranches();

  p_buses = new gridpack::dynamic_simulation::DSFullBus*[p_numBus];
  for (i=0; i<p_numBus; i++) {
    p_buses[i] = dynamic_cast<DSFullBus*>(p_network->getBus(i).get());
  }
  p_branches = new gridpack::dynamic_simulation::DSFullBranch*[p_numBranch];
  for (i=0; i<p_numBranch; i++) {
    p_branches[i] = dynamic_cast<DSFullBranch*>(p_network->getBranch(i).get());
  }
#ifdef USE_FNCS
  p_hash = new gridpack::hash_distr::HashDistribution<DSFullNetwork,
               gridpack::ComplexType,gridpack::ComplexType>
               (network);
  int dsize = sizeof(gridpack::dynamic_simulation::DSFullBus::voltage_data);
  p_busIO = new gridpack::serial_io::SerialBusIO<DSFullNetwork>(dsize,
                network);
#endif
}

/**
 * Basic destructor
 */
gridpack::dynamic_simulation::DSFullFactory::~DSFullFactory()
{
  delete [] p_buses;
  delete [] p_branches;
#ifdef USE_FNCS
  delete p_busIO;
  delete p_hash;
#endif
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::dynamic_simulation::DSFullFactory::setYBus(void)
{
  int i;

  // Invoke setYBus method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->setYBus();
  }

  // Invoke setYBus method on all branch objects
  for (i=0; i<p_numBranch; i++) {
    p_branches[i]->setYBus();
  }
}

/**
 * Get the updating factor for posfy11 stage ybus
 */
gridpack::ComplexType
gridpack::dynamic_simulation::DSFullFactory::setFactor(int sw2_2, int sw3_2)
{
  gridpack::ComplexType dummy(-999.0, -999.0);

  int i;

  // Invoke getPosfy11YbusUpdateFactor method on all branch objects
  for (i=0; i<p_numBranch; i++) {
    gridpack::ComplexType ret = p_branches[i]
      ->getPosfy11YbusUpdateFactor(sw2_2, sw3_2);
    if (ret != dummy) {
      return ret;
    }
  }

  return dummy;
}

/**
 * Apply an event to all branches in the system
 * @param event a struct describing a fault
 */
void gridpack::dynamic_simulation::DSFullFactory::setEvent(const
    gridpack::dynamic_simulation::DSFullBranch::Event &event)
{
  int i;
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->clearEvent();
  }
  for (i=0; i<p_numBranch; i++) {
    p_branches[i]->setEvent(event);
  }
}

/**
 * Check network to see if there is a process with no generators
 * @return true if all processors have at least on generator
 */
bool gridpack::dynamic_simulation::DSFullFactory::checkGen(void)
{
  int i, count;
  count = 0;
  for (i=0; i<p_numBus; i++) {
    if (p_network->getActiveBus(i)) {
      count += p_buses[i]->getNumGen();
    }
  }
  //printf("p[%d] number of generators: %d\n",p_network->communicator().rank(),count);
  int iok = 0;
  if (count > 0) iok = 1;
  int ok;
  boost::mpi::all_reduce(p_network->communicator(),iok,ok,std::plus<int>());
  bool ret = false;
  if (ok == p_network->communicator().size()) ret = true;
  return ret;
}

/**
 * Initialize init vectors for integration
 * @param ts time step
 */
void gridpack::dynamic_simulation::DSFullFactory::initDSVect(double ts)
{
  int i;

  // Invoke initDSVect method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->initDSVect(ts);
  }
}

/**
 * Update vectors in each integration time step (Predictor)
 */
void gridpack::dynamic_simulation::DSFullFactory::predictor_currentInjection(bool flag)
{
  int i;

  // Invoke method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->predictor_currentInjection(flag);
  }
}

/**
 * Update vectors in each integration time step (Predictor)
 */
void gridpack::dynamic_simulation::DSFullFactory::predictor(double t_inc, bool flag)
{
  int i;

  // Invoke updateDSVect method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->predictor(t_inc,flag);
  }
}

/**
 * Update vectors in each integration time step (Corrector)
 */
void gridpack::dynamic_simulation::DSFullFactory::corrector_currentInjection(bool flag)
{
  int i;

  // Invoke method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->corrector_currentInjection(flag);
  }
}

/**
 * Update vectors in each integration time step (Corrector)
 */
void gridpack::dynamic_simulation::DSFullFactory::corrector(double t_inc, bool flag)
{
  int i;

  // Invoke updateDSVect method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->corrector(t_inc,flag);
  }
}

/**
 * Set volt from volt_full
 */
void gridpack::dynamic_simulation::DSFullFactory::setVolt(bool flag)
{
  int i;

  // Invoke updateDSVect method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->setVolt(flag);
  }
}

/**
 * Add constant impedance load admittance to diagonal elements of
 * Y-matrix
 */
void gridpack::dynamic_simulation::DSFullFactory::addLoadAdmittance()
{
  int i;

  // Invoke addLoadAdmittance method on all bus objects
  for (i=0; i<p_numBus; i++) {
    p_buses[i]->addLoadAdmittance();
  }
}

bool gridpack::dynamic_simulation::DSFullFactory::securityCheck()
{
  int i;
  bool secure;

  double maxAngle = -999999.0;
  double minAngle = 999999.0;
  double angle;
  for (i = 0; i < p_numBus; i++) {
    int nGen = p_buses[i]->getNumGen();
    if (nGen > 0) {
      angle = p_buses[i]->getAngle();
      //printf("angle = %f\n", angle);
      if (angle > maxAngle) {
        maxAngle = angle;
        //printf("   max = %f\n", maxAngle);
      }
      if (angle < minAngle) {
        minAngle = angle;
        //printf("   min = %f\n", minAngle);
      }
      double pi = 4.0*atan(1.0);
      if ((maxAngle - minAngle) > 360 * pi / 180)
        secure = false;
      else
        secure = true; 
    }
  }
  //printf("max = %f, min = %f\n", maxAngle, minAngle);
  return secure;
}

#ifdef USE_FNCS
/**
 * Scatter load from FNCS framework to buses
 */
void gridpack::dynamic_simulation::DSFullFactory::scatterLoad()
{
  std::vector<int> keys;
  std::vector<gridpack::ComplexType> values;

  if (p_network->communicator().rank() == 0) {
    /**
     * Add FNCS code here that puts bus IDs in keys and corresponding LOAD_PL and
     * LOAD_QL as the real and imaginary parts of a single complex number in
     * values
     */
  }

  // distribute loads to all processors
  p_hash->distributeBusValues(keys, values);

  // assign loads to individual processors
  int nvals = keys.size();
  int i;
  for (i=0; i<nvals; i++) {
    p_buses[keys[i]]->setLoad(real(values[i]), imag(values[i]));
  }
}

/**
 * Gather voltage magnitude and phase angle to root node and transmit to
 * FNCS framework
 */
void gridpack::dynamic_simulation::DSFullFactory::gatherVoltage()
{
  // Get data from network buses
  std::vector<gridpack::dynamic_simulation::DSFullBus::voltage_data> data;
  p_busIO->gatherData(data);
  
  // Loop over all data items
  int i;
  if (p_network->communicator().rank() == 0) {
    int nvals = data.size();
    for (i=0; i<nvals; i++) {
      int id = data[i].busID;
      gridpack::ComplexType voltage = data[i].voltage;
      /**
       * Add FNCS code here that maps complex voltage to bus id
       */
    }
  }
}

#endif

} // namespace dynamic_simulation
} // namespace gridpack
