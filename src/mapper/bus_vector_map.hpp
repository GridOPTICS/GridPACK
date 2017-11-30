/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/* *****************************************************************************
 * bus_vector_map.hpp
 * gridpack
 * kglass, bjpalmer
 * Jul 22, 2013
 *
 * Documentation status
 * Revision status
 * File contents
 **************************************************************************** */

#ifndef BUSVECTORMAP_HPP_
#define BUSVECTORMAP_HPP_

#include <boost/smart_ptr/shared_ptr.hpp>
#include <ga.h>
#include "gridpack/parallel/parallel.hpp"
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/component/base_component.hpp>
#include <gridpack/network/base_network.hpp>
#include <gridpack/math/vector.hpp>

//#define DBG_CHECK

namespace gridpack {
namespace mapper {

template <class _network>
class BusVectorMap {
  public: 
/**
 * Initialize mapper for the given network and the current mode. Create global
 * arrays that contain offsets that will be used to create vector from the
 * network component objects
 * @param network network that will generate vector
 */
BusVectorMap(boost::shared_ptr<_network> network)
  : p_network(network)
{
  p_Offsets                        = NULL;
  p_ISize                          = NULL;
  int                     iSize    = 0;
  p_contributingBuses              = NULL;
  p_Indices                        = NULL;

  p_timer = NULL;
  p_timer = gridpack::utility::CoarseTimer::instance();

  p_GAgrp = network->communicator().getGroup();
  p_me = GA_Pgroup_nodeid(p_GAgrp);
  p_nNodes = GA_Pgroup_nnodes(p_GAgrp);


  p_nBuses = p_network->numBuses();

  contributions();

  setBusIndexArrays();
}

~BusVectorMap()
{
  if (p_Offsets != NULL) delete [] p_Offsets;
  if (p_ISize != NULL) delete [] p_ISize;
  if (p_contributingBuses != NULL) delete [] p_contributingBuses;
  if (p_Indices != NULL) delete [] p_Indices;
}

/**
 * Create a vector from the current bus state
 * @return return pointer to new vector
 */
boost::shared_ptr<gridpack::math::Vector> mapToVector(void)
{
  gridpack::parallel::Communicator comm = p_network->communicator();
  int t_new(0), t_bus(0), t_set(0);
  if (p_timer) t_new = p_timer->createCategory("Vector Map: New Vector");
  if (p_timer) p_timer->start(t_new);
  boost::shared_ptr<gridpack::math::Vector>
             Ret(new gridpack::math::Vector(comm, p_numValues));
  if (p_timer) p_timer->stop(t_new);
  if (p_timer) t_bus = p_timer->createCategory("Vector Map: Load Bus Data");
  if (p_timer) p_timer->start(t_bus);
  loadBusData(*Ret,true);
  if (p_timer) p_timer->stop(t_bus);
  if (p_timer) t_set = p_timer->createCategory("Vector Map: Set Vector");
  if (p_timer) p_timer->start(t_set);
  Ret->ready();
  if (p_timer) p_timer->stop(t_set);
  return Ret;
}

/**
 * Create a vector from the current bus state
 * @return return pointer to new vector
 */
boost::shared_ptr<gridpack::math::RealVector> mapToRealVector(void)
{
  gridpack::parallel::Communicator comm = p_network->communicator();
  int t_new(0), t_bus(0), t_set(0);
  if (p_timer) t_new = p_timer->createCategory("Vector Map: New Vector");
  if (p_timer) p_timer->start(t_new);
  boost::shared_ptr<gridpack::math::RealVector>
             Ret(new gridpack::math::RealVector(comm, p_numValues));
  if (p_timer) p_timer->stop(t_new);
  if (p_timer) t_bus = p_timer->createCategory("Vector Map: Load Bus Data");
  if (p_timer) p_timer->start(t_bus);
  loadRealBusData(*Ret,true);
  if (p_timer) p_timer->stop(t_bus);
  if (p_timer) t_set = p_timer->createCategory("Vector Map: Set Vector");
  if (p_timer) p_timer->start(t_set);
  Ret->ready();
  if (p_timer) p_timer->stop(t_set);
  return Ret;
}

/**
 * Create a vector from the current bus state and return a conventional pointer
 * to it. Used for Fortran interface
 * @return return pointer to new vector
 */
gridpack::math::Vector* intMapToVector(void)
{
  gridpack::parallel::Communicator comm = p_network->communicator();
  int t_new, t_bus, t_set;
  if (p_timer) t_new = p_timer->createCategory("Vector Map: New Vector");
  if (p_timer) p_timer->start(t_new);
  gridpack::math::Vector*
     Ret(new gridpack::math::Vector(comm, p_numValues));
  if (p_timer) p_timer->stop(t_new);
  if (p_timer) t_bus = p_timer->createCategory("Vector Map: Load Bus Data");
  if (p_timer) p_timer->start(t_bus);
  loadBusData(*Ret,true);
  if (p_timer) p_timer->stop(t_bus);
  if (p_timer) t_set = p_timer->createCategory("Vector Map: Set Vector");
  if (p_timer) p_timer->start(t_set);
  Ret->ready();
  if (p_timer) p_timer->stop(t_set);
  return Ret;
}

/**
 * Reset a vector from the current bus state (vector should be created with same
 * mapper)
 * @param vector existing vector that should be reset
 */
void mapToVector(gridpack::math::Vector &vector)
{
  int t_bus(0), t_set(0);
  if (p_timer) t_set = p_timer->createCategory("Vector Map: Set Vector");
  if (p_timer) p_timer->start(t_set);
  vector.zero();
  if (p_timer) p_timer->stop(t_set);
  if (p_timer) t_bus = p_timer->createCategory("Vector Map: Load Bus Data");
  if (p_timer) p_timer->start(t_bus);
  loadBusData(vector,false);
  if (p_timer) p_timer->stop(t_bus);
  if (p_timer) p_timer->start(t_set);
  vector.ready();
  if (p_timer) p_timer->stop(t_set);
}

/**
 * Reset a vector from the current bus state (vector should be created with same
 * mapper)
 * @param vector existing vector that should be reset
 */
void mapToRealVector(gridpack::math::RealVector &vector)
{
  int t_bus(0), t_set(0);
  if (p_timer) t_set = p_timer->createCategory("Vector Map: Set Vector");
  if (p_timer) p_timer->start(t_set);
  vector.zero();
  if (p_timer) p_timer->stop(t_set);
  if (p_timer) t_bus = p_timer->createCategory("Vector Map: Load Bus Data");
  if (p_timer) p_timer->start(t_bus);
  loadRealBusData(vector,false);
  if (p_timer) p_timer->stop(t_bus);
  if (p_timer) p_timer->start(t_set);
  vector.ready();
  if (p_timer) p_timer->stop(t_set);
}

/**
 * Reset a vector from the current bus state (vector should be created with same
 * mapper). 
 * @param vector existing vector that should be reset
 */
void mapToVector(boost::shared_ptr<gridpack::math::Vector> &vector)
{
  mapToVector(*vector);
}

/**
 * Reset a vector from the current bus state (vector should be created with same
 * mapper). 
 * @param vector existing vector that should be reset
 */
void mapToRealVector(boost::shared_ptr<gridpack::math::RealVector> &vector)
{
  mapToRealVector(*vector);
}

/**
 * Push data from vector onto buses. Vector must be created with the
 * mapToVector method using the same BusVectorMap
 * @param vector vector containing data to be pushed to buses
 */
void mapToBus(const gridpack::math::Vector &vector)
{
  ComplexType *values = new ComplexType[p_numValues];
  int t_get, t_unpack;
  int i, j, isize, offset;
  int one = 1;
  if (p_timer) t_get = p_timer->createCategory("mapToBus: get Data");
  if (p_timer) p_timer->start(t_get);
  vector.getElements(p_numValues, p_Indices, values);
  if (p_timer) p_timer->stop(t_get);
  ComplexType *vptr = values;
  if (p_timer) t_unpack = p_timer->createCategory("mapToBus: set Data");
  if (p_timer) p_timer->start(t_unpack);
  for (i=0; i<p_busContribution; i++) {
    p_contributingBuses[i]->setValues(vptr);
    isize = p_ISize[i];
    vptr += isize;
  }
  if (p_timer) p_timer->stop(t_unpack);
  delete [] values;
}

/**
 * Push data from vector onto buses. Vector must be created with the
 * mapToVector method using the same BusVectorMap
 * @param vector vector containing data to be pushed to buses
 */
void mapToBus(const gridpack::math::RealVector &vector)
{
  RealType *values = new RealType[p_numValues];
  int t_get, t_unpack;
  int i, j, isize, offset;
  int one = 1;
  if (p_timer) t_get = p_timer->createCategory("mapToBus: get Data");
  if (p_timer) p_timer->start(t_get);
  vector.getElements(p_numValues, p_Indices, values);
  if (p_timer) p_timer->stop(t_get);
  RealType *vptr = values;
  if (p_timer) t_unpack = p_timer->createCategory("mapToBus: set Data");
  if (p_timer) p_timer->start(t_unpack);
  for (i=0; i<p_busContribution; i++) {
    p_contributingBuses[i]->setValues(vptr);
    isize = p_ISize[i];
    vptr += isize;
  }
  if (p_timer) p_timer->stop(t_unpack);
  delete [] values;
}

/**
 * Push data from vector onto buses. Vector must be created with the
 * mapToVector method using the same BusVectorMap
 * @param vector vector containing data to be pushed to buses
 */
void mapToBus(boost::shared_ptr<gridpack::math::Vector> &vector)
{
  mapToBus(*vector);
}

void mapToBus(boost::shared_ptr<gridpack::math::RealVector> &vector)
{
  mapToBus(*vector);
}

private:
/**
 * Add block contributions from buses to vector
 * @param vector vector to which contributions are added
 * @param flag flag to distinguish new vector (true) from existing vector * (false)
 */
void loadBusData(gridpack::math::Vector &vector, bool flag)
{
  int i,idx,isize,icnt;
  // Add vector elements
  boost::shared_ptr<gridpack::component::BaseBusComponent> bus;
  int t_bus(0);
  if (p_timer) t_bus = p_timer->createCategory("loadBusData: Add Vector Elements");
  if (p_timer) p_timer->start(t_bus);
  int j;
  int t_pack, t_add;
  if (p_timer) t_pack = p_timer->createCategory("loadBusData: Fill Buffer");
  if (p_timer) p_timer->start(t_pack);
  ComplexType *vbuf = new ComplexType[p_numValues];
  ComplexType *vptr = vbuf;
  int *ibuf = new int[p_numValues];
  icnt = 0;
  for (i=0; i<p_busContribution; i++) {
    p_contributingBuses[i]->vectorValues(vptr);
    isize = p_ISize[i];
    idx = p_Offsets[i];
    for (j=0; j<isize; j++) {
      ibuf[icnt] = idx;
      idx++;
      icnt++;
    }
    vptr += isize;
  }
  if (p_timer) p_timer->stop(t_pack);
  if (p_timer) t_add = p_timer->createCategory("loadBusData: Add Elements");
  if (p_timer) p_timer->start(t_add);
  vector.addElements(p_numValues,ibuf,vbuf);
  if (p_timer) p_timer->stop(t_add);
  delete [] vbuf;
  delete [] ibuf;
  if (p_timer) p_timer->stop(t_bus);
}

/**
 * Add block contributions from buses to vector
 * @param vector vector to which contributions are added
 * @param flag flag to distinguish new vector (true) from existing vector * (false)
 */
void loadRealBusData(gridpack::math::RealVector &vector, bool flag)
{
  int i,idx,isize,icnt;
  // Add vector elements
  boost::shared_ptr<gridpack::component::BaseBusComponent> bus;
  int t_bus(0);
  if (p_timer) t_bus = p_timer->createCategory("loadBusData: Add Vector Elements");
  if (p_timer) p_timer->start(t_bus);
  int j;
  int t_pack, t_add;
  if (p_timer) t_pack = p_timer->createCategory("loadBusData: Fill Buffer");
  if (p_timer) p_timer->start(t_pack);
  RealType *vbuf = new RealType[p_numValues];
  RealType *vptr = vbuf;
  int *ibuf = new int[p_numValues];
  icnt = 0;
  for (i=0; i<p_busContribution; i++) {
    p_contributingBuses[i]->vectorValues(vptr);
    isize = p_ISize[i];
    idx = p_Offsets[i];
    for (j=0; j<isize; j++) {
      ibuf[icnt] = idx;
      idx++;
      icnt++;
    }
    vptr += isize;
  }
  if (p_timer) p_timer->stop(t_pack);
  if (p_timer) t_add = p_timer->createCategory("loadBusData: Add Elements");
  if (p_timer) p_timer->start(t_add);
  vector.addElements(p_numValues,ibuf,vbuf);
  if (p_timer) p_timer->stop(t_add);
  delete [] vbuf;
  delete [] ibuf;
  if (p_timer) p_timer->stop(t_bus);
}

/**
 * Add block contributions from buses to vector
 * @param vector vector to which contributions are added
 * @param flag flag to distinguish new vector (true) from existing vector * (false)
 */
void loadBusData(boost::shared_ptr<gridpack::math::Vector> &vector, bool flag)
{
  loadBusData(*vector, flag);
}

/**
 * Add block contributions from buses to vector
 * @param vector vector to which contributions are added
 * @param flag flag to distinguish new vector (true) from existing vector * (false)
 */
void loadRealBusData(boost::shared_ptr<gridpack::math::RealVector> &vector, bool flag)
{
  loadRealBusData(*vector, flag);
}

/**
 * Calculate how many buses contribute to vector
 */
void contributions(void)
{
  int i;
  // Get number of contributions from buses
  int isize;
  p_busContribution = 0;
  p_numValues = 0;
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      if (p_network->getBus(i)->vectorSize(&isize)) {
        p_busContribution++;
        p_numValues += isize;
      }
    }
  }
  p_contributingBuses 
    = new gridpack::component::BaseBusComponent*[p_busContribution];
  p_ISize = new int[p_busContribution];
  int icnt = 0;
  for (i=0; i<p_nBuses; i++) {
    if (p_network->getActiveBus(i)) {
      if (p_network->getBus(i)->vectorSize(&isize)) {
        p_contributingBuses[icnt] = p_network->getBus(i).get();
        p_ISize[icnt] = isize;
        icnt++;
      }
    }
  }
}

/**
 * Set up internal arrays of indices that are used in mapper
 */
void setBusIndexArrays(void)
{
  // Exchange number of values contributed by each process
  std::vector<int> nVals(p_nNodes);
  int i, j, nsize;
  for (i=0; i<p_nNodes; i++) {
    nVals[i] = 0;
  }
  nVals[p_me] = p_numValues;

  char cplus[2];
  strcpy(cplus,"+");
  GA_Pgroup_igop(p_GAgrp,&nVals[0],p_nNodes,cplus);

  //Evaluate starting index on each process;
  int offset = 0;
  for (i=0; i<p_me; i++) {
    offset += nVals[i];
  }
  
  // Loop over contributing buses and get indices
  p_Offsets = new int[p_busContribution];
  p_Indices = new int[p_numValues];
  int icnt = 0;
  int jcnt = offset;
  for (i=0; i<p_busContribution; i++) {
    nsize = p_ISize[i];
    p_Offsets[i] = jcnt;
    jcnt += nsize;
    for (j=0; j<nsize; j++) {
      p_Indices[icnt] = offset+icnt;
      icnt++;
    }
  }
}

    // GA information
int                         p_me;
int                         p_nNodes;

    // network information
boost::shared_ptr<_network> p_network;
int                         p_nBuses;

    // vector information
int                         p_busContribution;  // Number of buses contributing to vector
int                         p_numValues; // Number of values contributed to vector

int*                        p_Offsets;
int*                        p_ISize;
int*                        p_Indices;
gridpack::component::BaseBusComponent **p_contributingBuses;

    // global vector block size array
int                         p_GAgrp; // GA group

    // pointer to timer
gridpack::utility::CoarseTimer *p_timer;

};

} /* namespace mapper */
} /* namespace gridpack */

#endif //BUSVECTORMAP_HPP_
