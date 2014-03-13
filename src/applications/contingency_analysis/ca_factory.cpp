/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_factory.cpp
 * @author Yousu Chen 
 * @date   Feb 11, 2014 
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "ca_factory.hpp"
#include "gridpack/mapper/bus_vector_map.hpp"
#include "gridpack/mapper/full_map.hpp"


namespace gridpack {
namespace contingency_analysis {

// Powerflow factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
CAFactory::CAFactory(CAFactory::NetworkPtr network)
  : gridpack::factory::BaseFactory<CANetwork>(network)
{
  p_network = network;
}

/**
 * Basic destructor
 */
gridpack::contingency_analysis::CAFactory::~CAFactory()
{
}

/**
 * Reset voltages to initial values on all buses
 */
void gridpack::contingency_analysis::CAFactory::resetVoltage(void)
{
  int numBus = p_network->numBuses();
  int i;
  // Invoke resetVoltage method on all buses
  for (i=0; i<numBus; i++) {
    dynamic_cast<gridpack::contingency_analysis::CABus*>(p_network->getBus(i).get())
       ->resetVoltage();
  }
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::contingency_analysis::CAFactory::setYBus(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // Invoke setYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
    dynamic_cast<gridpack::contingency_analysis::CABranch*>(p_network->getBranch(i).get())
       ->setYBus();
  }

  // Invoke setYBus method on all bus objects
  for (i=0; i<numBus; i++) {
    dynamic_cast<gridpack::contingency_analysis::CABus*>(p_network->getBus(i).get())
       ->setYBus();
  }

}

/**
  * Make SBus vector 
  */
void gridpack::contingency_analysis::CAFactory::setSBus(void)
{
  int numBus = p_network->numBuses();
  int i;

  // Invoke setSBus method on all bus objects
  for (i=0; i<numBus; i++) {
    dynamic_cast<gridpack::contingency_analysis::CABus*>(p_network->getBus(i).get())
       ->setSBus();
  }
}

/**
  * Create the PQ 
  */
void gridpack::contingency_analysis::CAFactory::setPQ(void)
{
  int numBus = p_network->numBuses();
  int i;
  ComplexType values[2];

  for (i=0; i<numBus; i++) {
    dynamic_cast<gridpack::contingency_analysis::CABus*>(p_network->getBus(i).get())
       ->vectorValues(values);
  }
}

/**
 * Set contingency
 * @param contingency the contigency that is to be set
 */
void gridpack::contingency_analysis::CAFactory::setContingency(
    gridpack::contingency_analysis::Contingency contingency)
{
  int numBranch = p_network->numBranches();
  int i,j,idx1,idx2;
  int size = contingency.p_from.size();
  p_saveStatus.clear();
  gridpack::contingency_analysis::CABranch *branch;
  bool found = false;
  for (i=0; i<numBranch; i++) {
    branch = dynamic_cast<gridpack::contingency_analysis::CABranch*>(p_network->getBranch(i).get());
    idx1 = branch->getBus1OriginalIndex();
    idx2 = branch->getBus2OriginalIndex();
    // check branch indices against contingencies
    if (contingency.p_type == gridpack::contingency_analysis::Branch) {
      for (j = 0; j<size; j++) {
        if (contingency.p_from[j] == idx1 && contingency.p_to[j] == idx2) {
          // contingency matches branch indices. Find tag that matches contingency
          std::vector<bool> status = branch->getLineStatus();
          std::vector<std::string> tags = branch->getLineTags();
          int l;
          int lsize = status.size();
          for (l=0; l<lsize; l++) {
            if (tags[l] == contingency.p_ckt[j]) {
              p_saveStatus.push_back(status[j]);
              branch->setLineStatus(tags[l], false);
              found = true;
            }
          }
        }
      }
    }
  }
  int numBus = p_network->numBuses();
  gridpack::contingency_analysis::CABus *bus;
  int idx;
  size = contingency.p_busid.size();
  for (i=0; i<numBus; i++) {
    bus = dynamic_cast<gridpack::contingency_analysis::CABus*>(p_network->getBus(i).get());
    idx = bus->getOriginalIndex();
    // check bus indices against contingencies
    if (contingency.p_type == gridpack::contingency_analysis::Generator) {
      for (j = 0; j<size; j++) {
        if (contingency.p_busid[j] == idx) {
          // continency matches bus ID. Find generator ID that matches
          // contingency
          std::string genid = contingency.p_genid[j];

          std::vector<int> status = bus->getGenStatus();
          std::vector<std::string>  genids = bus->getGenerators();
          int l;
          int lsize = status.size();
          int sumIsPV = 0;
          for (l=0; l<lsize; l++) {
            if (genids[l] == genid ) {
              p_saveStatus.push_back(static_cast<bool>(status[l]));
              bus->setGenStatus(genids[l], false);
              found = true;
            } else {
              //at least one gen on the same bus is still on, still PV bus
              if (status[l] == 1 ) sumIsPV = sumIsPV + 1;
            }   
          }
          if (sumIsPV == 0) {
              bus->setIsPV(false);
          }
        }
      }
    }
  }
  if (!found && p_network->communicator().rank() == 0) {
    printf("No changes for contingency: %s. Check input file\n",
        contingency.p_name.c_str());
  }

}

/**
 * Clear contingency and set branch to its pre-contingency state
 */
void gridpack::contingency_analysis::CAFactory::clearContingency(
    gridpack::contingency_analysis::Contingency contingency)
{
  if (p_saveStatus.size() == 0) return;
  int numBranch = p_network->numBranches();
  int i,j,idx1,idx2;
  int size = contingency.p_from.size();
  int count = 0;
  gridpack::contingency_analysis::CABranch *branch;
  for (i=0; i<numBranch; i++) {
    branch = dynamic_cast<gridpack::contingency_analysis::CABranch*>(p_network->getBranch(i).get());
    idx1 = branch->getBus1OriginalIndex();
    idx2 = branch->getBus2OriginalIndex();
    // check branch indices against contingencies
    for (j = 0; j<size; j++) {
      if (contingency.p_from[j] == idx1 && contingency.p_to[j] == idx2) {
        // contingency matches branch indices. Find tag that matches contingency
        std::vector<bool> status = branch->getLineStatus();
        std::vector<std::string> tags = branch->getLineTags();
        int l;
        int lsize = status.size();
        for (l=0; l<lsize; l++) {
          if (tags[l] == contingency.p_ckt[j]) {
            branch->setLineStatus(tags[l], p_saveStatus[count]);
            count++;
          }
        }
      }
    }
  }
  int numBus = p_network->numBuses();
  gridpack::contingency_analysis::CABus *bus;
  int idx;
  size = contingency.p_busid.size();
  for (i=0; i<numBus; i++) {
    bus = dynamic_cast<gridpack::contingency_analysis::CABus*>(p_network->getBus(i).get());
    idx = bus->getOriginalIndex();

    // check bus indices against contingencies
    for (j = 0; j<size; j++) {
      if (contingency.p_type == gridpack::contingency_analysis::Generator) {
        if (contingency.p_busid[j] == idx) {
          // continency matches bus ID. Find generator ID that matches
          // contingency
          std::string genid = contingency.p_genid[j];
          std::vector<int> status = bus->getGenStatus();
          bus->resetIsPV();
          std::vector<std::string>  genids = bus->getGenerators();
          int l;
          int lsize = status.size();
          for (l=0; l<lsize; l++) {
            if (genids[l] == genid) {
              bus->setGenStatus(genids[l],static_cast<int>(p_saveStatus[count]));
              count++;
            }
          }
        }
      }
    }
  }
}

/**
 * Check for lone buses in the system. Do this by looking for buses that
 * have no branches attached to them or for whom all the branches attached
 * to the bus have all transmission elements with status false (the element
 * is off). Set status of bus to isolated so that it does not contribute to
 * powerflow matrix
 * @param stream optional stream pointer that can be used to print out IDs
 * of isolated buses
 * @return false if there is an isolated bus in the network
 */
bool gridpack::contingency_analysis::CAFactory::checkLoneBus(std::ofstream *stream)
{
  int numBus = p_network->numBuses();
  int i, j, k;
  bool bus_ok = true;
  char buf[128];
  p_saveIsolatedStatus.clear();
  for (i=0; i<numBus; i++) {
    if (!p_network->getActiveBus(i)) continue;
    gridpack::contingency_analysis::CABus *bus =
      dynamic_cast<gridpack::contingency_analysis::CABus*>
      (p_network->getBus(i).get());
    std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > branches;
    bus->getNeighborBranches(branches);
    int size = branches.size();
    bool ok = true;
    if (size == 0) {
      ok = false;
      sprintf(buf,"Lone bus %d found\n",bus->getOriginalIndex());
      if (stream != NULL) *stream << buf;
    }
    if (ok) {
      ok = false;
      for (j=0; j<size; j++) {
	bool branch_ok = false;
	std::vector<bool> status = dynamic_cast<gridpack::contingency_analysis::CABranch*>
	  (branches[j].get())->getLineStatus(); 
	int nlines = status.size();
	for (k=0; k<nlines; k++) {
	  if (status[k]) branch_ok = true;
	}
	if (branch_ok) ok = true;
      }
    }
    if (!ok) {
      sprintf(buf,"Lone bus %d found\n",bus->getOriginalIndex());
      p_saveIsolatedStatus.push_back(bus->isIsolated());
      bus->setIsolated(true);
      if (stream != NULL) *stream << buf;
    }
    if (!ok) bus_ok = false;
  }
  // Check whether bus_ok is true on all processors
  return checkTrue(bus_ok);
}

/**
 * Set lone buses back to their original status.
 */
void gridpack::contingency_analysis::CAFactory::clearLoneBus()
{
  if (p_saveIsolatedStatus.size() == 0) return;
  int numBus = p_network->numBuses();
  int i, j, k;
  int ncount = 0;
  for (i=0; i<numBus; i++) {
    if (!p_network->getActiveBus(i)) continue;
    gridpack::contingency_analysis::CABus *bus =
      dynamic_cast<gridpack::contingency_analysis::CABus*>
      (p_network->getBus(i).get());
    std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > branches;
    bus->getNeighborBranches(branches);
    int size = branches.size();
    bool ok = true;
    if (size == 0) {
      ok = false;
    }
    if (ok) {
      ok = false;
      for (j=0; j<size; j++) {
	bool branch_ok = false;
	std::vector<bool> status = dynamic_cast<gridpack::contingency_analysis::CABranch*>
	  (branches[j].get())->getLineStatus(); 
	int nlines = status.size();
	for (k=0; k<nlines; k++) {
	  if (status[k]) branch_ok = true;
	}
	if (branch_ok) ok = true;
      }
    }
    if (!ok) {
      bus->setIsolated(p_saveIsolatedStatus[ncount]);
      ncount++;
    }
  }
}

/**
 * Check to see if there any violations on the network
 * @param minV maximum voltage limit
 * @param maxV maximum voltage limit
 * @return true if no violations found
 */
void  gridpack::contingency_analysis::CAFactory::checkContingencies(
    double minV, double maxV, bool *bus_ok_r, bool *branch_ok_r)
{
  int numBus = p_network->numBuses();
  int i;
  bool bus_ok = true;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::contingency_analysis::CABus *bus =
        dynamic_cast<gridpack::contingency_analysis::CABus*>
        (p_network->getBus(i).get());
      bus->setVoltageLimits(minV, maxV);
      double V = bus->getVoltage();
      if (V < minV || V > maxV) bus_ok = false;
    }
  }
  int numBranch = p_network->numBranches();
  bool branch_ok = true;
  for (i=0; i<numBranch; i++) {
    if (p_network->getActiveBranch(i)) {
      gridpack::contingency_analysis::CABranch *branch =
	dynamic_cast<gridpack::contingency_analysis::CABranch*>
	(p_network->getBranch(i).get());
      gridpack::contingency_analysis::CABus *bus =
        dynamic_cast<gridpack::contingency_analysis::CABus*>
	(branch->getBus1().get());
      // Loop over all lines in the branch and choose the smallest rating value
      int nlines;
      p_network->getBranchData(i)->getValue(BRANCH_NUM_ELEMENTS,&nlines);
      std::vector<std::string> tags = branch->getLineTags();
      double rateA;
      for (int k = 0; k<nlines; k++) {
	if (p_network->getBranchData(i)->getValue(BRANCH_RATING_A,&rateA,k)) {
	  if (rateA > 0.0) {
	    gridpack::ComplexType s = branch->getComplexPower(tags[k]);
	    double p = real(s);
	    double q = imag(s);
	    double pq = sqrt(p*p+q*q);
            if (pq > rateA) branch_ok = false;
	  }
	}
      }
    }
  }
  *bus_ok_r = checkTrue(bus_ok);
  *branch_ok_r = checkTrue(branch_ok);
}


} // namespace contingency_analysis 
} // namespace gridpack
