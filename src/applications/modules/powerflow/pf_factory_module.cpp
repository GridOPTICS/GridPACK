/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_factory.cpp
 * @author Bruce Palmer
 * @date   2014-01-28 11:31:23 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "pf_factory_module.hpp"


namespace gridpack {
namespace powerflow {

// Powerflow factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
PFFactoryModule::PFFactoryModule(PFFactoryModule::NetworkPtr network)
  : gridpack::factory::BaseFactory<PFNetwork>(network)
{
  p_network = network;
}

/**
 * Basic destructor
 */
gridpack::powerflow::PFFactoryModule::~PFFactoryModule()
{
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::powerflow::PFFactoryModule::setYBus(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // Invoke setYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
    dynamic_cast<PFBranch*>(p_network->getBranch(i).get())->setYBus();
  }

  // Invoke setYBus method on all bus objects
  for (i=0; i<numBus; i++) {
    dynamic_cast<PFBus*>(p_network->getBus(i).get())->setYBus();
  }

}

/**
  * Make SBus vector 
  */
void gridpack::powerflow::PFFactoryModule::setSBus(void)
{
  int numBus = p_network->numBuses();
  int i;

  // Invoke setSBus method on all bus objects
  for (i=0; i<numBus; i++) {
    dynamic_cast<PFBus*>(p_network->getBus(i).get())->setSBus();
  }
}

/**
  * Create the PQ 
  */
void gridpack::powerflow::PFFactoryModule::setPQ(void)
{
  int numBus = p_network->numBuses();
  int i;
  ComplexType values[2];

  for (i=0; i<numBus; i++) {
    dynamic_cast<PFBus*>(p_network->getBus(i).get())->vectorValues(values);
  }
}

/**
 * Update pg of specified bus element based on their genID
 * @param name 
 * @param busID
 * @param genID
 * @param value
 */
//void gridpack::powerflow::PFFactoryModule::updatePg(std::string &name, int busID, std::string genID, double value)
void gridpack::powerflow::PFFactoryModule::updatePg(int busID, std::string genID, double value)
{
  int numBus = p_network->numBuses();
  int i;
  int genIndex=0;
  for (i=0; i<numBus; i++) {
//    dynamic_cast<PFBus*>(p_network->getBus(i).get())->setParam(name,busID, genID, value);
    dynamic_cast<PFBus*>(p_network->getBus(i).get())->setParam(GENERATOR_PG,
        busID, genID, value);
  }
}

void gridpack::powerflow::PFFactoryModule::updateQg(int busID, std::string genID, double value)
{
  int numBus = p_network->numBuses();
  int i;
  int genIndex=0;
  for (i=0; i<numBus; i++) {
//    dynamic_cast<PFBus*>(p_network->getBus(i).get())->setParam(name,busID, genID, value);
    dynamic_cast<PFBus*>(p_network->getBus(i).get())->setParam(GENERATOR_QG,
        busID, genID, value);
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
bool gridpack::powerflow::PFFactoryModule::checkLoneBus(std::ofstream *stream)
{
  int numBus = p_network->numBuses();
  int i, j, k;
  bool bus_ok = true;
  char buf[128];
  p_saveIsolatedStatus.clear();
  for (i=0; i<numBus; i++) {
    if (!p_network->getActiveBus(i)) continue;
    gridpack::powerflow::PFBus *bus =
      dynamic_cast<gridpack::powerflow::PFBus*>
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
        std::vector<bool> status =
          dynamic_cast<gridpack::powerflow::PFBranch*>
          (branches[j].get())->getLineStatus();
        int nlines = status.size();
        for (k=0; k<nlines; k++) {
          if (status[k]) branch_ok = true;
        }
        if (branch_ok) ok = true;
      }
    }
    if (!ok) {
      sprintf(buf,"\nLone bus %d found\n",bus->getOriginalIndex());
      p_saveIsolatedStatus.push_back(bus->isIsolated());
      bus->setIsolated(true);
      if (stream != NULL) *stream << buf;
    }
    if (!ok) bus_ok = false;
  }
  // Check whether bus_ok is true on all processors
  return checkTrue(!bus_ok);
}

/**
 * Set lone buses back to their original status.
 */
void gridpack::powerflow::PFFactoryModule::clearLoneBus()
{
  if (p_saveIsolatedStatus.size() == 0) return;
  int numBus = p_network->numBuses();
  int i, j, k;
  int ncount = 0;
  for (i=0; i<numBus; i++) {
    if (!p_network->getActiveBus(i)) continue;
    gridpack::powerflow::PFBus *bus =
      dynamic_cast<gridpack::powerflow::PFBus*>
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
        std::vector<bool> status =
          dynamic_cast<gridpack::powerflow::PFBranch*>
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
 * Check to see if there are any voltage violations in the network
 * @param minV maximum voltage limit
 * @param maxV maximum voltage limit
 * @return true if no violations found
 */
bool gridpack::powerflow::PFFactoryModule::checkVoltageViolations(
    double Vmin, double Vmax)
{
  int numBus = p_network->numBuses();
  int i;
  bool bus_ok = true;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      //bus->setVoltageLimits(Vmin, Vmax);
      if (!bus->getIgnore()) {
        double V = bus->getVoltage();
        if (V < Vmin || V > Vmax) bus_ok = false;
      }
    }
  }
  return checkTrue(bus_ok);
}

/**
 * Check to see if there are any voltage violations in the network
 * @param minV maximum voltage limit
 * @param maxV maximum voltage limit
 * @param area only check for voltage violations in this area
 * @return true if no violations found
 */
bool gridpack::powerflow::PFFactoryModule::checkVoltageViolations(
    int area, double Vmin, double Vmax)
{
  int numBus = p_network->numBuses();
  int i;
  bool bus_ok = true;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      if (!bus->getIgnore() && bus->getArea() == area) {
        double V = bus->getVoltage();
        if (V < Vmin || V > Vmax) bus_ok = false;
      }
    }
  }
  return checkTrue(bus_ok);
}

/**
 * Set "ignore" parameter on all buses with violations so that subsequent
 * checks are not counted as violations
 * @param minV maximum voltage limit
 * @param maxV maximum voltage limit
 */
void gridpack::powerflow::PFFactoryModule::ignoreVoltageViolations(double Vmin,
    double Vmax)
{
  int numBus = p_network->numBuses();
  int i;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      double V = bus->getVoltage();
      if (V < Vmin || V > Vmax) bus->setIgnore(true);
    }
  }
}


/**
 * Clear "ignore" parameter on all buses
 */
void gridpack::powerflow::PFFactoryModule::clearVoltageViolations()
{
  int numBus = p_network->numBuses();
  int i;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      bus->setIgnore(false);
    }
  }
}

/**
 * Check to see if there are any line overload violations in the
 * network
 * @return true if no violations found
 */
bool gridpack::powerflow::PFFactoryModule::checkLineOverloadViolations()
{
  int numBranch = p_network->numBranches();
  int i;
  bool branch_ok = true;
  for (i=0; i<numBranch; i++) {
    if (p_network->getActiveBranch(i)) {
      gridpack::powerflow::PFBranch *branch =
        dynamic_cast<gridpack::powerflow::PFBranch*>
        (p_network->getBranch(i).get());
      // Loop over all lines in the branch and choose the smallest rating value
      int nlines;
      p_network->getBranchData(i)->getValue(BRANCH_NUM_ELEMENTS,&nlines);
      std::vector<std::string> tags = branch->getLineTags();
      double rateA;
      for (int k = 0; k<nlines; k++) {
        if (!branch->getIgnore(tags[k])) {
          if (p_network->getBranchData(i)->getValue(BRANCH_RATING_A,&rateA,k)) {
            if (rateA > 0.0) {
              gridpack::ComplexType s = branch->getComplexPower(tags[k]);
              double pq = abs(s);
              if (pq > rateA) branch_ok = false;
            }
          }
        }
      }
    }
  }
  return checkTrue(branch_ok);
}

/**
 * Check to see if there are any line overload violations in the
 * network
 * @param area only check for voltage violations in this area
 * @return true if no violations found
 */
bool gridpack::powerflow::PFFactoryModule::checkLineOverloadViolations(int area)
{
  int numBranch = p_network->numBranches();
  int i;
  bool branch_ok = true;
  for (i=0; i<numBranch; i++) {
    if (p_network->getActiveBranch(i)) {
      gridpack::powerflow::PFBranch *branch =
        dynamic_cast<gridpack::powerflow::PFBranch*>
        (p_network->getBranch(i).get());
      // get buses at either end
      gridpack::powerflow::PFBus *bus1 =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (branch->getBus1().get());
      gridpack::powerflow::PFBus *bus2 =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (branch->getBus2().get());
      // Loop over all lines in the branch and choose the smallest rating value
      if (bus1->getArea() == area || bus2->getArea() == area) {
        int nlines;
        p_network->getBranchData(i)->getValue(BRANCH_NUM_ELEMENTS,&nlines);
        std::vector<std::string> tags = branch->getLineTags();
        double rateA;
        for (int k = 0; k<nlines; k++) {
          if (!branch->getIgnore(tags[k])) {
            if (p_network->getBranchData(i)->getValue(BRANCH_RATING_A,&rateA,k)) {
              if (rateA > 0.0) {
                gridpack::ComplexType s = branch->getComplexPower(tags[k]);
                double pq = abs(s);
                if (pq > rateA) branch_ok = false;
              }
            }
          }
        }
      }
    }
  }
  return checkTrue(branch_ok);
}

/**
 * Set "ignore" paramter on all lines with violations so that subsequent
 * checks are not counted as violations
 */
void gridpack::powerflow::PFFactoryModule::ignoreLineOverloadViolations()
{
  int numBranch = p_network->numBranches();
  int i;
  for (i=0; i<numBranch; i++) {
    if (p_network->getActiveBranch(i)) {
      gridpack::powerflow::PFBranch *branch =
        dynamic_cast<gridpack::powerflow::PFBranch*>
        (p_network->getBranch(i).get());
      // Loop over all lines in the branch and choose the smallest rating value
      int nlines;
      p_network->getBranchData(i)->getValue(BRANCH_NUM_ELEMENTS,&nlines);
      std::vector<std::string> tags = branch->getLineTags();
      double rateA;
      for (int k = 0; k<nlines; k++) {
        if (p_network->getBranchData(i)->getValue(BRANCH_RATING_A,&rateA,k)) {
          if (rateA > 0.0) {
            gridpack::ComplexType s = branch->getComplexPower(tags[k]);
            double pq = abs(s);
            if (pq > rateA) branch->setIgnore(tags[k],true);
          }
        }
      }
    }
  }
}

/**
 * Clear "ignore" parameter on all lines
 */
void gridpack::powerflow::PFFactoryModule::clearLineOverloadViolations()
{
  int numBranch = p_network->numBranches();
  int i;
  for (i=0; i<numBranch; i++) {
    if (p_network->getActiveBranch(i)) {
      gridpack::powerflow::PFBranch *branch =
        dynamic_cast<gridpack::powerflow::PFBranch*>
        (p_network->getBranch(i).get());
      // Loop over all lines in the branch and choose the smallest rating value
      int nlines;
      p_network->getBranchData(i)->getValue(BRANCH_NUM_ELEMENTS,&nlines);
      std::vector<std::string> tags = branch->getLineTags();
      double rateA;
      for (int k = 0; k<nlines; k++) {
        branch->setIgnore(tags[k],false);
      }
    }
  }
}

/**
 * Reinitialize voltages
 */
void gridpack::powerflow::PFFactoryModule::resetVoltages()
{
  int numBus = p_network->numBuses();
  int i;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      bus->resetVoltage();
    }
  }
}

} // namespace powerflow
} // namespace gridpack
