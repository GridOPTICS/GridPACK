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
#include "gridpack/parallel/global_vector.hpp"
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
  p_rateB = false;
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
      printf("%s",buf);
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
      printf("\nLone bus %d reset\n",bus->getOriginalIndex());
      bus->setIsolated(p_saveIsolatedStatus[ncount]);
      ncount++;
    }
  }
}

/**
 * Set voltage limits on all buses
 * @param Vmin lower bound on voltages
 * @param Vmax upper bound on voltages
 */
void gridpack::powerflow::PFFactoryModule::setVoltageLimits(double Vmin,
    double Vmax)
{
  int numBus = p_network->numBuses();
  int i;
  bool bus_ok = true;
  for (i=0; i<numBus; i++) {
    gridpack::powerflow::PFBus *bus =
      dynamic_cast<gridpack::powerflow::PFBus*>(p_network->getBus(i).get());
    bus->setVoltageLimits(Vmin,Vmax);
  }
}


/**
 * Check to see if there are any voltage violations in the network
 * @param minV maximum voltage limit
 * @param maxV maximum voltage limit
 * @return true if no violations found
 */
bool gridpack::powerflow::PFFactoryModule::checkVoltageViolations()
{
  int numBus = p_network->numBuses();
  int i;
  bool bus_ok = true;
  char buf[128];
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      if (!bus->getIgnore()) {
        if (!bus->checkVoltageViolation()) bus_ok = false;
      }
    }
  }
  return checkTrue(bus_ok);
}

/**
 * Check to see if there are any voltage violations in the network
 * @param area only check for voltage violations in this area
 * @return true if no violations found
 */
bool gridpack::powerflow::PFFactoryModule::checkVoltageViolations(
    int area)
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
        if (!bus->checkVoltageViolation()) {
          bus_ok = false;
          gridpack::powerflow::PFFactoryModule::Violation violation;
          int idx = bus->getOriginalIndex();
          violation.bus_violation = true;
          violation.line_violation = false;
          violation.bus1 = idx;
          violation.bus2 = -1;
          p_violations.push_back(violation);
        }
      }
    }
  }
  return checkTrue(bus_ok);
}

/**
 * Set "ignore" parameter on all buses with violations so that subsequent
 * checks are not counted as violations
 */
void gridpack::powerflow::PFFactoryModule::ignoreVoltageViolations()
{
  int numBus = p_network->numBuses();
  int i;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      if (bus->checkVoltageViolation()) bus->setIgnore(true);
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
      double rate;
      for (int k = 0; k<nlines; k++) {
        if (!branch->getIgnore(tags[k])) {
          bool foundRating=false;
          if (p_rateB) {
            if (p_network->getBranchData(i)->getValue(BRANCH_RATING_B,&rate,k)) {
              foundRating = true;
            } else {
              if (p_network->getBranchData(i)->getValue(BRANCH_RATING_A,&rate,k)) {
                foundRating = true;
              }
            }
          } else {
            if (p_network->getBranchData(i)->getValue(BRANCH_RATING_A,&rate,k)) {
              foundRating = true;
            }
          }
          if (foundRating) {
            if (rate > 0.0) {
              gridpack::ComplexType s = branch->getComplexPower(tags[k]);
              double pq = abs(s);
              if (pq > rate) {
                branch_ok = false;
                gridpack::powerflow::PFFactoryModule::Violation violation;
                violation.bus_violation = false;
                violation.line_violation = true;
                violation.bus1 = branch->getBus1OriginalIndex();
                violation.bus2 = branch->getBus2OriginalIndex();
                strncpy(violation.tag,tags[k].c_str(),2);
                violation.tag[2] = '\0';
                p_violations.push_back(violation);
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
        double rate;
        for (int k = 0; k<nlines; k++) {
          if (!branch->getIgnore(tags[k])) {
            bool foundRating=false;
            if (p_rateB) {
              if (p_network->getBranchData(i)->getValue(BRANCH_RATING_B,&rate,k)) {
                foundRating = true;
              } else {
                if (p_network->getBranchData(i)->getValue(BRANCH_RATING_A,&rate,k)) {
                  foundRating = true;
                }
              }
            } else {
              if (p_network->getBranchData(i)->getValue(BRANCH_RATING_A,&rate,k)) {
                foundRating = true;
              }
            }
            if (foundRating) {
              if (rate > 0.0) {
                gridpack::ComplexType s = branch->getComplexPower(tags[k]);
                double pq = abs(s);
                if (pq > rate) branch_ok = false;
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
 * Check to see if there are any line overload violations on
 * specific lines.
 * @param bus1 original index of "from" bus for branch
 * @param bus2 original index of "to" bus for branch
 * @param tags line IDs for individual lines
 * @param violations false if violation detected on branch, true otherwise
 * @return true if no violations found
 */
bool gridpack::powerflow::PFFactoryModule::checkLineOverloadViolations(
    std::vector<int> &bus1, std::vector<int> &bus2,
    std::vector<std::string> &tags, std::vector<bool> &violations)
{
  bool branch_ok = true;
  int nbranch = bus1.size();
  if (nbranch != bus2.size() || nbranch != tags.size()) {
    printf("checkLineOverloadViolations: number of entries"
        " in bus1 and bus2 or bus1 and tags not equal\n");
    return false;
  }
  int i;
  violations.clear();
  std::vector<int> failure;
  failure.resize(nbranch);
  for (i=0; i<nbranch; i++) {
    failure[i] = 0;
    std::vector<int> indices = p_network->getLocalBranchIndices(bus1[i],bus2[i]);
    int j;
    for (j=0; j<indices.size(); j++) {
      gridpack::powerflow::PFBranch *branch = p_network->getBranch(indices[j]).get();
      // Loop over all lines in the branch and choose the smallest rating value
      int nlines;
      p_network->getBranchData(indices[j])->getValue(BRANCH_NUM_ELEMENTS,&nlines);
      std::vector<std::string> alltags = branch->getLineTags();
      double rate;
      for (int k = 0; k<nlines; k++) {
        if (tags[i] == alltags[k] && !branch->getIgnore(tags[k])) {
          bool foundRating=false;
          if (p_rateB) {
            if (p_network->getBranchData(indices[j])->getValue(BRANCH_RATING_B,
                  &rate,k)) {
              foundRating = true;
            } else {
              if (p_network->getBranchData(indices[j])->getValue(BRANCH_RATING_A,
                    &rate,k)) {
                foundRating = true;
              }
            }
          } else {
            if (p_network->getBranchData(indices[j])->getValue(BRANCH_RATING_A,
                  &rate,k)) {
              foundRating = true;
            }
          }
          if (p_network->getBranchData(indices[j])->getValue(BRANCH_RATING_A,
                &rate,k)) {
            if (rate > 0.0) {
              gridpack::ComplexType s = branch->getComplexPower(tags[k]);
              double pq = abs(s);
              if (pq > rate) {
                failure[i] = 1;
              }
            }
          }
        }
      }
    }
  }
  p_network->communicator().sum(&failure[0],nbranch);
  for (i=0; i<nbranch; i++) {
    if (failure[i] == 0) {
      violations.push_back(true);
    } else {
      violations.push_back(false);
      branch_ok = false;
    }
  }
  return branch_ok;
}

/**
 * Set "ignore" parameter on all lines with violations so that subsequent
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
 * Check to see if there are any voltage violations in the network
 * @return true if no violations found
 */
bool gridpack::powerflow::PFFactoryModule::checkQlimViolations()
{
  int numBus = p_network->numBuses();
  int i;
  bool bus_ok = true;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      if (bus->chkQlim()) bus_ok = false;
    }
  }
  p_network->updateBuses();
  for (i=0; i<numBus; i++) {
    if (!p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      bus->pushIsPV();
    }
  }
  return checkTrue(bus_ok);
}

/**
 * Check to see if there are any voltage violations in the network
 * @param area only check for voltage violations in this area
 * @return true if no violations found
 */
bool gridpack::powerflow::PFFactoryModule::checkQlimViolations(int area)
{
  int numBus = p_network->numBuses();
  int i;
  bool bus_ok = true;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      if (bus->getArea() == area) {
        if (!bus->chkQlim()) bus_ok = false;
      }
    }
  }
  p_network->updateBuses();
  for (i=0; i<numBus; i++) {
    if (!p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      bus->pushIsPV();
    }
  }
  return checkTrue(bus_ok);
}

/**
 * Clear changes that were made for Q limit violations and reset
 * system to its original state
 */
void gridpack::powerflow::PFFactoryModule::clearQlimViolations()
{
  int numBus = p_network->numBuses();
  int i;
  bool bus_ok = true;
  for (i=0; i<numBus; i++) {
    if (p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      bus->clearQlim();
    }
  }
  p_network->updateBuses();
  for (i=0; i<numBus; i++) {
    if (!p_network->getActiveBus(i)) {
      gridpack::powerflow::PFBus *bus =
        dynamic_cast<gridpack::powerflow::PFBus*>
        (p_network->getBus(i).get());
      bus->pushIsPV();
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

/**
 * Scale generator real power. If zone less than 1 then scale all
 * generators in the area.
 * @param scale factor to scale real power generation
 * @param area index of area for scaling generation
 * @param zone index of zone for scaling generation
 */
void gridpack::powerflow::PFFactoryModule::scaleGeneratorRealPower(
    double scale, int area, int zone)
{
  int nbus = p_network->numBuses();
  int i, izone;
  for (i=0; i<nbus; i++) {
    gridpack::powerflow::PFBus *bus = p_network->getBus(i).get();
    if (zone > 0) {
      izone = bus->getZone();
    } else {
      izone = zone;
    }
    if (bus->getArea() == area && zone == izone) {
      std::vector<std::string> tags = bus->getGenerators();
      int j;
      for (j=0; j<tags.size(); j++) {
        bus->scaleGeneratorRealPower(tags[j],scale);
      }
    }
  }
}

/**
 * Scale load real power. If zone less than 1 then scale all
 * loads in the area.
 * @param scale factor to scale load real power
 * @param area index of area for scaling load
 * @param zone index of zone for scaling load
 */
void gridpack::powerflow::PFFactoryModule::scaleLoadPower(
    double scale, int area, int zone)
{
  int nbus = p_network->numBuses();
  int i, izone;
  for (i=0; i<nbus; i++) {
    gridpack::powerflow::PFBus *bus = p_network->getBus(i).get();
    if (zone > 0) {
      izone = bus->getZone();
    } else {
      izone = zone;
    }
    if (bus->getArea() == area && zone == izone) {
      std::vector<std::string> tags = bus->getLoads();
      int j;
      for (j=0; j<tags.size(); j++) {
        bus->scaleLoadPower(tags[j],scale);
      }
    }
  }
}

/**
 * Return the total real power load for all loads in the zone. If zone
 * less than 1, then return the total load for the area
 * @param area index of area
 * @param zone index of zone
 * @return total load
 */
double gridpack::powerflow::PFFactoryModule::getTotalLoadRealPower(int area, int zone)
{
  double ret = 0.0;
  int nbus = p_network->numBuses();
  int i, j, izone;
  for (i=0; i<nbus; i++) {
    gridpack::powerflow::PFBus *bus = p_network->getBus(i).get();
    if (zone > 0) {
      izone = bus->getZone();
    } else {
      izone = zone;
    }
    if (bus->getArea() == area && zone == izone) {
      std::vector<std::string> tags;
      std::vector<double> pl;
      std::vector<double> ql;
      std::vector<int> status;
      bus->getLoadPower(tags,pl,ql,status);
      for (j=0; j<tags.size(); j++) {
        if (static_cast<bool>(status[j])) {
          ret += pl[j];
        }
      }
    }
  }
  p_network->communicator().sum(&ret,1);
  return ret;
}

/**
 * Return the current real power generation and the maximum and minimum total
 * power generation for all generators in the zone. If zone is less than 1
 * then return values for all generators in the area
 * @param area index of area
 * @param zone index of zone
 * @param total total real power generation
 * @param pmin minimum allowable real power generation
 * @param pmax maximum available real power generation
 */
void gridpack::powerflow::PFFactoryModule::getGeneratorMargins(int area, int zone,
    double *total, double *pmin, double *pmax)
{
  *total = 0.0;
  *pmin = 0.0;
  *pmax = 0.0;
  int nbus = p_network->numBuses();
  int i, j, izone;
  for (i=0; i<nbus; i++) {
    gridpack::powerflow::PFBus *bus = p_network->getBus(i).get();
    if (zone > 0) {
      izone = bus->getZone();
    } else {
      izone = zone;
    }
    if (bus->getArea() == area && zone == izone) {
      std::vector<std::string> tags;
      std::vector<double> tcurrent;
      std::vector<double> tpmin;
      std::vector<double> tpmax;
      std::vector<int> status;
      bus->getGeneratorMargins(tags,tcurrent,tpmin,tpmax,status);
      for (j=0; j<tcurrent.size(); j++) {
        if (status[j] != 0) {
          *total += tcurrent[j];
          *pmin += tpmin[j];
          *pmax += tpmax[j];
        }
      }
    }
  }
}

/**
 * Reset power of loads and generators to original values
 */
void gridpack::powerflow::PFFactoryModule::resetPower()
{
  int nbus = p_network->numBuses();
  int i;
  for (i=0; i<nbus; i++) {
    p_network->getBus(i)->resetPower();
  }
}

/**
 * Set parameters for real time path rating diagnostics
 * @param src_area generation area
 * @param src_zone generation zone
 * @param load_area load area
 * @param load_zone load zone
 * @param gen_scale scale factor for generation
 * @param load_scale scale factor for loads
 */
void gridpack::powerflow::PFFactoryModule::setRTPRParams(
    int src_area, int src_zone, int load_area,
    int load_zone, double gen_scale, double load_scale)
{
  int nbus = p_network->numBuses();
  int i, j, izone;
  for (i=0; i<nbus; i++) {
    gridpack::powerflow::PFBus *bus = p_network->getBus(i).get();
    int tarea = bus->getArea();
    int tzone = bus->getZone();
    if (src_zone > 0) {
      izone = tzone;
    } else {
      izone = src_zone;
    }
    bus->setScale(1.0);
    if (tarea == src_area && src_zone == izone) {
      bus->setSource(true);
      bus->setScale(gen_scale);
    } else {
      bus->setSource(false);
    }
    if (load_zone > 0) {
      izone = tzone;
    } else {
      izone = load_zone;
    }
    if (tarea == load_area && load_zone == izone) {
      bus->setSink(true);
      bus->setScale(load_scale);
    } else {
      bus->setSink(false);
    }
  }
}

/**
 * Return vector describing all violations
 * @return violation vector
 */
std::vector<gridpack::powerflow::PFFactoryModule::Violation>
gridpack::powerflow::PFFactoryModule::getViolations()
{
  std::vector<gridpack::powerflow::PFFactoryModule::Violation> ret;
  gridpack::parallel::GlobalVector<gridpack::powerflow::PFFactoryModule::Violation>
    sumVec(p_network->communicator());
  int nproc = p_network->communicator().size();
  int me = p_network->communicator().rank();
  std::vector<int> sizes(nproc);
  int i;
  for (i=0; i<nproc; i++) sizes[i] = 0;
  sizes[me] = p_violations.size();

  p_network->communicator().sum(&sizes[0],nproc);
  int offset = 0;
  for (i=1; i<me; i++) offset += sizes[i];
  int total = 0;
  for (i=0; i<nproc; i++) total += sizes[i];
  if (total == 0) return ret;
  std::vector<int> idx;
  int last = offset+sizes[me];
  for (i=offset; i<last; i++) idx.push_back(i);
  sumVec.addElements(idx,p_violations);
  sumVec.upload();
  sumVec.getAllData(ret);
  return ret;
}

/**
 * Clear violation vector
 */
void gridpack::powerflow::PFFactoryModule::clearViolations()
{
  p_violations.clear();
}

/**
 * User rate B parameter for line overload violations
 * @param flag if true, use RATEB parameter
 */
void gridpack::powerflow::PFFactoryModule::useRateB(bool flag)
{
  if (flag) {
    p_rateB = true;
  } else {
    p_rateB = false;
  }
}

} // namespace powerflow
} // namespace gridpack
