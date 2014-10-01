/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   rg_components.cpp
 * @author Bruce Palmer
 * @date   2014-02-05 08:25:59 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "rg_components.hpp"
#include "gridpack/parser/dictionary.hpp"

//#define LARGE_MATRIX

/**
 *  Simple constructor
 */
gridpack::resistor_grid::RGBus::RGBus(void)
{
  p_lead = false;
  p_voltage = 0.0;
}

/**
 *  Simple destructor
 */
gridpack::resistor_grid::RGBus::~RGBus(void)
{
}


/**
 * Load values stored in DataCollection object into RGBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::resistor_grid::RGBus::load(const
         boost::shared_ptr<gridpack::component::DataCollection> &data)
{
   int type;
   data->getValue(BUS_TYPE,&type);
   if (type == 2) {
     p_lead = true;
     data->getValue(BUS_BASEKV,&p_voltage);
   }
}

/**
 * Is bus attached to external voltage
 * @return true if voltage is fixed
 */
bool gridpack::resistor_grid::RGBus::isLead() const
{
  return p_lead;
}

/**
 * Return value of voltage at bus
 * @return voltage
 */
double gridpack::resistor_grid::RGBus::voltage() const
{
  return p_voltage;
}

/**
 * Return size of matrix block on the diagonal contributed by component
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix
 *        element
 */
bool gridpack::resistor_grid::RGBus::matrixDiagSize(int *isize, int *jsize) const
{
  if (!p_lead) {
   *isize = 1;
   *jsize = 1;
   return true;
  } else {
    return false;
  }
}

/**
 * Return the values of for a diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute
 *        matrix element
 */
bool gridpack::resistor_grid::RGBus::matrixDiagValues(ComplexType *values)
{
  if (!p_lead) {
    gridpack::ComplexType ret(0.0,0.0);
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    int size = branches.size();
    int i;
    for (i=0; i<size; i++) {
      gridpack::resistor_grid::RGBranch *branch
        = dynamic_cast<gridpack::resistor_grid::RGBranch*>(branches[i].get());
      ret += branch->resistance();
    }
    values[0] = ret;
    return true;
  } else {
    return false;
  }
}

/**
 * Return size of vector block contributed by component
 * @param isize: number of vector elements
 * @return: false if network component does not contribute
 *        vector element
 */
bool gridpack::resistor_grid::RGBus::vectorSize(int *isize) const
{
  if (!p_lead) {
    *isize = 1;
    return true;
  } else {
    return false;
  }
}

/**
 * Return the values of the vector block
 * @param values: pointer to vector values
 * @return: false if network component does not contribute
 *        vector element
 */
bool gridpack::resistor_grid::RGBus::vectorValues(ComplexType *values)
{
  if (!p_lead) {
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    int size = branches.size();
    int i;
    gridpack::ComplexType ret(0.0,0.0);
    for (i=0; i<size; i++) {
      gridpack::resistor_grid::RGBranch *branch
        = dynamic_cast<gridpack::resistor_grid::RGBranch*>(branches[i].get());
      gridpack::resistor_grid::RGBus *bus1
        = dynamic_cast<gridpack::resistor_grid::RGBus*>(branch->getBus1().get());
      gridpack::resistor_grid::RGBus *bus2
        = dynamic_cast<gridpack::resistor_grid::RGBus*>(branch->getBus2().get());
      if (bus1 != this && bus1->isLead()) {
        ret += branch->resistance()*bus1->voltage();
      } else if (bus2 != this && bus2->isLead()) {
        ret += branch->resistance()*bus2->voltage();
      }
    }
    values[0] = ret;
    return true;
  } else {
    return false;
  }
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto buses 
 * @param values array containing voltage magnitude and angle
 */
void gridpack::resistor_grid::RGBus::setValues(gridpack::ComplexType *values)
{
  if (!p_lead) {
    p_voltage = real(values[0]);
  }
}

/**
 * Write output from buses to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool gridpack::resistor_grid::RGBus::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  sprintf(string,"Voltage on bus %d: %12.6f\n",getOriginalIndex(),p_voltage);
}

/**
 *  Simple constructor
 */
gridpack::resistor_grid::RGBranch::RGBranch(void)
{
  p_resistance = 0.0;
}

/**
 *  Simple destructor
 */
gridpack::resistor_grid::RGBranch::~RGBranch(void)
{
}

/**
 * Load values stored in DataCollection object into RGBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::resistor_grid::RGBranch::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  // Resistance is stored as a vector of values in DataCollection object
  data->getValue(BRANCH_R,&p_resistance,0);
}

/**
 * Return resistance of this branch
 * @return resistance
 */
double gridpack::resistor_grid::RGBranch::resistance(void) const
{
  return p_resistance;
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::resistor_grid::RGBranch::matrixForwardSize(int *isize, int *jsize) const
{
  gridpack::resistor_grid::RGBus *bus1
    = dynamic_cast<gridpack::resistor_grid::RGBus*>(getBus1().get());
  gridpack::resistor_grid::RGBus *bus2
    = dynamic_cast<gridpack::resistor_grid::RGBus*>(getBus2().get());
  if (!bus1->isLead() && !bus2->isLead()) {
    *isize = 1;
    *jsize = 1;
    return true;
  } else {
    return false;
  }
}

bool gridpack::resistor_grid::RGBranch::matrixReverseSize(int *isize, int *jsize) const
{
  gridpack::resistor_grid::RGBus *bus1
    = dynamic_cast<gridpack::resistor_grid::RGBus*>(getBus1().get());
  gridpack::resistor_grid::RGBus *bus2
    = dynamic_cast<gridpack::resistor_grid::RGBus*>(getBus2().get());
  if (!bus1->isLead() && !bus2->isLead()) {
    *isize = 1;
    *jsize = 1;
    return true;
  } else {
    return false;
  }
}

/**
 * Return the values of the forward/reverse matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::resistor_grid::RGBranch::matrixForwardValues(ComplexType *values)
{
  gridpack::resistor_grid::RGBus *bus1
    = dynamic_cast<gridpack::resistor_grid::RGBus*>(getBus1().get());
  gridpack::resistor_grid::RGBus *bus2
    = dynamic_cast<gridpack::resistor_grid::RGBus*>(getBus2().get());
  if (!bus1->isLead() && !bus2->isLead()) {
    values[0] = -p_resistance;
    return true;
  } else {
    return false;
  }
}

bool gridpack::resistor_grid::RGBranch::matrixReverseValues(ComplexType *values)
{
  gridpack::resistor_grid::RGBus *bus1
    = dynamic_cast<gridpack::resistor_grid::RGBus*>(getBus1().get());
  gridpack::resistor_grid::RGBus *bus2
    = dynamic_cast<gridpack::resistor_grid::RGBus*>(getBus2().get());
  if (!bus1->isLead() && !bus2->isLead()) {
    values[0] = -p_resistance;
    return true;
  } else {
    return false;
  }
}

/**
 * Write output from branches to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if branch is contributing string to output, false otherwise
 */
bool gridpack::resistor_grid::RGBranch::serialWrite(char *string, const int
    bufsize,  const char *signal)
{

  gridpack::resistor_grid::RGBus *bus1
    = dynamic_cast<gridpack::resistor_grid::RGBus*>(getBus1().get());
  gridpack::resistor_grid::RGBus *bus2
    = dynamic_cast<gridpack::resistor_grid::RGBus*>(getBus2().get());
  double v1 = bus1->voltage();
  double v2 = bus2->voltage();
  double icur = (v1 - v2)*p_resistance;
  sprintf(string,"Current on line from bus %d to %d is: %12.6f\n",
      bus1->getOriginalIndex(),bus2->getOriginalIndex(),icur);
  return true;
}
