/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_components.cpp
 * @author Yousu Chen
 * @date   2016-07-14 13:50:44 d3g096
 * 
 * @brief  
 * 
 * @update Yousu Chen
 *         Adding functions of bad data dection, chi-square testing 
 * @date   2025-03-05
 *         Adding more functions to handle measurements more efficiently
 * @date   2025-04-02
 *         Adding more functions to handle bad data detection and create more comprehensive outputs
 * @date   2025-04-20
 *         Added IIJ and IJI measurements
 *         Cleaned up debug information
 * @date   2025-11-19
 *         Implemented parallel version
 * @date   2025-11-25
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>
#include <chrono>
#include <ctime>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "se_components.hpp"

//#define LARGE_MATRIX

namespace {
  // Helper function to trim trailing whitespace from circuit IDs
  std::string trimCircuitID(const std::string& str) {
    size_t end = str.find_last_not_of(" \t\r\n");
    return (end == std::string::npos) ? "" : str.substr(0, end + 1);
  }

  // Helper function to compare circuit IDs ignoring trailing whitespace
  bool circuitIDsMatch(const std::string& ckt1, const std::string& ckt2) {
    return trimCircuitID(ckt1) == trimCircuitID(ckt2);
  }
}

/**
 *  Simple constructor
 */
gridpack::state_estimation::SEBus::SEBus(void)
{
  p_shunt_gs = 0.0;
  p_shunt_bs = 0.0;
  p_v = 0.0;
  p_a = 0.0;
  p_theta = 0.0;
  p_angle = 0.0;
  p_voltage = 0.0;
  p_pl = 0.0;
  p_ql = 0.0;
  p_sbase = 0.0;
  p_mode = YBus;
  p_rowJidx.clear();
  p_rowRidx.clear();
  p_colJidx.clear();
  p_colRidx.clear();
  setReferenceBus(false);
  // Default voltage limits 
  p_v_min = 0.9;
  p_v_max = 1.1;
  p_enforce_v_limits = false;
  
  // Initialize Jacobian optimization cache
  p_cache_valid = false;
  p_last_jacobian_update = -1;
  p_cached_neighbors.clear();
  
  // Initialize performance profiling counters
  p_cache_hits = 0;
  p_cache_misses = 0;
  p_total_cache_time = 0.0;
  p_total_jacobian_time = 0.0;
  p_jacobian_calls = 0;
}

/**
 *  Simple destructor
 */
gridpack::state_estimation::SEBus::~SEBus(void)
{
}

/**
 *  Return size of matrix block contributed by the component
 *  @param isize, jsize: number of rows and columns of matrix block
 *  @return: false if network component does not contribute matrix element
 */
bool gridpack::state_estimation::SEBus::matrixDiagSize(int *isize, int *jsize) const
{
  if (p_mode == YBus) {
    return YMBus::matrixDiagSize(isize,jsize);
  } else if (p_mode == Jacobian_H) {
    if (!isIsolated()) {
      *isize = 1;
      *jsize = 1;
      return true;
    } else {
      return false;
    }
  }
  return false;
}

/**
 * Return the values of the matrix block. The values are
 * returned in row-major order.
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::state_estimation::SEBus::matrixDiagValues(ComplexType *values)
{
  if (p_mode == YBus) {
    return YMBus::matrixDiagValues(values);
/*  } else if (p_mode == Jacobian_H) {
    std::vector<gridpack::state_estimation::Measurement>
    meas = p_meas; //p_meas supposed to be all measurements on this bus
    int nmeas = meas.size();
    int i;
    for (i=0; i<nmeas; i++ ) {
       if (meas[i].p_type == "VM") {
          values[0] = 0.0; 
          values[1] = 1.0; 
       } else if (meas[i].p_type == "PI") {
         std::vector<boost::shared_ptr<BaseComponent> > branches;
         getNeighborBranches(branches);
         int size = branches.size();
         int j;
         double ret1 = 0.0;
         double ret2 = 0.0;
         for (j=0; j<size; j++) {
           gridpack::state_estimation::SEBranch *branch
             = dynamic_cast<gridpack::state_estimation::SEBranch*>(branches[i].get());
          branch->getVTheta(this,&v,&theta);
          ret1 += p_v * v * (p_ybusr_frwd*sin(theta) + p_ybusi_frwd*cos(theta)) - p_v * p_v * p_ybusi;
          ret2 +=  v * (p_ybusr_frwd*cos(theta) + p_ybusi_frwd*sin(theta)) + p_v * p_ybusr;
         }
          values[0] = ret1;
          values[1] = ret2;
       //} // to add other bus measurements
       }
*/
  }
  return false;
}

/**
 * Return the size of the block that this component contributes to the
 * vector
 * @param size: size of vector block
 * @return: false if component does not contribute to vector
 */
bool gridpack::state_estimation::SEBus::vectorSize(int *size) const
{
  if (p_mode == Voltage) {
    if (!isIsolated()) {
      if (getReferenceBus()) {
//        return false;
        *size = 1;
      } else {
        *size = 2;
      }
      return true;
    } else {
      return false;
    }
  }
  return false;
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto buses 
 * @param values array containing voltage magnitude and angle
 */
void gridpack::state_estimation::SEBus::setValues(gridpack::ComplexType *values)
{
  double vt = p_v;
  double at = p_a;
  if (p_mode == Voltage) {
   if (!isIsolated()) {
    if (getReferenceBus()) {
      p_v += real(values[0]);
    } else {
      p_a += real(values[0]);
      p_v += real(values[1]);
    }
    *p_vAng_ptr = p_a;
    *p_vMag_ptr = p_v;
    // Apply voltage projection using the configurable limits
    if (p_enforce_v_limits) {
      if (p_v > p_v_max) {
        p_v = p_v_max;
        *p_vMag_ptr = p_v;
      } else if (p_v < p_v_min) {
        p_v = p_v_min;
        *p_vMag_ptr = p_v;
      }
    }
    
    // Always invalidate cache when state variables change to ensure
    // Jacobian/residual use current neighbor data
    invalidateNeighborCache();
   }
  }
}

/**
 * Return number of elements in vector coming from component
 * @return number of elements contributed from component
 */
int gridpack::state_estimation::SEBus::vectorNumElements() const
{
  if (p_mode == Jacobian_H || p_mode == Residual) {
    // When computing residuals or Jacobian, we need the measurement count
    return p_meas.size();
  } else if (p_mode == Voltage) {
    // Original code for Voltage mode
    if (!isIsolated()) {
      if (getReferenceBus()) {
        return 1;
      } else {
        return 2;
      }
    }
  }
  return 0;
}

/**
 * Set indices corresponding to the elements contributed by this component
 * @param ielem index of element contributed by this component (e.g.
 * if component contributes 3 elements then ielem is between 0 and 2)
 * @param idx vector index of element ielem
 */
void gridpack::state_estimation::SEBus::vectorSetElementIndex(int ielem, int idx)
{
  if (p_mode == Jacobian_H || p_mode == Residual) {
    if (ielem < p_vecZidx.size()) {
      p_vecZidx[ielem] = idx;
    } else {
      p_vecZidx.push_back(idx);
    }
  }
}

/**
 * Get list of element indices from component
 * @param idx list of indices that component maps onto
 */
void gridpack::state_estimation::SEBus::vectorGetElementIndices(int *idx)
{
  if (p_mode == Jacobian_H || p_mode == Residual) {
    int nsize = p_vecZidx.size();
    int i;
    for (i=0; i<nsize; i++) {
      idx[i] = p_vecZidx[i];
    }
  }
}

/**
 * Transfer vector values to component
 * @param values list of vector element values
 */
void gridpack::state_estimation::SEBus::vectorSetElementValues(ComplexType *values)
{
  //TODO: Is this function needed?
}

/**
 * Return the size of the buffer used in data exchanges on the network.
 * For this problem, the voltage magnitude and phase angle need to be exchanged
 * @return size of buffer
 */
int gridpack::state_estimation::SEBus::getXCBufSize(void)
{
  return 2*sizeof(double);
}

/**
 * Assign pointers for voltage magnitude and phase angle
 */
void gridpack::state_estimation::SEBus::setXCBuf(void *buf)
{
  p_vAng_ptr = static_cast<double*>(buf);
  p_vMag_ptr = p_vAng_ptr+1;
  // Note: we are assuming that the load function has been called BEFORE
  // the factory setExchange method, so p_a and p_v are set with their initial
  // values.
  *p_vAng_ptr = p_a;
  *p_vMag_ptr = p_v;
}

/**
 * Load values stored in DataCollection object into SEBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::state_estimation::SEBus::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  YMBus::load(data);

  p_sbase = 100.0;
  data->getValue(CASE_SBASE, &p_sbase);
  data->getValue(BUS_VOLTAGE_ANG, &p_angle);
  data->getValue(BUS_VOLTAGE_MAG, &p_voltage); 
  p_v = p_voltage;
  double pi = 4.0*atan(1.0);
  p_angle = p_angle*pi/180.0;
  p_a = p_angle;
  int itype;
  data->getValue(BUS_TYPE, &itype);
  if (itype == 3) {
    setReferenceBus(true);
  }
  p_shunt_gs = 0.0;
  p_shunt_bs = 0.0;
  data->getValue(BUS_SHUNT_GL, &p_shunt_gs,0);
  data->getValue(BUS_SHUNT_BL, &p_shunt_bs,0);
  // if BUS_TYPE = 2 then bus is a PV bus
  p_isPV = false;
  p_isSlack = false;
  
  // Set bus type flags based on bus type
  if (itype == 2) {
    p_isPV = true;   // PV (generator) bus
  } else if (itype == 3) {
    p_isSlack = true; // Slack (reference) bus
  }

  // added p_pg,p_qg,p_pl,p_ql,p_sbase;
  p_load = true;
  p_pl = 0.0;
  p_ql = 0.0;
  data->getValue(LOAD_PL, &p_pl);
  data->getValue(LOAD_QL, &p_ql);
  bool lgen;
  int i, ngen, gstatus;
  double pg, qg, vs, qmin, qmax;
  ngen = 0;
  if (data->getValue(GENERATOR_NUMBER, &ngen)) {
    double qtot = 0.0;
    for (i=0; i<ngen; i++) {
      lgen = true;
      lgen = lgen && data->getValue(GENERATOR_PG, &pg,i);
      lgen = lgen && data->getValue(GENERATOR_QG, &qg,i);
      lgen = lgen && data->getValue(GENERATOR_VS, &vs,i);
      lgen = lgen && data->getValue(GENERATOR_STAT, &gstatus,i);
      lgen = lgen && data->getValue(GENERATOR_QMIN, &qmin,i);
      lgen = lgen && data->getValue(GENERATOR_QMAX, &qmax,i);
      if (lgen) {
        p_pg.push_back(pg);
        p_qg.push_back(qg);
        p_gstatus.push_back(gstatus);
        p_qmin.push_back(qmin);
        p_qmax.push_back(qmax);
        if (gstatus == 1) {
          p_v = vs; //reset initial PV voltage to set voltage
          if (itype == 2) p_isPV = true;
        }
      }
    }
  }
}

/**
 * Set values of YBus matrix. These can then be used in subsequent
 * calculations
 */
void gridpack::state_estimation::SEBus::setYBus(void)
{
  YMBus::setYBus();
  gridpack::ComplexType ret;
  ret = YMBus::getYBus();
  p_ybusr = real(ret);
  p_ybusi = imag(ret);
}

/**
 * Get values of YBus matrix. These can then be used in subsequent
 * calculations
 */
gridpack::ComplexType gridpack::state_estimation::SEBus::getYBus(void)
{
  return YMBus::getYBus();
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::state_estimation::SEBus::setMode(int mode)
{
  // Save the mode internally
  p_mode = mode;
  
  // Handle specific setup for different modes
  if (mode == YBus) {
    YMBus::setMode(gridpack::ymatrix::YBus);
  } else if (mode == Jacobian_H) {
    // Setup for Jacobian matrix calculation
    // Nothing special needed here beyond setting p_mode
  } else if (mode == R_inv) {
    // Setup for measurement error covariance inverse
    // Nothing special needed here beyond setting p_mode
  } else if (mode == Voltage) {
    // Setup for voltage state handling
    // Nothing special needed here beyond setting p_mode
  } else if (mode == Residual) {
    // Setup for residual calculation
    // Make sure latest state values are used
    if (p_vAng_ptr != NULL && p_vMag_ptr != NULL) {
      p_a = *p_vAng_ptr;
      p_v = *p_vMag_ptr;
    }
    // Calculate power injections if needed for residuals
    calculatePowerInjections();
  }
}

/**
 * Calculate real and reactive power injections at this bus
 * Used for residual calculations
 */
void gridpack::state_estimation::SEBus::calculatePowerInjections(void)
{
  // Only calculate if we're not isolated
  if (!isIsolated()) {
    // Get the latest state variables
    p_a = *p_vAng_ptr;
    p_v = *p_vMag_ptr;
    
    // Initialize power injections
    p_Pinj = 0.0;
    p_Qinj = 0.0;
    
    // Sum contributions from all branches
    std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > branches;
    getNeighborBranches(branches);
    int nsize = branches.size();
    
    for (int j = 0; j < nsize; j++) {
      SEBranch *branch = dynamic_cast<SEBranch*>(branches[j].get());
      double p, q;
      branch->getPQ(this, &p, &q);
      p_Pinj += p;
      p_Qinj += q;
    }
    
    // Add shunt contributions
    p_Pinj += p_v * p_v * p_ybusr;
    p_Qinj += p_v * p_v * (-p_ybusi);
  }
}

/**
 * Return the value of the voltage magnitude on this bus
 * @return: voltage magnitude
 */
double gridpack::state_estimation::SEBus::getVoltage()
{
  p_v = *p_vMag_ptr;
  return *p_vMag_ptr;
}

/**
 * Return whether or not the bus is a PV bus (V held fixed in powerflow
 * equations)
 * @return true if bus is PV bus
 */
bool gridpack::state_estimation::SEBus::isPV(void)
{
  return p_isPV;
}

bool gridpack::state_estimation::SEBus::isSlack(void)
{
  return p_isSlack;
}

/**
 * Return whether or not a bus is isolated
 * @return true if bus is isolated
 */
bool gridpack::state_estimation::SEBus::isIsolated(void) const
{
  return YMBus::isIsolated();
}

/**
 * Return the value of the phase angle on this bus
 * @return: phase angle
 */
double gridpack::state_estimation::SEBus::getPhase()
{
  p_a = *p_vAng_ptr;
  return *p_vAng_ptr;
}

/**
 * Write output from buses to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool gridpack::state_estimation::SEBus::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  if (!isIsolated()) { 
  p_v = *p_vMag_ptr;
  p_a = *p_vAng_ptr;

  if (signal == NULL) {
    double pi = 4.0*atan(1.0);
    double angle = p_a*180.0/pi;
    snprintf(string, bufsize, "     %6d      %12.6f         %12.6f\n",
        getOriginalIndex(),angle,p_v);
  } else if (!strcmp(signal,"pq")) {
    gridpack::ComplexType v[2];
    vectorValues(v);
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    snprintf(string, bufsize, "     %6d      %12.6f         %12.6f      %2d\n",
        getOriginalIndex(),real(v[0]),real(v[1]),
        static_cast<int>(branches.size()));
  } else if (!strcmp(signal,"se")) {
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    int nsize = branches.size();
    double p,q,P, Q;
    P = 0.0;
    Q = 0.0;
    for (int j=0; j<nsize; j++) {
      gridpack::state_estimation::SEBranch *branch
        = dynamic_cast<gridpack::state_estimation::SEBranch*>(branches[j].get());
      branch->getPQ(this, &p, &q);
      P += p;
      Q += q;
    }
    P += p_v*p_v*p_ybusr;
    Q += p_v*p_v*(-p_ybusi);
    p_Pinj = P;
    p_Qinj = Q;
    if (p_meas.size()>0) {
      int nmeas = p_meas.size();
      char buf[128];
      int ilen = 0;
      std::string meas_type,type;
      for (int i=0; i<nmeas; i++) {
        meas_type = p_meas[i].p_type;
        if (meas_type.length() == 3) {
          type = meas_type;
        } else if (meas_type.length() == 2) {
          type = " ";
          type.append(meas_type);
        }
        double estimate;
        buf[0] = '\0';
        if (meas_type == "VM") {
          estimate = p_v;
          snprintf(buf, sizeof(buf), "    %s %8d    %16.5f  %16.5f   %16.5f    %8.4f\n",
              type.c_str(),getOriginalIndex(),p_meas[i].p_value, estimate,
              estimate-p_meas[i].p_value,p_meas[i].p_deviation);
        } else if (meas_type == "PI") {
          estimate = p_Pinj;
          snprintf(buf, sizeof(buf), "    %s %8d    %16.5f  %16.5f   %16.5f    %8.4f\n",
              type.c_str(),getOriginalIndex(),p_meas[i].p_value, estimate,
              estimate-p_meas[i].p_value,p_meas[i].p_deviation);
        } else if (meas_type == "QI") {
          estimate = p_Qinj;
          snprintf(buf, sizeof(buf), "    %s %8d    %16.5f  %16.5f   %16.5f    %8.4f\n",
              type.c_str(),getOriginalIndex(),p_meas[i].p_value, estimate,
              estimate-p_meas[i].p_value,p_meas[i].p_deviation);
        }
        int buflen = strlen(buf);
        if (buflen + ilen < bufsize) {
          snprintf(string, bufsize - ilen, "%s", buf);
          string += buflen;
          ilen += buflen;
        }
      }
      if (ilen == 0) return false;
      return true;
    } else {
      return false;
    }
  }
  return true;
  } else {
  return false;
  }
}

/**
 * Add a measurement to the bus
 * @param measurement a measurement struct that will be used to
 * assign
 * internal paramters
 */
void gridpack::state_estimation::SEBus::addMeasurement(
    gridpack::state_estimation::Measurement measurement)
{
  p_meas.push_back(measurement);
}

/**
 * Sort measurements so that they are in a consistent order
 */
void gridpack::state_estimation::SEBus::sortMeasurements(void)
{
  // Use a multi-map to create an ordered list of measurements
  std::multimap<std::string,gridpack::state_estimation::Measurement> map;
  int nsize = p_meas.size();
  int i;
  char sbuf[32];
  for (i=0; i<nsize; i++) {
    // Create a unique (hopefully) key for each measurement
    snprintf(sbuf, sizeof(sbuf), "%s%s%f", p_meas[i].p_type, p_meas[i].p_ckt, p_meas[i].p_value);
    std::string key = sbuf;
    map.insert(std::pair<std::string,
        gridpack::state_estimation::Measurement>(key,p_meas[i]));
  }
  p_meas.clear();
  std::multimap<std::string,gridpack::state_estimation::Measurement>::iterator it;
  it = map.begin();
  while (it != map.end()) {
    p_meas.push_back(it->second);
    it++;
  }
}

/**
 * Return the complex voltage on this bus
 * @return the complex voltage
 */
gridpack::ComplexType gridpack::state_estimation::SEBus::getComplexVoltage(void)
{
  p_a = *p_vAng_ptr;
  p_v = *p_vMag_ptr;
  gridpack::ComplexType ret(cos(p_a),sin(p_a));
  ret = ret*p_v;
  return ret;
}

/**
 * Configure buses with state estimation parameters. These can be
 * used in other methods
 */
void gridpack::state_estimation::SEBus::configureSE(void)
{
  // Calculate the number of matrix values associated with this bus
  int nmeas = p_meas.size(); // Suppose p_meas is the vector of all the measurements on this bus
  int ncnt = 0;
  int i, j, nsize;
  int busid = getOriginalIndex();
  for (i=0; i<nmeas; i++) {
//   if (!isIsolated()) {
    std::string type = p_meas[i].p_type;
    if (type == "VM" || type == "VA" || type == "VUL" || type == "VLL") {
      if (!getReferenceBus()) {
        ncnt += 2;
      } else {
        ncnt++;
      }
    } else if (type == "PI" || type == "QI") {
      std::vector<boost::shared_ptr<BaseComponent> > branch_nghbrs;
      getNeighborBranches(branch_nghbrs);
      nsize = branch_nghbrs.size();
      for (j=0; j<nsize; j++) {
        SEBranch *branch
          = dynamic_cast<SEBranch*>(branch_nghbrs[j].get());
        SEBus *bus = dynamic_cast<SEBus*>(branch->getBus1().get());
        if (bus == this) bus = dynamic_cast<SEBus*>(branch->getBus2().get());
        if (!bus->getReferenceBus()) {
          ncnt += 2;
        } else {
          ncnt++;
        }
      }
      if (!getReferenceBus()) { 
        ncnt += 2;
      } else {
        ncnt++;
      }
    }
   //}
  } 
  p_numElements = ncnt;
}

/**
 * Save state variables inside the component to a DataCollection object.
 * This can be used as a way of moving data in a way that is useful for
 * creating output or for copying state data from one network to another.
 * @param data data collection object into which new values are inserted
 */
void gridpack::state_estimation::SEBus::saveData(
    boost::shared_ptr<gridpack::component::DataCollection> data)
{
  double rval;
  int i;
  if (!data->setValue("BUS_SE_VMAG",*p_vMag_ptr)) {
    data->addValue("BUS_SE_VMAG",*p_vMag_ptr);
  }
  rval = *p_vAng_ptr;
  double pi = 4.0*atan(1.0);
  rval = 180.0*rval/pi;
  if (!data->setValue("BUS_SE_VANG",rval)) {
    data->addValue("BUS_SE_VANG",rval);
  }
  int ngen=p_qmin.size();
  double sum_qmin = 0.0;
  double sum_qmax = 0.0;
  for (i=0; i<ngen; i++) {
    sum_qmin += p_qmin[i];
    sum_qmax += p_qmax[i];
  }
  for (i=0; i<ngen; i++) {
    rval = p_pg[i];
    if (!data->setValue("GENERATOR_SE_PGEN",rval,i)) {
      data->addValue("GENERATOR_SE_PGEN",rval,i);
    }
    rval = p_qmin[i] + (p_Qinj-sum_qmin)*(p_qmax[i]-p_qmin[i])
      / (sum_qmax - sum_qmin);
    rval += p_ql/p_sbase;
    if (!data->setValue("GENERATOR_SE_QGEN",rval,i)) {
      data->addValue("GENERATOR_SE_QGEN",rval,i);
    }
  }
}

/**
 * Return number of rows in matrix from component
 * @return number of rows from component
 */
int gridpack::state_estimation::SEBus::matrixNumRows() const
{
  return p_meas.size();
}

/**
 * Return number of cols in matrix from component
 * @return number of cols from component
 */
int gridpack::state_estimation::SEBus::matrixNumCols() const
{
  if (p_mode == Jacobian_H) {
    // Check to see if this bus has measurements or is attached to anything that
    // has measurements
//   if (!isIsolated()) {
    bool meas = false;
    if (p_meas.size() > 0) meas = true;
    if (!meas) {
      std::vector<boost::shared_ptr<BaseComponent> > branch_nghbrs;
      getNeighborBranches(branch_nghbrs);
      std::vector<boost::shared_ptr<BaseComponent> > bus_nghbrs;
      getNeighborBuses(bus_nghbrs);
      int nsize = branch_nghbrs.size();
      int i;
      gridpack::state_estimation::SEBus *bus;
      gridpack::state_estimation::SEBranch *branch;
      for (i=0; i<nsize && !meas; i++) {
        bus = dynamic_cast<SEBus*>(bus_nghbrs[i].get());
        branch = dynamic_cast<SEBranch*>(branch_nghbrs[i].get());
        if (bus->matrixNumRows() > 0) meas = true;
        if (branch->matrixNumRows() > 0) meas = true;
      }
    }
    if (!meas) return 0;
    // Bus has measurements associated with it.
    if (!isIsolated()) {
      if (!getReferenceBus()) {
        return 2;
      } else {
        return 1;
      }
    } else {
       return 0;
    } 
//   } else {
//     return 0;
//   }
  } else if (p_mode == R_inv) {
    return p_meas.size();
  }
  return 0;
}

/**
 * Set row indices corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @param idx matrix index of row irow
 */
void gridpack::state_estimation::SEBus::matrixSetRowIndex(int irow, int idx)
{
  if (p_mode == Jacobian_H) {
    if (irow < p_rowJidx.size()) {
      p_rowJidx[irow] = idx;
    } else {
      p_rowJidx.push_back(idx);
    }
  } else if (p_mode == R_inv) {
    if (irow < p_rowRidx.size()) {
      p_rowRidx[irow] = idx;
    } else {
      p_rowRidx.push_back(idx);
    }
  }
}

/**
 * Set column indices corresponding to the columns contributed by this component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @param idx matrix index of column icol
 */
void gridpack::state_estimation::SEBus::matrixSetColIndex(int icol, int idx)
{
  if (p_mode == Jacobian_H) {
    if (icol < p_colJidx.size()) {
      p_colJidx[icol] = idx;
    } else {
      p_colJidx.push_back(idx);
    }
  } else if (p_mode == R_inv) {
    if (icol < p_colRidx.size()) {
      p_colRidx[icol] = idx;
    } else {
      p_colRidx.push_back(idx);
    }
  }
}

/**
 * Get the row index corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @return matrix index of row irow
 */
int gridpack::state_estimation::SEBus::matrixGetRowIndex(int idx)
{
  if (p_mode == Jacobian_H) {
    if (idx >= p_rowJidx.size())
      printf("violation in bus:matrixGetRowIndex bus: %d size: %d idx: %d\n",
          getOriginalIndex(),idx,static_cast<int>(p_rowJidx.size()));
    return p_rowJidx[idx];
  } else if (p_mode == R_inv) {
    return p_rowRidx[idx];
  }
  return -1;
}

/**
 * Get the column index corresponding to the columns contributed by this component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @return matrix index of column icol
 */
int gridpack::state_estimation::SEBus::matrixGetColIndex(int idx)
{
  if (p_mode == Jacobian_H) {
    if (idx >= p_colJidx.size())
      printf("violation in bus:matrixGetColIndex bus: %d size: %d idx: %d\n",
          getOriginalIndex(),idx,static_cast<int>(p_colJidx.size()));
    return p_colJidx[idx];
  } else if (p_mode == R_inv) {
    return p_colRidx[idx];
  }
  return -1;
}

/**
 * Return the number of matrix values contributed by this component
 * @return number of matrix values
 */
int gridpack::state_estimation::SEBus::matrixNumValues() const
{
  if (p_mode == Jacobian_H) {
    return p_numElements;
  } else if (p_mode == R_inv) {
    return p_meas.size();
  }
  return 0;
}

/**
 * Return values from a matrix block
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
*/
void gridpack::state_estimation::SEBus::matrixGetValues(int *nvals,ComplexType *values, int *rows, int *cols)
{
  // Start timing for Jacobian calculation
  auto jacobian_start = std::chrono::high_resolution_clock::now();
  
  p_v = *p_vMag_ptr;
  p_a = *p_vAng_ptr;
  if (p_mode == Jacobian_H) {
   if (!isIsolated()) {
    int nmeas = p_meas.size(); // Suppose p_meas is the vector of all the measurements on this bus
    int ncnt = 0;
    int i, j, im, jm, nsize;
    double v, theta, yfbusr,yfbusi;
    std::string ctk, type;
    for (i=0; i<nmeas; i++) {
      im = matrixGetRowIndex(i);
      ctk = p_meas[i].p_ckt;
      type = p_meas[i].p_type;
      if (type == "VM") {
        if (!getReferenceBus()) { 
          jm = matrixGetColIndex(0);
          values[ncnt] = gridpack::ComplexType(0.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
          jm = matrixGetColIndex(1);
          values[ncnt] = gridpack::ComplexType(1.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
        } else {
          jm = matrixGetColIndex(0);
          values[ncnt] = gridpack::ComplexType(1.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
        }
      } else if (type == "VUL") {  // ADD THIS SECTION FOR VIRTUAL UPPER LIMIT
         p_v = *p_vMag_ptr;
         double limit = p_meas[i].p_value;

         // Always add Jacobian elements to maintain consistent matrix structure
         // The residual will be zero when not violated, so this won't affect the solution
         if (!getReferenceBus()) {
           // Angle derivative (zero)
           jm = matrixGetColIndex(0);
           values[ncnt] = gridpack::ComplexType(0.0, 0.0);
           rows[ncnt] = im;
           cols[ncnt] = jm;
           ncnt++;

           // Voltage derivative (one)
           jm = matrixGetColIndex(1);
           values[ncnt] = gridpack::ComplexType(1.0, 0.0);
           rows[ncnt] = im;
           cols[ncnt] = jm;
           ncnt++;
         } else {
           // Reference bus case
           jm = matrixGetColIndex(0);
           values[ncnt] = gridpack::ComplexType(1.0, 0.0);
           rows[ncnt] = im;
           cols[ncnt] = jm;
           ncnt++;
         }
      } else if (type == "VLL") { // ADD THIS SECTION FOR VIRTUAL LOWER LIMIT
         p_v = *p_vMag_ptr;
         double limit = p_meas[i].p_value;

         // Always add Jacobian elements (maintain matrix structure)
         // The residual will be zero when not violated
         if (!getReferenceBus()) {
           // Angle derivative (zero)
           jm = matrixGetColIndex(0);
           values[ncnt] = gridpack::ComplexType(0.0, 0.0);
           rows[ncnt] = im;
           cols[ncnt] = jm;
           ncnt++;

           // Voltage derivative (one)
           jm = matrixGetColIndex(1);
           values[ncnt] = gridpack::ComplexType(1.0, 0.0);
           rows[ncnt] = im;
           cols[ncnt] = jm;
           ncnt++;
         } else {
           // Reference bus case
           jm = matrixGetColIndex(0);
           values[ncnt] = gridpack::ComplexType(1.0, 0.0);
           rows[ncnt] = im;
           cols[ncnt] = jm;
           ncnt++;
         }
      } else if (type == "PI") {
        // Use cached neighbor data to avoid repeated getNeighborBranches() calls
        if (!p_cache_valid) {
          auto cache_start = std::chrono::high_resolution_clock::now();
          cacheNeighborBranchData();
          auto cache_end = std::chrono::high_resolution_clock::now();
          p_total_cache_time += std::chrono::duration<double>(cache_end - cache_start).count();
          p_cache_misses++;
        } else {
          p_cache_hits++;
        }
        
        nsize = p_cached_neighbors.size();
        double ret1 = 0.0;
        double ret2 = 0.0;
        
        for (j=0; j<nsize; j++) {
          const NeighborBranchData& neighbor = p_cached_neighbors[j];
          
          if (neighbor.isValid) {
            // Use cached values
            v = neighbor.v;
            theta = neighbor.theta;
            yfbusr = neighbor.yfbusr;
            yfbusi = neighbor.yfbusi;
            SEBus* bus = neighbor.otherBus;
            
            // Calculate Jacobian elements using cached data
            ret1 += p_v * v * (-yfbusr*sin(theta) + yfbusi*cos(theta));
            if (!bus->getReferenceBus()) {
              values[ncnt] = gridpack::ComplexType(p_v*v*(yfbusr*sin(theta)-yfbusi*cos(theta)),0.0);
              jm = bus->matrixGetColIndex(0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus->matrixGetColIndex(1);
            } else {
              jm = bus->matrixGetColIndex(0);
            }
            ret2 += v * (yfbusr*cos(theta) + yfbusi*sin(theta));
            values[ncnt] = gridpack::ComplexType(p_v*(yfbusr*cos(theta)+yfbusi*sin(theta)),0.0);
            rows[ncnt] = im;
            cols[ncnt] = jm;
            ncnt++;
          }
        }
        if (!getReferenceBus()) {
          jm = matrixGetColIndex(0);
          values[ncnt] = gridpack::ComplexType(ret1,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
          jm = matrixGetColIndex(1);
        } else {
          jm = matrixGetColIndex(0);
        }
        ret2 += 2 * p_v * p_ybusr;
        values[ncnt] = gridpack::ComplexType(ret2,0.0); 
        rows[ncnt] = im;
        cols[ncnt] = jm;
        ncnt++;
      } else if (type == "QI") {
        // Use cached neighbor data to avoid repeated getNeighborBranches() calls
        if (!p_cache_valid) {
          auto cache_start = std::chrono::high_resolution_clock::now();
          cacheNeighborBranchData();
          auto cache_end = std::chrono::high_resolution_clock::now();
          p_total_cache_time += std::chrono::duration<double>(cache_end - cache_start).count();
          p_cache_misses++;
        } else {
          p_cache_hits++;
        }
        
        nsize = p_cached_neighbors.size();
        double ret1 = 0.0;
        double ret2 = 0.0;
        
        for (j=0; j<nsize; j++) {
          const NeighborBranchData& neighbor = p_cached_neighbors[j];
          
          if (neighbor.isValid) {
            // Use cached values
            v = neighbor.v;
            theta = neighbor.theta;
            yfbusr = neighbor.yfbusr;
            yfbusi = neighbor.yfbusi;
            SEBus* bus = neighbor.otherBus;
            
            // Calculate Jacobian elements using cached data
            ret1 += p_v * v * (yfbusr*cos(theta) + yfbusi*sin(theta));
            if (!bus->getReferenceBus()) {
              values[ncnt] = gridpack::ComplexType(p_v*v*(-yfbusr*cos(theta)-yfbusi*sin(theta)),0.0);
              jm = bus->matrixGetColIndex(0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus->matrixGetColIndex(1);
            } else {
              jm = bus->matrixGetColIndex(0);
            }
            ret2 += v * (yfbusr*sin(theta) - yfbusi*cos(theta));
            values[ncnt] = gridpack::ComplexType(p_v*v*(yfbusr*sin(theta)-yfbusi*cos(theta)),0.0);
            rows[ncnt] = im;
            cols[ncnt] = jm;
            ncnt++;
          }
        }
//        ret1 += p_v * p_v * p_ybusr;
//        ret2 += p_v * p_ybusi
//        ret1 -= p_v * p_v * p_ybusr;
        if (!getReferenceBus()) {
          jm = matrixGetColIndex(0);
          values[ncnt] = gridpack::ComplexType(ret1,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
          jm = matrixGetColIndex(1);
        } else {
          jm = matrixGetColIndex(0);
        }
        ret2 -= 2 * p_v * p_ybusi;
        values[ncnt] = gridpack::ComplexType(ret2,0.0); 
        rows[ncnt] = im;
        cols[ncnt] = jm;
        ncnt++;
      } else if (type == "VA") {
#if 0
        std::vector<boost::shared_ptr<BaseComponent> > branch_nghbrs;
        getNeighborBranches(branch_nghbrs);
        nsize = branch_nghbrs.size();
        for (j=0; j<nsize; j++) {
          jm = matrixGetColIndex(0);
          values[ncnt] = gridpack::ComplexType(0.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
          jm = matrixGetColIndex(1);
          values[ncnt] = gridpack::ComplexType(0.0,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
        }
#endif
        if (!getReferenceBus()) {
          jm = matrixGetColIndex(0);
          // Convert radians to degrees
          values[ncnt] = gridpack::ComplexType(180.0/M_PI,0.0); 
          rows[ncnt] = im;
          cols[ncnt] = jm;
          ncnt++;
          jm = matrixGetColIndex(1);
        } else {
          jm = matrixGetColIndex(0);
        }
        values[ncnt] = gridpack::ComplexType(0.0,0.0); 
        rows[ncnt] = im;
        cols[ncnt] = jm;
        ncnt++;
      }
    }
    *nvals = ncnt;
   }
  } else if (p_mode == R_inv) {
   if (!isIsolated()) {
    int nsize = p_meas.size();
    int i;
    for (i=0; i<nsize; i++) {
      if (p_meas[i].p_deviation != 0.0) {
        values[i] = 1.0/(p_meas[i].p_deviation*p_meas[i].p_deviation);
      } else {
        values[i] = 0.0;
      }
      rows[i] = matrixGetRowIndex(i);
      cols[i] = matrixGetColIndex(i);
    }
    *nvals = nsize;
   }
  }
  
  // End timing for Jacobian calculation
  if (p_mode == Jacobian_H) {
    auto jacobian_end = std::chrono::high_resolution_clock::now();
    p_total_jacobian_time += std::chrono::duration<double>(jacobian_end - jacobian_start).count();
    p_jacobian_calls++;
  }
}

/**
 * Return values from a vector
 * @param values: pointer to vector values (z-h(x))
 * @param idx: pointer to vector index 
*/
void gridpack::state_estimation::SEBus::vectorGetElementValues(ComplexType *values, int *idx)
{
  p_a = *p_vAng_ptr;
  p_v = *p_vMag_ptr;
  if (p_mode == Jacobian_H) {
   if (!isIsolated()) {
    int nmeas = p_meas.size(); // Suppose p_meas is the vector of all the measurements on this bus
    int ncnt = 0;
    int i, j, im, jm, nsize;
    double v, theta,yfbusr,yfbusi;

    vectorGetElementIndices(idx);
    for (i=0; i<nmeas; i++) {
       std::string type = p_meas[i].p_type;
       if (type == "VM") {
         int index = getGlobalIndex();
         values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-p_v),0.0);
         ncnt++;
       } else if (type == "PI") {
         // Use cached neighbor data to avoid repeated getNeighborBranches() calls
         if (!p_cache_valid) {
           auto cache_start = std::chrono::high_resolution_clock::now();
           cacheNeighborBranchData();
           auto cache_end = std::chrono::high_resolution_clock::now();
           p_total_cache_time += std::chrono::duration<double>(cache_end - cache_start).count();
           p_cache_misses++;
         } else {
           p_cache_hits++;
         }
         
         nsize = p_cached_neighbors.size();
         double ret = 0.0;
         
         for (j=0; j<nsize; j++) {
           const NeighborBranchData& neighbor = p_cached_neighbors[j];
           
           if (neighbor.isValid) {
             // Use cached values
             v = neighbor.v;
             theta = neighbor.theta;
             yfbusr = neighbor.yfbusr;
             yfbusi = neighbor.yfbusi;
             
             ret += v * (yfbusr*cos(theta) + yfbusi*sin(theta));
           }
         }
         ret += p_v * p_ybusr;
         ret *= p_v; 
         int index = getGlobalIndex();
         values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret),0.0);
         ncnt++;
       } else if (type == "QI") {
         // Use cached neighbor data to avoid repeated getNeighborBranches() calls
         if (!p_cache_valid) {
           auto cache_start = std::chrono::high_resolution_clock::now();
           cacheNeighborBranchData();
           auto cache_end = std::chrono::high_resolution_clock::now();
           p_total_cache_time += std::chrono::duration<double>(cache_end - cache_start).count();
           p_cache_misses++;
         } else {
           p_cache_hits++;
         }
         
         nsize = p_cached_neighbors.size();
         double ret = 0.0;
         
         for (j=0; j<nsize; j++) {
           const NeighborBranchData& neighbor = p_cached_neighbors[j];
           
           if (neighbor.isValid) {
             // Use cached values
             v = neighbor.v;
             theta = neighbor.theta;
             yfbusr = neighbor.yfbusr;
             yfbusi = neighbor.yfbusi;
             
             ret += v * (yfbusr*sin(theta) - yfbusi*cos(theta));
           }
         }
         ret -= p_v * p_ybusi;
         ret *= p_v; 
         int index = getGlobalIndex();
         values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret),0.0);
         ncnt++;

      } else if (type == "VA") {
         int index = getGlobalIndex();
         values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value - (p_a * 180.0/M_PI)),0.0);
         ncnt++;
      }
    } 
   }
  } else if (p_mode == Residual) {
    vectorGetElementIndices(idx);
    // Only proceed if not isolated and have measurements
    if (!isIsolated() && p_meas.size() > 0) {
      // Get latest state variables
      p_a = *p_vAng_ptr;
      p_v = *p_vMag_ptr;
      
      // Calculate power injections
      double Pinj = 0.0, Qinj = 0.0;
      std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > branches;
      getNeighborBranches(branches);
      
      for (int j = 0; j < branches.size(); j++) {
        SEBranch* branch = dynamic_cast<SEBranch*>(branches[j].get());
        double p, q;
        branch->getPQ(this, &p, &q);
        Pinj += p;
        Qinj += q;
      }
      
      // Add shunt contribution
      Pinj += p_v * p_v * p_ybusr;
      Qinj += p_v * p_v * (-p_ybusi);
      
      // Store for later use
      p_Pinj = Pinj;
      p_Qinj = Qinj;
      
      // Process each measurement
      for (int i = 0; i < p_meas.size(); i++) {
        std::string type = p_meas[i].p_type;
        double measured = p_meas[i].p_value;
        double estimated = 0.0;
        
        if (type == "VM") {
          estimated = p_v;
        } 
        // Handle virtual limit measurements
	else if (type == "VUL") {
          // Upper limit: only compute residual if limit is violated (consistent with Jacobian)
          if (p_v > measured) {
            estimated = p_v;  // Violation: residual = measured - p_v (negative, pushes voltage down)
          } else {
            estimated = measured;  // No violation: residual = 0
          }
        }
	else if (type == "VLL") {
          // Lower limit: only compute residual if limit is violated (consistent with Jacobian)
          if (p_v < measured) {
            estimated = p_v;  // Violation: residual = measured - p_v (positive, pushes voltage up)
          } else {
            estimated = measured;  // No violation: residual = 0
          }
        }
	 else if (type == "VA") { // angle to radians
          estimated = p_a * 180.0 / M_PI;
        } else if (type == "PI") {
          estimated = Pinj;
        } else if (type == "QI") {
          estimated = Qinj;
        }
        
        // Calculate residual
        double residual = measured - estimated;
        values[i] = ComplexType(residual, 0.0);
      }
    }
  } else if (p_mode == R_inv) {
    // Implement proper R_inv handling for branch measurements
    vectorGetElementIndices(idx);
    int nmeas = p_meas.size();

    if (nmeas > 0) {
      int bus_idx = getOriginalIndex();
      
      for (int i = 0; i < nmeas; i++) {
        // Set diagonal elements of R_inv matrix (inverse of variance)
        double sigma = p_meas[i].p_deviation;
        if (sigma > 0.0) {
          double weight = 1.0/(sigma*sigma);
          values[i] = ComplexType(weight, 0.0);
        } else {
          values[i] = ComplexType(0.0, 0.0);
        }
      }
    }
  }
}


/**
 *  Get shunt gs and bs
 */

void gridpack::state_estimation::SEBus::getShuntGsBs(double *gs, double *bs) 
{
  *gs = p_shunt_gs/p_sbase; 
  *bs = p_shunt_bs/p_sbase; 
//  *gs = p_shunt_gs;
//  *bs = p_shunt_bs;
}

/**
 * Check if residual index matches this bus and report bad data
 * @param idx: index to check
 * @param report: if true, print information about the bad data
 * @return: true if index found on this bus
 */
bool gridpack::state_estimation::SEBus::checkResidualIndex(int idx, bool report)
{
  int nsize = p_vecZidx.size();
  for (int i = 0; i < nsize; i++) {
    if (p_vecZidx[i] == idx) {
      if (report) {
        printf("Bad data detected on bus %d, measurement type: %s, value: %f, deviation: %f\n",
               getOriginalIndex(), p_meas[i].p_type, p_meas[i].p_value, p_meas[i].p_deviation);
      }
      return true;
    }
  }
  return false;
}

bool gridpack::state_estimation::SEBus::getResidualDetails(int idx, char* buffer)
{
  int nsize = p_vecZidx.size();
  for (int i = 0; i < nsize; i++) {
    if (p_vecZidx[i] == idx) {
      snprintf(buffer, 256, "Bus %d, Type: %s, Value: %f, Deviation: %f", 
              getOriginalIndex(), p_meas[i].p_type, 
              p_meas[i].p_value, p_meas[i].p_deviation);
      return true;
    }
  }
  return false;
}

/**
 * Return the values of the residual vector
 * @param values: pointer to vector values
 * @return: false if network component does not contribute vector element
 */

bool gridpack::state_estimation::SEBus::vectorValues(ComplexType *values)
{
  if (p_mode == Voltage) {
    // Existing code for Voltage mode
    return true;
  } else if (p_mode == Residual) {
    // Only proceed if we have measurements and aren't isolated
    if (!isIsolated() && p_meas.size() > 0) {
      // Get latest state variables
      p_a = *p_vAng_ptr;
      p_v = *p_vMag_ptr;
      
      // Calculate power injections
      double Pinj = 0.0, Qinj = 0.0;
      std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > branches;
      getNeighborBranches(branches);
      
      for (int j = 0; j < branches.size(); j++) {
        SEBranch* branch = dynamic_cast<SEBranch*>(branches[j].get());
        double p, q;
        branch->getPQ(this, &p, &q);
        Pinj += p;
        Qinj += q;
      }
      
      // Add shunt contribution
      Pinj += p_v * p_v * p_ybusr;
      Qinj += p_v * p_v * (-p_ybusi);
      
      // Store for later use
      p_Pinj = Pinj;
      p_Qinj = Qinj;

      // Process each measurement
      for (int i = 0; i < p_meas.size(); i++) {
        std::string type = p_meas[i].p_type;
        double measured = p_meas[i].p_value;
        double estimated = 0.0;
        
        if (type == "VM") {
          estimated = p_v;
        } else if (type == "VUL") {
          if (p_v > measured) {
            estimated = p_v;
          } else {
            estimated = measured;
          }
        } else if (type == "VLL") {
          if (p_v < measured) {
            estimated = p_v;
          } else {
            estimated = measured;
          }
        } else if (type == "VA") {
          estimated = p_a * 180.0 / M_PI;
        } else if (type == "PI") {
          estimated = Pinj;
        } else if (type == "QI") {
          estimated = Qinj;
        }
        
        // Calculate residual
        double residual = measured - estimated;
        values[i] = ComplexType(residual, 0.0);
      }

      return true;  // Critical to return true
    }
  }
  return false;
}


bool gridpack::state_estimation::SEBus::serialWritePreCheck(char *string, const int bufsize)
{
  if (!isIsolated()) {
    int nmeas = vectorNumElements();
    if (nmeas > 0) {
      char buf[128];
      int idx = getOriginalIndex();
      
      // Loop through all measurements on this bus
      for (int i = 0; i < p_meas.size(); i++) {
        std::string type = p_meas[i].p_type;
        double value = p_meas[i].p_value;
        double deviation = p_meas[i].p_deviation;
        
        // Check for suspicious values based on measurement type
        bool suspicious = false;
        std::string reason = "";
        
        if (type == "VM") {
          // Voltage magnitudes are typically in range 0.9-1.1
          if (value < 0.8 || value > 1.2) {
            suspicious = true;
            reason = "unusual voltage magnitude";
          }
        } else if (type == "PI" || type == "QI") {
          // Power injections should be reasonable for system size
          // For a 14-bus system, values above 5.0 are suspicious
          if (fabs(value) > 5.0) {
            suspicious = true;
            reason = "unusually large power value";
          }
        }
        
        if (suspicious) {
          snprintf(buf, sizeof(buf), "SUSPICIOUS MEASUREMENT: Bus %d, Type %s, Value %.4f - %s\n",
                  idx, type.c_str(), value, reason.c_str());
          printf("%s", buf);
          return true;
        }
      }
    }
  }
  return false;
}

/**
 * Set voltage magnitude limits
 * @param v_min minimum voltage limit
 * @param v_max maximum voltage limit
 * @param enforce whether to enforce limits
 */
void gridpack::state_estimation::SEBus::setVoltageLimits(double v_min, double v_max, bool enforce)
{
  p_v_min = v_min;
  p_v_max = v_max;
  p_enforce_v_limits = enforce;
}

// In SEBranch class
bool gridpack::state_estimation::SEBranch::getResidualDetails(int idx, char* buffer)
{
  int nsize = p_vecZidx.size();
  for (int i = 0; i < nsize; i++) {
    if (p_vecZidx[i] == idx) {
      boost::shared_ptr<gridpack::component::BaseComponent> b1 = getBus1();
      SEBus *bus1 = dynamic_cast<SEBus*>(b1.get());
      boost::shared_ptr<gridpack::component::BaseComponent> b2 = getBus2();
      SEBus *bus2 = dynamic_cast<SEBus*>(b2.get());
      snprintf(buffer, 256, "Branch from Bus %d to Bus %d, Circuit %s, Type: %s, Value: %f, Deviation: %f", 
              bus1->getOriginalIndex(), bus2->getOriginalIndex(), 
              p_meas[i].p_ckt, p_meas[i].p_type, 
              p_meas[i].p_value, p_meas[i].p_deviation);
      return true;
    }
  }
  return false;
}

/**
 * Cache neighbor branch data for efficient Jacobian computation
 */
void gridpack::state_estimation::SEBus::cacheNeighborBranchData() const
{
  // Performance optimization: Return early if cache is already valid
  if (p_cache_valid) {
    p_cache_hits++;
    return;
  }
  
  // Clear existing cache
  p_cached_neighbors.clear();
  
  // Get current state variables 
  double current_v = *p_vMag_ptr;
  double current_a = *p_vAng_ptr;
  
  // Get neighbor branches
  std::vector<boost::shared_ptr<gridpack::component::BaseComponent> > branch_nghbrs;
  getNeighborBranches(branch_nghbrs);
  
  // Cache data for each neighbor branch
  for (int j = 0; j < branch_nghbrs.size(); j++) {
    NeighborBranchData data;
    
    data.branch = dynamic_cast<SEBranch*>(branch_nghbrs[j].get());
    SEBus *bus = dynamic_cast<SEBus*>(data.branch->getBus1().get());
    SEBus *bus2 = dynamic_cast<SEBus*>(data.branch->getBus2().get());
    
    // Check if both buses are valid (not isolated)
    data.isValid = !bus->isIsolated() && !bus2->isIsolated();
    
    if (data.isValid) {
      // Determine direction and get other bus
      if (bus == this) {
        data.isForward = true;
        data.otherBus = dynamic_cast<SEBus*>(data.branch->getBus2().get());
        gridpack::ComplexType yfbus = data.branch->getForwardYBus();
        data.yfbusr = real(yfbus);
        data.yfbusi = imag(yfbus);
      } else {
        data.isForward = false;
        data.otherBus = bus;
        gridpack::ComplexType yfbus = data.branch->getReverseYBus();
        data.yfbusr = real(yfbus);
        data.yfbusi = imag(yfbus);
      }
      
      // Get voltage and angle data  
      data.branch->getVTheta(const_cast<SEBus*>(this), &data.v, &data.theta);
    } else {
      data.otherBus = nullptr;
      data.yfbusr = 0.0;
      data.yfbusi = 0.0;
      data.v = 0.0;
      data.theta = 0.0;
      data.isForward = true;
    }
    
    p_cached_neighbors.push_back(data);
  }
  
  p_cache_valid = true;
}

/**
 * Invalidate the neighbor branch cache
 */
void gridpack::state_estimation::SEBus::invalidateNeighborCache()
{
  p_cache_valid = false;
  p_cached_neighbors.clear();
}

/**
 * Reset performance profiling statistics
 */
void gridpack::state_estimation::SEBus::resetPerformanceCounters()
{
  p_cache_hits = 0;
  p_cache_misses = 0;
  p_total_cache_time = 0.0;
  p_total_jacobian_time = 0.0;
  p_jacobian_calls = 0;
}

/**
 * Get performance profiling statistics
 */
void gridpack::state_estimation::SEBus::getPerformanceStats(int& cache_hits, int& cache_misses, 
                                                           double& cache_time, double& jacobian_time, int& jacobian_calls) const
{
  cache_hits = p_cache_hits;
  cache_misses = p_cache_misses;
  cache_time = p_total_cache_time;
  jacobian_time = p_total_jacobian_time;
  jacobian_calls = p_jacobian_calls;
}


/**
 *  SEBranch constructor
 */
gridpack::state_estimation::SEBranch::SEBranch(void)
{
  p_reactance.clear();
  p_resistance.clear();
  p_tap_ratio.clear();
  p_phase_shift.clear();
  p_charging.clear();
  p_shunt_admt_g1.clear();
  p_shunt_admt_b1.clear();
  p_shunt_admt_g2.clear();
  p_shunt_admt_b2.clear();
  p_xform.clear();
  p_shunt.clear();
  p_branch_status.clear();
  p_elems = 0;
  p_theta = 0.0;
  p_sbase = 0.0;
  p_rowJidx.clear();
  p_rowRidx.clear();
  p_colJidx.clear();
  p_colRidx.clear();
  p_mode = YBus;
}

/**
 *  SEBranch destructor
 */
gridpack::state_estimation::SEBranch::~SEBranch(void)
{
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::state_estimation::SEBranch::matrixForwardSize(int *isize, int *jsize) const
{
  if (p_mode == YBus) {
    return YMBranch::matrixForwardSize(isize,jsize);
  }
  return false;
}
bool gridpack::state_estimation::SEBranch::matrixReverseSize(int *isize, int *jsize) const
{
  if (p_mode == YBus) {
    return YMBranch::matrixReverseSize(isize,jsize);
  }
  return false;
}

/**
 * Return the values of the off-diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::state_estimation::SEBranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == YBus) {
    return YMBranch::matrixForwardValues(values);
  }
  return false;
}

bool gridpack::state_estimation::SEBranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == YBus) {
    return YMBranch::matrixReverseValues(values);
  }
  return false;
}

/**
  * Calculate contributions to the admittance matrix from the branches
**/
void gridpack::state_estimation::SEBranch::setYBus(void)
{
  YMBranch::setYBus();
  gridpack::ComplexType ret;
  ret = YMBranch::getForwardYBus();
  p_ybusr_frwd = real(ret);
  p_ybusi_frwd = imag(ret);
  ret = YMBranch::getReverseYBus();
  p_ybusr_rvrs = real(ret);
  p_ybusi_rvrs = imag(ret);
  gridpack::state_estimation::SEBus *bus1 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
  double pi = 4.0*atan(1.0);
  bool ok = !bus1->isIsolated();
  ok = ok && !bus2->isIsolated();
  if (ok) {
    p_theta = (bus1->getPhase() - bus2->getPhase());
  }
}

/**
 * Load values stored in DataCollection object into SEBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::state_estimation::SEBranch::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  YMBranch::load(data);
  gridpack::state_estimation::SEBus *bus1 = 
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
  bool ok1 = !bus1->isIsolated();
  ok1 = ok1 && !bus2->isIsolated();
  if (ok1) {

  data->getValue(BRANCH_NUM_ELEMENTS, &p_elems);
  double rvar;
  int ivar;
  std::string svar;
  double pi = 4.0*atan(1.0);
  p_active = false;
  p_sbase = 100.0;
  data->getValue(CASE_SBASE, &p_sbase);
  int idx;
  for (idx = 0; idx<p_elems; idx++) {
    data->getValue(BRANCH_X, &rvar, idx);
    if (rvar <1.0e-5 && rvar >=0.0) rvar = 1.0e-5;
    if (rvar >-1.0e-5 && rvar <0.0) rvar =-1.0e-5;
    p_reactance.push_back(rvar);
    data->getValue(BRANCH_R, &rvar, idx);
    p_resistance.push_back(rvar);
    rvar = 0.0; 
    data->getValue(BRANCH_SHIFT, &rvar, idx);
    rvar = -rvar*pi/180.0; 
    p_phase_shift.push_back(rvar);
    rvar = 0.0;
    data->getValue(BRANCH_TAP, &rvar, idx);
    p_tap_ratio.push_back(rvar); 
    //ok = ok && data->getValue(BRANCH_CKT, &svar, idx);
    data->getValue(BRANCH_CKT, &svar, idx);
    p_tag.push_back(svar);
    if (rvar != 0.0) {
      p_xform.push_back(true);
    } else {
      p_xform.push_back(false);
    }
    ivar = 1;
    data->getValue(BRANCH_STATUS, &ivar, idx);
    p_branch_status.push_back(static_cast<bool>(ivar));
    if (ivar == 1) p_active = true;
    bool shunt = true;
    rvar = 0.0;
    shunt = shunt && data->getValue(BRANCH_B, &rvar, idx);
    p_charging.push_back(rvar);
    rvar = 0.0;
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_G1, &rvar, idx);
    p_shunt_admt_g1.push_back(rvar);
    rvar = 0.0;
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_B1, &rvar, idx);
    p_shunt_admt_b1.push_back(rvar);
    rvar = 0.0;
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_G2, &rvar, idx);
    p_shunt_admt_g2.push_back(rvar);
    rvar = 0.0;
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_B2, &rvar, idx);
    p_shunt_admt_b2.push_back(rvar);
    p_shunt.push_back(shunt);
  }
  }
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::state_estimation::SEBranch::setMode(int mode)
{
  // Save the mode internally
  p_mode = mode;
  
  // Handle specific setup for different modes
  if (mode == YBus) {
    YMBranch::setMode(gridpack::ymatrix::YBus);
  } else if (mode == Jacobian_H) {
    // Setup for Jacobian matrix calculation
    // Nothing special needed here beyond setting p_mode
  } else if (mode == R_inv) {
    // Setup for measurement error covariance inverse
    // Nothing special needed here beyond setting p_mode
  } else if (mode == Voltage) {
    // Setup for voltage state handling
    // Nothing special needed here beyond setting p_mode
  } else if (mode == Residual) {
    // Setup for residual calculation
    // Update branch calculations if needed
      boost::shared_ptr<gridpack::component::BaseComponent> b1 = getBus1();
      SEBus *bus1 = dynamic_cast<SEBus*>(b1.get());
      boost::shared_ptr<gridpack::component::BaseComponent> b2 = getBus2();
      SEBus *bus2 = dynamic_cast<SEBus*>(b2.get());
    
    // Only continue if both buses are not isolated
    if (!bus1->isIsolated() && !bus2->isIsolated()) {
      // Calculate current phase angle
      p_theta = bus1->getPhase() - bus2->getPhase();
    }
  }
}

/**
 * Return the complex admittance of the branch
 * @return: complex addmittance of branch
 */
gridpack::ComplexType gridpack::state_estimation::SEBranch::getAdmittance(void)
{
  int i;
  gridpack::ComplexType ret(0.0,0.0);
  for (i=0; i<p_elems; i++) {
    gridpack::ComplexType tmp(p_resistance[i], p_reactance[i]);
    if (!p_xform[i] && p_branch_status[i]) {
      tmp = -1.0/tmp;
    } else {
      tmp = gridpack::ComplexType(0.0,0.0);
    }
    ret += tmp;
  }
  return ret;
}

/**
 * Return transformer contribution from the branch to the calling
 * bus
 * @param bus: pointer to the bus making the call
 * @return: contribution to Y matrix from branch
 */
gridpack::ComplexType
gridpack::state_estimation::SEBranch::getTransformer(gridpack::state_estimation::SEBus *bus)
{
  int i;
  gridpack::ComplexType ret(0.0,0.0);
  for (i=0; i<p_elems; i++) {
    gridpack::ComplexType tmp(p_resistance[i],p_reactance[i]);
    gridpack::ComplexType tmpB(0.0,0.5*p_charging[i]);
    if (p_xform[i] && p_branch_status[i]) {
      tmp = -1.0/tmp;
      tmp = tmp - tmpB;
      gridpack::ComplexType a(cos(p_phase_shift[i]),sin(p_phase_shift[i]));
      a = p_tap_ratio[i]*a;
      if (bus == getBus1().get()) {
        tmp = tmp/(conj(a)*a);
      } else if (bus == getBus2().get()) {
        // tmp is unchanged
      }
    } else {
      tmp = gridpack::ComplexType(0.0,0.0);
    }
    ret += tmp;
  }
  return ret;
}

/**
 * Return the contribution to a bus from shunts
 * @param bus: pointer to the bus making the call
 * @return: contribution to Y matrix from shunts associated with branches
 */
gridpack::ComplexType
gridpack::state_estimation::SEBranch::getShunt(gridpack::state_estimation::SEBus *bus)
{
  double retr, reti;
  retr = 0.0;
  reti = 0.0;
  int i;
  for (i=0; i<p_elems; i++) {
    double tmpr, tmpi;
    if (p_shunt[i] && p_branch_status[i]) {
      tmpr = 0.0;
      tmpi = 0.0;
      if (!p_xform[i]) {
        tmpi = 0.5*p_charging[i];
        tmpr = 0.0;
      }
      // HACK: pointer comparison, maybe could handle this better
      if (bus == getBus1().get()) {
        tmpr += p_shunt_admt_g1[i];
        tmpi += p_shunt_admt_b1[i];
      } else if (bus == getBus2().get()) {
        tmpr += p_shunt_admt_g2[i];
        tmpi += p_shunt_admt_b2[i];
      } else {
        // TODO: Some kind of error
      }
    } else {
      tmpr = 0.0;
      tmpi = 0.0;
    }
    retr += tmpr;
    reti += tmpi;
  }
  return gridpack::ComplexType(retr,reti);
}

/**
 * Write output from branches to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if branch is contributing string to output, false otherwise
 */
bool gridpack::state_estimation::SEBranch::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  gridpack::ComplexType v1, v2, y, s;
  gridpack::state_estimation::SEBus *bus1 = 
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
  bool ok = !bus1->isIsolated();
  ok = ok && !bus2->isIsolated();
  if (ok) {
    v1 = bus1->getComplexVoltage();
    v2 = bus2->getComplexVoltage();
    y = gridpack::ComplexType(p_ybusr_frwd,p_ybusi_frwd);
    s = -v1*conj(y*(v1-v2));
    double p = real(s);
    double q = imag(s);
  //  double pi = 4.0*atan(1.0);
  //  double angle = p_a*180.0/pi;
    if (signal == NULL) {
      snprintf(string, bufsize, "     %6d      %6d      %12.6f         %12.6f\n",
          bus1->getOriginalIndex(),bus2->getOriginalIndex(),p,q);
    } else if (!strcmp(signal,"se")) {
      if (p_meas.size()>0) {
        int nmeas = p_meas.size();
        char buf[128];
        int ilen = 0;
        std::string meas_type,type, ckt;
        for (int i=0; i<nmeas; i++) {
          meas_type = p_meas[i].p_type;
          ckt = p_meas[i].p_ckt;
          if (meas_type.length() == 3) {
            type = meas_type;
          } else if (meas_type.length() == 2) {
            type = " ";
            type.append(meas_type);
          }
          double estimate;
          buf[0] = '\0';
          if (meas_type == "PIJ") {
            s = getComplexPower(p_meas[i].p_ckt);
            estimate = real(s)/p_sbase;
            snprintf(buf, sizeof(buf), "    %s  %8d  %8d   %s %16.5f  %16.5f   %16.5f    %8.4f\n",
              type.c_str(),bus1->getOriginalIndex(),bus2->getOriginalIndex(),ckt.c_str(),
              p_meas[i].p_value, estimate, estimate-p_meas[i].p_value,p_meas[i].p_deviation);
          } else if (meas_type == "QIJ") {
            s = getComplexPower(p_meas[i].p_ckt);
            estimate = imag(s)/p_sbase;
            snprintf(buf, sizeof(buf), "    %s  %8d  %8d   %s %16.5f  %16.5f   %16.5f    %8.4f\n",
              type.c_str(),bus1->getOriginalIndex(),bus2->getOriginalIndex(),ckt.c_str(),
              p_meas[i].p_value, estimate, estimate-p_meas[i].p_value,p_meas[i].p_deviation);
          } else if (meas_type == "IIJ") {
            s = getComplexCurrent(p_meas[i].p_ckt);
            estimate = abs(s);
            snprintf(buf, sizeof(buf), "    %s  %8d  %8d   %s %16.5f  %16.5f   %16.5f    %8.4f\n",
              type.c_str(),bus1->getOriginalIndex(),bus2->getOriginalIndex(),ckt.c_str(),
              p_meas[i].p_value, estimate, estimate-p_meas[i].p_value,p_meas[i].p_deviation);
          } else if (meas_type == "PJI") {
            s = getRvrsComplexPower(p_meas[i].p_ckt);
            estimate = real(s)/p_sbase;
            snprintf(buf, sizeof(buf), "    %s  %8d  %8d   %s %16.5f  %16.5f   %16.5f    %8.4f\n",
              type.c_str(),bus1->getOriginalIndex(),bus2->getOriginalIndex(),ckt.c_str(),
              p_meas[i].p_value, estimate, estimate-p_meas[i].p_value,p_meas[i].p_deviation);
          } else if (meas_type == "QJI") {
            s = getRvrsComplexPower(p_meas[i].p_ckt);
            estimate = imag(s)/p_sbase;
            snprintf(buf, sizeof(buf), "    %s  %8d  %8d   %s %16.5f  %16.5f   %16.5f    %8.4f\n",
                type.c_str(),bus1->getOriginalIndex(),bus2->getOriginalIndex(),ckt.c_str(),
                p_meas[i].p_value, estimate, estimate-p_meas[i].p_value,p_meas[i].p_deviation);
          } else if (meas_type == "IJI") {
            s = getRvrsComplexCurrent(p_meas[i].p_ckt);
            estimate = abs(s);
            snprintf(buf, sizeof(buf), "    %s  %8d  %8d   %s %16.5f  %16.5f   %16.5f    %8.4f\n",
              type.c_str(),bus1->getOriginalIndex(),bus2->getOriginalIndex(),ckt.c_str(),
              p_meas[i].p_value, estimate, estimate-p_meas[i].p_value,p_meas[i].p_deviation);
          }
          int buflen = strlen(buf);
          if (buflen + ilen < bufsize) {
            snprintf(string, bufsize - ilen, "%s", buf);
            string += buflen;
            ilen += buflen;
          }
        }
        if (ilen == 0) return false;
        return true;
      } else {
        return false;
      }
    }
    return true;
  } else {
    return false;
  }
}
/**
 * Add a measurement to the branch
 * @param measurement a measurement struct that will be used to
 * assign
 * internal paramters
 */
void gridpack::state_estimation::SEBranch::addMeasurement(
    gridpack::state_estimation::Measurement measurement)
{
  p_meas.push_back(measurement);
}

/**
 * Sort measurements so that they are in a consistent order
 */
void gridpack::state_estimation::SEBranch::sortMeasurements(void)
{
  // Use a multi-map to create an ordered list of measurements
  std::multimap<std::string,gridpack::state_estimation::Measurement> map;
  int nsize = p_meas.size();
  int i;
  char sbuf[32];
  for (i=0; i<nsize; i++) {
    // Create a unique (hopefully) key for each measurement
    snprintf(sbuf, sizeof(sbuf), "%s%s%f", p_meas[i].p_type, p_meas[i].p_ckt, p_meas[i].p_value);
    std::string key = sbuf;
    map.insert(std::pair<std::string,
        gridpack::state_estimation::Measurement>(key,p_meas[i]));
  }
  p_meas.clear();
  std::multimap<std::string,gridpack::state_estimation::Measurement>::iterator it;
  it = map.begin();
  while (it != map.end()) {
    p_meas.push_back(it->second);
    it++;
  }
}

/**
  * Return contribution to constraints
  * @param v: voltage at the other bus
  * @param theta: angle difference between two buses
  */
void gridpack::state_estimation::SEBranch::getVTheta(gridpack::state_estimation::SEBus *bus, double *v, double *theta)
{
  gridpack::state_estimation::SEBus *bus1 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());

  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
  bool ok = !bus1->isIsolated();
  ok = ok && !bus2->isIsolated();
  if (ok) {
    double v1 = bus1->getVoltage();
    double v2 = bus2->getVoltage();
    if (bus == bus1) {
       *v = v2;
       *theta = bus1->getPhase() - bus2->getPhase();  
    }  else if (bus == bus2) {
       *v = v1;
       *theta = bus2->getPhase() - bus1->getPhase();  
    }
 }
}
 
/**
  * Return contribution to constraints
  * @param v1, v2: voltages at buses
  * @param theta: angle difference between two buses
  */
void gridpack::state_estimation::SEBranch::getV1V2Theta(gridpack::state_estimation::SEBranch *branch, double *v1, double *v2, double *theta)
{
  gridpack::state_estimation::SEBus *bus1 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());

  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
  bool ok = !bus1->isIsolated();
  ok = ok && !bus2->isIsolated();
  if (ok) {
    *v1 = bus1->getVoltage();
    *v2 = bus2->getVoltage();
    *theta = bus1->getPhase() - bus2->getPhase();  
  }
}
 
/**
 * Configure branches with state estimation parameters. These can be
 * used in other methods
 */
void gridpack::state_estimation::SEBranch::configureSE(void)
{
  // Calculate the number of matrix elements associated with this branch
  int reference = 1; // TBD: to be read from XML
  gridpack::state_estimation::SEBus *bus1 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());

  int nmeas = p_meas.size(); // Suppose p_meas is the vector of all the measurements on this branch
  int ncnt = 0;
  int i, j, im, jm, nsize;
  bool ok = !bus1->isIsolated();
  ok = ok && !bus2->isIsolated();
  p_numElements = 0;
  if (ok) {
   for (i=0; i<nmeas; i++) {
    std::string type = p_meas[i].p_type;
    std::string ckt = p_meas[i].p_ckt;
    if (type == "PIJ" || type == "QIJ" || type == "IIJ" || type == "PJI" || type == "QJI" || type == "IJI") {
      int nsize = p_tag.size();
      for (j=0; j<nsize; j++) {
        if (circuitIDsMatch(p_tag[j], ckt) && p_active) {
          if (!bus1->getReferenceBus()) {
            ncnt += 2;
          } else {  // reference bus, only for dPIJ/DVI
            ncnt++;
          }
          if (!bus2->getReferenceBus()) {
            ncnt += 2;
          } else {  // reference bus, only for dPIJ/DVJ
            ncnt++;
          }
        } 
      }
    }
  p_numElements = ncnt;
   }
  }
}

/**
 * Return number of rows in matrix from component
 * @return number of rows from component
 */
int gridpack::state_estimation::SEBranch::matrixNumRows() const
{
  return p_meas.size();
}

/**
 * Return number of cols in matrix from component
 * @return number of cols from component
 */
int gridpack::state_estimation::SEBranch::matrixNumCols() const
{
  if (p_mode == Jacobian_H) {
    return 0;
  } else if (p_mode == R_inv) {
    return p_meas.size();
  }
  return 0;
}

/**
 * Set row indices corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @param idx matrix index of row irow
 */
void gridpack::state_estimation::SEBranch::matrixSetRowIndex(int irow, int idx)
{
  if (p_mode == Jacobian_H) {
    if (irow < p_rowJidx.size()) {
      p_rowJidx[irow] = idx;
    } else {
      p_rowJidx.push_back(idx);
    }
  } else if (p_mode == R_inv) {
    if (irow < p_rowRidx.size()) {
      p_rowRidx[irow] = idx;
    } else {
      p_rowRidx.push_back(idx);
    }
  }
}

/**
 * Set column indices corresponding to the columns contributed by this component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @param idx matrix index of column icol
 */
void gridpack::state_estimation::SEBranch::matrixSetColIndex(int icol, int idx)
{
  if (p_mode == Jacobian_H) {
    if (icol < p_colJidx.size()) {
      p_colJidx[icol] = idx;
    } else {
      p_colJidx.push_back(idx);
    }
  } else if (p_mode == R_inv) {
    if (icol < p_colRidx.size()) {
      p_colRidx[icol] = idx;
    } else {
      p_colRidx.push_back(idx);
    }
  }
}

/**
 * Get the row index corresponding to the rows contributed by this component
 * @param irow index of row contributed by this component (e.g. if component
 * contributes 3 rows then irow is between 0 and 2)
 * @return matrix index of row irow
 */
int gridpack::state_estimation::SEBranch::matrixGetRowIndex(int idx)
{
  if (p_mode == Jacobian_H) {
    return p_rowJidx[idx];
  } else if (p_mode == R_inv) {
    return p_rowRidx[idx];
  }
  return -1;
}

/**
 * Get the column index corresponding to the columns contributed by this component
 * @param icol index of column contributed by this component (e.g. if component
 * contributes 3 columns then icol is between 0 and 2)
 * @return matrix index of column icol
 */
int gridpack::state_estimation::SEBranch::matrixGetColIndex(int idx)
{
  if (p_mode == Jacobian_H) {
    if (idx >= p_colJidx.size())
      printf("violation in branch:matrixGetColIndex branch: %d %d size: %d idx: %d\n",
          dynamic_cast<SEBus*>(getBus1().get())->getOriginalIndex(),
          dynamic_cast<SEBus*>(getBus2().get())->getOriginalIndex(),idx,
          static_cast<int>(p_colJidx.size()));
    return p_colJidx[idx];
  } else if (p_mode == R_inv) {
    return p_colRidx[idx];
  }
  return -1;
}

/**
 * Return the number of matrix values contributed by this component
 * @return number of matrix values
 */
int gridpack::state_estimation::SEBranch::matrixNumValues() const
{
  if (p_mode == Jacobian_H) {
    return p_numElements;
  } else if (p_mode == R_inv) {
    return p_meas.size();
  }
  return 0;
}

/**
 * Return values from a matrix block
 * @param nvals: number of values inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
*/
void gridpack::state_estimation::SEBranch::matrixGetValues(int *nvals,ComplexType *values, int *rows, int *cols)
{
      boost::shared_ptr<gridpack::component::BaseComponent> b1 = getBus1();
      SEBus *bus1 = dynamic_cast<SEBus*>(b1.get());
      boost::shared_ptr<gridpack::component::BaseComponent> b2 = getBus2();
      SEBus *bus2 = dynamic_cast<SEBus*>(b2.get());
  bool ok = !bus1->isIsolated();
  ok = ok && !bus2->isIsolated();
  if (ok) {
  if (p_mode == Jacobian_H) {
    int nmeas = p_meas.size(); // Suppose p_meas is the vector of all the measurements on this branch
    int ncnt = 0;
    int i, j, im, jm, nsize;
    double v1, v2, theta;
    double t,gij,bij,shifter; 
    std::string ckt, type;

    v1 = bus1->getVoltage();
    v2 = bus2->getVoltage();
    double gs1, bs1, gs2, bs2;
    bus1->getShuntGsBs(&gs1,&bs1);
    bus2->getShuntGsBs(&gs2,&bs2);
    theta = bus1->getPhase() - bus2->getPhase();  

    for (i=0; i<nmeas; i++) {
      im = matrixGetRowIndex(i);
      ckt = p_meas[i].p_ckt;
      type = p_meas[i].p_type;
      bool found = false;
      if (type == "PIJ") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (circuitIDsMatch(p_tag[j], ckt) && p_active) {
            found = true;
            gridpack::ComplexType ret(p_resistance[j],p_reactance[j]);
            double B=0.5*p_charging[j];
            ret = 1.0/ret;
            gij=real(ret);
            bij=imag(ret);
            if (p_tap_ratio[j] != 0.0) {
              shifter=p_phase_shift[j];
              t = p_tap_ratio[j];
            } else {
              shifter=0.0;
              t=1.0;
            }

            if (!bus1->getReferenceBus()) {
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(1.0/t*v1*v2*(gij*sin(theta+shifter)
                    -bij*cos(theta+shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus1->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v2*(gij*cos(theta+shifter)
                  +bij*sin(theta+shifter)) +2*(gij+gs1)*v1/(t*t),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dPIJ/DVI
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v2*(gij*cos(theta+shifter)
                  +bij*sin(theta+shifter)) +2*(gij+gs1)*v1/(t*t),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
            if (!bus2->getReferenceBus()) {
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v1*v2*(gij*sin(theta+shifter)
                    -bij*cos(theta+shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus2->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v1*(gij*cos(theta+shifter)
                    +bij*sin(theta+shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dPIJ/DVJ
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v1*(gij*cos(theta+shifter)
                    +bij*sin(theta+shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
          } 
        }
      } else if (type == "QIJ") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (circuitIDsMatch(p_tag[j], ckt) && p_active) {
            found = true;
            gridpack::ComplexType ret(p_resistance[j],p_reactance[j]);
            double B=0.5*p_charging[j];
            ret = 1.0/ret;
            gij=real(ret);
            bij=imag(ret);
            if (p_tap_ratio[j] != 0.0) {
              shifter=p_phase_shift[j];
              t = p_tap_ratio[j];
            } else {
              shifter=0.0;
              t=1.0;
            }
            if (!bus1->getReferenceBus()) {
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v1*v2*(gij*cos(theta+shifter)
                     + bij*sin(theta+shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus1->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v2*(gij * sin(theta+shifter)
                    - bij*cos(theta+shifter))-2/(t*t)*(bij+B)*v1,0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dQIJ/DVI
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v2*(gij * sin(theta+shifter)
                    - bij*cos(theta+shifter))-2/(t*t)*(bij+B)*v1,0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
            if (!bus2->getReferenceBus()) {
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(1.0/t*v1*v2*(gij*cos(theta+shifter)
                    + bij*sin(theta+shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus2->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v1*(gij*sin(theta+shifter)
                    - bij*cos(theta+shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dQIJ/DVJ
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v1*(gij*sin(theta+shifter)
                    - bij*cos(theta+shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
          } 
        }
      } else if (type == "IIJ") { 
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (circuitIDsMatch(p_tag[j], ckt) && p_active) {
            found = true;
            gridpack::ComplexType ret(p_resistance[j],p_reactance[j]);
            double B=0.5*p_charging[j];
            ret = 1.0/ret;
            gij=real(ret);
            bij=imag(ret);
            if (p_tap_ratio[j] != 0.0) {
              shifter=p_phase_shift[j];
              t = p_tap_ratio[j];
            } else {
              shifter=0.0;
              t=1.0;
            }

            double delta1 = bus1->getPhase();
            double delta2 = bus2->getPhase();

            //double Iij = sqrt((gij*gij+bij*bij) *(v1*v1+v2*v2-2*v1*v2*cos(theta))); 
            // Define self Admittance (Bus 1 side) 

            gridpack::ComplexType Yij_series(gij, bij);
            gridpack::ComplexType Yij_prime(gij/(t*t), bij/(t*t) + B); // Y'_ij = Y_ij/t^2 + jB
            gridpack::ComplexType a_conj(t * cos(-shifter), t * sin(-shifter));
            gridpack::ComplexType Yij_mutual = Yij_series / a_conj;

            // V1 and V2 in rectangular form
            gridpack::ComplexType V1(v1 * cos(delta1), v1 * sin(delta1));
            gridpack::ComplexType V2(v2 * cos(delta2), v2 * sin(delta2));

            // Calculate complex current I_ij = Y'_ij * V1 - Y_ij * V2
            gridpack::ComplexType Iij = Yij_prime * V1 - Yij_mutual * V2;
            double Iij_mag = abs(Iij); 
      
            // Check for zero division (unlikely in state estimation)
            //if (Iij_mag < 1.0e-9) Iij_mag = 1.0e-9;

            gridpack::ComplexType V1_div_v1(cos(delta1), sin(delta1)); // V2 / |V1|
            gridpack::ComplexType V2_div_v2(cos(delta2), sin(delta2)); // V2 / |V2|

            if (!bus1->getReferenceBus()) {
              // dIij/d(delta_i) term (H_Iij, delta_i)
              gridpack::ComplexType dI_d_delta1 = Yij_prime * gridpack::ComplexType(0.0, 1.0) * V1;
              gridpack::ComplexType H_d_delta1 = (conj(Iij) * dI_d_delta1) / Iij_mag;
              jm = bus1->matrixGetColIndex(0);
	      if (Iij_mag < 1.0e-4) {
		values[ncnt] = gridpack::ComplexType(0.0,0.0);
	      } else {
                values[ncnt] = gridpack::ComplexType(real(H_d_delta1),0.0);
              } 

              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;

              // dIij/d(v_i) term (H_Iij, v_i)
              gridpack::ComplexType dI_d_v1 = Yij_prime * V1_div_v1;
              gridpack::ComplexType H_d_v1 = (conj(Iij) * dI_d_v1) / Iij_mag;
              jm = bus1->matrixGetColIndex(1);
	      if (Iij_mag < 1.0e-4) {
		values[ncnt] = gridpack::ComplexType(0.0,0.0);
	      } else {
                values[ncnt] = gridpack::ComplexType(real(H_d_v1),0.0);
              } 

              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dIIJ/DVI
              gridpack::ComplexType dI_d_v1 = Yij_prime * V1_div_v1;
              gridpack::ComplexType H_d_v1 = (conj(Iij) * dI_d_v1) / Iij_mag;
              jm = bus1->matrixGetColIndex(0);
	      if (Iij_mag < 1.0e-4) {
		values[ncnt] = gridpack::ComplexType(0.0,0.0);
	      } else {
                values[ncnt] = gridpack::ComplexType(real(H_d_v1),0.0);
              } 

              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }

            // dIij/d(delta_j) term (H_Iij, delta_j)
            if (!bus2->getReferenceBus()) {
              gridpack::ComplexType dI_d_delta2 = -Yij_mutual * gridpack::ComplexType(0.0, 1.0) * V2;
              gridpack::ComplexType H_d_delta2 = (conj(Iij) * dI_d_delta2) / Iij_mag;
              jm = bus2->matrixGetColIndex(0);
	      if (Iij_mag < 1.0e-4) {
		values[ncnt] = gridpack::ComplexType(0.0,0.0);
	      } else {
                values[ncnt] = gridpack::ComplexType(real(H_d_delta2),0.0);
              } 

              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;

              // dIij/d(v_j) term (H_Iij, v_j)
              gridpack::ComplexType dI_d_v2 = -Yij_mutual * V2_div_v2;
              gridpack::ComplexType H_d_v2 = (conj(Iij) * dI_d_v2) / Iij_mag;
              jm = bus2->matrixGetColIndex(1);
              if (Iij_mag < 1.0e-4) {
                values[ncnt] = gridpack::ComplexType(0.0,0.0);
              } else {
                values[ncnt] = gridpack::ComplexType(real(H_d_v2),0.0);
              }
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;

            } else {  // reference bus, only for dIIJ/DVJ
              // dIij/d(v_j) term (H_Iij, v_j)
              gridpack::ComplexType dI_d_v2 = -Yij_mutual * V2_div_v2;
              gridpack::ComplexType H_d_v2 = (conj(Iij) * dI_d_v2) / Iij_mag;
              jm = bus2->matrixGetColIndex(0);
              if (Iij_mag < 1.0e-4) {
                values[ncnt] = gridpack::ComplexType(0.0,0.0);
              } else {
                values[ncnt] = gridpack::ComplexType(real(H_d_v2),0.0);
              }
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
          } 
        }
      } else if (type == "PJI") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (circuitIDsMatch(p_tag[j], ckt) && p_active) {
            found = true;
            gridpack::ComplexType ret(p_resistance[j],p_reactance[j]);
            double B=0.5*p_charging[j];
            ret = 1.0/ret;
            gij=real(ret);
            bij=imag(ret);
            if (p_tap_ratio[j] != 0.0) {
              shifter=p_phase_shift[j];
              t = p_tap_ratio[j];
            } else {
              shifter=0.0;
              t=1.0;
            }

            if (!bus1->getReferenceBus()) {
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v1*v2*(gij*sin(-theta-shifter)
                    -bij*cos(-theta-shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus1->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v2*(gij*cos(-theta-shifter)
                  +bij*sin(-theta-shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dPIJ/DVI
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v2*(gij*cos(-theta-shifter)
                  +bij*sin(-theta-shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
            if (!bus2->getReferenceBus()) {
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(1.0/t*v1*v2*(gij*sin(-theta-shifter)
                    -bij*cos(-theta-shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus2->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v1*(gij*cos(-theta-shifter)
                    +bij*sin(-theta-shifter))+2*(gij+gs2)*v2,0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dPIJ/DVJ
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v1*(gij*cos(-theta-shifter)
                    +bij*sin(-theta-shifter))+2*(gij+gs2)*v2,0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
          } 
        }
      } else if (type == "QJI") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (circuitIDsMatch(p_tag[j], ckt) && p_active) {
            found = true;
            gridpack::ComplexType ret(p_resistance[j],p_reactance[j]);
            double B=0.5*p_charging[j];
            ret = 1.0/ret;
            gij=real(ret);
            bij=imag(ret);
            if (p_tap_ratio[j] != 0.0) {
              shifter=p_phase_shift[j];
              t = p_tap_ratio[j];
            } else {
              shifter=0.0;
              t=1.0;
            }
            if (!bus1->getReferenceBus()) {
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(1.0/t*v1*v2*(gij*cos(-theta-shifter)
                     + bij*sin(-theta-shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus1->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v2*(gij*sin(-theta-shifter)
                    - bij*cos(-theta-shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dQIJ/DVI
              jm = bus1->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v2*(gij*sin(-theta-shifter)
                    - bij*cos(-theta-shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
            if (!bus2->getReferenceBus()) {
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v1*v2*(gij*cos(-theta-shifter)
                    + bij*sin(-theta-shifter)),0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
              jm = bus2->matrixGetColIndex(1);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v1*(gij*sin(-theta-shifter)
                    - bij*cos(-theta-shifter))-2*(bij+B)*v2,0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            } else {  // reference bus, only for dQIJ/DVJ
              jm = bus2->matrixGetColIndex(0);
              values[ncnt] = gridpack::ComplexType(-1.0/t*v1*(gij*sin(-theta-shifter)
                    - bij*cos(-theta-shifter))-2*(bij+B)*v2,0.0);
              rows[ncnt] = im;
              cols[ncnt] = jm;
              ncnt++;
            }
          } 
        }
      } else if (type == "IJI") {  // DONE !
	int nsize = p_tag.size();
	for (j=0; j<nsize; j++) {
	  if (circuitIDsMatch(p_tag[j], ckt) && p_active) {
	    found = true;
	    gridpack::ComplexType ret(p_resistance[j],p_reactance[j]);
	    double B = 0.5*p_charging[j];
	    ret = 1.0/ret; // ret is series admittance Y_ij = Y_ji = G + jB
	    gij=real(ret);
	    bij=imag(ret);
	    if (p_tap_ratio[j] != 0.0) {
	      shifter=p_phase_shift[j];
	      t = p_tap_ratio[j];
	    } else {
	      shifter=0.0;
	      t=1.0;
	    }
	    
	    // Phase angle calculation
	    double delta1 = bus1->getPhase();
	    double delta2 = bus2->getPhase();
	    
	    // Define Admittances
	    gridpack::ComplexType Yji_series(gij, bij); // Series admittance Y_ji
	    gridpack::ComplexType Yji_prime(gij, bij + B); // Y'_ji = Y_ji + jB (since tap is on bus i side)
            gridpack::ComplexType a(t * cos(shifter), t * sin(shifter));
            gridpack::ComplexType Yji_mutual = Yji_series / a;

	    // V1 and V2 in rectangular form
	    gridpack::ComplexType V1(v1 * cos(delta1), v1 * sin(delta1));
	    gridpack::ComplexType V2(v2 * cos(delta2), v2 * sin(delta2));

	    // Calculate complex current I_ji = Y'_ji * V2 - Y_ji * V1
	    gridpack::ComplexType Iji = Yji_prime * V2 - Yji_mutual * V1;
	    double Iji_mag = abs(Iji); 
	    
	    // Check for zero division
	    //if (Iji_mag < 1.0e-9) Iji_mag = 1.0e-9; 
	    
	    // Calculate derivatives of I_ji and Jacobian elements
	    gridpack::ComplexType V1_div_v1(cos(delta1), sin(delta1)); // V1 / |V1|
	    gridpack::ComplexType V2_div_v2(cos(delta2), sin(delta2)); // V2 / |V2|

	    // dIji/d(delta_i) term (H_Iji, delta_i) - Bus 1
	    if (!bus1->getReferenceBus()) {
	      gridpack::ComplexType dI_d_delta1 = -Yji_mutual * gridpack::ComplexType(0.0, 1.0) * V1;
	      gridpack::ComplexType H_d_delta1 = (conj(Iji) * dI_d_delta1) / Iji_mag;
	      jm = bus1->matrixGetColIndex(0);

	      if (Iji_mag < 1.0e-4) {
		values[ncnt] = gridpack::ComplexType(0.0,0.0);
	      } else {
                values[ncnt] = gridpack::ComplexType(real(H_d_delta1),0.0);
              } 
	      rows[ncnt] = im;
	      cols[ncnt] = jm;
	      ncnt++;

	      // dIji/d(v_i) term (H_Iji, v_i)
	      gridpack::ComplexType dI_d_v1 = -Yji_mutual * V1_div_v1;
	      gridpack::ComplexType H_d_v1 = (conj(Iji) * dI_d_v1) / Iji_mag;
	      jm = bus1->matrixGetColIndex(1);
	      if (Iji_mag < 1.0e-4) {
	        values[ncnt] = gridpack::ComplexType(0.0,0.0);
	      } else {
	        values[ncnt] = gridpack::ComplexType(real(H_d_v1),0.0);
	      }
	      rows[ncnt] = im;
	      cols[ncnt] = jm;
	      ncnt++;
	    } else {  // reference bus, only for dIJI/DVI
	      // dIji/d(v_i) term (H_Iji, v_i)
	      gridpack::ComplexType dI_d_v1 = -Yji_mutual * V1_div_v1;
	      gridpack::ComplexType H_d_v1 = (conj(Iji) * dI_d_v1) / Iji_mag;
	      jm = bus1->matrixGetColIndex(0);
	      if (Iji_mag < 1.0e-4) {
		values[ncnt] = gridpack::ComplexType(0.0,0.0);
	      } else {
                values[ncnt] = gridpack::ComplexType(real(H_d_v1),0.0);
              } 
	      rows[ncnt] = im;
	      cols[ncnt] = jm;
	      ncnt++;
	    }

	    // dIji/d(delta_j) term (H_Iji, delta_j) - Bus 2
	    if (!bus2->getReferenceBus()) {
	      gridpack::ComplexType dI_d_delta2 = Yji_prime * gridpack::ComplexType(0.0, 1.0) * V2;
	      gridpack::ComplexType H_d_delta2 = (conj(Iji) * dI_d_delta2) / Iji_mag;
	      jm = bus2->matrixGetColIndex(0);
	      if (Iji_mag < 1.0e-4) {
	        values[ncnt] = gridpack::ComplexType(0.0,0.0);
	      } else {
	        values[ncnt] = gridpack::ComplexType(real(H_d_delta2),0.0);
	      }
	      rows[ncnt] = im;
	      cols[ncnt] = jm;
	      ncnt++;

	      // dIji/d(v_j) term (H_Iji, v_j)
	      gridpack::ComplexType dI_d_v2 = Yji_prime * V2_div_v2;
	      gridpack::ComplexType H_d_v2 = (conj(Iji) * dI_d_v2) / Iji_mag;
	      jm = bus2->matrixGetColIndex(1);
	      if (Iji_mag < 1.0e-4) {
		values[ncnt] = gridpack::ComplexType(0.0,0.0);
	      } else {
	        values[ncnt] = gridpack::ComplexType(real(H_d_v2),0.0);
              }
	      rows[ncnt] = im;
	      cols[ncnt] = jm;
	      ncnt++;
	    } else {  // reference bus, only for dIJI/DVJ
	      // dIji/d(v_j) term (H_Iji, v_j)
	      gridpack::ComplexType dI_d_v2 = Yji_prime * V2_div_v2;
	      gridpack::ComplexType H_d_v2 = (conj(Iji) * dI_d_v2) / Iji_mag;
	      jm = bus2->matrixGetColIndex(0);
	      if (Iji_mag < 1.0e-4) {
		values[ncnt] = gridpack::ComplexType(0.0,0.0);
	      } else {
	        values[ncnt] = gridpack::ComplexType(real(H_d_v2),0.0);
              }
	      rows[ncnt] = im;
	      cols[ncnt] = jm;
	      ncnt++;
	    }
	  }
	}
      }
      if (!found) {
        printf("No match found for branch measurement\n   type: %s\n"
               "   branch: %d %d\n   ckt_id: %s\n   value: %f\n"
               "   deviation: %f\n",p_meas[i].p_type,p_meas[i].p_fbusid,
               p_meas[i].p_tbusid,p_meas[i].p_ckt,p_meas[i].p_value,
               p_meas[i].p_deviation);

      }
    }
    *nvals = ncnt;
  } else if (p_mode == R_inv) {
    int nsize = p_meas.size();
    int i;
    for (i=0; i<nsize; i++) {
      if (p_meas[i].p_deviation != 0.0) {
        values[i] = 1.0/(p_meas[i].p_deviation*p_meas[i].p_deviation);
      } else {
        values[i] = 0.0;
      }
      rows[i] = matrixGetRowIndex(i);
      cols[i] = matrixGetColIndex(i);
    }
    *nvals = nsize;
  }
  }
}

/**
 * Return number of elements in vector coming from component
 * @return number of elements contributed from component
 */

int gridpack::state_estimation::SEBranch::vectorNumElements() const
{
  if (p_mode == Jacobian_H || p_mode == Residual) {
    // When computing residuals or Jacobian, we need the measurement count
    return p_meas.size();
  }
  return 0;
}

/**
 * Set indices corresponding to the elements contributed by this component
 * @param ielem index of element contributed by this component (e.g.
 * if component contributes 3 elements then ielem is between 0 and 2)
 * @param idx vector index of element ielem
 */
void gridpack::state_estimation::SEBranch::vectorSetElementIndex(int ielem, int idx)
{
  if (p_mode == Jacobian_H || p_mode == Residual) {
    if (ielem < p_vecZidx.size()) {
      p_vecZidx[ielem] = idx;
    } else {
      p_vecZidx.push_back(idx);
    }
  }
}

/**
 * Get list of element indices from component
 * @param idx list of indices that component maps onto
 */
void gridpack::state_estimation::SEBranch::vectorGetElementIndices(int *idx)
{
  if (p_mode == Jacobian_H || p_mode == Residual) {
    int nsize = p_vecZidx.size();
    int i;
    for (i=0; i<nsize; i++) {
      idx[i] = p_vecZidx[i];
    }
  }
}

/**
 * Transfer vector values to component
 * @param values list of vector element values
 */
void gridpack::state_estimation::SEBranch::vectorSetElementValues(ComplexType *values)
{
  //TODO: PLACE HOLDER
}


/**
 * Return values from a vector
 * @param values: pointer to vector values (z-h(x))
 * @param idx: pointer to vector index 
*/
void gridpack::state_estimation::SEBranch::vectorGetElementValues(ComplexType *values, int *idx)
{
  if (p_mode == Jacobian_H) {
    gridpack::state_estimation::SEBus *bus1 =
      dynamic_cast<gridpack::state_estimation::SEBus*>(this->getBus1().get());
    gridpack::state_estimation::SEBus *bus2 =
      dynamic_cast<gridpack::state_estimation::SEBus*>(this->getBus2().get());
    bool ok = !bus1->isIsolated();
    ok = ok && !bus2->isIsolated();
    if (ok) {
    int nmeas = p_meas.size(); // Suppose p_meas is the vector of all the measurements on this branch
    int ncnt = 0;
    int i, j, im, jm, nsize;
    double ret1=0.0;
    double ret2=0.0;
    double ret3=0.0;
    double v1=0.0;
    double v2=0.0;
    double theta=0.0;
    vectorGetElementIndices(idx);
    v1 = bus1->getVoltage();
    v2 = bus2->getVoltage();
    theta = bus1->getPhase() - bus2->getPhase();  
    for (i=0; i<nmeas; i++) {
      std::string type = p_meas[i].p_type;
      int idx1, idx2;
      double gij,bij,t,shifter;
      idx1 = bus1->getOriginalIndex();
      idx2 = bus2->getOriginalIndex();
      gridpack::ComplexType a(1.0,0.0);
//      printf("branch %d %d type: %s row: %d\n",idx1,idx2,type.c_str(),idx[i]);
      if (type == "PIJ") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (circuitIDsMatch(p_tag[j], p_meas[i].p_ckt)) {
            gridpack::ComplexType ret(p_resistance[j],p_reactance[j]);
            ret = 1.0/ret;
            gij=real(ret);
            bij=imag(ret);
            if (p_tap_ratio[j] != 0.0) {
              shifter=p_phase_shift[j];
              t = p_tap_ratio[j];
            } else {
              shifter=0.0;
              t=1.0;
            }
            ret1 =  v1*v1*gij/(t*t) - 1.0/t*v1*v2*(gij*cos(theta+shifter) + bij*sin(theta+shifter));
          }
        }
        values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret1),0.0);
        ncnt++;
      } else if (type == "QIJ") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (circuitIDsMatch(p_tag[j], p_meas[i].p_ckt)) {
            gridpack::ComplexType ret(p_resistance[j],p_reactance[j]);
            double B=0.5*p_charging[j];
            ret = 1.0/ret;
            gij=real(ret);
            bij=imag(ret);
            if (p_tap_ratio[j] != 0.0) {
              shifter=p_phase_shift[j];
              t = p_tap_ratio[j];
            } else {
              shifter=0.0;
              t=1.0;
            }
            ret2 = - v1*v1* (bij+B)/(t*t) - 1.0/t*v1*v2*(gij*sin(theta+shifter) - bij*cos(theta+shifter));
          }
        }
        values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret2),0.0);
        ncnt++;
      } else if (type == "IIJ") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (circuitIDsMatch(p_tag[j], p_meas[i].p_ckt)) {
            gridpack::ComplexType ret(p_resistance[j],p_reactance[j]);
            double B=0.5*p_charging[j];
            ret = 1.0/ret;
            gij=real(ret);
            bij=imag(ret);
            if (p_tap_ratio[j] != 0.0) {
              t = p_tap_ratio[j];
              shifter = p_phase_shift[j];
            } else {
              t=1.0;
              shifter = 0.0;
            }

            // Calculate I_ij using complex form matching Jacobian (line 2536)
            double delta1 = bus1->getPhase();
            double delta2 = bus2->getPhase();
            gridpack::ComplexType V1(v1 * cos(delta1), v1 * sin(delta1));
            gridpack::ComplexType V2(v2 * cos(delta2), v2 * sin(delta2));
            gridpack::ComplexType Yij_prime(gij/(t*t), bij/(t*t) + B);
            gridpack::ComplexType Y_series(gij, bij);
            gridpack::ComplexType a_conj(t * cos(-shifter), t * sin(-shifter));
            gridpack::ComplexType Yij_mutual = Y_series / a_conj;
            gridpack::ComplexType Iij = Yij_prime * V1 - Yij_mutual * V2;
            ret3 = abs(Iij);
          }
        }
        values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret3),0.0);
        ncnt++;
      } else if (type == "PJI") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (circuitIDsMatch(p_tag[j], p_meas[i].p_ckt)) {
            gridpack::ComplexType ret(p_resistance[j],p_reactance[j]);
            ret = 1.0/ret;
            gij=real(ret);
            bij=imag(ret);
            if (p_tap_ratio[j] != 0.0) {
              shifter=p_phase_shift[j];
              t = p_tap_ratio[j];
            } else {
              shifter=0.0;
              t=1.0;
            }
            ret1 =  v2*v2*gij - 1.0/t*v1*v2*(gij*cos(-theta-shifter) + bij*sin(-theta-shifter));
          }
        }
        values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret1),0.0);
        ncnt++;
      } else if (type == "QJI") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (circuitIDsMatch(p_tag[j], p_meas[i].p_ckt)) {
            gridpack::ComplexType ret(p_resistance[j],p_reactance[j]);
            double B=0.5*p_charging[j];
            ret = 1.0/ret;
            gij=real(ret);
            bij=imag(ret);
            if (p_tap_ratio[j] != 0.0) {
              shifter=p_phase_shift[j];
              t = p_tap_ratio[j];
            } else {
              shifter=0.0;
              t=1.0;
            }
            ret2 = - v2*v2* (bij+B) - 1.0/t*v1*v2*(gij*sin(-theta-shifter) - bij*cos(-theta-shifter));
          }
        }
        values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret2),0.0);
        ncnt++;
      } else if (type == "IJI") {
        int nsize = p_tag.size();
        for (j=0; j<nsize; j++) {
          if (circuitIDsMatch(p_tag[j], p_meas[i].p_ckt)) {
            gridpack::ComplexType ret(p_resistance[j],p_reactance[j]);
            double Bch=0.5*p_charging[j];
            ret = 1.0/ret;
            gij=real(ret);
            bij=imag(ret);
            if (p_tap_ratio[j] != 0.0) {
              t = p_tap_ratio[j];
              shifter = p_phase_shift[j];
            } else {
              t=1.0;
              shifter = 0.0;
            }

            // Calculate I_ji using complex form matching Jacobian (line 2759)
            double delta1 = bus1->getPhase();
            double delta2 = bus2->getPhase();
            gridpack::ComplexType V1(v1 * cos(delta1), v1 * sin(delta1));
            gridpack::ComplexType V2(v2 * cos(delta2), v2 * sin(delta2));
            gridpack::ComplexType Yji_prime(gij, bij + Bch);
            gridpack::ComplexType Y_series(gij, bij);
            gridpack::ComplexType a(t * cos(shifter), t * sin(shifter));
            gridpack::ComplexType Yji_mutual = Y_series / a;
            gridpack::ComplexType Iji = Yji_prime * V2 - Yji_mutual * V1;
            ret3 = abs(Iji);
          }
        }
        values[ncnt] = gridpack::ComplexType(static_cast<double>(p_meas[i].p_value-ret3),0.0);
        ncnt++;
      }
    } 
    }
  } else if (p_mode == Residual) {
    gridpack::state_estimation::SEBus *bus1 = 
      dynamic_cast<gridpack::state_estimation::SEBus*>(this->getBus1().get());
    gridpack::state_estimation::SEBus *bus2 = 
      dynamic_cast<gridpack::state_estimation::SEBus*>(this->getBus2().get());
    bool ok = !bus1->isIsolated() && !bus2->isIsolated();
    
    if (ok && p_meas.size() > 0) {
      for (int i = 0; i < p_meas.size(); i++) {
        std::string type = p_meas[i].p_type;
        std::string ckt = p_meas[i].p_ckt;
        double measured = p_meas[i].p_value;
        double estimated = 0.0;
        
        if (type == "PIJ") {
          gridpack::ComplexType s = getComplexPower(ckt);
          estimated = real(s)/p_sbase;
        } else if (type == "QIJ") {
          gridpack::ComplexType s = getComplexPower(ckt);
          estimated = imag(s)/p_sbase;
        } else if (type == "PJI") {
          gridpack::ComplexType s = getRvrsComplexPower(ckt);
          estimated = real(s)/p_sbase;
        } else if (type == "QJI") {
          gridpack::ComplexType s = getRvrsComplexPower(ckt);
          estimated = imag(s)/p_sbase;
        } else if (type == "IIJ") {
          gridpack::ComplexType s = getComplexCurrent(ckt);
          estimated = abs(s);
        } else if (type == "IJI") {
          gridpack::ComplexType s = getRvrsComplexCurrent(ckt);
          estimated = abs(s);
        }
        
        // Calculate residual
        double residual = measured - estimated;
        values[i] = ComplexType(residual, 0.0);
      }
    }
  } else if (p_mode == R_inv) {
    vectorGetElementIndices(idx);
    int nmeas = p_meas.size();

    if (nmeas > 0) {
      for (int i = 0; i < nmeas; i++) {
        double sigma = p_meas[i].p_deviation;
        if (sigma > 0.0) {
          double weight = 1.0/(sigma*sigma);
          values[i] = ComplexType(weight, 0.0);
          
        } else {
          values[i] = ComplexType(0.0, 0.0);
        }
      }
    }
  }
}
/**
 * Return the values of the vector block
 * @param values: pointer to vector values
 * @return: false if network component does not contribute vector element
 */
bool gridpack::state_estimation::SEBranch::vectorValues(ComplexType *values)
{
  // Space holder...
  return false;
}

/**
 * Return contribution to constraints
 * @param p: real part of constraint
 * @param q: imaginary part of constraint
 */
void gridpack::state_estimation::SEBranch::getPQ(gridpack::state_estimation::SEBus *bus, double *p, double *q)
{
  gridpack::state_estimation::SEBus *bus1 = 
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  double v1 = bus1->getVoltage();
  gridpack::state_estimation::SEBus *bus2 = 
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
  double v2 = bus2->getVoltage();
  double cs, sn;
  double ybusr, ybusi;
  p_theta = bus1->getPhase() - bus2->getPhase();
  if (bus == bus1) {
    cs = cos(p_theta);
    sn = sin(p_theta);
    ybusr = p_ybusr_frwd;
    ybusi = p_ybusi_frwd;
  } else if (bus == bus2) {
    cs = cos(-p_theta);
    sn = sin(-p_theta);
    ybusr = p_ybusr_rvrs;
    ybusi = p_ybusi_rvrs;
  } else {
    // TODO: Some kind of error
    return;
  }
  *p = v1*v2*(ybusr*cs+ybusi*sn);
  *q = v1*v2*(ybusr*sn-ybusi*cs);
}


/**
 * Return complex power for line element at to end
 * @param tag describing line element on branch
 * @return complex power
 */
gridpack::ComplexType gridpack::state_estimation::SEBranch::getRvrsComplexPower(
        std::string tag)
{
  gridpack::ComplexType vi, vj, Yii, Yij, s;
  s = ComplexType(0.0,0.0);
  gridpack::state_estimation::SEBus *bus1 = 
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  vj = bus1->getComplexVoltage();
  gridpack::state_estimation::SEBus *bus2 = 
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
  vi = bus2->getComplexVoltage();
  getRvrsLineElements(tag,&Yii,&Yij);
  s = vi*conj(Yii*vi+Yij*vj)*p_sbase;
  return s;
}


/**
 * Return complex power for line element
 * @param tag describing line element on branch
 * @return complex power
 */
gridpack::ComplexType gridpack::state_estimation::SEBranch::getComplexPower(
        std::string tag)
{
  gridpack::ComplexType vi, vj, Yii, Yij, s;
  s = ComplexType(0.0,0.0);
  gridpack::state_estimation::SEBus *bus1 = 
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  vi = bus1->getComplexVoltage();
  gridpack::state_estimation::SEBus *bus2 = 
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());
  vj = bus2->getComplexVoltage();
  getLineElements(tag,&Yii,&Yij);
  s = vi*conj(Yii*vi+Yij*vj)*p_sbase;
  return s;
}



/**
 * Return reverse complex current for line element (IJI)
 * Calculates I_ji (from Bus 2 to Bus 1) using explicit admittance calculation
 * to match the IJI Jacobian matrix exactly (see line 2759).
 * Formula: I_ji = Y'_ji * V_j - Y_ji * V_i
 * where Y'_ji = gij + j(bij + Bch)  [charging added]
 *       Y_ji = gij + jbij            [plain series admittance]
 * @param tag describing line element on branch
 * @return complex current I_ji
 */
gridpack::ComplexType gridpack::state_estimation::SEBranch::getRvrsComplexCurrent(
        std::string tag)
{
  gridpack::ComplexType vi, vj, Iji;
  Iji = ComplexType(0.0,0.0);

  gridpack::state_estimation::SEBus *bus1 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());

  if (bus1 && bus2) {
    // Find the line index corresponding to this tag
    int line_idx = -1;
    for (int j = 0; j < p_tag.size(); j++) {
      if (p_tag[j] == tag) {
        line_idx = j;
        break;
      }
    }

    if (line_idx >= 0) {
      // Calculate series admittance Y_series = gij + jbij
      gridpack::ComplexType ret(p_resistance[line_idx], p_reactance[line_idx]);
      double B = 0.5 * p_charging[line_idx];
      ret = 1.0 / ret;
      double gij = real(ret);
      double bij = imag(ret);
      gridpack::ComplexType Y_series(gij,bij);

      // Get tap ratio
      double t = (p_tap_ratio[line_idx] != 0.0) ? p_tap_ratio[line_idx] : 1.0;
      double shifter = p_phase_shift[line_idx]; 

      // Define admittances (tap assumed on Bus 1 side)
      // Self Admittance (Y'_ji) on non-tap side:
      gridpack::ComplexType Yji_prime(gij, bij + B);  // Self admittance with tap and charging

      gridpack::ComplexType a(t * cos(shifter), t * (sin(shifter)));
      gridpack::ComplexType Yji_mutual = Y_series / a; // Mutual admittance

      // Get complex voltages
      vi = bus1->getComplexVoltage();
      vj = bus2->getComplexVoltage();

      // Calculate I_ji 
      Iji = Yji_prime * vj - Yji_mutual * vi;
    }
  }
  return Iji;
}


/**
 * Return complex current for line element (IIJ)
 * Calculates I_ij (from Bus 1 to Bus 2) using explicit admittance calculation
 * to match the IIJ Jacobian matrix exactly (see line 2536).
 * Formula: I_ij = Y'_ij * V_i - Y_ij * V_j
 * where Y'_ij = gij/t² + j(bij/t² + B)  [tap on self admittance]
 *       Y_ij = gij + jbij                [plain series admittance]
 * @param tag describing line element on branch
 * @return complex current I_ij
 */
gridpack::ComplexType gridpack::state_estimation::SEBranch::getComplexCurrent(
        std::string tag)
{
  gridpack::ComplexType vi, vj, Iij;
  Iij = ComplexType(0.0,0.0);

  gridpack::state_estimation::SEBus *bus1 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus1().get());
  gridpack::state_estimation::SEBus *bus2 =
    dynamic_cast<gridpack::state_estimation::SEBus*>(getBus2().get());

  if (bus1 && bus2) {
    // Find the line index corresponding to this tag
    int line_idx = -1;
    for (int j = 0; j < p_tag.size(); j++) {
      if (p_tag[j] == tag) {
        line_idx = j;
        break;
      }
    }

    if (line_idx >= 0) {
      // Calculate series admittance Y_series = gij + jbij
      gridpack::ComplexType ret(p_resistance[line_idx], p_reactance[line_idx]);
      double B = 0.5 * p_charging[line_idx];
      ret = 1.0 / ret;
      double gij = real(ret);
      double bij = imag(ret);
      gridpack::ComplexType Y_series(gij, bij);

      double t = (p_tap_ratio[line_idx] != 0.0) ? p_tap_ratio[line_idx] : 1.0;
      double shifter = p_phase_shift[line_idx];

      // Define admittances (Tap assumed on Bus 1 side)
      gridpack::ComplexType Yij_prime(gij / (t * t), (bij / (t * t) + B));  // Self admittance with charging

      // Mutual Admittance correction for IIJ
      // For I_ij, the mutual term is Y_series / conj(a) = Y_series / (t * e^-j*phi)
      gridpack::ComplexType a_conj(t * cos(-shifter), t * sin(-shifter));
      gridpack::ComplexType Yij_mutual = Y_series / a_conj;

      // Get complex voltages
      vi = bus1->getComplexVoltage();
      vj = bus2->getComplexVoltage();

      // Calculate I_ij
      Iij = Yij_prime * vi - Yij_mutual * vj;
    }
  }
  return Iij;
}


/**
 * Check if residual index matches this branch and report bad data
 * @param idx: index to check
 * @param report: if true, print information about the bad data
 * @return: true if index found on this branch
 */
bool gridpack::state_estimation::SEBranch::checkResidualIndex(int idx, bool report)
{
  int nsize = p_vecZidx.size();
  for (int i = 0; i < nsize; i++) {
    if (p_vecZidx[i] == idx) {
      if (report) {
      boost::shared_ptr<gridpack::component::BaseComponent> b1 = getBus1();
      SEBus *bus1 = dynamic_cast<SEBus*>(b1.get());
      boost::shared_ptr<gridpack::component::BaseComponent> b2 = getBus2();
      SEBus *bus2 = dynamic_cast<SEBus*>(b2.get());
        printf("Bad data detected on branch from bus %d to bus %d, circuit %s, " 
               "measurement type: %s, value: %f, deviation: %f\n",
               bus1->getOriginalIndex(), bus2->getOriginalIndex(), 
               p_meas[i].p_ckt, p_meas[i].p_type, 
               p_meas[i].p_value, p_meas[i].p_deviation);
      }
      return true;
    }
  }
  return false;
}

/**
 * Adjust the weight of a measurement by modifying its deviation
 */
bool gridpack::state_estimation::SEBus::adjustMeasurementWeight(
    int idx, double factor, double& oldDeviation, double& newDeviation)
{
  if (isIsolated()) return false;
  
  // Safety check for TAMU500 issue
  if (p_vecZidx.size() != p_meas.size()) {
    printf("ERROR: SEBus::adjustMeasurementWeight - size mismatch: p_vecZidx.size()=%d, p_meas.size()=%d\n",
           (int)p_vecZidx.size(), (int)p_meas.size());
    return false;
  }
  
  // Find the measurement with this index in our vector
  for (int i = 0; i < p_vecZidx.size(); i++) {
    if (p_vecZidx[i] == idx) {
      if (i >= p_meas.size()) {
        printf("ERROR: SEBus::adjustMeasurementWeight - index %d out of bounds (p_meas.size()=%d)\n",
               i, (int)p_meas.size());
        return false;
      }
      // Found the measurement, adjust its deviation
      oldDeviation = p_meas[i].p_deviation;
      p_meas[i].p_deviation *= factor;
      newDeviation = p_meas[i].p_deviation;
      return true;
    }
  }
  return false;
}

/**
 * Adjust the weight of a measurement by modifying its deviation
 */
bool gridpack::state_estimation::SEBranch::adjustMeasurementWeight(
    int idx, double factor, double& oldDeviation, double& newDeviation)
{
  if (!p_active) return false;
  
  // Safety check for TAMU500 issue
  if (p_vecZidx.size() != p_meas.size()) {
    printf("ERROR: SEBranch::adjustMeasurementWeight - size mismatch: p_vecZidx.size()=%d, p_meas.size()=%d\n",
           (int)p_vecZidx.size(), (int)p_meas.size());
    return false;
  }
  
  for (int i = 0; i < p_vecZidx.size(); i++) {
    if (p_vecZidx[i] == idx) {
      if (i >= p_meas.size()) {
        printf("ERROR: SEBranch::adjustMeasurementWeight - index %d out of bounds (p_meas.size()=%d)\n",
               i, (int)p_meas.size());
        return false;
      }
      // Found the measurement, adjust its deviation
      oldDeviation = p_meas[i].p_deviation;
      p_meas[i].p_deviation *= factor;
      newDeviation = p_meas[i].p_deviation;
      return true;
    }
  }
  return false;
}


