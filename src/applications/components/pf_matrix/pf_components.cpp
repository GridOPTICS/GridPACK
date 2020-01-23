/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   pf_components.cpp
 * @author Bruce Palmer
 * @date   2019-12-03 07:45:40 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>
#include <cstring>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "pf_components.hpp"
#include "gridpack/parser/dictionary.hpp"

//#define LARGE_MATRIX

/**
 *  Simple constructor
 */
gridpack::powerflow::PFBus::PFBus(void)
{
  p_shunt_gs = 0.0;
  p_shunt_bs = 0.0;
  p_v = 0.0;
  p_a = 0.0;
  p_theta = 0.0;
  p_angle = 0.0;
  p_voltage = 0.0;
  /*p_pl = 0.0;
  p_ql = 0.0;
  p_ip = 0.0;
  p_iq = 0.0;
  p_yp = 0.0;
  p_yq = 0.0;*/
  p_sbase = 0.0;
  p_mode = YBus;
  setReferenceBus(false);
  p_ngen = 0;
  p_data = NULL;
  p_ignore = false;
  p_vMag_ptr = NULL;
  p_vAng_ptr = NULL;
  p_PV_ptr = NULL;
}

/**
 *  Simple destructor
 */
gridpack::powerflow::PFBus::~PFBus(void)
{
}

/**
 *  Return size of matrix block contributed by the component
 *  @param isize, jsize: number of rows and columns of matrix block
 *  @return: false if network component does not contribute matrix element
 */
bool gridpack::powerflow::PFBus::matrixDiagSize(int *isize, int *jsize) const
{
  if (p_mode == Jacobian) {
    if (!isIsolated()) {
#ifdef LARGE_MATRIX
      *isize = 2;
      *jsize = 2;
      return true;
#else
      if (getReferenceBus()) {
        return false;
      } else if (p_isPV) {
        *isize = 1;
        *jsize = 1;
        return true;
      } else {
        *isize = 2;
        *jsize = 2;
        return true;
      }
#endif
    } else {
      return false;
    }
  } else if (p_mode == YBus) {
    return YMBus::matrixDiagSize(isize,jsize);
  }
  return true;
}

/**
 * Return the values of the matrix block. The values are
 * returned in row-major order.
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::powerflow::PFBus::matrixDiagValues(ComplexType *values)
{
  if (p_mode == YBus) {
    return YMBus::matrixDiagValues(values);
  } else if (p_mode == Jacobian) {
    double rvals[4];
    int nvals = diagonalJacobianValues(rvals);
    for (int i=0; i<nvals; i++) values[i] = rvals[i];
    if (nvals == 0) {
      return false;
    } else  {
      return true;
    }
  }
  return false;
}

bool gridpack::powerflow::PFBus::matrixDiagValues(RealType *values)
{
  if (p_mode == Jacobian) {
    int nvals = diagonalJacobianValues(values);
    if (nvals == 0) {
      return false;
    } else  {
      return true;
    }
  }
  return false;
}

/**
 * Return the size of the block that this component contributes to the
 * vector
 * @param size: size of vector block
 * @return: false if component does not contribute to vector
 */
bool gridpack::powerflow::PFBus::vectorSize(int *size) const
{
  if (p_mode == RHS || p_mode == State) {
    if (!isIsolated()) {
#ifdef LARGE_MATRIX
      *size = 2;
      return true;
#else
      if (getReferenceBus()) {
        return false;
      } else if (p_isPV) {
        *size = 1;
      } else {
        *size = 2;
      }
      return true;
#endif
    } else {
      return false;
    }
  } else if (p_mode == S_Cal){
    *size = 1;
  } else {
    *size = 2;
  }
  return true;
}

/**
 * Return the values of the vector block
 * @param values: pointer to vector values
 * @return: false if network component does not contribute
 *        vector element
 */
bool gridpack::powerflow::PFBus::vectorValues(ComplexType *values)
{
  if (p_mode == S_Cal)  {
    double retr = p_v * cos(p_a);
    double reti = p_v * sin(p_a);
    gridpack::ComplexType ret(retr, reti);
    values[0] = ret;
    return true;
  } else if (p_mode == State) {
    values[0] = p_v;
    values[1] = p_a;
    return true;
  } else if (p_mode == RHS) {
    double rvals[2];
    int nvals = rhsValues(rvals);
    for (int i=0; i<nvals; i++) values[i] = rvals[i];
    if (nvals == 0) {
      return false;
    } else {
      return true;
    }
  }
  return false;
}

bool gridpack::powerflow::PFBus::vectorValues(RealType *values)
{
  if (p_mode == State) {
    values[0] = p_v;
    values[1] = p_a;
    return true;
  }
  if (p_mode == RHS) {
    int nvals = rhsValues(values);
    if (nvals == 0) {
      return false;
    } else {
      return true;
    }
  }
  return false;
}

/**
 * Check QLIM
 * @return false: violations exist
 * @return true:  no violations
 * 
 */
bool gridpack::powerflow::PFBus::chkQlim(void)
{
  if (p_isPV) {
    double qmax, qmin, ppl;
    qmax = 0.0;
    qmin = 0.0;
    ppl = 0.0;
    for (int i=0; i<p_gstatus.size(); i++) {
      if (p_gstatus[i] == 1) {
        qmax += p_qmax[i];
        qmin += p_qmin[i];
        ppl  += p_pg[i];
      }
    }
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    int size = branches.size();
    double P, Q, p, q;
    int ngen=p_pFac.size();
    double pl =0.0;
    double ql =0.0;
    for (int i=0; i<p_lstatus.size(); i++) {
      if (p_lstatus[i] == 1) {
        pl += p_pl[i];
        ql += p_ql[i];
      }
    }
    p_save2isPV = p_isPV;
    double pval = p_Pinj*p_sbase+pl;
    double qval = p_Qinj*p_sbase+ql;

//  If qval exceeds the total generator Q capacity, perform PV->PQ
//
    if (qval > qmax ) {
      printf("\nWarning: Gen(s) at bus %d exceeds the QMAX %8.3f vs %8.3f, converted to PQ bus\n", getOriginalIndex(),qval, qmax);  
      ql = ql-qmax;
      p_save2isPV = p_isPV;
      p_isPV = false;
      *p_PV_ptr = false;
      pl -= ppl;
    //p_gstatus.clear();
      for (int i=0; i<p_gstatus.size(); i++) {
        p_gstatus_save.push_back(p_gstatus[i]);
        p_gstatus[i] = 0;
        p_qg[i] = p_qmax[i];
      }
      for (int i=0; i<p_lstatus.size(); i++) {
        if (p_lstatus[i] == 1) {
          p_pl[i] = pl;
          p_ql[i] = ql;
        }
      }
      if (p_PV_ptr) *p_PV_ptr = p_isPV;
      return true;
    } else if (qval < qmin) {
      printf("\nWarning: Gen(s) at bus %d exceeds the QMIN %8.3f vs %8.3f, converted to PQ bus\n", getOriginalIndex(),qval, qmin);  
      ql = ql-qmin;
      p_save2isPV = p_isPV;
      p_isPV = false;
      pl -= ppl;
    //  p_gstatus.clear();
      for (int i=0; i<p_gstatus.size(); i++) {
        p_gstatus_save.push_back(p_gstatus[i]);
        p_gstatus[i] = 0;
        p_qg[i] = p_qmin[i];
      }
      for (int i=0; i<p_lstatus.size(); i++) {
        if (p_lstatus[i] == 1) {
          p_pl[i] = pl;
          p_ql[i] = ql;
        }
      }

      if (p_PV_ptr) *p_PV_ptr = p_isPV;
      return true;
    } else {
       if (p_PV_ptr) *p_PV_ptr = p_isPV;
       return false;
    }
  } else {
    if (p_PV_ptr) *p_PV_ptr = p_isPV;
    return false;
  } 
  if (p_PV_ptr) *p_PV_ptr = p_isPV;
  return false;
}

/**
 * Clear changes that were made for Q limit violations and reset
 * bus to its original state
 */
void gridpack::powerflow::PFBus::clearQlim()
{
  int size = p_gstatus_save.size();
  p_gstatus.clear();
  int i;
  for (i=0; i<size; i++) {
    p_gstatus.push_back(p_gstatus_save[i]);
  }
  p_isPV = p_save2isPV;
  if (p_PV_ptr) *p_PV_ptr = p_isPV;
}


/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto buses 
 * @param values array containing voltage magnitude and angle
 */
void gridpack::powerflow::PFBus::setValues(gridpack::ComplexType *values)
{
  double vt = p_v;
  double at = p_a;
  p_a -= real(values[0]);
#ifdef LARGE_MATRIX
  p_v -= real(values[1]);
#else
  if (!p_isPV) {
    p_v -= real(values[1]);
  }
#endif
  *p_vMag_ptr = p_v;
  double pi = 4.0*atan(1.0);
  if (p_a >= 0.0) {
    *p_vAng_ptr = fmod(p_a+pi,2.0*pi)-pi;
  } else {
    *p_vAng_ptr = fmod(p_a-pi,2.0*pi)+pi;
  }
}

void gridpack::powerflow::PFBus::setValues(gridpack::RealType *values)
{
  double vt = p_v;
  double at = p_a;
  p_a -= values[0];
#ifdef LARGE_MATRIX
  p_v -= real(values[1]);
#else
  if (!p_isPV) {
    p_v -= values[1];
  }
#endif
  *p_vMag_ptr = p_v;
  double pi = 4.0*atan(1.0);
  if (p_a >= 0.0) {
    *p_vAng_ptr = fmod(p_a+pi,2.0*pi)-pi;
  } else {
    *p_vAng_ptr = fmod(p_a-pi,2.0*pi)+pi;
  }
}

/**
 * Return the size of the buffer used in data exchanges on the network.
 * For this problem, the voltage magnitude and phase angle need to be exchanged
 * @return size of buffer
 */
int gridpack::powerflow::PFBus::getXCBufSize(void)
{
  return (2*sizeof(double)+sizeof(bool));
}

/**
 * Assign pointers for voltage magnitude and phase angle
 */
void gridpack::powerflow::PFBus::setXCBuf(void *buf)
{
  p_vAng_ptr = static_cast<double*>(buf);
  p_vMag_ptr = p_vAng_ptr+1;
  void *ptr = static_cast<void*>(p_vMag_ptr+1);
  p_PV_ptr = static_cast<bool*>(ptr);
  // Note: we are assuming that the load function has been called BEFORE
  // the factory setExchange method, so p_a and p_v are set with their initial
  // values.
  *p_vMag_ptr = p_v;
  double pi = 4.0*atan(1.0);
  if (p_a >= 0.0) {
    *p_vAng_ptr = fmod(p_a+pi,2.0*pi)-pi;
  } else {
    *p_vAng_ptr = fmod(p_a-pi,2.0*pi)+pi;
  }
  *p_PV_ptr = p_isPV;
  
}

/**
 * Load values stored in DataCollection object into PFBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::powerflow::PFBus::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  p_data = data.get();
  YMBus::load(data);

  // This routine may be called more than once, so clear all vectors
  p_pg.clear();
  p_qg.clear();
  p_pFac.clear();
  p_pFac_orig.clear();
  p_gstatus.clear();
  p_gstatus_save.clear();
  p_qmin.clear();
  p_qmax.clear();
  p_qmin_orig.clear();
  p_qmax_orig.clear();
  p_gid.clear();
  p_pt.clear();
  p_pb.clear();
  p_pl.clear();
  p_ql.clear();
  p_lstatus.clear();
  p_lid.clear();

  bool ok = data->getValue(CASE_SBASE, &p_sbase);
  data->getValue(BUS_VOLTAGE_ANG, &p_angle);
  data->getValue(BUS_VOLTAGE_MAG, &p_voltage); 
  p_v = p_voltage;
  double pi = 4.0*atan(1.0);
  p_angle = p_angle*pi/180.0;
  p_a = p_angle;
  data->getValue(BUS_TYPE, &p_type);
  if (p_type == 3) {
    setReferenceBus(true);
  }
  if (isIsolated()) {
    p_original_isolated = true;
  } else {
    p_original_isolated = false;
  }
  p_area = 1;
  data->getValue(BUS_AREA, &p_area);
  p_zone = 1;
  data->getValue(BUS_ZONE, &p_zone);

  // if BUS_TYPE = 2, and gstatus is 1, then bus is a PV bus
  p_isPV = false;

  // added p_pg,p_qg,p_pl,p_ql,p_sbase;

  bool lgen;
  int i, gstatus;
  double pg, qg, vs,qmax,qmin;
  int ngen = 0;
  p_ngen = 0;
  if (data->getValue(GENERATOR_NUMBER, &ngen)) {
    double qtot = 0.0;
    for (i=0; i<ngen; i++) {
      lgen = true;
      lgen = lgen && data->getValue(GENERATOR_PG, &pg,i);
      lgen = lgen && data->getValue(GENERATOR_QG, &qg,i);
      lgen = lgen && data->getValue(GENERATOR_VS, &vs,i);
      lgen = lgen && data->getValue(GENERATOR_STAT, &gstatus,i);
      lgen = lgen && data->getValue(GENERATOR_QMAX, &qmax,i);
      lgen = lgen && data->getValue(GENERATOR_QMIN, &qmin,i);
      double pt = 0.0;
      double pb = 0.0;
      ok =  data->getValue(GENERATOR_PMAX,&pt,i);
      ok =  data->getValue(GENERATOR_PMIN,&pb,i);
      if (lgen) {
        p_pg.push_back(pg);
        p_savePg.push_back(pg);
        p_qg.push_back(qg);
        p_gstatus.push_back(gstatus);
        if (gstatus == 0) {
          qmax = 0.0;
          qmin = 0.0;
        }
        p_pFac.push_back(qmax);
        p_qmax.push_back(qmax);
        qtot += qmax;
        p_qmin.push_back(qmin);
        p_pFac_orig.push_back(qmax);
        p_qmax_orig.push_back(qmax);
        p_qmin_orig.push_back(qmin);
        p_pt.push_back(pt);
        p_pb.push_back(pb);
        if (gstatus == 1) {
          p_v = vs; //reset initial PV voltage to set voltage
          if (p_type == 2) p_isPV = true;
        }
        std::string id("-1");
        data->getValue(GENERATOR_ID,&id,i);
        p_gid.push_back(id);
        p_ngen++;
      }
    }
    if (qtot != 0.0 && p_ngen > 1) {
      for (i=0; i<p_ngen; i++) {
        p_pFac[i] = p_pFac[i]/qtot;
        p_pFac_orig[i] = p_pFac[i];
      }
    } else {
      p_pFac[0] = 1.0;
      p_pFac_orig[0] = p_pFac[0];
    }
  }
  p_saveisPV = p_isPV;

// Add load
  int lstatus;
  double pl,ql,ip,iq,yp,yq;
  p_load = true;
  p_load = p_load && data->getValue(LOAD_PL, &pl,0);
  p_load = p_load && data->getValue(LOAD_QL, &ql,0);
  p_load = p_load && data->getValue(LOAD_IP, &ip,0);
  p_load = p_load && data->getValue(LOAD_IQ, &iq,0);
  p_load = p_load && data->getValue(LOAD_YP, &yp,0);
  p_load = p_load && data->getValue(LOAD_YQ, &yq,0);
  int nld = 0;
  p_nload = 0;
  if (data->getValue(LOAD_NUMBER, &nld)) {
    for (i=0; i<nld; i++) {
      p_load = true;
      p_load = p_load && data->getValue(LOAD_PL, &pl,i);
      p_load = p_load && data->getValue(LOAD_QL, &ql,i);
      p_load = p_load && data->getValue(LOAD_STATUS, &lstatus,i);
      if (p_load) {
        p_pl.push_back(pl);
        p_savePl.push_back(pl);
        p_ql.push_back(ql);
        p_saveQl.push_back(ql);
        p_lstatus.push_back(lstatus);
        std::string id("-1");
        data->getValue(LOAD_ID,&id,i);
        p_lid.push_back(id);
        p_nload++;
      }
    }
  }
  // If this is being called a second time, then update pointers
  if (p_vMag_ptr) *p_vMag_ptr = p_v;
  if (p_vAng_ptr) {
    double pi = 4.0*atan(1.0);
    if (p_a >= 0.0) {
      *p_vAng_ptr = fmod(p_a+pi,2.0*pi)-pi;
    } else {
      *p_vAng_ptr = fmod(p_a-pi,2.0*pi)+pi;
    }
  }
}

/**
 * Set values of YBus matrix. These can then be used in subsequent
 * calculations
 */
void gridpack::powerflow::PFBus::setYBus(void)
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
gridpack::ComplexType gridpack::powerflow::PFBus::getYBus(void)
{
  return YMBus::getYBus();
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::powerflow::PFBus::setMode(int mode)
{
  if (mode == YBus) {
    YMBus::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Reset voltage and phase angle to initial values
 */
void gridpack::powerflow::PFBus::resetVoltage(void)
{
  p_v = p_voltage;
  p_a = p_angle;
  if (p_vMag_ptr) *p_vMag_ptr = p_v;
  if (p_vAng_ptr) {
    double pi = 4.0*atan(1.0);
    if (p_a >= 0.0) {
      *p_vAng_ptr = fmod(p_a+pi,2.0*pi)-pi;
    } else {
      *p_vAng_ptr = fmod(p_a-pi,2.0*pi)+pi;
    }
  }
}

/**
 * Set voltage limits on bus
 * @param vmin lower value of voltage
 * @param vmax upper value of voltage
 */
void gridpack::powerflow::PFBus::setVoltageLimits(double vmin, double vmax)
{
  p_vmin = vmin;
  p_vmax = vmax;
}

/**
 * Check voltage for violation
 * @return false if there is a voltage violation
 */
bool gridpack::powerflow::PFBus::checkVoltageViolation()
{
  bool ret = true;
  if (*p_vMag_ptr > p_vmax || *p_vMag_ptr < p_vmin) ret = false;
  return ret;
}

/**
 * Return the value of the voltage magnitude on this bus
 * @return: voltage magnitude
 */
double gridpack::powerflow::PFBus::getVoltage()
{
  return *p_vMag_ptr;
}

/**
 * Return whether or not the bus is a PV bus (V held fixed in powerflow
 * equations)
 * @return true if bus is PV bus
 */
bool gridpack::powerflow::PFBus::isPV(void)
{
  return p_isPV;
}

/**
 * Return the value of the phase angle on this bus
 * @return: phase angle
 */
double gridpack::powerflow::PFBus::getPhase()
{
  return *p_vAng_ptr;
}

/**
 * Get generator status
 * @param gen_id generator ID
 * @return vector of generator statuses
 */
bool gridpack::powerflow::PFBus::getGenStatus(std::string gen_id)
{
  int i;
  int gsize = p_gstatus.size();
  for (i=0; i<gsize; i++) {
    if (gen_id == p_gid[i]) {
      return p_gstatus[i];
    }
  }
  return false;
}

/**
 * Get list of generator IDs
 * @return vector of generator IDs
 */
std::vector<std::string> gridpack::powerflow::PFBus::getGenerators()
{
  return p_gid;
}

/**
 * Get list of load IDs
 * @return vector of generator IDs
 */
std::vector<std::string> gridpack::powerflow::PFBus::getLoads()
{
  return p_lid;
}

/**
 * Set generator status
 * @param gen_id generator ID
 * @param status generator status
 */
void gridpack::powerflow::PFBus::setGenStatus(std::string gen_id, bool status)
{
  int i;
  int gsize = p_gstatus.size();
  for (i=0; i<gsize; i++) {
    if (gen_id == p_gid[i]) {
      p_gstatus[i] = status;
      if (status == 0) {
        p_pFac[i] = 0.0;
        p_qmax[i] = 0.0;
        p_qmin[i] = 0.0;
      } else {
        p_pFac[i] = p_pFac_orig[i];
        p_qmax[i] = p_qmax_orig[i];
        p_qmin[i] = p_qmin_orig[i];
      }
      return;
    }
  }
}

/**
 * Set isPV status
 * @param status isPV status
 */
void gridpack::powerflow::PFBus::setIsPV(int status)
{
  p_saveisPV = p_isPV;
  p_isPV = status;
  if (p_PV_ptr) *p_PV_ptr = status;
  p_v = p_voltage;
}

/**
 * Reset isPV status
 */
void gridpack::powerflow::PFBus::resetIsPV()
{
  p_isPV = p_saveisPV;
  if (p_PV_ptr) *p_PV_ptr = p_saveisPV;
}

/**
 * setSBus
 * BUS = (CG*(GEN(ON,PG) + J*GEN(ON,QG)-(PD+J*QD))/BASEMVA
 */
void gridpack::powerflow::PFBus::setSBus(void)
{
  int i;
  double pg, qg, pl, ql;
  pg = 0.0;
  qg = 0.0;
  pl = 0.0;
  ql = 0.0;
  bool usegen = false;
  for (i=0; i<p_gstatus.size(); i++) {
    if (p_gstatus[i] == 1) {
      pg += p_pg[i];
      qg += p_qg[i];
      usegen = true;
    }
  }
  for (i=0; i<p_lstatus.size(); i++) {
    if (p_lstatus[i] == 1) {
      pl += p_pl[i];
      ql += p_ql[i];
    }
  }
  if (p_gstatus.size() > 0 && usegen) {
    gridpack::ComplexType sBus((pg - pl) / p_sbase, (qg - ql) / p_sbase);
    p_P0 = real(sBus);
    p_Q0 = imag(sBus);
  } else {
    gridpack::ComplexType sBus((- pl) / p_sbase, (- ql) / p_sbase);
    p_P0 = real(sBus);
    p_Q0 = imag(sBus);
  } 
}

/**
 ** Update pg of specified bus element based on their genID
 ** @param busID
 ** @param genID
 ** @param value
 **/
/*
void gridpack::powerflow::PFBus::updatePg(int busID, std::string genID, double value)
{
  if (getOriginalIndex() == busID) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        if (p_gid[i] == genID) {
          p_pg[i] += value;
        }
      }
    }
  }
}
*/

/**
 * Write output from buses to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool gridpack::powerflow::PFBus::serialWrite(char *string, const int bufsize,
                                             const char *signal)
{
  if (signal == NULL) {
    double pi = 4.0*atan(1.0);
    double angle = p_a*180.0/pi;
    if (!isIsolated()) {
      sprintf(string, "     %6d      %12.6f         %12.6f\n",
            getOriginalIndex(),angle,p_v);
    } else {
      return false;
      /*
      sprintf(string, "     %6d      %12.6f         %12.6f\n",
      getOriginalIndex(),0.0,0.0);
      */
    }
  } else if (!strcmp(signal,"vr_str")) {
    double pi = 4.0*atan(1.0);
    double angle = p_a*180.0/pi;
    int use_vmag = 1;
    if (p_saveisPV || p_original_isolated) use_vmag = 0;
    int changed = 0;
    if (p_isPV != p_saveisPV) changed = 1;
    sprintf(string, "%6d %20.12e %20.12e %d %d\n",
        getOriginalIndex(),angle,p_v,use_vmag,changed);
  } else if (!strcmp(signal,"vfail_str")) {
    int use_vmag = 1;
    if (p_saveisPV || p_original_isolated) use_vmag = 0;
    int changed = 0;
    sprintf(string, "%6d %20.12e %20.12e %d %d\n",
        getOriginalIndex(),0.0,0.0,use_vmag,changed);
  } else if (!strcmp(signal,"ca")) {
    double pi = 4.0*atan(1.0);
    double angle = p_a*180.0/pi;
    bool found = false;
    if ((p_v > p_vmax) || (p_v < p_vmin) ) {
      found = true;
      if (!isIsolated()) {
        sprintf(string, "     %6d      %12.6f         %12.6f\n",
            getOriginalIndex(),angle,p_v);
      } else {
        sprintf(string, "     %6d      %12.6f         %12.6f\n",
            getOriginalIndex(),0.0,0.0);
      }
    }
    return found;
  } else if (!strcmp(signal,"pq")) {
    gridpack::ComplexType v[2];
    vectorValues(v);
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    if (!isIsolated()) {
      sprintf(string, "     %6d      %12.6f      %12.6f      %2d\n",
            getOriginalIndex(),real(v[0]),real(v[1]),
            static_cast<int>(branches.size()));
    } else {
      sprintf(string, "     %6d      %12.6f      %12.6f      %2d\n",
            getOriginalIndex(),0.0, 0.0,
            static_cast<int>(branches.size()));
    }
  } else if (!strcmp(signal,"record")) {
    char sbuf[128];
    char *cptr = string;
    int slen = 0;
    int nld = p_nload;
    int i;
    double pl =0.0;
    double ql =0.0;
    for (i=0; i<nld; i++) {
      if (p_lstatus[i]) pl += p_pl[i];
      if (p_lstatus[i]) ql += p_ql[i];
    }
    sprintf(sbuf,"%8d, %4d, %16.8f, %16.8f,",getOriginalIndex(),p_type,
           pl/p_sbase,ql/p_sbase);
    int len = strlen(sbuf);
    if (len<=bufsize) {
      sprintf(cptr,"%s",sbuf);
      slen += len;
      cptr += len;
    }
    double pgen = 0.0;
    double qgen = 0.0;
    double qmin = 0.0;
    double qmax = 0.0;
    int ngen = p_ngen;
    for (i=0; i<ngen; i++) {
      if (p_gstatus[i]) pgen += p_pg[i];
      if (p_gstatus[i]) qgen += p_qg[i];
      if (p_gstatus[i]) qmin += p_qmin[i];
      if (p_gstatus[i]) qmax += p_qmax[i];
    }
    sprintf(sbuf," %16.8f, %16.8f, %16.8f, %16.8f,",
            pgen/p_sbase,qgen/p_sbase,qmax/p_sbase,qmin/p_sbase);
    len = strlen(sbuf);
    if (slen+len<=bufsize) {
      sprintf(cptr,"%s",sbuf);
      slen += len;
      cptr += len;
    }
    int gstatus = 0;
    double pt = 0.0;
    double pb = 0.0;
    for (i=0; i<ngen; i++) {
      if (p_gstatus[i]) gstatus = 1;
      if (p_gstatus[i]) pt += p_pt[i];
      if (p_gstatus[i]) pb += p_pb[i];
    }
    sprintf(sbuf," %16.8f, %16.8f, %1d,",pt,pb,gstatus);
    len = strlen(sbuf);
    if (slen+len<=bufsize) {
      sprintf(cptr,"%s",sbuf);
      slen += len;
      cptr += len;
    }
    double gl, bl;
    YMBus::getShuntValues(&bl, &gl);
    int area;
    p_data->getValue(BUS_AREA,&area);
    sprintf(sbuf," %16.8f, %16.8f, %8d,",gl,bl,area);
    len = strlen(sbuf);
    if (slen+len<=bufsize) {
      sprintf(cptr,"%s",sbuf);
      slen += len;
      cptr += len;
    }
    double zero = 0.0;
    int nzone;
    double basekv;
    p_data->getValue(BUS_ZONE,&nzone);
    p_data->getValue(BUS_BASEKV,&basekv);
    double pi = 4.0*atan(1.0);
    double angle = p_a*180.0/pi;
    if (!isIsolated()) {
      sprintf(sbuf," %16.8f, %16.8f, %16.8f, %8d, %4.2f, %4.2f\n",p_v,angle,basekv,nzone,p_vmax,p_vmin);
    } else {
      sprintf(sbuf," %16.8f, %16.8f, %16.8f, %8d, %4.2f, %4.2f\n",0.0,0.0,basekv,nzone,p_vmax,p_vmin);
    }
    len = strlen(sbuf);
    if (slen+len<=bufsize) {
      sprintf(cptr,"%s",sbuf);
      slen += len;
      cptr += len;
    }
  } else if (!strcmp(signal,"power") || !strcmp(signal,"gen_str")) {
    char sbuf[128];
    char *cptr = string;
    int i, len, slen = 0;
    int ngen=p_pFac.size();
    // Evalate p_Pinj and p_Qinj if bus is reference bus. This is skipped when
    // evaluating matrix elements.
#ifndef LARGE_MATRIX
    if (getReferenceBus() || isIsolated()) {
      std::vector<boost::shared_ptr<BaseComponent> > branches;
      getNeighborBranches(branches);
      int size = branches.size();
      double P, Q, p, q;
      P = 0.0;
      Q = 0.0;
      for (i=0; i<size; i++) {
        gridpack::powerflow::PFBranch *branch
          = dynamic_cast<gridpack::powerflow::PFBranch*>(branches[i].get());
        branch->getPQ(this, &p, &q);
        P += p;
        Q += q;
      }
      // Also add bus i's own Pi, Qi
      P += p_v*p_v*p_ybusr;
      Q += p_v*p_v*(-p_ybusi);
      p_Pinj = P;
      p_Qinj = Q;
    }
#endif
    double pl =0.0;
    double ql =0.0;
    for (i=0; i<p_pl.size(); i++) {
      if (p_lstatus[i] == 1) {
        pl += p_pl[i];
        ql += p_ql[i];
      }
    }
    for (i=0; i<ngen; i++) {
      double pval = p_pFac[i]*(p_Pinj+pl/p_sbase);
      double qval = p_pFac[i]*(p_Qinj+ql/p_sbase);
      if (!strcmp(signal,"power")) {
        sprintf(sbuf, "     %6d      %s   %12.6f      %12.6f\n",
            getOriginalIndex(),p_gid[i].c_str(),pval,qval);
      } else {
        sprintf(sbuf, "%6d %s %20.12e %20.12e\n",
            getOriginalIndex(),p_gid[i].c_str(),pval,qval);
      }
      len = strlen(sbuf);
      if (slen+len<=bufsize) {
        sprintf(cptr,"%s",sbuf);
        slen += len;
        cptr += len;
      }
    }
    if (slen>0) {
      return true;
    } else {
      return false;
    }
  } else if (!strcmp(signal,"pfail_str")) {
    char sbuf[128];
    char *cptr = string;
    int i, len, slen = 0;
    int ngen=p_pFac.size();
    for (i=0; i<ngen; i++) {
      sprintf(sbuf, "%6d %s %20.12e %20.12e\n",
            getOriginalIndex(),p_gid[i].c_str(),0.0,0.0);
      len = strlen(sbuf);
      if (slen+len<=bufsize) {
        sprintf(cptr,"%s",sbuf);
        slen += len;
        cptr += len;
      }
    }
    if (slen>0) {
      return true;
    } else {
      return false;
    }
  } else if (!strcmp(signal,"src_gen")) {
    if (p_source) {
      char sbuf[128];
      char *cptr = string;
      int i, len, slen = 0;
      std::string status;
      for (i=0; i<p_ngen; i++) {
        if (p_gstatus[i]) {
          status = "  active";
        } else {
          status = "inactive";
        }
        sprintf(sbuf,"%8d %s %s %4d %4d %14.4f %14.4f %14.4f %14.4f\n",
              getOriginalIndex(),
              p_gid[i].c_str(),status.c_str(),p_area,p_zone,p_savePg[i],
              p_pg[i],p_pb[i],p_pt[i]);
        len = strlen(sbuf);
        if (slen+len <= bufsize) {
          sprintf(cptr,"%s",sbuf);
          slen += len;
          cptr += len;
        }
      }
      if (slen>0) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  } else if (!strcmp(signal,"sink_load")) {
    if (p_sink) {
      char sbuf[128];
      char *cptr = string;
      int i, len, slen = 0;
      std::string status;
      for (i=0; i<p_nload; i++) {
        if (p_lstatus[i]) {
          status = "  active";
        } else {
          status = "inactive";
        }
        sprintf(sbuf,"%8d %s %s %4d %4d %14.4f %14.4f %14.4f %14.4f\n",
              getOriginalIndex(),
              p_lid[i].c_str(),status.c_str(),p_area,p_zone,p_savePl[i],
              p_pl[i],p_saveQl[i],p_ql[i]);
        len = strlen(sbuf);
        if (slen+len <= bufsize) {
          sprintf(cptr,"%s",sbuf);
          slen += len;
          cptr += len;
        }
      }
      if (slen>0) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }
  return true;
}

/**
 * Return the complex voltage on this bus
 * @return the complex voltage
 */
gridpack::ComplexType gridpack::powerflow::PFBus::getComplexVoltage(void)
{
  p_a = *p_vAng_ptr;
  p_v =  *p_vMag_ptr;
  gridpack::ComplexType ret(cos(p_a),sin(p_a));
  ret = ret*p_v;
  return ret;
}

/**
 * Save state variables inside the component to a DataCollection object.
 * This can be used as a way of moving data in a way that is useful for
 * creating output or for copying state data from one network to another.
 * @param data data collection object into which new values are inserted
 */
void gridpack::powerflow::PFBus::saveData(
    boost::shared_ptr<gridpack::component::DataCollection> data)
{
  double rval;
  int i;
  if (!data->setValue("BUS_PF_VMAG",*p_vMag_ptr)) {
    data->addValue("BUS_PF_VMAG",*p_vMag_ptr);
  }
  rval = *p_vAng_ptr;
  double pi = 4.0*atan(1.0);
  rval = 180.0*rval/pi;
  if (!data->setValue("BUS_PF_VANG",rval)) {
    data->addValue("BUS_PF_VANG",rval);
  }
  if (!data->setValue("BUS_TYPE",p_type)) {
    data->addValue("BUS_TYPE",p_type);
  }
  int ngen=p_pFac.size();
  // Evalate p_Pinj and p_Qinj if bus is reference bus. This is skipped when
  // evaluating matrix elements.
#ifndef LARGE_MATRIX
  if (getReferenceBus() || isIsolated()) {
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    int size = branches.size();
    double P, Q, p, q;
    P = 0.0;
    Q = 0.0;
    for (i=0; i<size; i++) {
      gridpack::powerflow::PFBranch *branch
        = dynamic_cast<gridpack::powerflow::PFBranch*>(branches[i].get());
      branch->getPQ(this, &p, &q);
      P += p;
      Q += q;
    }
    // Also add bus i's own Pi, Qi
    P += p_v*p_v*p_ybusr;
    Q += p_v*p_v*(-p_ybusi);
    p_Pinj = P;
    p_Qinj = Q;
  }
#endif
  double pl=0.0;
  double ql=0.0;
  for (i=0; i<p_pl.size(); i++) {
    if (p_lstatus[i] == 1) {
      pl += p_pl[i];
      ql += p_ql[i];
    }
  }
  for (i=0; i<ngen; i++) {
    rval = p_pFac[i]*(p_Pinj+pl/p_sbase);
    if (!data->setValue("GENERATOR_PF_PGEN",rval,i)) {
      data->addValue("GENERATOR_PF_PGEN",rval,i);
    }
    rval = p_pFac[i]*(p_Qinj+ql/p_sbase);
    if (!data->setValue("GENERATOR_PF_QGEN",rval,i)) {
      data->addValue("GENERATOR_PF_QGEN",rval,i);
    }
  }
}

/**
 * Modify parameters inside the bus module. This is designed to be
 * extensible
 * @param name character string describing parameter to be modified
 * @param busID generator bus number
 * @param genID specified genID
 * @param value new value of parameter
 */
void gridpack::powerflow::PFBus::setParam(std::string name, int busID, 
    std::string genID, double value) 
{
  if (getOriginalIndex() == busID) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        if (p_gid[i] == genID) {
          if (name == GENERATOR_PG) {
            p_pg[i] = value;
          } else if (name == GENERATOR_QG) {
            p_qg[i] = value;
          }
        }
      }
    }
  }
}

/**
 * Modify parameters inside the bus module. This is designed to be
 * extensible
 * @param name character string describing parameter to be modified
 * @param busID generator bus number
 * @param genID specified genID
 * @param value new value of parameter
 */
void gridpack::powerflow::PFBus::setParam(int busID, 
    std::string genID, double value) 
{
//  if (name == GENERATOR_PG) {
   if (getOriginalIndex() == busID) {
     if (p_ngen > 0) {
       for (int i = 0; i < p_ngen; i++) {
         if (p_gid[i] == genID) {
           p_pg[i] += value;
         }
       }
     }
   }
   // if (idx >= 0 && idx<p_pg.size()) {
   //   p_pg[idx] = value;
   // }
//  } else if (name == GENERATOR_QG) {
   // if (idx >= 0 && idx<p_qg.size()) {
   //   p_qg[idx] = value;
   // }
//  }
}

/**
 * Access parameters inside the bus module. This is designed to be
 * extensible
 * @param name character string describing parameter to be accessed
 * @param value value of parameter
 * @param idx index (if necessary) of variable to be accessed
 */
void gridpack::powerflow::PFBus::getParam(std::string &name,
    double *value, int idx)
{
  if (name == GENERATOR_PG) {
    if (idx >= 0 && idx<p_pg.size()) {
      *value = p_pg[idx];
    }
  } else if (name == GENERATOR_QG) {
    if (idx >= 0 && idx<p_qg.size()) {
      *value = p_qg[idx];
    }
  }
}

void gridpack::powerflow::PFBus::getParam(std::string &name,
    int *value, int idx)
{
  if (name == GENERATOR_NUMBER) {
    *value = p_pg.size();
  }
}

/**
 * Get index of internal bus element based on character string identifier
 * @param name character string describing element
 * @param tag character string specifying bus element
 * @return index of element
 */
int gridpack::powerflow::PFBus::getElementIndex(std::string &name, std::string &tag)
{
  if (name == "GENERATOR") {
    int i;
    int nsize = static_cast<int>(p_gid.size());
    for (i=0; i<nsize; i++) {
      if (tag == p_gid[i]) {
        return i;
      }
    }
  }
  return -1;
}

/**
 * Set parameter to ignore voltage violations
 * @param flag value of ignore parameter
 */
void gridpack::powerflow::PFBus::setIgnore(bool flag)
{
  p_ignore = flag;
}

/**
 * Get parameter to ignore voltage violations
 * @return value of ignore parameter
 */
bool gridpack::powerflow::PFBus::getIgnore()
{
  return p_ignore;
}

/**
 * Get area parameter for bus
 * @return bus area index
 */
int gridpack::powerflow::PFBus::getArea()
{
  return p_area;
}

/**
 * Get zone parameter for bus
 * @return bus zone index
 */
int gridpack::powerflow::PFBus::getZone()
{
  return p_zone;
}

/**
 * Evaluate diagonal block of Jacobian for power flow calculation and
 * return result as an array of real values
 * @param rvals values of Jacobian block
 * @return number of values returned
 */
int gridpack::powerflow::PFBus::diagonalJacobianValues(double *rvals)
{
  if (!isIsolated()) {
#ifdef LARGE_MATRIX
    if (!getReferenceBus()) {
      rvals[0] = -p_Qinj - p_ybusi * p_v *p_v; 
      rvals[1] = p_Pinj - p_ybusr * p_v *p_v; 
      rvals[2] = p_Pinj / p_v + p_ybusr * p_v; 
      rvals[3] = p_Qinj / p_v - p_ybusi * p_v; 
      // Fix up matrix elements if bus is PV bus
      if (p_isPV) {
        rvals[1] = 0.0;
        rvals[2] = 0.0;
        rvals[3] = 1.0;
      }
      return 4;
    } else {
      rvals[0] = 1.0;
      rvals[1] = 0.0;
      rvals[2] = 0.0;
      rvals[3] = 1.0;
      return 4;
    }
#else
    if (!getReferenceBus() && !p_isPV) {
      rvals[0] = -p_Qinj - p_ybusi * p_v *p_v; 
      rvals[1] = p_Pinj - p_ybusr * p_v *p_v; 
      rvals[2] = p_Pinj / p_v + p_ybusr * p_v; 
      rvals[3] = p_Qinj / p_v - p_ybusi * p_v; 
      // Fix up matrix elements if bus is PV bus
      return 4;
    } else if (!getReferenceBus() && p_isPV) {
      rvals[0] = -p_Qinj - p_ybusi * p_v *p_v; 
      return 1;
    } else {
      return 0;
    }
#endif
  } else {
    return 0;
  }
}

/**
 * Push p_isPV values from exchange buffer to p_isPV variable
 */
void gridpack::powerflow::PFBus::pushIsPV()
{
  if (p_PV_ptr) p_isPV = *p_PV_ptr;
}


/**
 * Evaluate RHS values for powerflow equation and return result as
 * an array of real values
 * @param rvals values of Jacobian block
 * @return number of values returned
 */
int gridpack::powerflow::PFBus::rhsValues(double *rvals)
{
  if (!isIsolated()) {
    if (!getReferenceBus()) {
      int nvals;
      std::vector<boost::shared_ptr<BaseComponent> > branches;
      getNeighborBranches(branches);
      int size = branches.size();
      int i;
      double P, Q, p, q;
      P = 0.0;
      Q = 0.0;
      for (i=0; i<size; i++) {
        gridpack::powerflow::PFBranch *branch
          = dynamic_cast<gridpack::powerflow::PFBranch*>(branches[i].get());
        branch->getPQ(this, &p, &q);
        P += p;
        Q += q;
      }
      // Also add bus i's own Pi, Qi
      P += p_v*p_v*p_ybusr;
      Q += p_v*p_v*(-p_ybusi);
      p_Pinj = P;
      p_Qinj = Q;
      P -= p_P0;
      Q -= p_Q0;
      rvals[0] = P;
#ifdef LARGE_MATRIX
      if (!p_isPV) {
        rvals[1] = Q;
      } else {
        rvals[1] = 0.0;
      }
      nvals = 2;
#else
      nvals = 1;
      if (!p_isPV) {
        rvals[1] = Q;
        nvals = 2;
      }
#endif
      return nvals;
    } else {
#ifdef LARGE_MATRIX
      std::vector<boost::shared_ptr<BaseComponent> > branches;
      getNeighborBranches(branches);
      int size = branches.size();
      int i;
      double P, Q, p, q;
      P = 0.0;
      Q = 0.0;
      for (i=0; i<size; i++) {
        gridpack::powerflow::PFBranch *branch
          = dynamic_cast<gridpack::powerflow::PFBranch*>(branches[i].get());
        branch->getPQ(this, &p, &q);
        P += p;
        Q += q;
      }
      // Also add bus i's own Pi, Qi
      P += p_v*p_v*p_ybusr;
      Q += p_v*p_v*(-p_ybusi);
      p_Pinj = P;
      p_Qinj = Q;
      rvals[0] = 0.0;
      rvals[1] = 0.0;
      return 2;
#else
      return 0;
#endif
    }
  } else {
    return false;
  }
}

/**
 * Get vector containing generator participation
 * @return vector of generator participation factors
 */
std::vector<double> gridpack::powerflow::PFBus::getGeneratorParticipation()
{
  return p_pFac;
}

/**
 * Set value of real power on individual generators
 * @param tag generator ID
 * @param value new value of real power
 * @param data data collection object associated with bus
 */
void gridpack::powerflow::PFBus::setGeneratorRealPower(
    std::string tag, double value, gridpack::component::DataCollection *data)
{
  int i, idx;
  idx = -1;
  for (i=0; i<p_ngen; i++) {
    if (p_gid[i] == tag) {
      idx = i;
      break;
    }
  }
  if (idx != -1) {
    if (!data->setValue(GENERATOR_PG,value,idx)) {
      data->addValue(GENERATOR_PG,value,idx);
    }
  } else {
    printf("No generator found for tag: (%s)\n",tag.c_str());
  }
}

/**
 * Scale value of real power on all generators
 * @param character ID for generator
 * @param value scale factor for real power
 */
void gridpack::powerflow::PFBus::scaleGeneratorRealPower(std::string tag,
    double value)
{
  int i;
  for (i=0; i<p_ngen; i++) {
    if (p_gid[i] == tag && p_gstatus[i]) {
      if (value > 0.0) {
        double excess = p_pt[i]-p_pg[i];
        if (excess < 0.0) {
          printf("bus: %d generator: %s excess (pt): %f (pg): %f\n",
              getOriginalIndex(),tag.c_str(),p_pt[i],p_pg[i]);
        }
        p_pg[i] += value*excess;
      } else {
        double slack = p_pg[i]-p_pb[i];
        p_pg[i] += value*slack;
      }
      break;
    }
  }
}

/**
 * Set value of real power on individual generators
 * @param tag generator ID
 * @param value new value of real power
 * @param data data collection object associated with bus
 */
void gridpack::powerflow::PFBus::setLoadRealPower(
    std::string tag, double value, gridpack::component::DataCollection *data)
{
  int i, idx;
  idx = -1;
  for (i=0; i<p_nload; i++) {
    if (p_lid[i] == tag) {
      idx = i;
      break;
    }
  }
  if (idx != -1) {
    if (!data->setValue(LOAD_PL,value,idx)) {
      data->addValue(LOAD_PL,value,idx);
    }
  } else {
    printf("No load found for tag: (%s)\n",tag.c_str());
  }
}

/**
 * Scale value of real and reactive power on loads
 * @param character ID for load
 * @param value scale factor for real power
 */
void gridpack::powerflow::PFBus::scaleLoadPower(std::string tag, double value)
{
  int i;
  for (i=0; i<p_nload; i++) {
    if (p_lid[i] == tag && p_lstatus[i] == 1) {
      p_pl[i] = value*p_pl[i];
      p_ql[i] = value*p_ql[i];
      break;
    }
  }
}

/**
 * Reset power for generators and loads back to original values
 */
void gridpack::powerflow::PFBus::resetPower()
{
  resetVoltage();
  int i;
  for (i=0; i<p_ngen; i++) {
    p_pg[i] = p_savePg[i];
  }
  for (i=0; i<p_nload; i++) {
    p_pl[i] = p_savePl[i];
    p_ql[i] = p_saveQl[i];
  }
}

/**
 * Get available margin for generator
 * @param tag character ID for generator
 * @param current initial generation
 * @param pmin minimum allowable generation
 * @param pmax maximum allowable generation
 * @param status current status of generator
 */
void gridpack::powerflow::PFBus::getGeneratorMargins(
    std::vector<std::string> &tag, std::vector<double> &current,
    std::vector<double> &pmin, std::vector<double> &pmax,
    std::vector<int> &status)
{
  tag.clear();
  current.clear();
  pmin.clear();
  pmax.clear();
  status.clear();
  int i;
  for (i=0; i<p_ngen; i++) {
    tag.push_back(p_gid[i]);
    current.push_back(p_savePg[i]);
    pmin.push_back(p_pb[i]);
    pmax.push_back(p_pt[i]);
    status.push_back(p_gstatus[i]);
  }
}

/**
 * Get current value of loads
 * @param tag character ID for load
 * @param pl initial value of load real power
 * @param ql initial value of load reactive power
 * @param status current status of load
 */
void gridpack::powerflow::PFBus::getLoadPower(
    std::vector<std::string> &tag, std::vector<double> &pl,
    std::vector<double> &ql, std::vector<int> &status)
{
  tag.clear();
  pl.clear();
  ql.clear();
  status.clear();
  int i;
  for (i=0; i<p_nload; i++) {
    tag.push_back(p_lid[i]);
    pl.push_back(p_savePl[i]);
    ql.push_back(p_saveQl[i]);
    status.push_back(p_lstatus[i]);
  }
}

/**
 * Label bus as a source for real time path rating
 * @param flag identify bus as source
 */
void gridpack::powerflow::PFBus::setSource(bool flag)
{
  p_source = flag;
}

/**
 * Label bus as a sink for real time path rating
 * @param flag identify bus as sink
 */
void gridpack::powerflow::PFBus::setSink(bool flag)
{
  p_sink = flag;
}

/**
 * Store scale factor
 * @param scale factor for scaling generation or loads
 */
void gridpack::powerflow::PFBus::setScale(double scale)
{
  p_rtpr_scale = scale;
}

/**
 *  Simple constructor
 */
gridpack::powerflow::PFBranch::PFBranch(void)
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
  p_elems = 0;
  p_theta = 0.0;
  p_sbase = 0.0;
  p_mode = YBus;
}

/**
 *  Simple destructor
 */
gridpack::powerflow::PFBranch::~PFBranch(void)
{
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::powerflow::PFBranch::matrixForwardSize(int *isize, int *jsize) const
{
  if (p_mode == Jacobian) {
    gridpack::powerflow::PFBus *bus1
      = dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
    gridpack::powerflow::PFBus *bus2
      = dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
    bool ok = !bus1->getReferenceBus();
    ok = ok && !bus2->getReferenceBus();
    ok = ok && !bus1->isIsolated();
    ok = ok && !bus2->isIsolated();
    ok = ok && (p_active);
    if (ok) {
#ifdef LARGE_MATRIX
      *isize = 2;
      *jsize = 2;
      return true;
#else
      bool bus1PV = bus1->isPV();
      bool bus2PV = bus2->isPV();
      if (bus1PV && bus2PV) {
        *isize = 1;
        *jsize = 1;
        return true;
      } else if (bus1PV) {
        *isize = 1;
        *jsize = 2;
        return true;
      } else if (bus2PV) {
        *isize = 2;
        *jsize = 1;
        return true;
      } else {
        *isize = 2;
        *jsize = 2;
        return true;
      }
#endif
    } else {
      return false;
    }
  } else if (p_mode == YBus) {
    return YMBranch::matrixForwardSize(isize,jsize);
  }
  return false;
}
bool gridpack::powerflow::PFBranch::matrixReverseSize(int *isize, int *jsize) const
{
  if (p_mode == Jacobian) {
    gridpack::powerflow::PFBus *bus1
      = dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
    gridpack::powerflow::PFBus *bus2
      = dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
    bool ok = !bus1->getReferenceBus();
    ok = ok && !bus2->getReferenceBus();
    ok = ok && !bus1->isIsolated();
    ok = ok && !bus2->isIsolated();
    ok = ok && (p_active);
    if (ok) {
#ifdef LARGE_MATRIX
      *isize = 2;
      *jsize = 2;
      return true;
#else
      bool bus1PV = bus1->isPV();
      bool bus2PV = bus2->isPV();
      if (bus1PV && bus2PV) {
        *isize = 1;
        *jsize = 1;
        return true;
      } else if (bus1PV) {
        *isize = 2;
        *jsize = 1;
        return true;
      } else if (bus2PV) {
        *isize = 1;
        *jsize = 2;
        return true;
      } else {
        *isize = 2;
        *jsize = 2;
        return true;
      }
#endif
    } else {
      return false;
    }
  } else if (p_mode == YBus) {
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
bool gridpack::powerflow::PFBranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == Jacobian) {
    double rvals[4];
    int nvals = forwardJacobianValues(rvals);
    for (int i=0; i<nvals; i++) values[i] = rvals[i];
    if (nvals == 0) {
      return false;
    } else {
      return true;
    }
  } else if (p_mode == YBus) {
    return YMBranch::matrixForwardValues(values);
  }
  return false;
}

bool gridpack::powerflow::PFBranch::matrixForwardValues(RealType *values)
{
  if (p_mode == Jacobian) {
    int nvals = forwardJacobianValues(values);
    if (nvals == 0) {
      return false;
    } else {
      return true;
    }
  }
  return false;
}

bool gridpack::powerflow::PFBranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == Jacobian) {
    double rvals[4];
    int nvals = reverseJacobianValues(rvals);
    for (int i=0; i<nvals; i++) values[i] = rvals[i];
    if (nvals == 0) {
      return false;
    } else {
      return true;
    }
  } else if (p_mode == YBus) {
    return YMBranch::matrixForwardValues(values);
  }
  return false;
}

bool gridpack::powerflow::PFBranch::matrixReverseValues(RealType *values)
{
  if (p_mode == Jacobian) {
    int nvals = reverseJacobianValues(values);
    if (nvals == 0) {
      return false;
    } else {
      return true;
    }
  }
  return false;
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::powerflow::PFBranch::setYBus(void)
{
  YMBranch::setYBus();
  gridpack::ComplexType ret;
  ret = YMBranch::getForwardYBus();
  p_ybusr_frwd = real(ret);
  p_ybusi_frwd = imag(ret);
  ret = YMBranch::getReverseYBus();
  p_ybusr_rvrs = real(ret);
  p_ybusi_rvrs = imag(ret);
  // Not really a contribution to the admittance matrix but might as well
  // calculate phase angle difference between buses at each end of branch
  gridpack::powerflow::PFBus *bus1 =
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  gridpack::powerflow::PFBus *bus2 =
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
  double pi = 4.0*atan(1.0);
  p_theta = (bus1->getPhase() - bus2->getPhase());
}

/**
 * Load values stored in DataCollection object into PFBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::powerflow::PFBranch::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  YMBranch::load(data);

  // This routine may be called more than once so clear all vectors
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
  p_rateA.clear();
  p_rateB.clear();
  p_rateC.clear();
  p_branch_status.clear();
  p_ckt.clear();
  p_ignore.clear();

  bool ok = true;
  data->getValue(BRANCH_NUM_ELEMENTS, &p_elems);
  double rvar;
  int ivar;
  double pi = 4.0*atan(1.0);
  p_active = false;
  ok = data->getValue(CASE_SBASE, &p_sbase);
  int idx;
  for (idx = 0; idx<p_elems; idx++) {
    bool xform = true;
    xform = xform && data->getValue(BRANCH_X, &rvar, idx);
    if (rvar <1.0e-5 && rvar >=0.0) rvar = 1.0e-5;
    if (rvar >-1.0e-5 && rvar <0.0) rvar =-1.0e-5;
    p_reactance.push_back(rvar);
    xform = xform && data->getValue(BRANCH_R, &rvar, idx);
    p_resistance.push_back(rvar);
    ok = data->getValue(BRANCH_SHIFT, &rvar, idx);
    rvar = -rvar*pi/180.0; 
    p_phase_shift.push_back(rvar);
    ok = data->getValue(BRANCH_TAP, &rvar, idx);
    p_tap_ratio.push_back(rvar); 
    if (rvar != 0.0) {
      p_xform.push_back(xform);
    } else {
      p_xform.push_back(false);
    }
    ivar = 1;
    ok = data->getValue(BRANCH_STATUS, &ivar, idx);
    p_branch_status.push_back(static_cast<bool>(ivar));
    if (ivar == 1 && ok) p_active = true;
    std::string tag;
    ok = data->getValue(BRANCH_CKT, &tag, idx);
    p_ckt.push_back(tag);
    bool shunt = true;
    shunt = shunt && data->getValue(BRANCH_B, &rvar, idx);
    p_charging.push_back(rvar);
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_G1, &rvar, idx);
    p_shunt_admt_g1.push_back(rvar);
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_B1, &rvar, idx);
    p_shunt_admt_b1.push_back(rvar);
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_G2, &rvar, idx);
    p_shunt_admt_g2.push_back(rvar);
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_B2, &rvar, idx);
    p_shunt_admt_b2.push_back(rvar);
    p_shunt.push_back(shunt);
    bool rate = true;
    rate = rate && data->getValue(BRANCH_RATING_A,&rvar,idx);
    p_rateA.push_back(rvar);
    rate = rate && data->getValue(BRANCH_RATING_B,&rvar,idx);
    p_rateB.push_back(rvar);
    rate = rate && data->getValue(BRANCH_RATING_C,&rvar,idx);
    p_rateC.push_back(rvar);
    p_ignore.push_back(false);
  }
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::powerflow::PFBranch::setMode(int mode)
{
  if (mode == YBus) {
    YMBranch::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Return the contribution to the Jacobian for the powerflow equations from
 * a branch
 * @param bus: pointer to the bus making the call
 * @param values: an array of 4 doubles that holds return metrix elements
 */
void gridpack::powerflow::PFBranch::getJacobian(gridpack::powerflow::PFBus *bus, double *values)
{
  double v;
  double cs, sn;
  double ybusr, ybusi;
  if (bus == getBus1().get()) {
    gridpack::powerflow::PFBus *bus2 =
      dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
    v = bus2->getVoltage();
    cs = cos(p_theta);
    sn = sin(p_theta);
    ybusr = p_ybusr_frwd;
    ybusi = p_ybusi_frwd;
  } else if (bus == getBus2().get()) {
    gridpack::powerflow::PFBus *bus1 =
      dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
    v = bus1->getVoltage();
    cs = cos(-p_theta);
    sn = sin(-p_theta);
    ybusr = p_ybusr_rvrs;
    ybusi = p_ybusi_rvrs;
  } else {
    // TODO: Some kind of error
    return;
  }
  values[0] = v*(ybusr*sn - ybusi*cs);
  values[1] = -v*(ybusr*cs + ybusi*sn);
  values[2] = (ybusr*cs + ybusi*sn);
  values[3] = (ybusr*sn - ybusi*cs);
}

/**
 * Return contribution to constraints
 * @param p: real part of constraint
 * @param q: imaginary part of constraint
 */
void gridpack::powerflow::PFBranch::getPQ(gridpack::powerflow::PFBus *bus, double *p, double *q)
{
  gridpack::powerflow::PFBus *bus1 = 
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  double v1 = bus1->getVoltage();
  gridpack::powerflow::PFBus *bus2 =
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
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
 * Return complex power for line element
 * @param tag describing line element on branch
 * @return complex power
 */
gridpack::ComplexType gridpack::powerflow::PFBranch::getComplexPower(
        std::string tag)
{
  gridpack::ComplexType vi, vj, Yii, Yij, s;
  s = ComplexType(0.0,0.0);
  gridpack::powerflow::PFBus *bus1 = 
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  vi = bus1->getComplexVoltage();
  gridpack::powerflow::PFBus *bus2 =
    dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
  vj = bus2->getComplexVoltage();
  getLineElements(tag,&Yii,&Yij);
  s = vi*conj(Yii*vi+Yij*vj)*p_sbase;
  return s;
}

/**
 * Write output from branches to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if branch is contributing string to output, false otherwise
 */
bool gridpack::powerflow::PFBranch::serialWrite(char *string, const int bufsize,
                                                const char *signal)
{
  char buf[128];
  gridpack::powerflow::PFBus *bus1
    = dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  gridpack::powerflow::PFBus *bus2
    = dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
  bool ok = p_active;
  if (ok) {


  if (signal == NULL || !strcmp(signal,"flow_str")) {
    bool rating = false;
    if (signal != NULL) rating = !strcmp(signal,"flow_str");
    gridpack::ComplexType s;
    std::vector<std::string> tags = getLineTags();
    int i;
    int ilen = 0;
    for (i=0; i<p_elems; i++) {
      s = getComplexPower(tags[i]);
      double p = real(s);
      double q = imag(s);
      if (!p_branch_status[i]) p = 0.0;
      if (!p_branch_status[i]) q = 0.0;
      if (bus1->isIsolated() || bus2->isIsolated()) p=0.0;
      if (bus1->isIsolated() || bus2->isIsolated()) q=0.0;
      double perf = 0.0;
      int viol = 0;
      if (p_rateA[i] > 0.0) {
        perf = abs(s)/p_rateA[i];
        if (perf > 1.0) viol = 1;
        perf = perf*perf;
      }
      if (rating) {
        sprintf(buf, "%6d %6d %s %20.12e %20.12e %20.12e %20.12e %1d\n",
            getBus1OriginalIndex(),getBus2OriginalIndex(),tags[i].c_str(),
            p,q,perf,p_rateA[i],viol);
      } else {
        sprintf(buf, "     %6d      %6d     %s   %12.6f         %12.6f\n",
            getBus1OriginalIndex(),getBus2OriginalIndex(),tags[i].c_str(),
            p,q);
      }
      ilen += strlen(buf);
      if (ilen<bufsize) sprintf(string,"%s",buf);
      string += strlen(buf);
    } 
    return true;
  } else if (!strcmp(signal,"fail_str")) {
    std::vector<std::string> tags = getLineTags();
    int i;
    int ilen = 0;
    for (i=0; i<p_elems; i++) {
        sprintf(buf, "     %6d      %6d     %s   %12.6f         %12.6f %12.6f %12.6f %1d\n",
            getBus1OriginalIndex(),getBus2OriginalIndex(),tags[i].c_str(),
            0.0,0.0,0.0,0.0,0);
      ilen += strlen(buf);
      if (ilen<bufsize) sprintf(string,"%s",buf);
      string += strlen(buf);
    }
    return true;
  } else if (!strcmp(signal,"flow")) {
    gridpack::ComplexType s;
    std::vector<std::string> tags = getLineTags();
    int i;
    bool found = false;
    int ilen = 0;
    for (i=0; i<p_elems; i++) {
      s = getComplexPower(tags[i]);
      double p = real(s);
      double q = imag(s);
      if (!p_branch_status[i]) p = 0.0;
      if (!p_branch_status[i]) q = 0.0;
      if (bus1->isIsolated() || bus2->isIsolated()) p=0.0;
      if (bus1->isIsolated() || bus2->isIsolated()) q=0.0;
      double S = sqrt(p*p+q*q);
      if (S > p_rateA[i] && p_rateA[i] != 0.0){
        sprintf(buf, "     %6d      %6d        %s  %12.6f         %12.6f     %8.2f     %8.2f%s\n",
    	  getBus1OriginalIndex(),getBus2OriginalIndex(),tags[i].c_str(),
          p,q,p_rateA[i],S/p_rateA[i]*100,"%");
        ilen += strlen(buf);
        if (ilen<bufsize) sprintf(string,"%s",buf);
        string += strlen(buf);
        found = true;
      }
    }
    return found;
  } else if (!strcmp(signal,"record")) {
    char *cptr = string;
    double pi = 4.0*atan(1.0);
    int slen = 0;
    int i, idx, jdx;
    for (i = 0; i<p_elems; i++) {
      idx = getBus1OriginalIndex();
      jdx = getBus2OriginalIndex();
      sprintf(buf,"%8d, %8d, %2s,",idx,jdx,p_ckt[i].c_str());
      int len = strlen(buf);
      if (len<=bufsize) {
        sprintf(cptr,"%s",buf);
        slen += len;
        cptr += len;
      }
//      double yi = 0.0;
//      double yj = 0.0;
      sprintf(buf," %12.6f, %12.6f, %12.6f, 0.0, 0.0, %8.2f, %8.2f, %8.2f,",
         p_resistance[i],p_reactance[i],p_charging[i],p_rateA[i],p_rateB[i],p_rateC[i]);
      len = strlen(buf);
      if (slen+len<=bufsize) {
        sprintf(cptr,"%s",buf);
        slen += len;
        cptr += len;
      }
      double rval = -180.0*p_phase_shift[i]/pi;
      idx = 0;
      if (p_branch_status[i]) idx = 1;
      sprintf(buf," %12.4f, %12.4f, %1d\n",p_tap_ratio[i],rval,idx);
      if (slen+len<=bufsize) {
        sprintf(cptr,"%s",buf);
        slen += len;
        cptr += len;
      }
    }
    return true;
  }
  return false;
  } // ok loop
  return false;
}

/**
 * Get the status of the branch element
 * @param tag character string identifying branch element
 * @return status of branch element
 */
bool gridpack::powerflow::PFBranch::getBranchStatus(std::string tag)
{
  int i;
  int bsize = p_branch_status.size();
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      return p_branch_status[i];
    }
  }
  return false;
}

/**
 * Set the status of the branch element
 * @param tag character string identifying branch element
 * @param status status of branch element
 */
void gridpack::powerflow::PFBranch::setBranchStatus(std::string tag, bool status)
{
  int i;
  int bsize = p_branch_status.size();
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      p_branch_status[i] = status;
      YMBranch::setLineStatus(tag,status);
      return;
    }
  }
}

/**
 * get branch rating value
 * @param tag transmission element ID
 * @return branch rating value
 */
double gridpack::powerflow::PFBranch::getBranchRatingA(std::string tag)
{
  int i;
  int bsize = p_ckt.size();
  double ret = 0.0;
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      return p_rateA[i];
    }
  }
  return ret;
}

/**
 * get branch rating B value
 * @param tag transmission element ID
 * @return branch rating value
 */
double gridpack::powerflow::PFBranch::getBranchRatingB(std::string tag)
{
  int i;
  int bsize = p_ckt.size();
  double ret = 0.0;
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      return p_rateB[i];
    }
  }
  return ret;
}


/**
 * get branch rating C value
 * @param tag transmission element ID
 * @return branch rating value
 */
double gridpack::powerflow::PFBranch::getBranchRatingC(std::string tag)
{
  int i;
  int bsize = p_ckt.size();
  double ret = 0.0;
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      return p_rateC[i];
    }
  }
  return ret;
}

/**
 * Get list of line IDs
 * @return list of line identifiers
 */
std::vector<std::string> gridpack::powerflow::PFBranch::getLineIDs()
{
  return p_ckt;
}

/**
 * Set parameter to ignore voltage violations
 * @param tag identifier of line element
 * @param flag value of ignore parameter
 */
void gridpack::powerflow::PFBranch::setIgnore(std::string tag, bool flag)
{
  int i;
  int bsize = p_ckt.size();
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      p_ignore[i] = flag;
    }
  }
}

/**
 * Get parameter to ignore voltage violations
 * @param tag identifier of line element
 * @return value of ignore parameter
 */
bool gridpack::powerflow::PFBranch::getIgnore(std::string tag)
{
  int i;
  int bsize = p_ckt.size();
  bool ret = false;
  for (i=0; i<bsize; i++) {
    if (tag == p_ckt[i]) {
      return p_ignore[i];
    }
  }
  return ret;
}

/**
 * Evaluate off-diagonal block of Jacobian for power flow calculation
 * and return result as an array of real values
 * @param rvals values of Jacobian block
 * @return number of values returned
 */
int gridpack::powerflow::PFBranch::forwardJacobianValues(double *rvals)
{
  gridpack::powerflow::PFBus *bus1
    = dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  gridpack::powerflow::PFBus *bus2
    = dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
  bool ok = !bus1->getReferenceBus();
  ok = ok && !bus2->getReferenceBus();
  ok = ok && !bus1->isIsolated();
  ok = ok && !bus2->isIsolated();
  ok = ok && (p_active);
  int nvals;
  if (ok) {
    double t11, t12, t21, t22;
    double cs = cos(p_theta);
    double sn = sin(p_theta);
    bool bus1PV = bus1->isPV();
    bool bus2PV = bus2->isPV();
#ifdef LARGE_MATRIX
    rvals[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
    rvals[1] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
    rvals[2] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
    rvals[3] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
    rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
    rvals[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
    rvals[2] *= bus1->getVoltage();
    rvals[3] *= bus1->getVoltage();
    // fix up matrix if one or both buses at the end of the branch is a PV bus
    if (bus1PV && bus2PV) {
      rvals[1] = 0.0;
      rvals[2] = 0.0;
      rvals[3] = 0.0;
    } else if (bus1PV) {
      rvals[1] = 0.0;
      rvals[3] = 0.0;
    } else if (bus2PV) {
      rvals[2] = 0.0;
      rvals[3] = 0.0;
    }
    nvals = 4;
#else
    if (bus1PV && bus2PV) {
      rvals[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      nvals = 1;
    } else if (bus1PV) {
      rvals[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
      rvals[1] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[1] *= bus1->getVoltage();
      nvals = 2;
    } else if (bus2PV) {
      rvals[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
      rvals[1] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      nvals = 2;
    } else {
      rvals[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
      rvals[1] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
      rvals[2] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
      rvals[3] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[2] *= bus1->getVoltage();
      rvals[3] *= bus1->getVoltage();
      nvals = 4;
    }  
#endif
    return nvals;
  } else {
    return 0;
  }
}

int gridpack::powerflow::PFBranch::reverseJacobianValues(double *rvals)
{
  gridpack::powerflow::PFBus *bus1
    = dynamic_cast<gridpack::powerflow::PFBus*>(getBus1().get());
  gridpack::powerflow::PFBus *bus2
    = dynamic_cast<gridpack::powerflow::PFBus*>(getBus2().get());
  bool ok = !bus1->getReferenceBus();
  ok = ok && !bus2->getReferenceBus();
  ok = ok && !bus1->isIsolated();
  ok = ok && !bus2->isIsolated();
  ok = ok && (p_active);
  int nvals;
  if (ok) {
    double t11, t12, t21, t22;
    double cs = cos(-p_theta);
    double sn = sin(-p_theta);
    bool bus1PV = bus1->isPV();
    bool bus2PV = bus2->isPV();
#ifdef LARGE_MATRIX
    rvals[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
    rvals[1] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
    rvals[2] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
    rvals[3] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
    rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
    rvals[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
    rvals[2] *= bus2->getVoltage();
    rvals[3] *= bus2->getVoltage();
    // fix up matrix if one or both buses at the end of the branch is a PV bus
    if (bus1PV && bus2PV) {
      rvals[1] = 0.0;
      rvals[2] = 0.0;
      rvals[3] = 0.0;
    } else if (bus1PV) {
      rvals[2] = 0.0;
      rvals[3] = 0.0;
    } else if (bus2PV) {
      rvals[1] = 0.0;
      rvals[3] = 0.0;
    }
    nvals = 4;
#else
    if (bus1PV && bus2PV) {
      rvals[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      nvals = 1;
    } else if (bus1PV) {
      rvals[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
      rvals[1] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      nvals = 2;
    } else if (bus2PV) {
      rvals[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
      rvals[1] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[1] *= bus2->getVoltage();
      nvals = 2;
    } else {
      rvals[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
      rvals[1] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
      rvals[2] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
      rvals[3] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
      rvals[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      rvals[2] *= bus2->getVoltage();
      rvals[3] *= bus2->getVoltage();
      nvals = 4;
    } 
#endif
    return nvals;
  } else {
    return 0;
  }
}
