/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ymatrix_components.cpp
 * @author Bruce Palmer
 * @date   2014-01-28 11:30:59 d3g096
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
#include "ymatrix_components.hpp"
#include "gridpack/parser/dictionary.hpp"

//#define LARGE_MATRIX

/**
 *  Simple constructor
 */
gridpack::ymatrix::YMBus::YMBus(void)
{
  p_shunt_gs = 0.0;
  p_shunt_bs = 0.0;
  p_mode = YBus;
  setReferenceBus(false);
}

/**
 *  Simple destructor
 */
gridpack::ymatrix::YMBus::~YMBus(void)
{
}

/**
 *  Return size of matrix block contributed by the component
 *  @param isize, jsize: number of rows and columns of matrix block
 *  @return: false if network component does not contribute matrix element
 */
bool gridpack::ymatrix::YMBus::matrixDiagSize(int *isize, int *jsize) const
{
  if (p_mode == YBus) {
    *isize = 1;
    *jsize = 1;
    return true;
  } else {
    return false;
  }
}

/**
 * Return the values of the matrix block. The values are
 * returned in row-major order.
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::ymatrix::YMBus::matrixDiagValues(ComplexType *values)
{
  if (p_mode == YBus) {
    gridpack::ComplexType ret(p_ybusr,p_ybusi);
    values[0] = ret;
    return true;
  } else {
    return false;
  }
}

/**
 * Set values of YBus matrix. These can then be used in subsequent
 * calculations
 */
void gridpack::ymatrix::YMBus::setYBus(void)
{
  gridpack::ComplexType ret(0.0,0.0);
  std::vector<boost::shared_ptr<BaseComponent> > branches;
  getNeighborBranches(branches);
  int size = branches.size();
  int i;
  // HACK: Need to cast pointer, is there a better way?
  for (i=0; i<size; i++) {
    gridpack::ymatrix::YMBranch *branch
      = dynamic_cast<gridpack::ymatrix::YMBranch*>(branches[i].get());
    ret -= branch->getAdmittance();
    ret -= branch->getTransformer(this);
    ret += branch->getShunt(this);
  }
  if (p_shunt) {
    gridpack::ComplexType shunt(p_shunt_gs,p_shunt_bs);
    ret += shunt;
  }
  p_ybusr = real(ret);
  p_ybusi = imag(ret);
}

gridpack::ComplexType gridpack::ymatrix::YMBus::getYBus(void)
{
  gridpack::ComplexType ret(p_ybusr,p_ybusi);
  return ret;
}
  
/**
 * Load values stored in DataCollection object into YMBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::ymatrix::YMBus::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  double sbase;
  data->getValue(CASE_SBASE, &sbase);
  p_shunt = true;
  p_shunt = p_shunt && data->getValue(BUS_SHUNT_GL, &p_shunt_gs);
  p_shunt = p_shunt && data->getValue(BUS_SHUNT_BL, &p_shunt_bs);
  p_shunt_gs /= sbase;
  p_shunt_bs /= sbase;
  // Check to see if bus is reference bus
  int itype;
  data->getValue(BUS_TYPE, &itype);
  if (itype == 3) {
    setReferenceBus(true);
  }
}


/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::ymatrix::YMBus::setMode(int mode)
{
  p_mode = mode;
}

/**
 *  Simple constructor
 */
gridpack::ymatrix::YMBranch::YMBranch(void)
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
  p_mode = YBus;
}

/**
 *  Simple destructor
 */
gridpack::ymatrix::YMBranch::~YMBranch(void)
{
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::ymatrix::YMBranch::matrixForwardSize(int *isize, int *jsize) const
{
  if (p_mode == YBus) {
    if (p_active) {
      *isize = 1;
      *jsize = 1;
      return true;
    } else {
      return false;
    }
  }
}
bool gridpack::ymatrix::YMBranch::matrixReverseSize(int *isize, int *jsize) const
{
  if (p_mode == YBus) {
    if (p_active == 1) {
      *isize = 1;
      *jsize = 1;
      return true;
    } else {
      return false;
    }
  }
}

/**
 * Return the values of the off-diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::ymatrix::YMBranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == YBus) {
    if (p_active == 1) {
      values[0] = gridpack::ComplexType(p_ybusr_frwd,p_ybusi_frwd);
      return true;
    } else {
      return false;
    }
  }
}

bool gridpack::ymatrix::YMBranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == YBus) {
    if (p_active) {
      values[0] = gridpack::ComplexType(p_ybusr_rvrs,p_ybusi_rvrs);
      return true;
    } else {
      return false;
    }
  }
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::ymatrix::YMBranch::setYBus(void)
{
  int i;
  p_ybusr_frwd = 0.0;
  p_ybusi_frwd = 0.0;
  p_ybusr_rvrs = 0.0;
  p_ybusi_rvrs = 0.0;
  for (i=0; i<p_elems; i++) {
    gridpack::ComplexType ret(p_resistance[i],p_reactance[i]);
    ret = -1.0/ret;
    gridpack::ComplexType a(cos(p_phase_shift[i]),sin(p_phase_shift[i]));
    a = p_tap_ratio[i]*a;
    if (p_branch_status[i] == 1) {
      if (p_xform[i]) {
        p_ybusr_frwd += real(ret/conj(a));
        p_ybusi_frwd += imag(ret/conj(a));
        p_ybusr_rvrs += real(ret/a);
        p_ybusi_rvrs += imag(ret/a);
      } else {
        p_ybusr_frwd += real(ret);
        p_ybusi_frwd += imag(ret);
        p_ybusr_rvrs += real(ret);
        p_ybusi_rvrs += imag(ret);
      }
    }
  }
}

/**
 * Load values stored in DataCollection object into YMBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::ymatrix::YMBranch::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  bool ok = true;
  data->getValue(BRANCH_NUM_ELEMENTS, &p_elems);
  double rvar;
  int ivar;
  double pi = 4.0*atan(1.0);
  p_active = false;
  int idx;
  for (idx = 0; idx<p_elems; idx++) {
    bool xform = true;
    xform = xform && data->getValue(BRANCH_X, &rvar, idx);
    p_reactance.push_back(rvar);
    xform = xform && data->getValue(BRANCH_R, &rvar, idx);
    p_resistance.push_back(rvar);
    ok = ok && data->getValue(BRANCH_SHIFT, &rvar, idx);
    rvar = -rvar*pi/180.0; 
    p_phase_shift.push_back(rvar);
    ok = ok && data->getValue(BRANCH_TAP, &rvar, idx);
    p_tap_ratio.push_back(rvar); 
    if (rvar != 0.0) {
      p_xform.push_back(xform);
    } else {
      p_xform.push_back(false);
    }
    ivar = 1;
    data->getValue(BRANCH_STATUS, &ivar, idx);
    p_branch_status.push_back(ivar);
    if (ivar == 1) p_active = true;
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
  }
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::ymatrix::YMBranch::setMode(int mode)
{
  p_mode = mode;
}

/**
 * Return the complex admittance of the branch
 * @return: complex addmittance of branch
 */
gridpack::ComplexType gridpack::ymatrix::YMBranch::getAdmittance(void)
{
  int i;
  gridpack::ComplexType ret(0.0,0.0);
  for (i=0; i<p_elems; i++) {
    gridpack::ComplexType tmp(p_resistance[i], p_reactance[i]);
    if (!p_xform[i] && p_branch_status[i] == 1) {
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
gridpack::ymatrix::YMBranch::getTransformer(gridpack::ymatrix::YMBus *bus)
{
  int i;
  gridpack::ComplexType ret(0.0,0.0);
  for (i=0; i<p_elems; i++) {
    gridpack::ComplexType tmp(p_resistance[i],p_reactance[i]);
    gridpack::ComplexType tmpB(0.0,0.5*p_charging[i]);
    if (p_xform[i] && p_branch_status[i] == 1) {
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
gridpack::ymatrix::YMBranch::getShunt(gridpack::ymatrix::YMBus *bus)
{
  double retr, reti;
  retr = 0.0;
  reti = 0.0;
  int i;
  for (i=0; i<p_elems; i++) {
    double tmpr, tmpi;
    if (p_shunt[i] && p_branch_status[i] == 1) {
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
