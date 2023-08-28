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
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "ymatrix_components.hpp"
#include "gridpack/parser/dictionary.hpp"

//#define PRINT_DEBUG
//#define LARGE_MATRIX

/**
 *  Simple constructor
 */
gridpack::ymatrix::YMBus::YMBus(void)
{
  p_shunt_gs = 0.0;
  p_shunt_bs = 0.0;
  p_mode = YBus;
  p_isolated = false;
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
  if (p_mode == YBus && !p_isolated) {
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
  if (p_mode == YBus && !p_isolated) {
    gridpack::ComplexType ret(p_ybusr,p_ybusi);
    values[0] = ret;
#ifdef PRINT_DEBUG
    printf("%d %d %f %f\n",getOriginalIndex(),getOriginalIndex(),
        p_ybusr,p_ybusi);
#endif
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

/**
 * Modify diagonal values of matrix.
 * @param rval real part of diagonal matrix element
 * @param ival imaginary part of diagonal matrix element
 */
void gridpack::ymatrix::YMBus::setYBusDiag(double rval, double ival)
{
  p_ybusr = rval;
  p_ybusi = ival;
}

/**
 * Get values of YBus matrix. These can then be used in subsequent
 * calculations
 */
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
  double shunt_binit;
  p_shunt = p_shunt && data->getValue(BUS_SHUNT_GL, &p_shunt_gs,0);
  p_shunt = p_shunt && data->getValue(BUS_SHUNT_BL, &p_shunt_bs,0);
  bool binit = data->getValue(SHUNT_BINIT, &shunt_binit);
  if (binit) p_shunt = true;

  p_shunt_gs /= sbase;
  p_shunt_bs /= sbase;
  // update shunt based on shunt table 
  if (binit) {
    shunt_binit /= sbase; 
    p_shunt_bs = shunt_binit;
  }
  // Check to see if bus is reference bus
  int itype;
  data->getValue(BUS_TYPE, &itype);
  if (itype == 3) {
    setReferenceBus(true);
  } else if (itype == 4) {
    p_isolated = true;
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
 * Return whether or not a bus is isolated
 * @return true if bus is isolated
 */
bool gridpack::ymatrix::YMBus::isIsolated(void) const
{
  return p_isolated;
}

/**
 * Change isolated status of bus
 * @param flag true if bus is isolated
 */
void gridpack::ymatrix::YMBus::setIsolated(const bool flag)
{
  p_isolated = flag;
}

/**
 * Get shunt values
 * @param gl shunt GL value
 * @param bl shunt BL value
 */
void gridpack::ymatrix::YMBus::getShuntValues(double *bl,
    double *gl) const
{
  *bl = p_shunt_bs;
  *gl = p_shunt_gs;
}

/**
 * Set shunt values
 * @param gs shunt GL value
 * @param bs shunt BL value
 */
void gridpack::ymatrix::YMBus::addShuntValues(double gs,double bs)
{
  p_shunt = true;
  p_shunt_bs += bs;
  p_shunt_gs += gs;
}


/**
 * Set internal parameters inside the Y-bus component
 * @param name character string describing component to be modified
 * @param value of parameter to be modified
 * @param idx index (if necessary) of variable to be modified
 */
void gridpack::ymatrix::YMBus::setParam(std::string name, double value,
    int idx)
{
  if (name==BUS_SHUNT_BL) {
    p_shunt_bs = value;
  }
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
  p_switched.clear();
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
    gridpack::ymatrix::YMBus *bus1
      = dynamic_cast<gridpack::ymatrix::YMBus*>(getBus1().get());
    gridpack::ymatrix::YMBus *bus2
      = dynamic_cast<gridpack::ymatrix::YMBus*>(getBus2().get());
    bool ok = !bus1->isIsolated();
    ok = ok && !bus2->isIsolated();
    if (p_active && ok) {
      *isize = 1;
      *jsize = 1;
      return true;
    } else {
      return false;
    }
  }
  return false;
}
bool gridpack::ymatrix::YMBranch::matrixReverseSize(int *isize, int *jsize) const
{
  if (p_mode == YBus) {
    gridpack::ymatrix::YMBus *bus1
      = dynamic_cast<gridpack::ymatrix::YMBus*>(getBus1().get());
    gridpack::ymatrix::YMBus *bus2
      = dynamic_cast<gridpack::ymatrix::YMBus*>(getBus2().get());
    bool ok = !bus1->isIsolated();
    ok = ok && !bus2->isIsolated();
    if (p_active && ok) {
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
 * Return the values of the off-diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::ymatrix::YMBranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == YBus) {
    gridpack::ymatrix::YMBus *bus1
      = dynamic_cast<gridpack::ymatrix::YMBus*>(getBus1().get());
    gridpack::ymatrix::YMBus *bus2
      = dynamic_cast<gridpack::ymatrix::YMBus*>(getBus2().get());
    bool ok = !bus1->isIsolated();
    ok = ok && !bus2->isIsolated();
    if (p_active && ok) {
      values[0] = gridpack::ComplexType(p_ybusr_frwd,p_ybusi_frwd);
#ifdef PRINT_DEBUG
      printf("%d %d %f %f\n",bus1->getOriginalIndex(),
          bus2->getOriginalIndex(),p_ybusr_frwd,p_ybusi_frwd);
#endif
      return true;
    } else {
      return false;
    }
  }
  return false;
}

bool gridpack::ymatrix::YMBranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == YBus) {
    gridpack::ymatrix::YMBus *bus1
      = dynamic_cast<gridpack::ymatrix::YMBus*>(getBus1().get());
    gridpack::ymatrix::YMBus *bus2
      = dynamic_cast<gridpack::ymatrix::YMBus*>(getBus2().get());
    bool ok = !bus1->isIsolated();
    ok = ok && !bus2->isIsolated();
    if (p_active && ok) {
      values[0] = gridpack::ComplexType(p_ybusr_rvrs,p_ybusi_rvrs);
#ifdef PRINT_DEBUG
      printf("%d %d %f %f\n",bus1->getOriginalIndex(),
          bus2->getOriginalIndex(),p_ybusr_rvrs,p_ybusi_rvrs);
#endif
      return true;
    } else {
      return false;
    }
  }
  return false;
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::ymatrix::YMBranch::setYBus(void)
{
  int i;
  p_ybusr_frwd = 0.0;
  p_ybusi_frwd = 0.0;
  p_ybusr_rvrs = 0.0;
  p_ybusi_rvrs = 0.0;
#ifdef USE_ACOPF
  p_yffr.clear();
  p_yffi.clear();
  p_yttr.clear();
  p_ytti.clear();
  p_yftr.clear();
  p_yfti.clear();
  p_ytfr.clear();
  p_ytfi.clear();
#endif
  for (i=0; i<p_elems; i++) {
    gridpack::ComplexType ret(p_resistance[i],p_reactance[i]);
    ret = -1.0/ret;
    gridpack::ComplexType a(cos(p_phase_shift[i]),sin(p_phase_shift[i]));
    if (p_xform[i]) a = p_tap_ratio[i]*a;
    if (p_switched[i]) a = conj(a);
    if (p_branch_status[i]) {
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
#ifdef USE_ACOPF
      gridpack::ComplexType j(0.0,1.0);
      ret = -ret;
      gridpack::ComplexType yff,ytt,yft,ytf;
      ytt = ret + 0.5*p_charging[i]*j;
      yff = ytt/(a*conj(a));
      yft = -ret/conj(a);
      ytf = -ret/a;
      p_yffr.push_back(real(yff));
      p_yffi.push_back(imag(yff));
      p_yttr.push_back(real(ytt));
      p_ytti.push_back(imag(ytt));
      p_yftr.push_back(real(yft));
      p_yfti.push_back(imag(yft));
      p_ytfr.push_back(real(ytf));
      p_ytfi.push_back(imag(ytf));
#endif
    }
  }
}

/**
 * Get values of YBus matrix. These can then be used in subsequent
 * calculations
 * @return forward y-matrix element for branch
 */
gridpack::ComplexType gridpack::ymatrix::YMBranch::getForwardYBus(void)
{
  return gridpack::ComplexType(p_ybusr_frwd,p_ybusi_frwd);
}

/**
 * Get values of YBus matrix. These can then be used in subsequent
 * calculations
 * @return reverse y-matrix element for branch
 */
gridpack::ComplexType gridpack::ymatrix::YMBranch::getReverseYBus(void)
{
  return gridpack::ComplexType(p_ybusr_rvrs,p_ybusi_rvrs);
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
  p_branch_status.clear();
  p_switched.clear();
  p_tag.clear();
  bool ok = true;
  data->getValue(BRANCH_NUM_ELEMENTS, &p_elems);
  double rvar;
  int ivar;
  std::string svar;
  bool lvar;
  double pi = 4.0*atan(1.0);
  p_active = false;
  int idx;
  for (idx = 0; idx<p_elems; idx++) {
    bool xform = true;
    xform = xform && data->getValue(BRANCH_X, &rvar, idx);
    //if (rvar <1.0e-5 && rvar >=0.0) rvar = 1.0e-5;
//    if (rvar >-1.0e-5 && rvar >=0.0) rvar =-1.0e-5;
    p_reactance.push_back(rvar);
    xform = xform && data->getValue(BRANCH_R, &rvar, idx);
    p_resistance.push_back(rvar);
    rvar = 0.0;
    ok = data->getValue(BRANCH_SHIFT, &rvar, idx);
    rvar = rvar*pi/180.0; 
    p_phase_shift.push_back(rvar);
    rvar = 0.0;
    ok = data->getValue(BRANCH_TAP, &rvar, idx);
    p_tap_ratio.push_back(rvar); 
    ok = data->getValue(BRANCH_CKT, &svar, idx);
    p_tag.push_back(svar);
    if (rvar != 0.0) {
      p_xform.push_back(xform);
    } else {
      p_xform.push_back(false);
    }
    ivar = 1;
    data->getValue(BRANCH_STATUS, &ivar, idx);
    p_branch_status.push_back(static_cast<bool>(ivar));
    if (ivar == 1) p_active = true;
    ok = data->getValue(BRANCH_SWITCHED, &lvar, idx);
    if (!ok) lvar = false;
    p_switched.push_back(lvar);
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
gridpack::ymatrix::YMBranch::getTransformer(gridpack::ymatrix::YMBus *bus)
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
      if ((!p_switched[i] && bus == getBus1().get()) ||
          (p_switched[i] && bus == getBus2().get())) {
        tmp = tmp/(conj(a)*a);
      } else {
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
 * Return contributions to Y-matrix from a specific transmission element
 * @param tag character string for transmission element
 * @param Yii contribution from "from" bus
 * @param Yij contribution from line element
 */
void gridpack::ymatrix::YMBranch::getLineElements(const std::string tag,
   gridpack::ComplexType *Yii, gridpack::ComplexType *Yij)// , gridpack::ComplexType *yii)
{
  gridpack::ComplexType zero = gridpack::ComplexType(0.0,0.0);
  gridpack::ComplexType flow = zero;
  *Yii = zero;
  *Yij = zero;
  int i, idx;
  idx = -1;
  // find line element corresponding to tag
  for (i=0; i<p_elems; i++) {
    if (tag == p_tag[i]) {
      idx = i;
      break;
    }
  }
  if (idx >= 0) {
    gridpack::ComplexType yij,aij,bij;
    yij = gridpack::ComplexType(p_resistance[idx],p_reactance[idx]);
    bij = gridpack::ComplexType(0.0,p_charging[idx]);
    //bij = 0.0;
    if (yij != zero) yij = -1.0/yij;
    if (p_xform[idx]) {
      // evaluate flow for transformer
      aij = gridpack::ComplexType(cos(p_phase_shift[idx]),sin(p_phase_shift[idx]));
      aij = p_tap_ratio[idx]*aij;
      if (aij != zero) {
        if (!p_switched[idx]) {
          *Yij = yij/conj(aij);
          *Yii = -(yij-0.5*bij);
          *Yii = (*Yii)/(aij*conj(aij));
        } else {
          *Yij = yij/aij;
          *Yii = -(yij-0.5*bij);
        }
      }
    } else {
      // evaluate flow for regular line
      *Yij = yij;
      *Yii = -((*Yij)-0.5*bij);
    }
  }
}

/**
 * Return contributions to Y-matrix from a specific transmission element
 * @param tag character string for transmission element
 * @param Yii contribution from "from" bus
 * @param Yij contribution from line element
 */
void gridpack::ymatrix::YMBranch::getRvrsLineElements(const std::string tag,
   gridpack::ComplexType *Yii, gridpack::ComplexType *Yij)// , gridpack::ComplexType *yii)
{
  gridpack::ComplexType zero = gridpack::ComplexType(0.0,0.0);
  gridpack::ComplexType flow = zero;
  *Yii = zero;
  *Yij = zero;
  int i, idx;
  idx = -1;
  // find line element corresponding to tag
  for (i=0; i<p_elems; i++) {
    if (tag == p_tag[i]) {
      idx = i;
      break;
    }
  }
  if (idx >= 0) {
    gridpack::ComplexType yij,aij,bij;
    yij = gridpack::ComplexType(p_resistance[idx],p_reactance[idx]);
    bij = gridpack::ComplexType(0.0,p_charging[idx]);
    //bij = 0.0;
    if (yij != zero) yij = -1.0/yij;
    if (p_xform[idx]) {
      // evaluate flow for transformer
      aij = gridpack::ComplexType(cos(p_phase_shift[idx]),sin(p_phase_shift[idx]));
      aij = p_tap_ratio[idx]*aij;
      if (aij != zero) {
        if (p_switched[idx]) {
          *Yij = yij/conj(aij);
          *Yii = -(yij-0.5*bij);
          *Yii = (*Yii)/(aij*conj(aij));
        } else {
          *Yij = yij/aij;
          *Yii = -(yij-0.5*bij);
        }
      }
    } else {
      // evaluate flow for regular line
      *Yij = yij;
      *Yii = -((*Yij)-0.5*bij);
    }
  }
}

/**
 * Return status of all transmission elements
 * @return vector containing status of transmission elements
 */
std::vector<bool> gridpack::ymatrix::YMBranch::getLineStatus()
{
  return p_branch_status;
}

/**
 * Return status of a transmission element based on its tag name
 * @param tag name of transmission element
 * @return status of that transmission element
 */
bool gridpack::ymatrix::YMBranch::getLineStatus(std::string tag)
{
  int i;
  bool found = false;
  for (i=0; i<p_elems; i++) {
    if (tag == p_tag[i]) {
      return p_branch_status[i];
    }
  }
  return found;
}


/**
 * Return tags of all transmission elements
 * @return vector containging tag of transmission elements
 */
std::vector<std::string> gridpack::ymatrix::YMBranch::getLineTags()
{
  return p_tag;
}

/**
 * Set the status of a transmission element based on its tag name
 * @param tag name of transmission element
 * @param status that transmission element should be set to
 * @return false if no transmission element with that name exists
 */
bool gridpack::ymatrix::YMBranch::setLineStatus(std::string tag,
    bool status)
{
  int i;
  bool found = false;
  for (i=0; i<p_elems; i++) {
    if (tag == p_tag[i]) {
      p_branch_status[i] = status;
      found = true;
      break;
    }
  }
  return found;
}

/**
 * Get branch susceptance (charging)
 * @param tag string identifier for transmission element
 * @return value of susceptance
 */
double gridpack::ymatrix::YMBranch::getSusceptance(std::string tag)
{
  int i;
  for (i=0; i<p_elems; i++) {
    if (tag == p_tag[i]) {
      return p_charging[i];
    }
  }
  return 0.0;
}

/**
 * Set internal parameters inside the Y-branch component
 * @param name character string describing component to be modified
 * @param value of parameter to be modified
 * @param idx index (if necessary) of variable to be modified
 */
void gridpack::ymatrix::YMBranch::setParam(std::string name, double value,
    int idx)
{
  if (name==BRANCH_X) {
    p_reactance[idx] = value;
  } else if (name==BRANCH_R) {
    p_resistance[idx] = value;
  } else if (name==BRANCH_TAP) {
    p_tap_ratio[idx] = value;
  } else if (name==BRANCH_SHIFT) {
    p_phase_shift[idx] = value;
  }
}

#ifdef USE_ACOPF
/**
 * Return components from individual transmission elements
 * @param yffr list of real parts of Yff
 * @param yffr list of imaginary parts of Yff
 * @param yttr list of real parts of Ytt
 * @param yttr list of imaginary parts of Ytt
 * @param yftr list of real parts of Yft
 * @param yftr list of imaginary parts of Yft
 * @param ytfr list of real parts of Ytf
 * @param ytfr list of imaginary parts of Ytf
 * @param switched flag on whether line is switched or not
 */
void gridpack::ymatrix::YMBranch::getYElements(
    std::vector<double> &yffr, std::vector<double> &yffi,
    std::vector<double> &yttr, std::vector<double> &ytti,
    std::vector<double> &yftr, std::vector<double> &yfti,
    std::vector<double> &ytfr, std::vector<double> &ytfi,
    std::vector<bool> &switched) {
  yffr = p_yffr;
  yffi = p_yffi;
  yttr = p_yttr;
  ytti = p_ytti;
  yftr = p_yftr;
  yfti = p_yfti;
  ytfr = p_ytfr;
  ytfi = p_ytfi;
  switched = p_switched;
}
#endif
