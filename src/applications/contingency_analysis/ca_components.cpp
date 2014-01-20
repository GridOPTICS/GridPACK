/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ca_components.cpp
 * @author Yousu Chen
 * @date   January 20, 2014
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
#include "gridpack/applications/contingency_analysis/ca_components.hpp"
#include "gridpack/parser/dictionary.hpp"

/**
 *  Simple constructor
 */
gridpack::contingency_analysis::CABus::CABus(void)
{
  p_shunt_gs = 0.0;
  p_shunt_bs = 0.0;
  p_mode = YBUS;
  setReferenceBus(false);
  p_ngen = 0;
}

/**
 *  Simple destructor
 */
gridpack::contingency_analysis::CABus::~CABus(void)
{
}

/**
 *  Return size of matrix block contributed by the component
 *  @param isize, jsize: number of rows and columns of matrix block
 *  @return: false if network component does not contribute matrix element
 */
bool gridpack::contingency_analysis::CABus::matrixDiagSize(int *isize, int *jsize) const
{
  if (p_mode == YBUS) {
    *isize = 1;
    *jsize = 1;
  } else {
    *isize = 1;
    *jsize = 1;
  }
  return true;
}

/**
 * Return the values of the matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::contingency_analysis::CABus::matrixDiagValues(ComplexType *values)
{
  if (p_mode == YBUS) {
    gridpack::ComplexType ret(p_ybusr,p_ybusi);
    values[0] = ret;
    return true;
  } 
}

/**
 * Return the size of the block that this component contributes to the
 * vector
 * @param size: size of vector block
 * @return: false if component does not contribute to vector
 */
bool gridpack::contingency_analysis::CABus::vectorSize(int *size) const
{
}

/**
 * Return the values of the vector block
 * @param values: pointer to vector values
 * @return: false if network component does not contribute
 *        vector element
 */
bool gridpack::contingency_analysis::CABus::vectorValues(ComplexType *values)
{
}

void gridpack::contingency_analysis::CABus::setValues(ComplexType *values)
{
}

void gridpack::contingency_analysis::CABus::setYBus(void)
{
  gridpack::ComplexType ret(0.0,0.0);
  std::vector<boost::shared_ptr<BaseComponent> > branches;
  getNeighborBranches(branches);
  int size = branches.size();
  int i;
  // HACK: Need to cast pointer, is there a better way?
  for (i=0; i<size; i++) {
    gridpack::contingency_analysis::CABranch *branch
      = dynamic_cast<gridpack::contingency_analysis::CABranch*>(branches[i].get());
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
 * Load values stored in DataCollection object into CABus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::contingency_analysis::CABus::load(
  const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  p_sbase = 100.0;

  data->getValue(BUS_VOLTAGE_ANG, &p_angle);
  data->getValue(BUS_VOLTAGE_MAG, &p_voltage);

  double pi = 4.0*atan(1.0);
  p_angle = p_angle*pi/180.0; 
 
  p_shunt = true;
  p_shunt = p_shunt && data->getValue(BUS_SHUNT_GL, &p_shunt_gs);
  p_shunt = p_shunt && data->getValue(BUS_SHUNT_BL, &p_shunt_bs);
  p_shunt_gs /= p_sbase;
  p_shunt_bs /= p_sbase; 

  // Check to see if bus is reference bus
  data->getValue(BUS_TYPE, &p_type);
  if (p_type == 3) {
    setReferenceBus(true);
  }
 
  p_load = true;
  p_load = p_load && data->getValue(LOAD_PL, &p_pl);
  p_load = p_load && data->getValue(LOAD_QL, &p_ql);

  p_pl /= p_sbase;
  p_ql /= p_sbase;

  bool lgen;
  int i, gstatus;
  double pg, qg, mva, r, dstr, dtr;
  double h, d0;
  if (data->getValue(GENERATOR_NUMBER, &p_ngen)) {
    for (i=0; i<p_ngen; i++) {
      lgen = true;
      lgen = lgen && data->getValue(GENERATOR_PG, &pg,i);
      lgen = lgen && data->getValue(GENERATOR_QG, &qg,i);
      lgen = lgen && data->getValue(GENERATOR_STAT, &gstatus,i);
      pg /= p_sbase;
      qg /= p_sbase;
      lgen = lgen && data->getValue(GENERATOR_MBASE, &mva, i); 
      if (lgen) {
        p_pg.push_back(pg);
        p_qg.push_back(qg);
        p_gstatus.push_back(gstatus);
        p_mva.push_back(mva);
      }
    }
  }

}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::contingency_analysis::CABus::setMode(int mode)
{
  p_mode = mode;
}

/**
 * Return the value of the voltage magnitude on this bus
 * @return: voltage magnitude
 */
double gridpack::contingency_analysis::CABus::getVoltage(void)
{
}

/**
 * Return the value of the phase angle on this bus
 * @return: phase angle
 */
double gridpack::contingency_analysis::CABus::getPhase(void)
{
}

/**
 *  Simple constructor
 */
gridpack::contingency_analysis::CABranch::CABranch(void)
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
  p_mode = YBUS;
}

/**
 *  Simple destructor
 */
gridpack::contingency_analysis::CABranch::~CABranch(void)
{
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::contingency_analysis::CABranch::matrixForwardSize(int *isize, int *jsize) const
{
  if (p_mode == YBUS) {
    if (p_active) {
      *isize = 1;
      *jsize = 1;
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}
bool gridpack::contingency_analysis::CABranch::matrixReverseSize(int *isize, int *jsize) const
{
  if (p_mode == YBUS) {
    if (p_active) {
      *isize = 1;
      *jsize = 1;
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

/**
 * Return the values of the off-diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::contingency_analysis::CABranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == YBUS) {
    if (p_active) {
      values[0] = gridpack::ComplexType(p_ybusr_frwd,p_ybusi_frwd); 
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}
bool gridpack::contingency_analysis::CABranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == YBUS) {
    if (p_active) {
      values[0] = gridpack::ComplexType(p_ybusr_rvrs,p_ybusi_rvrs);
      return true;
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::contingency_analysis::CABranch::setYBus(void)
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
  // Not really a contribution to the admittance matrix but might as well
  // calculate phase angle difference between buses at each end of branch
  gridpack::contingency_analysis::CABus *bus1 =
    dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus1().get());
  gridpack::contingency_analysis::CABus *bus2 =
    dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus2().get());
  double pi = 4.0*atan(1.0);
  p_theta = (bus1->getPhase() - bus2->getPhase());
}

/**
 * Load values stored in DataCollection object into CABranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::contingency_analysis::CABranch::load(
  const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
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
void gridpack::contingency_analysis::CABranch::setMode(int mode)
{
  p_mode = mode;
}

/**
 * Return the complex admittance of the branch
 * @return: complex addmittance of branch
 */
gridpack::ComplexType gridpack::contingency_analysis::CABranch::getAdmittance(void)
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
gridpack::contingency_analysis::CABranch::getTransformer(gridpack::contingency_analysis::CABus *bus)
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
gridpack::contingency_analysis::CABranch::getShunt(gridpack::contingency_analysis::CABus *bus)
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

}

