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
#include "ca_components.hpp"
#include "gridpack/parser/dictionary.hpp"

//#define LARGE_MATRIX

/**
 *  Simple constructor
 */
gridpack::contingency_analysis::CABus::CABus(void)
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
  setReferenceBus(false);
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
  if (p_mode == Jacobian) {
    if (!YMBus::isIsolated()) {
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
bool gridpack::contingency_analysis::CABus::matrixDiagValues(ComplexType *values)
{
  if (p_mode == YBus) {
    return YMBus::matrixDiagValues(values);
  } else if (p_mode == Jacobian) {
    if (!YMBus::isIsolated()) {
#ifdef LARGE_MATRIX
      if (!getReferenceBus()) {
        values[0] = -p_Qinj - p_ybusi * p_v *p_v; 
        values[1] = p_Pinj - p_ybusr * p_v *p_v; 
        values[2] = p_Pinj / p_v + p_ybusr * p_v; 
        values[3] = p_Qinj / p_v - p_ybusi * p_v; 
        // Fix up matrix elements if bus is PV bus
        if (p_isPV) {
          values[1] = 0.0;
          values[2] = 0.0;
          values[3] = 1.0;
        }
        return true;
      } else {
        values[0] = 1.0;
        values[1] = 0.0;
        values[2] = 0.0;
        values[3] = 1.0;
        return true;
      }
#else
      if (!getReferenceBus() && !p_isPV) {
        values[0] = -p_Qinj - p_ybusi * p_v *p_v; 
        values[1] = p_Pinj - p_ybusr * p_v *p_v; 
        values[2] = p_Pinj / p_v + p_ybusr * p_v; 
        values[3] = p_Qinj / p_v - p_ybusi * p_v; 
        // Fix up matrix elements if bus is PV bus
        return true;
      } else if (!getReferenceBus() && p_isPV) {
        values[0] = -p_Qinj - p_ybusi * p_v *p_v; 
        return true;
      } else {
        return false;
      }
#endif
    } else {
      return false;
    }
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
  if (p_mode == RHS || p_mode == State) {
    if (!YMBus::isIsolated()) {
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
  } else if (p_mode == S_Cal ){
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
bool gridpack::contingency_analysis::CABus::vectorValues(ComplexType *values)
{
  if (p_mode == S_Cal)  {
    double retr = p_v * cos(p_a);
    double reti = p_v * sin(p_a);
    gridpack::ComplexType ret(retr, reti);
    values[0] = ret;
    //printf ("retr = %f, reti = %f, p_v = %f, p_a = %f\n", retr, reti, p_v, p_a);
    return true;
  }
  if (p_mode == State) {
    values[0] = p_v;
    values[1] = p_a;
    return true;
  }
  if (p_mode == RHS) {
    if (!YMBus::isIsolated()) {
      if (!getReferenceBus()) {
        std::vector<boost::shared_ptr<BaseComponent> > branches;
        getNeighborBranches(branches);
        int size = branches.size();
        int i;
        double P, Q, p, q;
        P = 0.0;
        Q = 0.0;
        for (i=0; i<size; i++) {
          gridpack::contingency_analysis::CABranch *branch
            = dynamic_cast<gridpack::contingency_analysis::CABranch*>(branches[i].get());
          branch->getPQ(this, &p, &q);
          P += p;
          Q += q;
          //printf("i=%d:p=%f, q=%f, P=%f, Q=%f\n", i,p,q,P,Q);
        }
        //printf("p_P0=%f,p_Q0=%f\n\n", p_P0,p_Q0);
        // Also add bus i's own Pi, Qi
        P += p_v*p_v*p_ybusr;
        Q += p_v*p_v*(-p_ybusi);
        p_Pinj = P;
        p_Qinj = Q;
        //printf("p = %f, q = %f\n", p_voltage*p_voltage*p_ybusr, p_voltage*p_voltage*(-p_ybusi));
        P -= p_P0;
        Q -= p_Q0;
        values[0] = P;
#ifdef LARGE_MATRIX
        if (!p_isPV) {
          values[1] = Q;
        } else {
          values[1] = 0.0;
        }
#else
        if (!p_isPV) {
          values[1] = Q;
        }
#endif
        return true;
      } else {
#ifdef LARGE_MATRIX
        values[0] = 0.0;
        values[1] = 0.0;
        return true;
#else
        return false;
#endif
      }
    } else {
      return false;
    }
  }
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto buses 
 * @param values array containing voltage magnitude and angle
 */
void gridpack::contingency_analysis::CABus::setValues(gridpack::ComplexType *values)
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
  *p_vAng_ptr = p_a;
  *p_vMag_ptr = p_v;
//  printf("at: %12.6f vt: %12.6f da: %12.6f dv: %12.6f  p_a: %12.6f p_v: %12.6f\n",
//      at,vt,real(values[0]),real(values[1]),p_a,p_v);
}

/**
 * Return the size of the buffer used in data exchanges on the network.
 * For this problem, the voltage magnitude and phase angle need to be exchanged
 * @return size of buffer
 */
int gridpack::contingency_analysis::CABus::getXCBufSize(void)
{
  return 2*sizeof(double);
}

/**
 * Assign pointers for voltage magnitude and phase angle
 */
void gridpack::contingency_analysis::CABus::setXCBuf(void *buf)
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
 * Load values stored in DataCollection object into CABus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::contingency_analysis::CABus::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  YMBus::load(data);

  bool ok = data->getValue(CASE_SBASE, &p_sbase);
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

  // if BUS_TYPE = 2 then bus is a PV bus
  p_isPV = false;
  // if (itype == 2) p_isPV = true;

  // added p_pg,p_qg,p_pl,p_ql,p_sbase;
  p_load = true;
  p_load = p_load && data->getValue(LOAD_PL, &p_pl);
  p_load = p_load && data->getValue(LOAD_QL, &p_ql);
  //printf("p_pl=%f,p_ql=%f\n",p_pl,p_ql);
  bool lgen;
  int i, ngen, gstatus;
  double pg, qg, vs;
  ngen = 0;
  if (data->getValue(GENERATOR_NUMBER, &ngen)) {
    for (i=0; i<ngen; i++) {
      lgen = true;
      lgen = lgen && data->getValue(GENERATOR_PG, &pg,i);
      lgen = lgen && data->getValue(GENERATOR_QG, &qg,i);
      lgen = lgen && data->getValue(GENERATOR_VS, &vs,i);
      lgen = lgen && data->getValue(GENERATOR_STAT, &gstatus,i);
      if (lgen) {
        p_pg.push_back(pg);
        p_qg.push_back(qg);
        p_gstatus.push_back(gstatus);
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
void gridpack::contingency_analysis::CABus::setYBus(void)
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
gridpack::ComplexType gridpack::contingency_analysis::CABus::getYBus(void)
{
  return YMBus::getYBus();
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::contingency_analysis::CABus::setMode(int mode)
{
  if (mode == YBus) {
    YMBus::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Return the value of the voltage magnitude on this bus
 * @return: voltage magnitude
 */
double gridpack::contingency_analysis::CABus::getVoltage()
{
  return *p_vMag_ptr;
}

/**
 * Return whether or not the bus is a PV bus (V held fixed in powerflow
 * equations)
 * @return true if bus is PV bus
 */
bool gridpack::contingency_analysis::CABus::isPV(void)
{
  return p_isPV;
}

/**
 * Return whether or not a bus is isolated
 * @return true if bus is isolated
 */
bool gridpack::contingency_analysis::CABus::isIsolated(void) const
{
  return YMBus::isIsolated();
}

/**
 * Return the value of the phase angle on this bus
 * @return: phase angle
 */
double gridpack::contingency_analysis::CABus::getPhase()
{
  return *p_vAng_ptr;
}

/**
 * setGBus
 */
void gridpack::contingency_analysis::CABus::setGBus(void)
{
  //if (p_gstatus == 1) 
}

/**
 * setSBus
 BUS = (CG*(GEN(ON,PG) + J*GEN(ON,QG)-(PD+J*QD))/BASEMVA
 */
void gridpack::contingency_analysis::CABus::setSBus(void)
{
  // need to update later to consider multiple generators located at the same bus 
  // Chen 8_27_2013 (DONE, 9/29/2013)
#if 1
  // TODO: Need to fix this so that is works for more than 1 generator per bus
  int i;
  double pg, qg;
  pg = 0.0;
  qg = 0.0;
  bool usegen = false;
  for (i=0; i<p_gstatus.size(); i++) {
    if (p_gstatus[i] == 1) {
      pg += p_pg[i];
      qg += p_qg[i];
      usegen = true;
    }
  }
  if (p_gstatus.size() > 0 && usegen) {
    gridpack::ComplexType sBus((pg - p_pl) / p_sbase, (qg - p_ql) / p_sbase);
    //p_sbusr = real(sBus);
    //p_sbusr = real(sBus);
    p_P0 = real(sBus);
    p_Q0 = imag(sBus);
    //printf("p_P0=%f, p_Q0=%f\n",p_P0,p_Q0);
  } else {
    gridpack::ComplexType sBus((- p_pl) / p_sbase, (- p_ql) / p_sbase);
    p_P0 = real(sBus);
    p_Q0 = imag(sBus);
  } 
  //printf("p_P0=%f, p_Q0=%f\n",p_P0,p_Q0);
#endif
}

/**
 * Write output from buses to standard out
 * @param string (output) string with information to be printed out
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool gridpack::contingency_analysis::CABus::serialWrite(char *string, const char *signal)
{
  if (signal == NULL) {
    double pi = 4.0*atan(1.0);
    double angle = p_a*180.0/pi;
    sprintf(string, "     %6d      %12.6f         %12.6f\n",
        getOriginalIndex(),angle,p_v);
  } else if (!strcmp(signal,"pq")) {
    gridpack::ComplexType v[2];
    vectorValues(v);
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    sprintf(string, "     %6d      %12.6f         %12.6f      %2d\n",
        getOriginalIndex(),real(v[0]),real(v[1]),branches.size());
  }
  return true;
}

/**
 * Return the complex voltage on this bus
 * @return the complex voltage
 */
gridpack::ComplexType gridpack::contingency_analysis::CABus::getComplexVoltage(void)
{
  gridpack::ComplexType ret(cos(p_a),sin(p_a));
  ret = ret*p_v;
  return ret;
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
  p_mode = YBus;
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
  if (p_mode == Jacobian) {
    gridpack::contingency_analysis::CABus *bus1
      = dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus1().get());
    gridpack::contingency_analysis::CABus *bus2
      = dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus2().get());
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
/*    *isize = 2;
      *jsize = 2;
      return true; */
    }
  } else if (p_mode == YBus) {
    return YMBranch::matrixForwardSize(isize,jsize);
  }
}
bool gridpack::contingency_analysis::CABranch::matrixReverseSize(int *isize, int *jsize) const
{
  if (p_mode == Jacobian) {
    gridpack::contingency_analysis::CABus *bus1
      = dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus1().get());
    gridpack::contingency_analysis::CABus *bus2
      = dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus2().get());
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
}

/**
 * Return the values of the off-diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::contingency_analysis::CABranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == Jacobian) {
    gridpack::contingency_analysis::CABus *bus1
      = dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus1().get());
    gridpack::contingency_analysis::CABus *bus2
      = dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus2().get());
    bool ok = !bus1->getReferenceBus();
    ok = ok && !bus2->getReferenceBus();
    ok = ok && !bus1->isIsolated();
    ok = ok && !bus2->isIsolated();
    ok = ok && (p_active);
    if (ok) {
      double t11, t12, t21, t22;
      double cs = cos(p_theta);
      double sn = sin(p_theta);
      bool bus1PV = bus1->isPV();
      bool bus2PV = bus2->isPV();
#ifdef LARGE_MATRIX
      values[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
      values[1] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
      values[2] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
      values[3] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
      values[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      values[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      values[2] *= bus1->getVoltage();
      values[3] *= bus1->getVoltage();
      // fix up matrix if one or both buses at the end of the branch is a PV bus
      if (bus1PV && bus2PV) {
        values[1] = 0.0;
        values[2] = 0.0;
        values[3] = 0.0;
      } else if (bus1PV) {
        values[1] = 0.0;
        values[3] = 0.0;
      } else if (bus2PV) {
        values[2] = 0.0;
        values[3] = 0.0;
      }
#else
      if (bus1PV && bus2PV) {
        values[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
        values[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      } else if (bus1PV) {
        values[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
        values[1] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
        values[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
        values[1] *= bus1->getVoltage();
      } else if (bus2PV) {
        values[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
        values[1] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
        values[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
        values[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      } else {
        values[0] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
        values[1] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
        values[2] = (p_ybusr_frwd*cs + p_ybusi_frwd*sn);
        values[3] = (p_ybusr_frwd*sn - p_ybusi_frwd*cs);
        values[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
        values[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
        values[2] *= bus1->getVoltage();
        values[3] *= bus1->getVoltage();
      }  
#endif
      return true;
    } else {
      return false;
    }
  } else if (p_mode == YBus) {
//    values[0] = p_ybusr_frwd;
//    values[1] = p_ybusi_frwd;
//    values[2] = -p_ybusi_frwd;
//    values[3] = p_ybusr_frwd;
    return YMBranch::matrixForwardValues(values);
  }
}

bool gridpack::contingency_analysis::CABranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == Jacobian) {
    gridpack::contingency_analysis::CABus *bus1
      = dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus1().get());
    gridpack::contingency_analysis::CABus *bus2
      = dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus2().get());
    bool ok = !bus1->getReferenceBus();
    ok = ok && !bus2->getReferenceBus();
    ok = ok && !bus1->isIsolated();
    ok = ok && !bus2->isIsolated();
    ok = ok && (p_active == 1);
    if (ok) {
      double t11, t12, t21, t22;
      double cs = cos(-p_theta);
      double sn = sin(-p_theta);
      bool bus1PV = bus1->isPV();
      bool bus2PV = bus2->isPV();
#ifdef LARGE_MATRIX
      values[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
      values[1] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
      values[2] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
      values[3] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
      values[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      values[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      values[2] *= bus2->getVoltage();
      values[3] *= bus2->getVoltage();
      // fix up matrix if one or both buses at the end of the branch is a PV bus
      if (bus1PV && bus2PV) {
        values[1] = 0.0;
        values[2] = 0.0;
        values[3] = 0.0;
      } else if (bus1PV) {
        values[2] = 0.0;
        values[3] = 0.0;
      } else if (bus2PV) {
        values[1] = 0.0;
        values[3] = 0.0;
      }
#else
      if (bus1PV && bus2PV) {
        values[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
        values[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      } else if (bus1PV) {
        values[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
        values[1] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
        values[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
        values[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      } else if (bus2PV) {
        values[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
        values[1] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
        values[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
        values[1] *= bus2->getVoltage();
      } else {
        values[0] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
        values[1] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
        values[2] = (p_ybusr_rvrs*cs + p_ybusi_rvrs*sn);
        values[3] = (p_ybusr_rvrs*sn - p_ybusi_rvrs*cs);
        values[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
        values[1] *= -((bus1->getVoltage())*(bus2->getVoltage()));
        values[2] *= bus2->getVoltage();
        values[3] *= bus2->getVoltage();
      } 
#endif
      return true;
    } else {
      return false;
    }
  } else if (p_mode == YBus) {
  //  values[0] = p_ybusr_rvrs;
  //  values[1] = p_ybusi_rvrs;
  //  values[2] = -p_ybusi_rvrs;
  //  values[3] = p_ybusr_rvrs;
    return YMBranch::matrixForwardValues(values);
  }
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::contingency_analysis::CABranch::setYBus(void)
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
  gridpack::contingency_analysis::CABus *bus1 =
    dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus1().get());
  gridpack::contingency_analysis::CABus *bus2 =
    dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus2().get());
//  if (p_xform) {
//    printf ("from %d-> to %d: p_phase_shift = %f, a = %f+%fi\n", bus1->getOriginalIndex(), bus2->getOriginalIndex(), p_phase_shift, real(a), imag(a) );
//  }
  //p_theta = bus1->getPhase() - bus2->getPhase();
  double pi = 4.0*atan(1.0);
  p_theta = (bus1->getPhase() - bus2->getPhase());
  //printf("p_phase_shift: %12.6f\n",p_phase_shift);
  //printf("p_theta: %12.6f\n",p_theta);
  //printf("p_tap_ratio: %12.6f\n",p_tap_ratio);

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
  YMBranch::load(data);

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
  if (mode == YBus) {
    YMBranch::setMode(gridpack::ymatrix::YBus);
  }
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

/**
 * Return the contribution to the Jacobian for the powerflow equations from
 * a branch
 * @param bus: pointer to the bus making the call
 * @param values: an array of 4 doubles that holds return metrix elements
 */
void gridpack::contingency_analysis::CABranch::getJacobian(gridpack::contingency_analysis::CABus *bus, double *values)
{
  double v;
  double cs, sn;
  double ybusr, ybusi;
  if (bus == getBus1().get()) {
    gridpack::contingency_analysis::CABus *bus2 =
      dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus2().get());
    v = bus2->getVoltage();
    cs = cos(p_theta);
    sn = sin(p_theta);
    ybusr = p_ybusr_frwd;
    ybusi = p_ybusi_frwd;
  } else if (bus == getBus2().get()) {
    gridpack::contingency_analysis::CABus *bus1 =
      dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus1().get());
    v = bus1->getVoltage();
    cs = cos(-p_theta);
    sn = sin(-p_theta);
    ybusr = p_ybusr_rvrs;
    ybusi = p_ybusi_rvrs;
  } else {
    // TODO: Some kind of error
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
void gridpack::contingency_analysis::CABranch::getPQ(gridpack::contingency_analysis::CABus *bus, double *p, double *q)
{
  gridpack::contingency_analysis::CABus *bus1 = 
    dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus1().get());
  double v1 = bus1->getVoltage();
  gridpack::contingency_analysis::CABus *bus2 =
    dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus2().get());
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
  }
  //*p = -v1*v2*(p_ybusr*cs-p_ybusi*sn);
  //*q = v1*v2*(p_ybusr*sn+p_ybusi*cs);
  *p = v1*v2*(ybusr*cs+ybusi*sn);
  *q = v1*v2*(ybusr*sn-ybusi*cs);
//  printf("v1=%f, v2=%f, cs=%f, sn=%f, p_ybusr=%f, p_ybusi=%f\n", v1,v2,cs,sn,p_ybusr,p_ybusi);
//  printf("*p=%f,*q=%f\n",*p,*q);
}

/**
 * Write output from branches to standard out
 * @param string (output) string with information to be printed out
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if branch is contributing string to output, false otherwise
 */
bool gridpack::contingency_analysis::CABranch::serialWrite(char *string, const char *signal)
{
  gridpack::ComplexType v1, v2, y, s;
  gridpack::contingency_analysis::CABus *bus1 = 
    dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus1().get());
  v1 = bus1->getComplexVoltage();
  gridpack::contingency_analysis::CABus *bus2 =
    dynamic_cast<gridpack::contingency_analysis::CABus*>(getBus2().get());
  v2 = bus2->getComplexVoltage();
  y = gridpack::ComplexType(p_ybusr_frwd,p_ybusi_frwd);
  s = v1*conj(y*(v1-v2));
  double p = real(s)*p_sbase;
  double q = imag(s)*p_sbase;
  sprintf(string, "     %6d      %6d      %12.6f         %12.6f\n",
      bus1->getOriginalIndex(),bus2->getOriginalIndex(),p,q);
  return true;
}
