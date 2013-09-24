// -------------------------------------------------------------
/**
 * @file   dynsim_components.cpp
 * @author Shuangshuang Jin 
 * @date   September 19, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/applications/dynamic_simulation/dynsim_components.hpp"
#include "gridpack/parser/dictionary.hpp"

/**
 *  Simple constructor
 */
gridpack::dynsim::DynSimBus::DynSimBus(void)
{
  p_shunt_gs = 0.0;
  p_shunt_bs = 0.0;
  p_mode = YBUS;
  setReferenceBus(false);
  //ngen = 0;
  p_ngen = 0;
}

/**
 *  Simple destructor
 */
gridpack::dynsim::DynSimBus::~DynSimBus(void)
{
}

/**
 *  Return size of matrix block contributed by the component
 *  @param isize, jsize: number of rows and columns of matrix block
 *  @return: false if network component does not contribute matrix element
 */
bool gridpack::dynsim::DynSimBus::matrixDiagSize(int *isize, int *jsize) const
{
  if (p_mode == JACOBIAN && getReferenceBus()) {
    *isize = 0;
    *jsize = 0;
    return false;
  } else if (p_mode == JACOBIAN) {
    *isize = 2;
    *jsize = 2;
  } else if (p_mode == GENERATOR) {
  } else if (p_mode == PERM) {
    if (p_ngen > 0) {
      *isize = 1;
      *jsize = p_ngen;
    } else {
      return false;
    }
  } else if (p_mode == YA) {
      *isize = 2*p_ngen;
      *jsize = 2*p_ngen;
    } else {
      return false;
    }
  }
  return true;
}

/**
 * Return the values of the matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::dynsim::DynSimBus::matrixDiagValues(void *values)
{
  if (p_mode == YL) {
    if (getReferenceBus()) { // getReferenceBus() == true: bus type is 3  
      p_pl = p_pl - p_pg[0];
      p_ql = p_ql - p_qg[0];
    }
    values[0] = p_ybusr + p_pl / (p_voltage * p_voltage); 
    values[1] = p_ybusi + (-p_ql) / (p_voltage * p_voltage);
    values[2] = -values[1];
    values[3] = values[0];
    return true;
  } else if (p_mode == PERM) {
    if (p_ngen > 0) {
      for (i = 0; i < p_ngen; i++) {
        values[i] = 1;
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == YA) {
    for (i = 0; i < p_ngen*p_ngen; i++) {
      int ip = i % p_ngen;
      int jp = (i - ip) / p_ngen;
      if (ip == jp) {
        double ra = p_r[ip] * p_sbase / p_mva[ip];
        double xd;
        if (p_dstr[ip] == 0) {
          xd = p_dtr[ip] * p_sbase / p_mva[ip];
        gridpack::ComplexType Y_a(ra, xd);
        Y_a = 1.0 / Y_a;
        values[4*i] = real(Y_a);
        values[4*i+1] = imag(Y_a);
        values[4*i+2] = -imag(Y_a);
        values[4*i+3] = real(Y_a);
      } else {
        values[4*i] = 0;
        values[4*i+1] = 0;
        values[4*i+2] = 0;
        values[4*i+3] = 0;
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == YB) {
    if (p_isGen) {

    } 
  }
 
  /*if (p_mode == YBUS) {
    gridpack::ComplexType ret(p_ybusr,p_ybusi);
    values[0] = p_ybusr;
    values[1] = p_ybusi;
    values[2] = -p_ybusi;
    values[3] = p_ybusr;
    return true;
  } else if (p_mode == JACOBIAN) {
    if (!getReferenceBus()) {
      double branch_values[4];
      // TODO: More stuff here
      std::vector<boost::shared_ptr<BaseComponent> > branches;
      getNeighborBranches(branches);
      int size = branches.size();
      int i;
      values[0] = 0.0;
      values[1] = 0.0;
      values[2] = 2.0*p_v*p_ybusr;
      values[3] = -2.0*p_v*p_ybusi;
      // HACK: Need to cast pointer, is there a better way?
      for (i=0; i<size; i++) {
	(dynamic_cast<gridpack::dynsim::DynSimBranch*>(branches[i].get()))->
          getJacobian(this, branch_values);
        values[0] -= p_v*branch_values[0];
        values[1] += p_v*branch_values[1];
        values[2] += p_v*branch_values[2];
        values[3] += p_v*branch_values[3];
      }
    } else {
      return false;
    }
  } else if (p_mode == GENERATOR) {
  }*/
}

/**
 * Return the size of the block that this component contributes to the
 * vector
 * @param size: size of vector block
 * @return: false if component does not contribute to vector
 */
bool gridpack::dynsim::DynSimBus::vectorSize(int *size) const
{
  if (p_mode == JACOBIAN && getReferenceBus()) {
    *size = 0;
    return false;
  if (p_mode == JACOBIAN) {
    *size = 2;
  } else if (p_mode == GENERATOR) {
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
bool gridpack::dynsim::DynSimBus::vectorValues(void *values)
{
  if (p_mode == JACOBIAN) {
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    int size = branches.size();
    int i;
    double P, Q, p, q;
    P = 0.0;
    Q = 0.0;
    for (i=0; i<size; i++) {
      gridpack::dynsim::DynSimBranch *branch
        = dynamic_cast<gridpack::dynsim::DynSimBranch*>(branches[i].get());
      branch->getPQ(this, &p, &q);
      P += p;
      Q += q;
    }
    P -= p_P0;
    Q -= p_Q0;
    values[0] = P;
    values[1] = Q;
    return true;
  } else if (p_mode == GENERATOR) {
  }
}

void gridpack::dynsim::DynSimBus::setYBus(void)
{
  gridpack::ComplexType ret(0.0,0.0);
  std::vector<boost::shared_ptr<BaseComponent> > branches;
  getNeighborBranches(branches);
  int size = branches.size();
  int i;
  // HACK: Need to cast pointer, is there a better way?
  for (i=0; i<size; i++) {
    gridpack::dynsim::DynSimBranch *branch
      = dynamic_cast<gridpack::dynsim::DynSimBranch*>(branches[i].get());
    ret -= branch->getAdmittance();
    ret += branch->getTransformer(this);
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
 * Load values stored in DataCollection object into DynSimBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::dynsim::DynSimBus::load(
  const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  data->getValue(BUS_VOLTAGE_MAG, &p_voltage);
  
  p_load = true;
  p_load = p_load && data->getValue(LOAD_PL, &p_pl);
  p_load = p_load && data->getValue(LOAD_QL, &p_ql);

  // Check to see if bus is reference bus
  int itype;
  data->getValue(BUS_TYPE, &itype);
  if (itype == 3) {
    setReferenceBus(true);
  }
 
  if (itype == 3 || itype == 2) {
    p_isGen = true;
  }

  bool lgen;
  int i, gstatus;
  double pg, qg, mva, r, dstr, dtr;
  if (data->getValue(GENERATOR_NUMBER, &p_ngen)) {
    for (i=0; i<p_ngen; i++) {
      //ngen++;
      lgen = true;
      lgen = lgen && data->getValue(GENERATOR_PG, &pg,i);
      lgen = lgen && data->getValue(GENERATOR_QG, &qg,i);
      lgen = lgen && data->getValue(GENERATOR_STAT, &gstatus,i);
      //p_gen = p_gen && data->getValue(GENERATOR_MBASE, &p_mva, 0); // mva?
      //p_gen = p_gen && data->getValue(GENERATOR_R, &p_r, 0); // r?
      //p_gen = p_gen && data->getValue(GENERATOR_DSTR, &p_dstr, 0); // dstr?
      //p_gen = p_gen && data->getValue(GENERATOR_DTR, &p_dtr, 0); // dtr
      //printf("ng=%d,pg=%f,qg=%f\n",i,pg,qg);
      if (lgen) {
        p_pg.push_back(pg);
        p_qg.push_back(qg);
        p_gstatus.push_back(gstatus);
        p_mva.push_back(mva);
        p_r.push_back(r);
        p_dstr.push_back(dstr);
        p_dtr.push_back(dtr);
      }
    }
  }

  p_sbase = 100.0;

}

/**
 * Return the value of the voltage magnitude on this bus
 * @return: voltage magnitude
 */
double gridpack::dynsim::DynSimBus::getVoltage(void)
{
}

/**
 * Return the value of the phase angle on this bus
 * @return: phase angle
 */
double gridpack::dynsim::DynSimBus::getPhase(void)
{
}

/**
 * Return whether or not the bus is a generator
 * @return: true if bus is geneartor
 */
 bool gridpack::dynsim::DynSimBus::isGen(void)
 {
   return p_isGen;
 }


/**
 *  Simple constructor
 */
gridpack::dynsim::DynSimBranch::DynSimBranch(void)
{
  p_reactance = 0.0;
  p_resistance = 0.0;
  p_tap_ratio = 1.0;
  p_phase_shift = 0.0;
  p_charging = 0.0;
  p_shunt_admt_g1 = 0.0;
  p_shunt_admt_b1 = 0.0;
  p_shunt_admt_g2 = 0.0;
  p_shunt_admt_b2 = 0.0;
  p_mode = YBUS;
}

/**
 *  Simple destructor
 */
gridpack::dynsim::DynSimBranch::~DynSimBranch(void)
{
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::dynsim::DynSimBranch::matrixForwardSize(int *isize, int *jsize) const
{
  if (p_mode == JACOBIAN) {
    gridpack::dynsim::DynSimBus *bus1 
       = dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus1().get());
    gridpack::dynsim::DynSimBus *bus2 
       = dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus2().get());
    bool ok = !bus1->getReferenceBus();
    ok = ok & !bus2->getReferenceBus();
    if (ok) {
      *isize = 2;
      *jsize = 2;
      return true;
    } else {
      *isize = 0;
      *jsize = 0;
      return false;
    }
  } else if (p_mode == YBUS) {
    *isize = 2;
    *jsize = 2;
    return true;
  } else if (p_mode == GENERATOR) {
  }
}
bool gridpack::dynsim::DynSimBranch::matrixReverseSize(int *isize, int *jsize) const
{
  if (p_mode == JACOBIAN) {
    gridpack::dynsim::DynSimBus *bus1 
       = dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus1().get());
    gridpack::dynsim::DynSimBus *bus2 
       = dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus2().get());
    bool ok = !bus1->getReferenceBus();
    ok = ok & !bus2->getReferenceBus();
    if (ok) {
      *isize = 2;
      *jsize = 2;
      return true;
    } else {
      *isize = 0;
      *jsize = 0;
      return false;
    }
  } else if (p_mode == YBUS) {
    *isize = 2;
    *jsize = 2;
    return true;
  } else if (p_mode == GENERATOR) {
  }
}

/**
 * Return the values of the off-diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::dynsim::DynSimBranch::matrixForwardValues(void *values)
{
  if (p_mode == JACOBIAN) {
    gridpack::dynsim::DynSimBus *bus1 
       = dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus1().get());
    gridpack::dynsim::DynSimBus *bus2 
       = dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus2().get());
    bool ok = !bus1->getReferenceBus();
    ok = ok & !bus2->getReferenceBus();
    if (ok) {
      double t11, t12, t21, t22;
      double cs = cos(p_theta);
      double sn = sin(p_theta);
      values[0] = (p_ybusr*sn - p_ybusi*cs);
      values[1] = (p_ybusr*cs + p_ybusi*sn);
      values[2] = (p_ybusr*cs + p_ybusi*sn);
      values[3] = (p_ybusr*sn - p_ybusi*cs);
      values[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      values[1] *= bus1->getVoltage();
      values[2] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      values[3] *= bus1->getVoltage();
      return true;
    } else {
      return false;
    }
  } else if (p_mode == YBUS) {
    values[0] = p_ybusr;
    values[1] = p_ybusi;
    values[2] = -p_ybusi;
    values[3] = p_ybusr;
    return true;
  } else if (p_mode == GENERATOR) {
  }
}
bool gridpack::dynsim::DynSimBranch::matrixReverseValues(void *values)
{
  if (p_mode == JACOBIAN) {
    gridpack::dynsim::DynSimBus *bus1 
       = dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus1().get());
    gridpack::dynsim::DynSimBus *bus2 
       = dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus2().get());
    bool ok = !bus1->getReferenceBus();
    ok = ok & !bus2->getReferenceBus();
    if (ok) {
      double t11, t12, t21, t22;
      double cs = cos(-p_theta);
      double sn = sin(-p_theta);
      values[0] = (p_ybusr*sn - p_ybusi*cs);
      values[1] = (p_ybusr*cs + p_ybusi*sn);
      values[2] = (p_ybusr*cs + p_ybusi*sn);
      values[3] = (p_ybusr*sn - p_ybusi*cs);
      values[0] *= ((bus1->getVoltage())*(bus2->getVoltage()));
      values[1] *= bus2->getVoltage();
      values[2] *= -((bus1->getVoltage())*(bus2->getVoltage()));
      values[3] *= bus2->getVoltage();
      return true;
    } else {
      return false;
    }
  } else if (p_mode == YBUS) {
    values[0] = p_ybusr;
    values[1] = p_ybusi;
    values[2] = -p_ybusi;
    values[3] = p_ybusr;
    return true;
  } else if (p_mode == GENERATOR) {
  }
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::dynsim::DynSimBranch::setYBus(void)
{
  gridpack::ComplexType ret(p_resistance,p_reactance);
  ret = -1.0/ret;
  gridpack::ComplexType a(cos(p_phase_shift),sin(p_phase_shift));
  a = p_tap_ratio*a;
  ret = ret - ret/conj(a);
  p_ybusr = real(ret);
  p_ybusi = imag(ret);
  // Not really a contribution to the admittance matrix but might as well
  // calculate phase angle difference between buses at each end of branch
  gridpack::dynsim::DynSimBus *bus1 
     = dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus1().get());
  gridpack::dynsim::DynSimBus *bus2 
     = dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus2().get());
  p_theta = bus1->getPhase() - bus2->getPhase();
   
}

/**
 * Load values stored in DataCollection object into DynSimBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::dynsim::DynSimBranch::load(
  const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  bool ok = true;
  ok = ok && data->getValue(BRANCH_REACTANCE, &p_reactance);
  ok = ok && data->getValue(BRANCH_RESISTANCE, &p_resistance);
  p_xform = true;
  p_xform = p_xform && data->getValue(BRANCH_TAP_RATIO, &p_tap_ratio);
  p_xform = p_xform && data->getValue(BRANCH_PHASE_SHIFT, &p_phase_shift);
  p_shunt = true;
  p_shunt = p_shunt && data->getValue(BRANCH_CHARGING, &p_charging);
  p_shunt = p_shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_G1, &p_shunt_admt_g1);
  p_shunt = p_shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_B1, &p_shunt_admt_b1);
  p_shunt = p_shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_G2, &p_shunt_admt_g2);
  p_shunt = p_shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_B2, &p_shunt_admt_b2);
}

/**
 * Return the complex admittance of the branch
 * @return: complex addmittance of branch
 */
/*gridpack::ComplexType gridpack::dynsim::DynSimBranch::getAdmittance(void)
{
  gridpack::ComplexType ret(p_resistance, p_reactance);
  return -1.0/ret;
}*/

/**
 * Return transformer contribution from the branch to the calling
 * bus
 * @param bus: pointer to the bus making the call
 * @return: contribution to Y matrix from branch
 */
/*gridpack::ComplexType
gridpack::dynsim::DynSimBranch::getTransformer(gridpack::dynsim::DynSimBus *bus)
{
  if (p_xform) {
    gridpack::ComplexType ret(p_resistance,p_reactance);
    ret = -1.0/ret;
    // HACK: pointer comparison, maybe could handle this better
    if (bus == getBus1().get()) {
      ret = ret/(p_tap_ratio*p_tap_ratio);
    } else if (bus == getBus2().get()) {
      // No further action required
    } else {
      // TODO: Some kind of error
    }
    return ret;
  } else {
    gridpack::ComplexType ret(0.0,0.0);
    return ret;
  }
}*/

/**
 * Return the contribution to a bus from shunts
 * @param bus: pointer to the bus making the call
 * @return: contribution to Y matrix from shunts associated with branches
 */
/*gridpack::ComplexType
gridpack::dynsim::DynSimBranch::getShunt(gridpack::dynsim::DynSimBus *bus)
{
  double retr, reti;
  if (p_shunt) {
    retr = 0.5*p_charging;
    reti = 0.0;
    // HACK: pointer comparison, maybe could handle this better
    if (bus == getBus1().get()) {
      retr += p_shunt_admt_g1;
      reti += p_shunt_admt_b1;
    } else if (bus == getBus2().get()) {
      retr += p_shunt_admt_g2;
      reti += p_shunt_admt_b2;
    } else {
      // TODO: Some kind of error
    }
  } else {
    retr = 0.0;
    reti = 0.0;
  }
  return gridpack::ComplexType(retr,reti);
}*/

/**
 * Return the contribution to the Jacobian for the dynsim equations from
 * a branch
 * @param bus: pointer to the bus making the call
 * @param values: an array of 4 doubles that holds return metrix elements
 */
/*void gridpack::dynsim::DynSimBranch::getJacobian(gridpack::dynsim::DynSimBus *bus, double *values)
{
  double v;
  double cs, sn;
  if (bus == getBus1().get()) {
    boost::shared_ptr<gridpack::dynsim::DynSimBus>
      bus2(dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus2().get()));
    v = bus2->getVoltage();
    cs = cos(p_theta);
    sn = sin(p_theta);
  } else if (bus == getBus2().get()) {
    boost::shared_ptr<gridpack::dynsim::DynSimBus>
      bus1(dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus1().get()));
    v = bus1->getVoltage();
    cs = cos(-p_theta);
    sn = sin(-p_theta);
  } else {
    // TODO: Some kind of error
  }
  values[0] = v*(p_ybusr*sn - p_ybusi*cs);
  values[1] = v*(p_ybusr*cs + p_ybusi*sn);
  values[2] = (p_ybusr*cs + p_ybusi*sn);
  values[3] = (p_ybusr*sn - p_ybusi*cs);
}*/

/**
 * Return contribution to constraints
 * @param p: real part of constraint
 * @param q: imaginary part of constraint
 */
/*void gridpack::dynsim::DynSimBranch::getPQ(gridpack::dynsim::DynSimBus *bus, double *p, double *q)
{
  double cs, sn;
  if (bus == getBus1().get()) {
    cs = cos(p_theta);
    sn = sin(p_theta);
  } else if (bus == getBus2().get()) {
    cs = cos(-p_theta);
    sn = sin(-p_theta);
  } else {
    // TODO: Some kind of error
  }
    boost::shared_ptr<gridpack::dynsim::DynSimBus>
      bus1(dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus1().get()));
    double v1 = bus1->getVoltage();
    boost::shared_ptr<gridpack::dynsim::DynSimBus>
      bus2(dynamic_cast<gridpack::dynsim::DynSimBus*>(getBus2().get()));
    double v2 = bus2->getVoltage();
    *p = v1*v2*(p_ybusr*cs+p_ybusi*sn);
    *q = v1*v2*(p_ybusr*sn-p_ybusi*cs);
}*/
