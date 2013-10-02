// -------------------------------------------------------------
/**
 * @file   ds_components.cpp
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
#include "gridpack/applications/dynamic_simulation/ds_components.hpp"
#include "gridpack/parser/dictionary.hpp"

/**
 *  Simple constructor
 */
gridpack::dynamic_simulation::DSBus::DSBus(void)
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
gridpack::dynamic_simulation::DSBus::~DSBus(void)
{
}

/**
 *  Return size of matrix block contributed by the component
 *  @param isize, jsize: number of rows and columns of matrix block
 *  @return: false if network component does not contribute matrix element
 */
bool gridpack::dynamic_simulation::DSBus::matrixDiagSize(int *isize, int *jsize) const
{
  /*if (p_mode == JACOBIAN && getReferenceBus()) {
    *isize = 0;
    *jsize = 0;
    return false;
  } else if (p_mode == JACOBIAN) {
    *isize = 2;
    *jsize = 2;
  } else if (p_mode == GENERATOR) {
  } else if (p_mode == PERM) {*/
  if (p_mode == YBUS || p_mode == YL) {
    *isize = 2;
    *jsize = 2;
  } else if (p_mode == PERM) {
    if (p_ngen > 0) {
      *isize = 1;
      *jsize = p_ngen;
    } else {
      return false;
    }
  } else if (p_mode == YA) {
    if (p_ngen > 0) {
      *isize = 2*p_ngen;
      *jsize = 2*p_ngen;
    } else {
      return false;
    }
  } else if (p_mode == PMatrix) {
    if (p_ngen > 0) {
      *isize = 2;
      *jsize = 2;
    } else {
      *isize = 2;
      *jsize = 0;
    } 
  } else {
    *isize = 2;
    *jsize = 2;
  }
  return true;
}

/**
 * Return the values of the matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::dynamic_simulation::DSBus::matrixDiagValues(ComplexType *values)
{
  if (p_mode == YBUS) {
    gridpack::ComplexType ret(p_ybusr,p_ybusi);
    values[0] = p_ybusr;
    values[1] = p_ybusi;
    values[2] = -p_ybusi;
    values[3] = p_ybusr;
    return true;
  } else if (p_mode == YL) {
    //if (p_type == 1) { 
      //p_pl = p_pl;// - p_pg[0];
      //p_ql = p_ql;// - p_qg[0];
    //}
    values[0] = p_ybusr + p_pl / (p_voltage * p_voltage); 
    values[1] = p_ybusi + (-p_ql) / (p_voltage * p_voltage);
    values[2] = -values[1];
    values[3] = values[0];
    return true;
  } else if (p_mode == PERM) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        values[i] = 1;
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == YA) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen*p_ngen; i++) {
        int ip = i % p_ngen;
        int jp = (i - ip) / p_ngen;
        int ngen2 = 2*p_ngen;
        int ii1 = 2*jp*ngen2+2*ip;
        int ii2 = 2*jp*ngen2+2*ip+1;
        int ii3 = (2*jp+1)*ngen2+2*ip;
        int ii4 = (2*jp+1)*ngen2+2*ip+1;
        if (ip == jp) {
          double ra = p_r[ip] * p_sbase / p_mva[ip];
          double xd;
          if (p_dstr[ip] == 0) { 
            xd = p_dtr[ip] * p_sbase / p_mva[ip];
          }
          gridpack::ComplexType Y_a(ra, xd);
          Y_a = 1.0 / Y_a;
          values[ii1] = real(Y_a);
          values[ii2] = imag(Y_a);
          values[ii3] = -imag(Y_a);
          values[ii4] = real(Y_a);
        } else {
          values[ii1] = 0;
          values[ii2] = 0;
          values[ii3] = 0;
          values[ii4] = 0;
        }
        return true;
      } 
    } else {
      return false;
    }
  } else if (p_mode == PMatrix) {
    if (p_ngen > 0) {
      values[0] = -1.0;
      values[0] = 0.0;
      values[0] = 0.0;
      values[0] = -1.0;
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
bool gridpack::dynamic_simulation::DSBus::vectorSize(int *size) const
{
  /*if (p_mode == JACOBIAN && getReferenceBus()) {
    *size = 0;
    return false;
  if (p_mode == JACOBIAN) {
    *size = 2;
  } else if (p_mode == GENERATOR) {
    *size = 2;
  }*/
  *size = 2;
  return true;
}

/**
 * Return the values of the vector block
 * @param values: pointer to vector values
 * @return: false if network component does not contribute
 *        vector element
 */
bool gridpack::dynamic_simulation::DSBus::vectorValues(ComplexType *values)
{
  /*if (p_mode == JACOBIAN) {
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    getNeighborBranches(branches);
    int size = branches.size();
    int i;
    double P, Q, p, q;
    P = 0.0;
    Q = 0.0;
    for (i=0; i<size; i++) {
      gridpack::dynamic_simulation::DSBranch *branch
        = dynamic_cast<gridpack::dynamic_simulation::DSBranch*>(branches[i].get());
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
  }*/
}

void gridpack::dynamic_simulation::DSBus::setYBus(void)
{
  gridpack::ComplexType ret(0.0,0.0);
  std::vector<boost::shared_ptr<BaseComponent> > branches;
  getNeighborBranches(branches);
  int size = branches.size();
  int i;
  // HACK: Need to cast pointer, is there a better way?
  for (i=0; i<size; i++) {
    gridpack::dynamic_simulation::DSBranch *branch
      = dynamic_cast<gridpack::dynamic_simulation::DSBranch*>(branches[i].get());
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
 * Load values stored in DataCollection object into DSBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::dynamic_simulation::DSBus::load(
  const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  p_sbase = 100.0;

  data->getValue(BUS_VOLTAGE_ANG, &p_angle);
  data->getValue(BUS_VOLTAGE_MAG, &p_voltage);

  double pi = 4.9*atan(1.0);
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
  if (data->getValue(GENERATOR_NUMBER, &p_ngen)) {
    for (i=0; i<p_ngen; i++) {
      lgen = true;
      lgen = lgen && data->getValue(GENERATOR_PG, &pg,i);
      lgen = lgen && data->getValue(GENERATOR_QG, &qg,i);
      lgen = lgen && data->getValue(GENERATOR_STAT, &gstatus,i);
      pg /= p_sbase;
      qg /= p_sbase;

      lgen = lgen && data->getValue(GENERATOR_MBASE, &mva, i); 
      lgen = lgen && data->getValue(GENERATOR_RESISTANCE, &r, i); // r
      lgen = lgen && data->getValue(GENERATOR_SUBTRANSIENT_REACTANCE, &dstr,i); // dstr
      lgen = lgen && data->getValue(GENERATOR_TRANSIENT_REACTANCE, &dtr,i); // dtr
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


}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::dynamic_simulation::DSBus::setMode(int mode)
{
  p_mode = mode;
}

/**
 * Return the value of the voltage magnitude on this bus
 * @return: voltage magnitude
 */
double gridpack::dynamic_simulation::DSBus::getVoltage(void)
{
}

/**
 * Return the value of the phase angle on this bus
 * @return: phase angle
 */
double gridpack::dynamic_simulation::DSBus::getPhase(void)
{
}


/**
 *  Simple constructor
 */
gridpack::dynamic_simulation::DSBranch::DSBranch(void)
{
  p_reactance = 0.0;
  p_resistance = 0.0;
  //p_tap_ratio = 1.0;
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
gridpack::dynamic_simulation::DSBranch::~DSBranch(void)
{
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::dynamic_simulation::DSBranch::matrixForwardSize(int *isize, int *jsize) const
{
  /*if (p_mode == JACOBIAN) {
    gridpack::dynamic_simulation::DSBus *bus1 
       = dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus1().get());
    gridpack::dynamic_simulation::DSBus *bus2 
       = dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus2().get());
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
  }*/
  if (p_mode == YBUS) {
    *isize = 2;
    *jsize = 2;
    return true;
  }
}
bool gridpack::dynamic_simulation::DSBranch::matrixReverseSize(int *isize, int *jsize) const
{
  /*if (p_mode == JACOBIAN) {
    gridpack::dynamic_simulation::DSBus *bus1 
       = dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus1().get());
    gridpack::dynamic_simulation::DSBus *bus2 
       = dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus2().get());
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
  } else if (p_mode == YBUS) {*/
  if (p_mode == YBUS || p_mode == YL) {
    *isize = 2;
    *jsize = 2;
    return true;
  } else {
    *isize = 2;
    *jsize = 2;
    return true;
  } 
}

/**
 * Return the values of the off-diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::dynamic_simulation::DSBranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == YBUS || p_mode == YL) {
    values[0] = p_ybusr_frwd;
    values[1] = p_ybusi_frwd;
    values[2] = -p_ybusi_frwd;
    values[3] = p_ybusr_frwd;
    return true;
  /*} else if (p_mode == PERM) {
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
    return true;*/
  } else {
    return true;
  }
  /*if (p_mode == JACOBIAN) {
    gridpack::dynamic_simulation::DSBus *bus1 
       = dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus1().get());
    gridpack::dynamic_simulation::DSBus *bus2 
       = dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus2().get());
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
    values[0] = p_ybusr_frwd;
    values[1] = p_ybusi_frwd;
    values[2] = -p_ybusi_frwd;
    values[3] = p_ybusr_frwd;
    return true;
  } else if (p_mode == GENERATOR) {
  }*/
}
bool gridpack::dynamic_simulation::DSBranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == YBUS || p_mode == YL) {
    values[0] = p_ybusr_rvrs;
    values[1] = p_ybusi_rvrs;
    values[2] = -p_ybusi_rvrs;
    values[3] = p_ybusr_rvrs;
    return true;
  /*} else if (p_mode == PERM) {
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
    return true;*/
  } else {
    return true;
  }
  /*if (p_mode == JACOBIAN) {
    gridpack::dynamic_simulation::DSBus *bus1 
       = dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus1().get());
    gridpack::dynamic_simulation::DSBus *bus2 
       = dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus2().get());
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
    values[0] = p_ybusr_rvrs;
    values[1] = p_ybusi_rvrs;
    values[2] = -p_ybusi_rvrs;
    values[3] = p_ybusr_rvrs;
    return true;
  } else if (p_mode == GENERATOR) {
  }*/
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::dynamic_simulation::DSBranch::setYBus(void)
{
  gridpack::ComplexType ret(p_resistance,p_reactance);
  ret = -1.0/ret;
  gridpack::ComplexType a(cos(p_phase_shift),sin(p_phase_shift));
  a = p_tap_ratio*a;
  if (p_xform) {
    p_ybusr_frwd = real(ret/conj(a));
    p_ybusi_frwd = imag(ret/conj(a));
    p_ybusr_rvrs = real(ret/a);
    p_ybusi_rvrs = imag(ret/a);
  } else {
    p_ybusr_frwd = real(ret);
    p_ybusi_frwd = imag(ret);
    p_ybusr_rvrs = real(ret);
    p_ybusi_rvrs = imag(ret);
  }
  // Not really a contribution to the admittance matrix but might as well
  // calculate phase angle difference between buses at each end of branch
  gridpack::dynamic_simulation::DSBus *bus1 = 
    dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus1().get());
  gridpack::dynamic_simulation::DSBus *bus2 =  
    dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus2().get());
  //p_theta = bus1->getPhase() - bus2->getPhase();
  double pi = 4.0*atan(1.0);
  p_theta = (bus1->getPhase() - bus2->getPhase());
  //printf("p_phase_shift: %12.6f\n",p_phase_shift);
  //printf("p_theta: %12.6f\n",p_theta);
  //printf("p_tap_ratio: %12.6f\n",p_tap_ratio);
   
}

/**
 * Load values stored in DataCollection object into DSBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::dynamic_simulation::DSBranch::load(
  const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  bool ok = true;
  ok = ok && data->getValue(BRANCH_X, &p_reactance);
  ok = ok && data->getValue(BRANCH_R, &p_resistance);
  ok = ok && data->getValue(BRANCH_SHIFT, &p_phase_shift);
  double temp;
  ok = ok && data->getValue(BRANCH_TAP, &temp);
  if (temp != 0.0) {
    p_tap_ratio = temp;
    p_xform = true;
    p_xform = p_xform && data->getValue(BRANCH_X, &p_reactance);
    p_xform = p_xform && data->getValue(BRANCH_R, &p_resistance);
  } else {
    p_xform = false;
  }
  p_shunt = true;
  p_shunt = p_shunt && data->getValue(BRANCH_B, &p_charging);
  p_shunt = p_shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_G1, &p_shunt_admt_g1);
  p_shunt = p_shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_B1, &p_shunt_admt_b1);
  p_shunt = p_shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_G2, &p_shunt_admt_g2);
  p_shunt = p_shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_B2, &p_shunt_admt_b2);
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::dynamic_simulation::DSBranch::setMode(int mode)
{
  p_mode = mode;
}

/**
 * Return the complex admittance of the branch
 * @return: complex addmittance of branch
 */
gridpack::ComplexType gridpack::dynamic_simulation::DSBranch::getAdmittance(void)
{
  gridpack::ComplexType ret(p_resistance, p_reactance);
  if (!p_xform) {
    ret = -1.0/ret;
  } else {
    ret = gridpack::ComplexType(0.0,0.0);
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
gridpack::dynamic_simulation::DSBranch::getTransformer(gridpack::dynamic_simulation::DSBus *bus)
{
  gridpack::ComplexType ret(p_resistance,p_reactance);
  if (p_xform) {
    ret = -1.0/ret;
    gridpack::ComplexType a(cos(p_phase_shift),sin(p_phase_shift));
    a = p_tap_ratio*a;
    if (bus == getBus1().get()) {
      ret = ret/(conj(a)*a);
    } else if (bus == getBus2().get()) {
      // ret is unchanged
    } 
  } else {
    ret = gridpack::ComplexType(0.0,0.0);
  }
  return ret;
}

/**
 * Return the contribution to a bus from shunts
 * @param bus: pointer to the bus making the call
 * @return: contribution to Y matrix from shunts associated with branches
 */
gridpack::ComplexType
gridpack::dynamic_simulation::DSBranch::getShunt(gridpack::dynamic_simulation::DSBus *bus)
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
}

