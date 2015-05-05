/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_components.cpp
 * @author Shuangshuang Jin 
 * @date   2015-04-30 14:48:06 d3g096
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
#include "ds_components.hpp"
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
  p_from_flag = false;
  p_to_flag = false;
  p_branch = NULL;
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
  //if (p_mode == YBUS || p_mode == YL || p_mode = FY || p_mode = POSFY) {
  if (YMBus::isIsolated()) return false;
  if (p_mode == YBUS || p_mode == YL || p_mode == updateYbus || p_mode == onFY ||
      p_mode == posFY) {
    //*isize = 1;
    //*jsize = 1;
    return YMBus::matrixDiagSize(isize,jsize);
  } else if (p_mode == PERM) {
    if (p_ngen > 0) {
      *isize = 1;
      *jsize = p_ngen;
    } else {
      return false;
    }
  } else if (p_mode == PERMTrans) {
    if (p_ngen > 0) {
      *isize = p_ngen;
      *jsize = 1;
    } else {
      return false;
    }
  } else if (p_mode == YA) {
    if (p_ngen > 0) {
      *isize = p_ngen;
      *jsize = p_ngen;
    } else {
      return false;
    }
  } else if (p_mode == YB) {
    // YB is the negative of the transpose of YC
    if (p_ngen > 0) {
      *isize = p_ngen;
      *jsize = 1;
    } else {
      *isize = 0;
      *jsize = 1;
    }
  } else if (p_mode == YC) {
    if (p_ngen > 0) {
      *isize = 1;
      *jsize = p_ngen;
    } else {
      *isize = 1;
      *jsize = 0;
    }
  } else if (p_mode == permYMOD) {
    if (p_ngen > 0) {
      *isize = 1;
      *jsize = 1;
    } else {
      return false;
    }
  } else if (p_mode == PMatrix) {
    if (p_ngen > 0) {
      *isize = 1;
      *jsize = 1;
      //*jsize = p_ngen;
    } else {
      *isize = 1;
      *jsize = 0;
    } 
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
bool gridpack::dynamic_simulation::DSBus::matrixDiagValues(ComplexType *values)
{
  if (p_mode == YBUS) {
    return YMBus::matrixDiagValues(values);
    //gridpack::ComplexType ret(p_ybusr,p_ybusi);
    //values[0] = ret;
    //return true;
  } else if (p_mode == onFY) {
    if (p_from_flag) {
      gridpack::ComplexType ret(0.0, 1.0e7);
      values[0] = ret;
      return true;
    } else {
      return false;
    }
  } else if (p_mode == posFY) {
    if (p_from_flag || p_to_flag) {
      values[0] = dynamic_cast<gridpack::dynamic_simulation::DSBranch*>(p_branch)
          ->getUpdateFactor();
      return true;
    } else {
      return false;
    }
  } else if (p_mode == YL) {
    //printf("p_ybusr=%f, p_ybusi=%f\n", p_ybusr, p_ybusi);
    gridpack::ComplexType ret(p_ybusr, p_ybusi);
    values[0] = ret;
    return true;
  } else if (p_mode == PERM) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        values[i] = 1.0;
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == PERMTrans) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        values[i] = 1.0;
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
        int ii = jp*p_ngen + ip;
        if (ip == jp) {
          double ra = p_r[ip] * p_sbase / p_mva[ip];
          double xd = 0.0;
          if (p_dstr[ip] == 0) { 
            xd = p_dtr[ip] * p_sbase / p_mva[ip];
          }
          gridpack::ComplexType Y_a(ra, xd);
          Y_a = 1.0 / Y_a;
          values[ii] = Y_a;
        } else {
          values[ii] = 0.0;
        }
      } 
      return true;
    } else {
      return false;
    }
  } else if (p_mode == YC) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        double ra = p_r[i] * p_sbase / p_mva[i];
        double xd;
        if (p_dstr[i] == 0) { 
          xd = p_dtr[i] * p_sbase / p_mva[i];
        }
        gridpack::ComplexType Y_a(ra, xd);
        Y_a = 1.0 / Y_a;
        values[i] = Y_a;
      } 
      return true;
    } else {
      return false;
    }
  } else if (p_mode == YB) {
    // YB is the negative of the transpose of YC
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        double ra = p_r[i] * p_sbase / p_mva[i];
        double xd;
        if (p_dstr[i] == 0) { 
          xd = p_dtr[i] * p_sbase / p_mva[i];
        }
        gridpack::ComplexType Y_a(ra, xd);
        Y_a = 1.0 / Y_a;
        values[i] = -Y_a;
      } 
      return true;
    } else {
      return false;
    }
  } else if (p_mode == permYMOD) {
    if (p_ngen > 0) {
      values[0] = 0.0;
      for (int ip = 0; ip < p_ngen; ip++) {
        double ra = p_r[ip] * p_sbase / p_mva[ip];
        double xd;
        if (p_dstr[ip] == 0) { 
          xd = p_dtr[ip] * p_sbase / p_mva[ip];
        }
        gridpack::ComplexType Y_a(ra, xd);
        Y_a = 1.0 / Y_a;
        values[0] += Y_a;
      } 
      return true;
    } else {
      return false;
    }
  } else if (p_mode == PMatrix) {
    if (p_ngen > 0) {
      /*for (int i = 0; i < p_ngen; i++) {
        values[i] = -1.0;
      }*/
      values[0] = 1.0; 
      return true;
    } else {
      return false;
    } 
  } else if (p_mode == updateYbus) {
    if (p_ngen > 0) {
      gridpack::ComplexType permYmod = 0.0;
      for (int ip = 0; ip < p_ngen; ip++) {
        double ra = p_r[ip] * p_sbase / p_mva[ip];
        double xd;
        if (p_dstr[ip] == 0) { 
          xd = p_dtr[ip] * p_sbase / p_mva[ip];
        }
        gridpack::ComplexType Y_a(ra, xd);
        Y_a = 1.0 / Y_a;
        permYmod += Y_a;
      } 
      double ur = p_ybusr + real(permYmod); 
      double ui = p_ybusi + imag(permYmod); 
      gridpack::ComplexType u(ur, ui);
      values[0] = u;
      return true;
    } else {
      gridpack::ComplexType u(p_ybusr, p_ybusi);
      values[0] = u;
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
  if (p_mode == updateYbus) {
    if (p_ngen > 0) {
      *size = 1;
    } else {
      *size = 0;
      return false;
    }
  } else if (p_mode == DAE_init) {
    if (p_ngen > 0) {
      //*size = 4;
    } else {
      return false;
    }
  } else if (p_mode == init_pelect || p_mode == init_eprime ||
      p_mode == init_mac_ang || p_mode == init_mac_spd ||
      p_mode == init_eqprime || p_mode == init_pmech ||
      p_mode == init_mva || p_mode == init_d0 || p_mode == init_h) {
    if (p_ngen > 0) {
      //*size = 1;
      *size = p_ngen;
     } else {
      return false;
     }
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
  if (p_mode == updateYbus) {
    printf("returning Value on %d\n",getOriginalIndex());
  } else if (p_mode == DAE_init) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        p_mva[i] = p_sbase / p_mva[i]; 
        double eterm = p_voltage;
        double pelect = p_pg[i];
        double qelect = p_qg[i];
        double currr = sqrt(pelect * pelect + qelect * qelect) / eterm * p_mva[i];
        double phi = atan2(qelect, pelect);  
        double pi = 4.0*atan(1.0);
        double vi = p_angle;
        gridpack::ComplexType v(0.0, vi);
        v = eterm * exp(v);
        double curri = p_angle - phi;
        gridpack::ComplexType curr(0.0, curri);
        curr = currr * exp(curr);
        gridpack::ComplexType jay(0.0, 1.0);
        gridpack::ComplexType eprime(0.0, 0.0);
        eprime = v + jay * p_dtr[i] * curr;
        double mac_ang = atan2(imag(eprime), real(eprime));
        double mac_spd = 1.0; 
        double eqprime = abs(eprime);
        double pmech = pelect;
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == init_pelect) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        values[i] = p_pg[i];
        p_pelect.push_back(values[i]);
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == init_eprime) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        p_mva[i] = p_sbase / p_mva[i]; 
        double eterm = p_voltage;
        double pelect = p_pg[i];
        double qelect = p_qg[i];
        double currr = sqrt(pelect * pelect + qelect * qelect) / eterm * p_mva[i];
        double phi = atan2(qelect, pelect);  
        double pi = 4.0*atan(1.0);
        double vi = p_angle;
        gridpack::ComplexType v(0.0, vi);
        v = eterm * exp(v);
        double curri = p_angle - phi;
        gridpack::ComplexType curr(0.0, curri);
        curr = currr * exp(curr);
        gridpack::ComplexType jay(0.0, 1.0);
        values[i] = v + jay * p_dtr[i] * curr;
        p_eprime.push_back(values[i]);
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == init_mac_ang) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        //values[0] = atan2(imag(p_eprime[i]), real(p_eprime[i]));
        values[i] = atan2(imag(p_eprime[i]), real(p_eprime[i]));
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == init_mac_spd) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        //values[0] = 1.0; 
        values[i] = 1.0; 
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == init_eqprime) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        //values[0] = abs(p_eprime[i]); 
        values[i] = abs(p_eprime[i]); 
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == init_pmech) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        //values[0] = abs(p_pelect[i]); 
        values[i] = abs(p_pelect[i]); 
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == init_mva) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        //values[0] = p_mva[i]; 
        values[i] = p_mva[i]; 
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == init_d0) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        //values[0] = p_d0[i]; 
        values[i] = p_d0[i]; 
      }
      return true;
    } else {
      return false;
    }
  } else if (p_mode == init_h) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        //values[0] = p_h[i]; 
        values[i] = p_h[i]; 
      }
      return true;
    } else {
      return false;
    }
  }
  return false;
}

void gridpack::dynamic_simulation::DSBus::setValues(ComplexType *values)
{
  int i;
  if (p_mode == updateYbus) {
    if (p_ngen > 0) {
      p_permYmod = values[0];
      //printf("p_permYmod = %f+%fi\n", getOriginalIndex(), real(p_permYmod), imag(p_permYmod));
    }
  } else if (p_mode == init_mac_ang) {
    p_mac_ang_final.clear();
    for (i=0; i<p_ngen; i++) {
      p_mac_ang_final.push_back(values[i]);
    }
  } else if (p_mode == init_mac_spd) {
    p_mac_spd_final.clear();
    for (i=0; i<p_ngen; i++) {
      p_mac_spd_final.push_back(values[i]);
    }
  } else if (p_mode == init_pmech) {
    p_mech_final.clear();
    for (i=0; i<p_ngen; i++) {
      p_mech_final.push_back(values[i]);
    }
  } else if (p_mode == init_pelect) {
    p_elect_final.clear();
    for (i=0; i<p_ngen; i++) {
      p_elect_final.push_back(values[i]);
    }
  }
}

void gridpack::dynamic_simulation::DSBus::setYBus(void)
{
  /*gridpack::ComplexType ret(0.0,0.0);
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
  }*/
  YMBus::setYBus();
  gridpack::ComplexType ret;
  ret = YMBus::getYBus();
  p_ybusr = real(ret);
  p_ybusi = imag(ret);
  p_ybusr = p_ybusr+p_pl/(p_voltage*p_voltage);
  p_ybusi = p_ybusi+(-p_ql)/(p_voltage*p_voltage);
}

/**
 * Get values of YBus matrix. These can then be used in subsequent
 * calculations
 */
gridpack::ComplexType gridpack::dynamic_simulation::DSBus::getYBus(void)
{
  return YMBus::getYBus();
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
  YMBus::load(data);

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
  p_ngen = 0;
  if (data->getValue(GENERATOR_NUMBER, &p_ngen)) {
    for (i=0; i<p_ngen; i++) {
      lgen = true;
      lgen = lgen && data->getValue(GENERATOR_PG, &pg,i);
      lgen = lgen && data->getValue(GENERATOR_QG, &qg,i);
      lgen = lgen && data->getValue(GENERATOR_STAT, &gstatus,i);
      pg /= p_sbase;
      qg /= p_sbase;

      lgen = lgen && data->getValue(GENERATOR_MBASE, &mva, i); 
      if (!data->getValue(GENERATOR_RESISTANCE, &r, i)) r=0.0; // r
      if (!data->getValue(GENERATOR_SUBTRANSIENT_REACTANCE, &dstr,i)) dstr = 0.0; // dstr
      gridpack::ComplexType zsrc;
      //lgen = lgen && data->getValue(GENERATOR_TRANSIENT_REACTANCE, &dtr,i); // dtr
      lgen = lgen && data->getValue(GENERATOR_ZSOURCE, &zsrc,i); // dtr
      dtr = imag(zsrc);
      // SJin: need to be added to parser
      lgen = lgen && data->getValue(GENERATOR_INERTIA_CONSTANT_H, &h, i); // h
      lgen = lgen && data->getValue(GENERATOR_DAMPING_COEFFICIENT_0, &d0, i); // d0
      if (lgen) {
        p_pg.push_back(pg);
        p_qg.push_back(qg);
        p_gstatus.push_back(gstatus);

        p_mva.push_back(mva);
        p_r.push_back(r);
        p_dstr.push_back(dstr);
        p_dtr.push_back(dtr);

        p_h.push_back(h);
        p_d0.push_back(d0);
        std::string id("-1");
        bool ok = data->getValue(GENERATOR_ID,&id,i);
        p_genid.push_back(id);
        p_watch.push_back(false);
      }
    }
  }

  /*// assume switch info is set up here instead of reading from the input file
  sw2_2 = 5; // 6-1
  sw3_2 = 6; // 7-1*/
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::dynamic_simulation::DSBus::setMode(int mode)
{
  if (mode == YBUS || mode == YL) {
    YMBus::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Return the value of the voltage magnitude on this bus
 * @return: voltage magnitude
 */
double gridpack::dynamic_simulation::DSBus::getVoltage(void)
{
  return p_voltage;
}

/**
 * Return the value of the phase angle on this bus
 * @return: phase angle
 */
double gridpack::dynamic_simulation::DSBus::getPhase(void)
{
  return p_angle;
}

/**
 * Return whether or not a bus is isolated
 * @return true if bus is isolated
 */
bool gridpack::dynamic_simulation::DSBus::isIsolated(void) const
{
  return YMBus::isIsolated();
}

/**
 * Return the number of generators on this bus
 * @return number of generators on bus
 */
int gridpack::dynamic_simulation::DSBus::getNumGen(void)
{
  return p_ngen;
}

void gridpack::dynamic_simulation::DSBus::setIFunc(void)
{
  if (p_ngen > 0) {
  } 
}

void gridpack::dynamic_simulation::DSBus::setIJaco(void)
{
}

/**
 * Check to see if a fault event applies to this bus and set an internal
 * flag marking the bus as the "from" or "to" bus for the event
 * @param from_idx index of "from" bus for fault event
 * @param to_idx index of "to" bus for fault event
 */
void gridpack::dynamic_simulation::DSBus::setEvent(int from_idx, int to_idx,
  gridpack::component::BaseBranchComponent* branch_ptr)
{
  if (from_idx == getOriginalIndex()) {
    p_from_flag = true;
  } else {
    p_from_flag = false;
  }
  if (to_idx == getOriginalIndex()) {
    p_to_flag = true;
  } else {
    p_to_flag = false;
  }
  if (p_to_flag || p_from_flag) {
    p_branch = branch_ptr;
  } else {
    p_branch = NULL;
  }
}

/**
 * Clear fault event from bus
 */
void gridpack::dynamic_simulation::DSBus::clearEvent()
{
  p_from_flag = false;
  p_to_flag = false;
  p_branch = NULL;
}

/**
 * Write output from buses to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool gridpack::dynamic_simulation::DSBus::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  if (signal == NULL) {
    if (p_ngen == 0) return false;
    int i;
    char buf[128];
    char *ptr = string;
    int idx = getOriginalIndex();
    int len = 0;
    for (i=0; i<p_ngen; i++) {
      sprintf(buf,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f\n",
          idx,p_genid[i].c_str(),real(p_mac_ang_final[i]),real(p_mac_spd_final[i]),
          real(p_mech_final[i]),real(p_elect_final[i]));
      int slen = strlen(buf);
      len += slen;
      if (len < bufsize) sprintf(ptr,"%s",buf);
      ptr += slen;
    }
    return true;
  } else if (!strcmp(signal,"watch_header")) {
    if (p_ngen == 0) return false;
    int i;
    char buf[128];
    char *ptr = string;
    int idx = getOriginalIndex();
    int len = 0;
    for (i=0; i<p_ngen; i++) {
      std::string tstr;
      if ((p_genid[i])[0] == ' ') {
        tstr = (p_genid[i])[1];
      } else {
        tstr = p_genid[i];
      }
      if (p_watch[i]) {
        sprintf(buf,", %d_%s_angle, %d_%s_speed",
            idx,tstr.c_str(),idx,tstr.c_str());
        int slen = strlen(buf);
        if (len+slen < bufsize) sprintf(ptr,"%s",buf);
        len += slen;
        ptr += slen;
      }
    }
    if (len > 0) return true;
  } else if (!strcmp(signal,"watch")) {
    if (p_ngen == 0) return false;
    int i;
    char buf[128];
    char *ptr = string;
    int len = 0;
    for (i=0; i<p_ngen; i++) {
      if (p_watch[i]) {
        sprintf(buf,", %f, %f",real(p_mac_ang_final[i]),
            real(p_mac_spd_final[i]));
        int slen = strlen(buf);
        if (len + slen < bufsize) sprintf(ptr,"%s",buf);
        len += slen;
        ptr += slen;
      }
    }
    if (len > 0) return true;
  }
  return false;
}

/**
 * Set an internal parameter that specifies that the rotor speed and angle
 * for the generator corresponding to the string tag are to be printed to
 * output
 * @param tag 2-character identifier of generator
 * @param flag set to true to monitor generator
 */
void gridpack::dynamic_simulation::DSBus::setWatch(std::string tag, bool flag)
{
  int i;
  for (i=0; i<p_genid.size(); i++) {
    if (tag == p_genid[i]) {
      p_watch[i] = flag;
      break;
    }
  }
}

/**
 *  Simple constructor
 */
gridpack::dynamic_simulation::DSBranch::DSBranch(void)
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
  p_event = false;
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
  if (p_mode == YBUS || p_mode == YL || p_mode == updateYbus ||
      p_mode == onFY || p_mode == posFY) {
    /*if (p_active) {
      *isize = 1;
      *jsize = 1;
      return true;
    } else {
      return false;
    }*/
    return YMBranch::matrixForwardSize(isize,jsize);
  } else {
    return false;
  }
  return false;
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
  if (p_mode == YBUS || p_mode == YL || p_mode == updateYbus ||
      p_mode == onFY || p_mode == posFY) {
    /*if (p_active) {
      *isize = 1;
      *jsize = 1;
      return true;
    } else {
      return false;
    }*/
    return YMBranch::matrixReverseSize(isize,jsize);
  } else {
    return false;
  }
  return false;
}

/**
 * Return the values of the off-diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::dynamic_simulation::DSBranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == YBUS || p_mode == YL || p_mode == updateYbus) {
    /*if (p_active) {
      values[0] = gridpack::ComplexType(p_ybusr_frwd,p_ybusi_frwd); 
      return true;
    } else {
      return false;
    }*/
    return YMBranch::matrixForwardValues(values);
  /*} else if (p_mode == PERM) {
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
    return true;*/
  } else if (p_mode == posFY) {
    if (p_event) {
      values[0] = -getUpdateFactor();
      return true;
    } else {
      return false;
    }
  } else {
    return false;
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
  return false;
}
bool gridpack::dynamic_simulation::DSBranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == YBUS || p_mode == YL || p_mode == updateYbus) {
    /*if (p_active) {
      values[0] = gridpack::ComplexType(p_ybusr_rvrs,p_ybusi_rvrs);
      return true;
      return true;
    } else {
      return false;
    }*/
    return YMBranch::matrixForwardValues(values);
  } else if (p_mode == posFY) {
    if (p_event) {
      values[0] = -getUpdateFactor();
      return true;
    } else {
      return false;
    }
  /*} else if (p_mode == PERM) {
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
    return true;*/
  } else {
    return false;
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
  return false;
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::dynamic_simulation::DSBranch::setYBus(void)
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
  gridpack::dynamic_simulation::DSBus *bus1 =
    dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus1().get());
  gridpack::dynamic_simulation::DSBus *bus2 =
    dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus2().get());
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
 * Load values stored in DataCollection object into DSBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::dynamic_simulation::DSBranch::load(
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
void gridpack::dynamic_simulation::DSBranch::setMode(int mode)
{
  if (mode == YBUS || mode == YL) {
    YMBranch::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Return the complex admittance of the branch
 * @return: complex addmittance of branch
 */
gridpack::ComplexType gridpack::dynamic_simulation::DSBranch::getAdmittance(void)
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
gridpack::dynamic_simulation::DSBranch::getTransformer(gridpack::dynamic_simulation::DSBus *bus)
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
gridpack::dynamic_simulation::DSBranch::getShunt(gridpack::dynamic_simulation::DSBus *bus)
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
      } else if (bus == getBus2().get()) {        tmpr += p_shunt_admt_g2[i];        tmpi += p_shunt_admt_b2[i];
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
 * Return the updating factor that will be applied to the ybus matrix at
 * the clear fault phase
 * @return: value of update factor
 */
gridpack::ComplexType 
gridpack::dynamic_simulation::DSBranch::getPosfy11YbusUpdateFactor(int sw2_2, int sw3_2)
{ 
  double retr, reti;
  int i;
  gridpack::dynamic_simulation::DSBus *bus1 =
    dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus1().get());
  gridpack::dynamic_simulation::DSBus *bus2 =
    dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus2().get());
  if (bus1->getOriginalIndex() == sw2_2+1 && bus2->getOriginalIndex() == sw3_2+1) {
    for (i=0; i<p_elems; i++) {
      gridpack::ComplexType myValue(p_resistance[i], p_reactance[i]);
      myValue = 1.0 / myValue;
      //printf("myValue = %f+%fi\n", real(myValue), imag(myValue));
      //printf("%f %f\n", p_resistance, p_reactance);
      retr = real(myValue);
      reti = imag(myValue);
      return gridpack::ComplexType(retr, reti);
    }
  } else {
    return gridpack::ComplexType(-999.0, -999.0); // return a dummy value
  }
  return gridpack::ComplexType(0.0,0.0);
/*  double retr, reti;
  gridpack::dynamic_simulation::DSBus *bus1 =
    dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus1().get());
  gridpack::dynamic_simulation::DSBus *bus2 =
    dynamic_cast<gridpack::dynamic_simulation::DSBus*>(getBus2().get());
  if (bus1->getOriginalIndex() == sw2_2+1 && bus2->getOriginalIndex() == sw3_2+1) {
    gridpack::ComplexType myValue(p_resistance, p_reactance);
    myValue = 1.0 / myValue;
    //printf("myValue = %f+%fi\n", real(myValue), imag(myValue));
    //printf("%f %f\n", p_resistance, p_reactance);
    retr = real(myValue);
    reti = imag(myValue);
    return gridpack::ComplexType(retr, reti);
  } else {
    return gridpack::ComplexType(-999.0, -999.0); // return a dummy value
  }
*/
}
gridpack::ComplexType 
gridpack::dynamic_simulation::DSBranch::getUpdateFactor()
{ 
  int i;
  gridpack::ComplexType ret(0.0,0.0);
  for (i=0; i<p_elems; i++) {
    gridpack::ComplexType tmp(p_resistance[i], p_reactance[i]);
    tmp = -1.0 / tmp;
    ret += tmp;
  }
  return ret;
}

/**
 * Check to see if an event applies to this branch and set appropriate internal
 * parameters
 * @param event a struct containing parameters that describe a fault event in
 * a dyanamic simulation
 */
void gridpack::dynamic_simulation::DSBranch::setEvent(const Event &event)
{
  int idx1 = getBus1OriginalIndex();
  int idx2 = getBus2OriginalIndex();
  // Check to see if event refers to this bus
  if (idx1 == event.from_idx && idx2 == event.to_idx) {
    p_event = true;
  } else {
    p_event = false;
  }
  if (p_event) {
    dynamic_cast<gridpack::dynamic_simulation::DSBus*>
      (getBus1().get())->setEvent(idx1,idx2,this);
    dynamic_cast<gridpack::dynamic_simulation::DSBus*>
      (getBus2().get())->setEvent(idx1,idx2,this);
  }
}
