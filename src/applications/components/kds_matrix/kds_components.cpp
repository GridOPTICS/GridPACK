/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   KalmanDS_components.cpp
 * @author Da Meng and Yousu Chen
 * @date   2018-07-17 11:36:24 d3g096
 * 
 * @brief  
 * 
 * Modified by Xinya Li, July 2015
 * Ensemble Kalman Filter to Dynamic Simulation
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>
#include <stdio.h>
#include <string.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parallel/random.hpp"
#include "gridpack/applications/components/kds_matrix/kds_components.hpp"
#include "gridpack/parser/dictionary.hpp"

//#define LARGE_MATRIX

/**
 *  Simple constructor
 */
gridpack::kalman_filter::KalmanBus::KalmanBus(void)
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
  p_mag_series = NULL;
  p_ang_series = NULL;
  p_delta1 = NULL;
  p_omega1 = NULL;
  p_delta_diff = NULL;
  p_omega_diff = NULL;
  p_V1 = NULL;
  p_V2 = NULL;
  p_V3 = NULL;
  p_delta1_dt = NULL;
  p_omega1_dt = NULL;
  p_delta2 = NULL;
  p_omega2 = NULL;
  p_delta2_dt = NULL;
  p_omega2_dt = NULL;
  p_delta3 = NULL;
  p_omega3 = NULL;
  p_delta_t = 0.0;
  p_t_len = 0;
  setReferenceBus(false);
  p_from_flag = false;
  p_to_flag = false;
  p_branch = NULL;
}

/**
 *  Simple destructor
 */
gridpack::kalman_filter::KalmanBus::~KalmanBus(void)
{
  if (p_mag_series != NULL) delete [] p_mag_series;
  if (p_ang_series != NULL) delete [] p_ang_series;
  if (p_ngen > 0) {
    if (p_delta1 != NULL) {
      int i;
      for (i=0; i<p_ngen; i++) {
        delete [] p_delta1[i];
        delete [] p_omega1[i];
        delete [] p_delta_diff[i];
        delete [] p_omega_diff[i];
      }
      delete [] p_delta1;
      delete [] p_omega1;
      delete [] p_delta_diff;
      delete [] p_omega_diff;
    }
    if (p_V1 != NULL) {
      int i;
      for (i=0; i<p_ngen; i++) {
        delete [] p_V1[i];
      }
      delete [] p_V1;
    }
    if (p_V2 != NULL) {
      int i;
      for (i=0; i<p_ngen; i++) {
        delete [] p_V2[i];
      }
      delete [] p_V2;
    }
    if (p_V3 != NULL) {
      delete [] p_V3;
    }
    if (p_delta1_dt != NULL) {
      int i;
      for (i=0; i<p_ngen; i++) {
        delete [] p_delta1_dt[i];
      }
      delete [] p_delta1_dt;
    }
    if (p_omega1_dt != NULL) {
      int i;
      for (i=0; i<p_ngen; i++) {
        delete [] p_omega1_dt[i];
      }
      delete [] p_omega1_dt;
    }
    if (p_delta2_dt != NULL) {
      int i;
      for (i=0; i<p_ngen; i++) {
        delete [] p_delta2_dt[i];
      }
      delete [] p_delta2_dt;
    }
    if (p_omega2_dt != NULL) {
      int i;
      for (i=0; i<p_ngen; i++) {
        delete [] p_omega2_dt[i];
      }
      delete [] p_omega2_dt;
    }
  }
}

/**
 *  Return size of matrix block contributed by the component
 *  @param isize, jsize: number of rows and columns of matrix block
 *  @return: false if network component does not contribute matrix element
 */
bool gridpack::kalman_filter::KalmanBus::matrixDiagSize(int *isize, int *jsize) const
{
  if (p_mode == YBus || p_mode == RefreshY || p_mode == onFY || p_mode == posFY) {
    return YMBus::matrixDiagSize(isize,jsize);
  } else if (p_mode == YC) {
    if (p_ngen > 0) {
      *isize = 1;
      *jsize = p_ngen;
    } else {
      *isize = 1;
      *jsize = 0;
    }
  }
  return true;
}

/**
 * Return the values of the matrix block. The values are
 * returned in row-major order.
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::kalman_filter::KalmanBus::matrixDiagValues(ComplexType *values)
{
  if (p_mode == YBus) {
    return YMBus::matrixDiagValues(values);
  } else if (p_mode == onFY){
    if (p_from_flag) {
      gridpack::ComplexType ret(1.0e6,0.0);
      values[0] = ret;
      return true;
    } else {
      return false;
    }
  } else if (p_mode == posFY) {
    if (p_from_flag || p_to_flag) {
      values[0] = dynamic_cast<gridpack::kalman_filter::KalmanBranch*>(p_branch)
          ->getUpdateFactor();
      return true;
    } else {
      return false;
    }
  } else if (p_mode == YC) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
        double ra = p_r[i] * p_sbase / p_mva[i];
        double xd = 0.0;
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
  } else if (p_mode == RefreshY) {
    double ur = p_ybusr;
    double ui = p_ybusi;
    if (p_ngen > 0) {
      gridpack::ComplexType permYmod = 0.0;
      for (int ip = 0; ip < p_ngen; ip++) {
        double ra = p_r[ip] * p_sbase / p_mva[ip];
        double xd=0.0;
        if (p_dstr[ip] == 0) {
          xd = p_dtr[ip] * p_sbase / p_mva[ip];
        }
        gridpack::ComplexType Y_a(ra, xd);
        Y_a = 1.0 / Y_a;
        permYmod += Y_a;
      }
      ur += real(permYmod);
      ui += imag(permYmod);
    }
    ur += p_pl/(p_v*p_v*p_sbase);
    ui += (-p_ql)/(p_v*p_v*p_sbase);
    gridpack::ComplexType u(ur, ui);
    values[0] = u;
    return true;
  } else {
    return false;
  }
}

/**
 * Set the number of ensemble samples on the buses
 * @param nsize number of ensembles
 */
void gridpack::kalman_filter::KalmanBus::setEnsembleSize(int nsize)
{
  p_nEnsemble = nsize;
}

/**
 * Set the distribution width
 * @param sigma width of guassian distribution
 */
void gridpack::kalman_filter::KalmanBus::setGaussianWidth(double sigma,
    double noise)
{
  p_sigma = sigma;
  p_noise = noise;
}

/**
 * Return number of rows and columns in matrix from component
 * Number of columns must be the same for all components
 * @return size of block contributed by component
 */
void gridpack::kalman_filter::KalmanBus::slabSize(int *rows,
    int *cols) const
{
  if (p_mode == EnsembleX || p_mode == Perturbation || p_mode == X_INC || p_mode == X_Update) {
    *rows = 2*p_ngen;
    *cols = p_nEnsemble;
  } else if (p_mode == V3) {
    *rows = 1;
    *cols = p_nEnsemble;
  } else if (p_mode == HX || p_mode == HA || p_mode == Measurements) {
    *rows = 2;
    *cols = p_nEnsemble;
  } else if (p_mode == V1 || p_mode == V2) {
//    if (p_ngen>0) {
      *rows = 1;
//      *rows = p_ngen;
//    } else {
//      *rows = 0;
//    }
    *cols = p_nEnsemble;
  } else if (p_mode == E_Ensemble1 || p_mode == E_Ensemble2 ||
      p_mode == E_Ensemble3) {
    *rows = p_ngen;
    *cols = p_nEnsemble;
  } else {
    *rows = 0;
    *cols = 0;
  }
}

/**
 * Set indices corresponding to the rows contributed by this
 * component
 * @param irow index of row contributed by this component
 * @param idx row index of row irow
 */
void gridpack::kalman_filter::KalmanBus::slabSetRowIndex(int irow,
    int idx)
{
  if (p_mode == EnsembleX || p_mode == Perturbation || p_mode == X_INC || p_mode == X_Update) {
    p_X_idx.push_back(idx);
  } else if (p_mode == V3) {
    p_V3_idx.push_back(idx);
  } else if (p_mode == HX || p_mode == HA || p_mode == Measurements) {
    p_HX_idx.push_back(idx);
  } else if (p_mode == E_Ensemble1 || p_mode == E_Ensemble2 ||
      p_mode == E_Ensemble3) {
    p_E_idx.push_back(idx);
  } else if (p_mode == V1 || p_mode == V2) {
    p_V_idx.push_back(idx);
  }
}

/**
 * Get list of row indices from component
 * @param idx list of row indices that component maps onto
 */
void gridpack::kalman_filter::KalmanBus::slabGetRowIndices(int *idx)
{
  if (p_mode == EnsembleX || p_mode == Perturbation || p_mode == X_INC || p_mode == X_Update) {
    int i;
    for (i=0; i<p_X_idx.size(); i++) {
      idx[i] = p_X_idx[i];
    }
  } else if (p_mode == V3) {
    int i;
    for (i=0; i<p_V3_idx.size(); i++) {
      idx[i] = p_V3_idx[i];
    }
  } else if (p_mode == HX || p_mode == HA || p_mode == Measurements) {
    int i;
    for (i=0; i<p_HX_idx.size(); i++) {
      idx[i] = p_HX_idx[i];
    }
  } else if (p_mode == E_Ensemble1 || p_mode == E_Ensemble2 ||
      p_mode == E_Ensemble3) {
    int i;
    for (i=0; i<p_E_idx.size(); i++) {
      idx[i] = p_E_idx[i];
    }
  } else if (p_mode == V1 || p_mode == V2) {
    int i;
    for (i=0; i<p_V_idx.size(); i++) {
      idx[i] = p_V_idx[i];
    }
  }
}

/**
 * Get a list of row values contributed by this component and their
 * indices
 * @param values list of values for rows
 * @param idx indices for the matrix rows
 */
void gridpack::kalman_filter::KalmanBus::slabGetValues(
    std::vector<ComplexType*> &values, int *idx)
{
  int i, k;
  gridpack::random::Random random;
  if (p_mode == EnsembleX) {
    if (p_ngen > 0) {
      for (i=0; i<p_ngen; i++) {
        idx[2*i] = p_X_idx[2*i];
        idx[2*i+1] = p_X_idx[2*i+1];
        for (k=0; k<p_nEnsemble; k++) {
          (values[2*i])[k] = (p_delta1[i])[k];
        }
        for (k=0; k<p_nEnsemble; k++) {
          (values[2*i+1])[k] = (p_omega1[i])[k];
        }
      }
    }
  } else if (p_mode == Perturbation) {
    if (p_ngen > 0) {
      for (i=0; i<p_ngen; i++) {
        idx[2*i] = p_X_idx[2*i];
        idx[2*i+1] = p_X_idx[2*i+1];
        for (k=0; k<p_nEnsemble; k++) {
          (values[2*i])[k] = (p_delta_diff[i])[k];
        }
        for (k=0; k<p_nEnsemble; k++) {
          (values[2*i+1])[k] = (p_omega_diff[i])[k];
        }
      }
    }
  } else if (p_mode == Measurements) {
    idx[0] = p_HX_idx[0];
    idx[1] = p_HX_idx[1];
    double d1 = p_ang_series[p_currentStep];
    double d2 = p_mag_series[p_currentStep];
    for (k=0; k<p_nEnsemble; k++) {
      (values[0])[k] = d1+p_noise*random.grand();
    }
    for (k=0; k<p_nEnsemble; k++) {
      (values[1])[k] = d2+p_noise*random.grand();
    }
  } else if (p_mode == E_Ensemble1) {
    if (p_ngen > 0) {
      for (i=0; i<p_ngen; i++) {
        idx[i] = p_E_idx[i];
        for (k=0; k<p_nEnsemble; k++) {
          (values[i])[k] = p_Ebase[i]*gridpack::ComplexType(cos((p_delta1[i])[k]),
              sin((p_delta1[i])[k]));
        }
      }
    }
  } else if (p_mode == E_Ensemble2) {
    if (p_ngen > 0) {
      for (i=0; i<p_ngen; i++) {
        idx[i] = p_E_idx[i];
        for (k=0; k<p_nEnsemble; k++) {
          (values[i])[k] = p_Ebase[i]*gridpack::ComplexType(cos((p_delta2[i])[k]),
              sin((p_delta2[i])[k]));
        }
      }
    }
  } else if (p_mode == E_Ensemble3) {
    for (i=0; i<p_ngen; i++) {
      idx[i] = p_E_idx[i];
      for (k=0; k<p_nEnsemble; k++) {
        (values[i])[k] = p_Ebase[i]*gridpack::ComplexType(cos((p_delta3[i])[k]),
            sin((p_delta3[i])[k]));
      }
    }
  } else if (p_mode == HX) {
    idx[0] = p_HX_idx[0];
    idx[1] = p_HX_idx[1];
    for (k=0; k<p_nEnsemble; k++) {
      (values[0])[k] = arg(p_V3[k]);
      (values[1])[k] = abs(p_V3[k]);
    }
  } else if (p_mode == HA) {
    // Get averages
    double ang = 0.0;
    double mag = 0.0;
    double n = static_cast<double>(p_nEnsemble);
    for (k=0; k<p_nEnsemble; k++) {
      ang += arg(p_V3[k]);
      mag += abs(p_V3[k]);
    }
    ang = ang/n;
    mag = mag/n;
    idx[0] = p_HX_idx[0];
    idx[1] = p_HX_idx[1];
    for (k=0; k<p_nEnsemble; k++) {
      (values[0])[k] = arg(p_V3[k])-ang;
      (values[1])[k] = abs(p_V3[k])-mag;
    }
  }  
}

/**
 * Transfer slab values to component
 * @param values list of slab values
 */
void gridpack::kalman_filter::KalmanBus::slabSetValues(ComplexType **values)
{
  if (p_mode == V1) {
    int i, k;
    if (p_ngen > 0) {
      for (i=0; i<p_ngen; i++) {
        for (k=0; k<p_nEnsemble; k++) {
          (p_V1[i])[k] = (values[i])[k];
        }
      }
    } 
  } else if (p_mode == V2) {
    int i, k;
    if (p_ngen > 0) {
      for (i=0; i<p_ngen; i++) {
        for (k=0; k<p_nEnsemble; k++) {
          (p_V2[i])[k] = (values[i])[k];
        }
      }
    }
  } else if (p_mode == V3) {
    int k;
    for (k=0; k<p_nEnsemble; k++) {
      p_V3[k] = (values[0])[k];
    }
  } else if (p_mode == X_INC) {
    int i, k;
    if (p_ngen > 0) {
      for (i=0; i<p_ngen; i++) {
        for (k=0; k<p_nEnsemble; k++) {
          (p_delta1[i])[k] = (p_delta3[i])[k] +  real((values[2*i])[k]);
        }
        for (k=0; k<p_nEnsemble; k++) {
          (p_omega1[i])[k] = (p_omega3[i])[k] + real((values[2*i+1])[k]);
        }
      }
    }
  } else if (p_mode == X_Update) {
    int i, k;
    if (p_ngen > 0) {
      for (i=0; i<p_ngen; i++) {
        for (k=0; k<p_nEnsemble; k++) {
          (p_delta1[i])[k] = (p_delta3[i])[k];
        }
        for (k=0; k<p_nEnsemble; k++) {
          (p_omega1[i])[k] = (p_omega3[i])[k];
        }
      }
    } 
  }
}

/**
 * Load values stored in DataCollection object into KalmanBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::kalman_filter::KalmanBus::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  YMBus::load(data);

  p_sbase = 100.0;
  bool ok = data->getValue(CASE_SBASE, &p_sbase);
  if (!data->getValue("BUS_PF_VANG", &p_angle))
    data->getValue("BUS_VOLTAGE_ANG",&p_angle);
  if (!data->getValue("BUS_PF_VMAG", &p_voltage))
    data->getValue("BUS_VOLTAGE_MAG", &p_voltage); 
//  data->getValue(BUS_VOLTAGE_ANG, &p_angle);
//  data->getValue(BUS_VOLTAGE_MAG, &p_voltage); 
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
  bool lgen, dgen;
  int i, ngen, gstatus;
  double pg, qg, vs, r, dstr, dtr;
  double h, d0;
  ngen = 0;
  if (data->getValue(GENERATOR_NUMBER, &ngen)) {
    for (i=0; i<ngen; i++) {
      lgen = true;
      lgen = lgen && data->getValue(GENERATOR_PG, &pg,i);
      lgen = lgen && data->getValue(GENERATOR_QG, &qg,i);
      lgen = lgen && data->getValue(GENERATOR_VS, &vs,i);
      lgen = lgen && data->getValue(GENERATOR_STAT, &gstatus,i);
      pg /= p_sbase;
      qg /= p_sbase;

      double genMVA;
      data->getValue(GENERATOR_MBASE, &genMVA, i);
      if (!data->getValue(GENERATOR_RESISTANCE, &r, i)) r=0.0; // r
      if (!data->getValue(GENERATOR_SUBTRANSIENT_REACTANCE, &dstr,i)) dstr = 0.0; // dstr

      gridpack::ComplexType zsource;
      data->getValue(GENERATOR_ZSOURCE, &zsource, i);
      double xd = imag(zsource);
      double pgen=0.0;
      if (!data->getValue("GENERATOR_PF_PGEN",&pgen,i)) {
        data->getValue(GENERATOR_PG,&pgen,i);
      }
      double qgen=0.0;
      if (!data->getValue("GENERATOR_PF_QGEN",&qgen,i)) {
        data->getValue(GENERATOR_QG,&qgen,i);
      }
      gridpack::ComplexType s = gridpack::ComplexType(pgen,qgen);
      gridpack::ComplexType ii = gridpack::ComplexType(0.0,1.0);
      gridpack::ComplexType vterm = p_v*(cos(p_a)+ii*sin(p_a));
      gridpack::ComplexType iterm = conj(s/vterm)*p_sbase/genMVA;
      gridpack::ComplexType Egen0 = vterm+iterm*ii*xd;
      std::string genID;
      data->getValue("GENERATOR_ID",&genID,i);
      double Emag = abs(Egen0);
      double delta = arg(Egen0);
      dgen = true;
      double dm=0.0;
      double hmach = 0.0;
      dgen = dgen && data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&dm,i);
      dgen = dgen && data->getValue(GENERATOR_INERTIA_CONSTANT_H,&hmach,i);
      double pm = pgen*p_sbase/genMVA;
      if (lgen && dgen) {
        p_omega_0.push_back(1.0);
        p_omega_n.push_back(1.0);
        p_delta_0.push_back(delta);
        p_delta_n.push_back(delta);
        p_Dm.push_back(dm);
        p_Pm.push_back(pm);
//        p_pg.push_back(pg);
//        p_qg.push_back(qg);
        p_pg.push_back(pgen);
        p_qg.push_back(qgen);
        p_hmach.push_back(hmach);
        p_gstatus.push_back(gstatus);
        p_mva.push_back(genMVA);
        p_r.push_back(r);
        p_dstr.push_back(dstr);
        p_dtr.push_back(xd);
        p_Ebase.push_back(Emag);
        p_genID.push_back(genID);
        if (gstatus == 1) {
          p_v = vs; //reset initial PV voltage to set voltage
          if (itype == 2) p_isPV = true;
        }
      } else if (lgen) {
        p_pl -= pgen*p_sbase;
        p_ql -= qgen*p_sbase;          
      }
    }
  }
  p_ngen = p_pg.size();

}

/**
 * Set values of YBus matrix. These can then be used in subsequent
 * calculations
 */
void gridpack::kalman_filter::KalmanBus::setYBus(void)
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
gridpack::ComplexType gridpack::kalman_filter::KalmanBus::getYBus(void)
{
  return YMBus::getYBus();
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::kalman_filter::KalmanBus::setMode(int mode)
{
  if (mode == YBus) {
    YMBus::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Return whether or not a bus is isolated
 * @return true if bus is isolated
 */
bool gridpack::kalman_filter::KalmanBus::isIsolated(void) const
{
  return YMBus::isIsolated();
}

/**
 * Check to see if a fault event applies to this bus and set an internal
 * flag marking the bus as the "from" or "to" bus for the event
 * @param from_idx index of "from" bus for fault event
 * @param to_idx index of "to" bus for fault event
 */
void gridpack::kalman_filter::KalmanBus::setEvent(int from_idx, int to_idx,
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
void gridpack::kalman_filter::KalmanBus::clearEvent()
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
bool gridpack::kalman_filter::KalmanBus::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  char tmp[128];
  if (!strcmp(signal,"xnew")) {
    if (p_ngen == 0) return false;
    double pi = 4.0*atan(1.0);
    int i, k;
    char *ptr = string;
    int slen = 0;
    for (i=0; i<p_ngen; i++) {
      double delta = 0.0;
      double omega = 0.0;
      for (k=0; k<p_nEnsemble; k++) {
        delta += (p_delta1[i])[k];
        omega += (p_omega1[i])[k];
      }
      delta = delta/static_cast<double>(p_nEnsemble);
      omega = omega/static_cast<double>(p_nEnsemble);
      sprintf(tmp,"  %8d      %s  %16.8f  %16.8f\n",getOriginalIndex(),
          p_genID[i].c_str(),delta,omega);
      slen += strlen(tmp);
      if (slen < bufsize) {
        sprintf(ptr,"%s",tmp);
        ptr += strlen(tmp);
      }
    }
  } else if (!strcmp(signal,"delta")) {
    if (p_ngen == 0) return false;
    double pi = 4.0*atan(1.0);
    int i, k;
    char *ptr = string;
    int slen = 0;
    for (i=0; i<p_ngen; i++) {
      double delta = 0.0;
      for (k=0; k<p_nEnsemble; k++) {
        delta += (p_delta1[i])[k];
      }
      delta = delta/static_cast<double>(p_nEnsemble);
      sprintf(tmp," %16.8e",delta);
      slen += strlen(tmp);
      if (slen < bufsize) {
        sprintf(ptr,"%s",tmp);
        ptr += strlen(tmp);
      }
    }
  } else if (!strcmp(signal,"omega")) {
    if (p_ngen == 0) return false;
    double pi = 4.0*atan(1.0);
    int i, k;
    char *ptr = string;
    int slen = 0;
    for (i=0; i<p_ngen; i++) {
      double omega = 0.0;
      for (k=0; k<p_nEnsemble; k++) {
        omega += (p_omega1[i])[k];
      }
      omega = omega/static_cast<double>(p_nEnsemble);
      sprintf(tmp," %16.8e",omega);
      slen += strlen(tmp);
      if (slen < bufsize) {
        sprintf(ptr,"%s",tmp);
        ptr += strlen(tmp);
      }
    }
  } else {
    return false;
  }
  return true;
}

/**
 * Set voltage angle time series data
 * @param ang array containing time series
 */
void gridpack::kalman_filter::KalmanBus::setVAngSeries(double *ang)
{
  p_ang_series = ang;
}

/**
 * Set voltage magnitude time series data
 * @param mag array containing time series
 */
void gridpack::kalman_filter::KalmanBus::setVMagSeries(double *mag)
{
  p_mag_series = mag;
}

/**
 * Set time increment and number of timesteps
 * @param delta_t time step increment
 * @param nsteps number of steps in time series
 */
void gridpack::kalman_filter::KalmanBus::setTimeSteps(double delta_t, int nsteps)
{
  p_delta_t = delta_t;
  p_t_len = nsteps;
}

void gridpack::kalman_filter::KalmanBus::printT()
{
  printf("Bus %d      (mag)      (ang)\n",getOriginalIndex());
  int i;
  for (i=0; i<p_t_len; i++) {
    printf("        %16.8f  %16.8f\n",p_mag_series[i],p_ang_series[i]);
  }
}

/**
 * Set current step
 * @param istep current step
 */
void gridpack::kalman_filter::KalmanBus::setCurrentTimeStep(int istep)
{
  p_currentStep = istep;
}

/**
 * Create the ensemble of rotor phases and rotor speeds 
 */
void gridpack::kalman_filter::KalmanBus::createEnsemble()
{
  int i, k;
  gridpack::random::Random random;
  if (p_V3 == NULL) {
    p_V3 = new gridpack::ComplexType[p_nEnsemble];
  }
  if (p_ngen > 0) {
    if (p_delta1 == NULL) {
      p_delta1 = new double*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_delta1[i] = new double[p_nEnsemble];
      }
    }
    if (p_omega1 == NULL) {
      p_omega1 = new double*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_omega1[i] = new double[p_nEnsemble];
      }
    }
    if (p_delta_diff == NULL) {
      p_delta_diff = new double*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_delta_diff[i] = new double[p_nEnsemble];
      }
    }
    if (p_omega_diff == NULL) {
      p_omega_diff = new double*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_omega_diff[i] = new double[p_nEnsemble];
      }
    }
    if (p_V1 == NULL) {
      p_V1 = new gridpack::ComplexType*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_V1[i] = new gridpack::ComplexType[p_nEnsemble];
      }
    }
    if (p_V2 == NULL) {
      p_V2 = new gridpack::ComplexType*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_V2[i] = new gridpack::ComplexType[p_nEnsemble];
      }
    }
    for (i=0; i<p_ngen; i++) {
      double angT = 0.0;
      double magT = 0.0;
      for (k=0; k<p_nEnsemble; k++) {
        (p_delta1[i])[k] = p_delta_0[i]*(1.0+p_sigma*random.grand());
      }
      for (k=0; k<p_nEnsemble; k++) {
        (p_omega1[i])[k] = p_omega_0[i]*(1.0+p_sigma*random.grand());
      }
    }
  }
}


/**
 * Evaluate dX_dt_1 and then get X2
 */
void gridpack::kalman_filter::KalmanBus::evaluateX2()
{
  if (p_ngen > 0) {
    int i, k;
    // Create new buffers, if necessary
    if (p_delta1_dt == NULL) {
      p_delta1_dt = new double*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_delta1_dt[i] = new double[p_nEnsemble];
      }
    }
    if (p_omega1_dt == NULL) {
      p_omega1_dt = new double*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_omega1_dt[i] = new double[p_nEnsemble];
      }
    }
    if (p_delta2 == NULL) {
      p_delta2 = new double*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_delta2[i] = new double[p_nEnsemble];
      }
    }
    if (p_omega2 == NULL) {
      p_omega2 = new double*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_omega2[i] = new double[p_nEnsemble];
      }
    }
//    printf("p[%d] bus[%d] (X2) Got to 1\n",GA_Nodeid(),getOriginalIndex());
    double pi = 4.0*atan(1.0);
    double vmag, vang;
    // Evaluate time derivatives of delta1, omega1
    for (i=0; i<p_ngen; i++) {
//        printf("Pm: %f Dm: %f Emag: %f Xd: %f Hm: %f\n",
//            p_Pm[i],p_Dm[i],p_Ebase[i],p_dtr[i],p_hmach[i]);
//        printf("omega1: ");
//      printf("X1_dt:");
      for (k=0; k<p_nEnsemble; k++) {
        if (p_dtr[i] == 0.0 || p_hmach[i] == 0.0) {
          printf("Zero value found for Xd: %f or Hmach: %f\n",
              p_dtr[i],p_hmach[i]);
        }
        vmag = abs((p_V1[i])[k]);
        vang = arg((p_V1[i])[k]);
//        printf("vmag,vang (%f, %f)\n",vmag,vang);
        (p_delta1_dt[i])[k] = 120.0*pi*((p_omega1[i])[k]-1.0);
        (p_omega1_dt[i])[k] = (p_Pm[i]-p_Dm[i]*((p_omega1[i])[k]
              - 1.0) - p_Ebase[i]*vmag*sin((p_delta1[i])[k]-vang)/p_dtr[i])
          /(2.0*p_hmach[i]);
//        printf("p_1_dt (%f,%f)\n",(p_delta1_dt[i])[k],(p_omega1_dt[i])[k]);
      }
//      printf("\n");
    }
//    printf("p[%d] bus[%d] (X2) Got to 2\n",GA_Nodeid(),getOriginalIndex());
    // Evaluate delta2, omega2
    for (i=0; i<p_ngen; i++) {
//        printf("delta2: ");
      for (k=0; k<p_nEnsemble; k++) {
        (p_delta2[i])[k] = (p_delta1[i])[k] + (p_delta1_dt[i])[k]*p_delta_t;
//        printf("p_1_d (%f,%f)\n",(p_delta1[i])[k],(p_delta2[i])[k]);
      }
//      printf("\n");
//        printf("omega2: ");
      for (k=0; k<p_nEnsemble; k++) {
        (p_omega2[i])[k] = (p_omega1[i])[k] + (p_omega1_dt[i])[k]*p_delta_t;
//        printf("p_1_o (%f,%f)\n",(p_omega1[i])[k],(p_omega2[i])[k]);
      }
//      printf("\n");
    }
//    printf("p[%d] bus[%d] (X2) Got to 3\n",GA_Nodeid(),getOriginalIndex());
  }
}

/**
 * Evaluate dX_dt_2 and then get X3
 */
void gridpack::kalman_filter::KalmanBus::evaluateX3()
{
  if (p_ngen > 0) {
    int i, k;
    // Create new buffers, if necessary
    if (p_delta2_dt == NULL) {
      p_delta2_dt = new double*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_delta2_dt[i] = new double[p_nEnsemble];
      }
    }
    if (p_omega2_dt == NULL) {
      p_omega2_dt = new double*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_omega2_dt[i] = new double[p_nEnsemble];
      }
    }
    if (p_delta3 == NULL) {
      p_delta3 = new double*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_delta3[i] = new double[p_nEnsemble];
      }
    }
    if (p_omega3 == NULL) {
      p_omega3 = new double*[p_ngen];
      for (i=0; i<p_ngen; i++) {
        p_omega3[i] = new double[p_nEnsemble];
      }
    }
    double pi = 4.0*atan(1.0);
    double vmag, vang;
    // Evaluate time derivatives of delta3, omega3
    for (i=0; i<p_ngen; i++) {
//      printf("X2_dt:");
      for (k=0; k<p_nEnsemble; k++) {
        if (p_dtr[i] == 0.0 || p_hmach[i] == 0.0) {
          printf("Zero value found for Xd: %f or Hmach: %f\n",
              p_dtr[i],p_hmach[i]);
        }
        vmag = abs((p_V2[i])[k]);
        vang = arg((p_V2[i])[k]);
        (p_delta2_dt[i])[k] = 120.0*pi*((p_omega2[i])[k]-1.0);
        (p_omega2_dt[i])[k] = (p_Pm[i]-p_Dm[i]*((p_omega2[i])[k]
              - 1.0) - p_Ebase[i]*vmag*sin((p_delta2[i])[k]-vang)/p_dtr[i])
          /(2.0*p_hmach[i]);
//        printf("p_2_dt (%f,%f)\n",(p_delta2_dt[i])[k],(p_omega2_dt[i])[k]);
      }
//      printf("\n");
    }
    // Evaluate delta3, omega3
    double ddelta_dt, domega_dt;
    for (i=0; i<p_ngen; i++) {
      double angT = 0.0;
      double magT = 0.0;      
      for (k=0; k<p_nEnsemble; k++) {
        ddelta_dt = 0.5*((p_delta1_dt[i])[k]+(p_delta2_dt[i])[k]);
        domega_dt = 0.5*((p_omega1_dt[i])[k]+(p_omega2_dt[i])[k]);
        (p_delta3[i])[k] = (p_delta1[i])[k] + ddelta_dt*p_delta_t;        
        (p_omega3[i])[k] = (p_omega1[i])[k] + domega_dt*p_delta_t;
        angT += (p_delta3[i])[k];
        magT += (p_omega3[i])[k];
      }
      angT = angT/static_cast<double>(p_nEnsemble);
      magT = magT/static_cast<double>(p_nEnsemble);
      for (k=0; k<p_nEnsemble; k++) {
        (p_delta_diff[i])[k] = (p_delta3[i])[k]-angT;
        (p_omega_diff[i])[k] = (p_omega3[i])[k]-magT;
      }
    }
  }
}

/**
 * The number of generators on this bus
 * @return number of generators
 */
int gridpack::kalman_filter::KalmanBus::numGenerators()
{
  return p_ngen;
}

/**
 *  Simple constructor
 */
gridpack::kalman_filter::KalmanBranch::KalmanBranch(void)
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
 *  Simple destructor
 */
gridpack::kalman_filter::KalmanBranch::~KalmanBranch(void)
{
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::kalman_filter::KalmanBranch::matrixForwardSize(int *isize, int *jsize) const
{
  if (p_mode == YBus || p_mode == RefreshY || p_mode == onFY || p_mode == posFY) {
    return YMBranch::matrixForwardSize(isize,jsize);
  } else {
    return false;
  }
}
bool gridpack::kalman_filter::KalmanBranch::matrixReverseSize(int *isize, int *jsize) const
{
  if (p_mode == YBus || p_mode == RefreshY || p_mode == onFY || p_mode == posFY) {
    return YMBranch::matrixReverseSize(isize,jsize);
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
bool gridpack::kalman_filter::KalmanBranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == YBus || p_mode == RefreshY) {
    return YMBranch::matrixForwardValues(values);
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
  return false;
}

bool gridpack::kalman_filter::KalmanBranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == YBus || p_mode == RefreshY) {
    return YMBranch::matrixReverseValues(values);
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
  return false;
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::kalman_filter::KalmanBranch::setYBus(void)
{
  YMBranch::setYBus();
  gridpack::ComplexType ret;
  ret = YMBranch::getForwardYBus();
  p_ybusr_frwd = real(ret);
  p_ybusi_frwd = imag(ret);
  ret = YMBranch::getReverseYBus();
  p_ybusr_rvrs = real(ret);
  p_ybusi_rvrs = imag(ret);
}

/**
 * Load values stored in DataCollection object into KalmanBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::kalman_filter::KalmanBranch::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  YMBranch::load(data);

  bool ok = true;
  data->getValue(BRANCH_NUM_ELEMENTS, &p_elems);
  double rvar;
  int ivar;
  std::string svar;
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
    ok = ok && data->getValue(BRANCH_CKT, &svar, idx);
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
void gridpack::kalman_filter::KalmanBranch::setMode(int mode)
{
  if (mode == YBus) {
    YMBranch::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Write output from branches to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if branch is contributing string to output, false otherwise
 */
bool gridpack::kalman_filter::KalmanBranch::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  return false;
}


/**
 * Return the updating factor that will be applied to the ybus matrix at
 * the clear fault phase
 * @return: value of update factor
 */
gridpack::ComplexType
gridpack::kalman_filter::KalmanBranch::getPosfy11YbusUpdateFactor(int sw2_2, int sw3_2)
{
  double retr, reti;
  int i;
  gridpack::kalman_filter::KalmanBus *bus1 =
    dynamic_cast<gridpack::kalman_filter::KalmanBus*>(getBus1().get());
  gridpack::kalman_filter::KalmanBus *bus2 =
    dynamic_cast<gridpack::kalman_filter::KalmanBus*>(getBus2().get());
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
}
gridpack::ComplexType
gridpack::kalman_filter::KalmanBranch::getUpdateFactor()
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
void gridpack::kalman_filter::KalmanBranch::setEvent(const Event &event)
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
    dynamic_cast<gridpack::kalman_filter::KalmanBus*>
      (getBus1().get())->setEvent(idx1,idx2,this);
    dynamic_cast<gridpack::kalman_filter::KalmanBus*>
      (getBus2().get())->setEvent(idx1,idx2,this);
  }
}

