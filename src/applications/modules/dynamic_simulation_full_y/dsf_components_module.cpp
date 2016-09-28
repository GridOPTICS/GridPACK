/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsf_components_module.cpp
 * @author Shuangshuang Jin 
 * @date   2016-07-14 14:27:38 d3g096
 * @date   2014-03-06 15:22:00 d3m956
 * @last modified date   2015-05-13 12:01:00 d3m956
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/include/gridpack.hpp"
#include "dsf_components_module.hpp"

/**
 *  Simple constructor
 */
gridpack::dynamic_simulation::DSFullBus::DSFullBus(void)
{
  p_shunt_gs = 0.0;
  p_shunt_bs = 0.0;
  p_mode = YBUS;
  setReferenceBus(false);
  p_ngen = 0;
  p_from_flag = false;
  p_to_flag = false;
  p_branch = NULL;
  p_isolated = false;
}

/**
 *  Simple destructor
 */
gridpack::dynamic_simulation::DSFullBus::~DSFullBus(void)
{
}

/**
 *  Return size of matrix block contributed by the component
 *  @param isize, jsize: number of rows and columns of matrix block
 *  @return: false if network component does not contribute matrix element
 */
bool gridpack::dynamic_simulation::DSFullBus::matrixDiagSize(int *isize, int *jsize) const
{
  if (YMBus::isIsolated()) return false;
  if (p_mode == YBUS || p_mode == YL || p_mode == PG || p_mode == onFY || p_mode == posFY
  || p_mode == jxd) {
    return YMBus::matrixDiagSize(isize,jsize);
  }  else {
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
bool gridpack::dynamic_simulation::DSFullBus::matrixDiagValues(ComplexType *values)
{
  if (p_mode == YBUS) {
    //bool status = YMBus::matrixDiagValues(values);
    //if (status) printf("idx: %d Real: %f Imag: %f\n",getOriginalIndex(),
    // real(values[0]), imag(values[0]));
    //return status;
    return YMBus::matrixDiagValues(values);
  } else if (p_mode == YL) {
    //printf("bus %d: p_pl = %f, p_ql = %f, p_voltage = %f\n", getOriginalIndex(), p_pl, p_ql, p_voltage);
    //printf("p_ybusr = %f, p_ybusi = %f\n", p_ybusr, p_ybusi);
    p_ybusr = p_ybusr+p_pl/(p_voltage*p_voltage);
    p_ybusi = p_ybusi+(-p_ql)/(p_voltage*p_voltage);
    /* TBD: p_ybusr = p_ybusr+(-p_pg)/(p_voltage*p_voltage);
    p_ybusi = p_ybusi+p_qg/(p_voltage*p_voltage);*/
    //printf(".. p_ybusr = %f, p_ybusi = %f\n", p_ybusr, p_ybusi);
    gridpack::ComplexType ret(p_ybusr, p_ybusi);
    values[0] = ret;
    return true;
  } else if (p_mode == PG) {
    /* TBD: p_ybusr = p_ybusr+(-p_pg)/(p_voltage*p_voltage);
    p_ybusi = p_ybusi+p_qg/(p_voltage*p_voltage);*/
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
         //printf("!!!!!!!%f, %f\n", p_ybusr, p_ybusi);
         if (p_pg[i] < 0) {
           printf("================\n");
           p_ybusr = p_ybusr+(-p_pg[i])/(p_voltage*p_voltage);
           p_ybusi = p_ybusi+p_qg[i]/(p_voltage*p_voltage);
           gridpack::ComplexType ret(p_ybusr, p_ybusi);
           values[0] = ret;
         } else {
           gridpack::ComplexType u(p_ybusr, p_ybusi);
           values[0] = u;
         }
      }
    } else {
      //printf("!!!!!!!%f, %f\n", p_ybusr, p_ybusi);
      gridpack::ComplexType u(p_ybusr, p_ybusi);
      values[0] = u;
    }
    return true;
  } else if (p_mode == jxd) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
#if 0
        double ra = p_r[i] * p_sbase / p_mva[i];
        double xd;
        if (p_dstr[i] == 0) {
          xd = p_dtr[i] * p_sbase / p_mva[i];
        }
        gridpack::ComplexType Y_a(ra, xd);
        Y_a = 1.0 / Y_a;
#else
        //printf("Bus %d here 1\n", getOriginalIndex());
        gridpack::ComplexType Y_a
          = p_generators[i]->NortonImpedence();
        //printf("here 2 real(Y_a): %f imag(Y_a): %f\n",real(Y_a),imag(Y_a));
#endif
        p_ybusr = p_ybusr + real(Y_a);
        p_ybusi = p_ybusi + imag(Y_a);
        gridpack::ComplexType ret(p_ybusr, p_ybusi);
        //printf("here 3 p_ybusr: %f p_ybusi: %f\n",p_ybusr,p_ybusi);
        values[0] = ret;
      }
      //return true;
    } else {
      gridpack::ComplexType u(p_ybusr, p_ybusi);
      values[0] = u;
      //return true;
    }
    //printf("idx: %d Real: %f Imag: %f\n",getOriginalIndex(),
      //real(values[0]), imag(values[0]));
    return true;
  } else if (p_mode == onFY) {
    if (p_from_flag) {
      gridpack::ComplexType ret(0.0, -1.0e12);
      values[0] = ret;
      return true;
    } else {
      return false;
    }
  } else if (p_mode == posFY) {
    if (p_from_flag || p_to_flag) {
      values[0] = dynamic_cast<gridpack::dynamic_simulation::DSFullBranch*>(p_branch)->getUpdateFactor();
      return true;
    } else {
      return false;
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
bool gridpack::dynamic_simulation::DSFullBus::vectorSize(int *size) const
{
  if (!p_isolated) {
    if (p_mode == make_INorton_full) {
      *size = 1;
    }else {
      *size = 2;
    }
    return true;
  } else {
    return false;
  }
}

/**
 * Return the values of the vector block
 * @param values: pointer to vector values
 * @return: false if network component does not contribute
 *        vector element
 */
bool gridpack::dynamic_simulation::DSFullBus::vectorValues(ComplexType *values)
{
  if (!p_isolated) {
    if (p_mode == make_INorton_full) {
      values[0] = 0;
      if (p_ngen > 0) {
        for (int i = 0; i < p_ngen; i++) {
#if 0
          values[0] += p_INorton[i];
#else
          values[0] += p_generators[i]->INorton();
#endif
        }
      }
      //printf("bus id = %d\n, values[0] = %f, %f\n", getOriginalIndex(), real(values[0]), imag(values[0])); 
      return true;
    } else {
      return false;
    }
  } 
  return false;
}

void gridpack::dynamic_simulation::DSFullBus::setValues(ComplexType *values)
{
  int i;
  if (p_mode == make_INorton_full) {
    p_volt_full = values[0];
  }
}

/**
 * Set values of YBus matrix. These can then be used in subsequent
 * calculations
 */
void gridpack::dynamic_simulation::DSFullBus::setYBus(void)
{
  YMBus::setYBus();
  gridpack::ComplexType ret;
  ret = YMBus::getYBus();
  p_ybusr = real(ret);
  p_ybusi = imag(ret);
}

/**
 * Set initial values of vectors for integration. 
 * These can then be used in subsequent calculations
 * @param ts time step
 */
void gridpack::dynamic_simulation::DSFullBus::initDSVect(double ts)
{
  if (p_ngen > 0) {
    for (int i = 0; i < p_ngen; i++) {
#if 0
      // mva
      p_mva[i] = p_sbase / p_mva[i];
      // d0
      p_d0[i] = p_d0[i] / p_mva[i];
      // h
      p_h[i] = p_h[i] / p_mva[i];
      // dtr
      p_dtr[i] = p_dtr[i] * p_mva[i];
      // pelect
      p_pelect.push_back(p_pg[i]);
      // volt
      double eterm = p_voltage;
      double vi = p_angle;
      gridpack::ComplexType v(0.0, vi);
      v = eterm * exp(v);
      p_volt.push_back(v);
      // eprime_s0
      double pelect = p_pg[i];
      double qelect = p_qg[i];
      double currr = sqrt(pelect * pelect + qelect * qelect) / eterm;
      double phi = atan2(qelect, pelect);
      double pi = 4.0*atan(1.0);
      double curri = p_angle - phi;
      gridpack::ComplexType curr(0.0, curri);
      curr = currr * exp(curr);
      gridpack::ComplexType jay(0.0, 1.0);
      gridpack::ComplexType temp = v + jay * (p_dtr[i] * p_mva[i]) * curr;
      p_eprime_s0.push_back(temp);
      // mac_ang_s0
      temp = atan2(imag(p_eprime_s0[i]), real(p_eprime_s0[i]));
      p_mac_ang_s0.push_back(temp);
      // mac_spd_s0
      p_mac_spd_s0.push_back(1.0);
      // eqprime
      p_eqprime.push_back(abs(p_eprime_s0[i]));
      // pmech
      p_pmech.push_back(abs(p_pelect[i]));
      //printf("%f+%f\n", real(p_pmech[i]), imag(p_pmech[i]));
      // Allocate and initialize other vectors 
      p_mac_ang_s1.push_back(0.0);
      p_mac_spd_s1.push_back(0.0);
      p_dmac_ang_s0.push_back(0.0);
      p_dmac_spd_s0.push_back(0.0);
      p_dmac_ang_s1.push_back(0.0);
      p_dmac_spd_s1.push_back(0.0);
      p_eprime_s1.push_back(0.0);
      p_INorton.push_back(0.0);
#else
      ///printf("\ngen %d:********************************************************************************\n", getOriginalIndex());
      p_generators[i]->init(p_voltage,p_angle, ts);
#endif
    }
  } 
}

/**
 * Update values for vectors in each integration time step (Predictor)
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::DSFullBus::predictor_currentInjection(bool flag)
{
  if (p_ngen == 0) return;
  int i;
  for (i = 0; i < p_ngen; i++) {
    p_generators[i]->predictor_currentInjection(flag);
  }
}

/**
 * Update values for vectors in each integration time step (Predictor)
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::DSFullBus::predictor(double t_inc, bool flag)
{
  if (p_ngen == 0) return;
#if 0
  gridpack::ComplexType jay, curr, pelect;
  double sysFreq = 60.0;
  double pi = 4.0*atan(1.0);
  const double basrad = 2.0*pi*sysFreq;
  int i;
  // Reset values of machine parameters after first time step
  if (!flag) {
    for (i = 0; i < p_ngen; i++) {
      p_mac_ang_s0[i] = p_mac_ang_s1[i];
      p_mac_spd_s0[i] = p_mac_spd_s1[i];
      p_eprime_s0[i] = p_eprime_s1[i];
    }
  }
  // Evaluate updated values of machine parameters for integration
  jay = gridpack::ComplexType(0.0,1.0);
  p_INorton.clear();
  for (i = 0; i < p_ngen; i++) {
    // --------- CALL mac_em11(k,S_Steps) to calculate
    // terminal curr: curr = (eprime - volt) / GEN_dtr) -----------
    // curr:
    //printf("eprime_s0 = %f + %fi\n", real(p_eprime_s0[i]), imag(p_eprime_s0[i]));
    curr = p_eprime_s0[i];
    curr -= p_volt[i];
    curr = curr/(jay*p_dtr[i]);
    //printf("curr = %f + %fi\n", real(curr), imag(curr));

    // CALL mac_em2(k, S_Steps):
    // pelect:
    imag(curr) = -imag(curr);
    pelect = real(p_eprime_s0[i] * curr);
    //printf("pelect = %f\n", pelect);
    // dmac_ang:
    p_dmac_ang_s0[i] = (p_mac_spd_s0[i] - 1.0) * basrad;
    //printf("dmac_ang_s0 = %f\n", p_dmac_ang_s0[i]);
    // dmac_spd:
    p_dmac_spd_s0[i] = (p_pmech[i] - pelect - p_d0[i] * (p_mac_spd_s0[i] - 1.0)) 
                        / (2.0 * p_h[i]);
    //printf("dmac_spd_s0 = %f\n", p_dmac_spd_s0[i]);

    p_mac_ang_s1[i] = p_mac_ang_s0[i] + p_dmac_ang_s0[i] * t_inc;
    p_mac_spd_s1[i] = p_mac_spd_s0[i] + p_dmac_spd_s0[i] * t_inc;
    //printf("mac_ang_s1 = %f\n", p_mac_ang_s1[i]);
    //printf("mac_spd_s1 = %f\n", p_mac_spd_s1[i]);

    // CALL mac_em1(k, S_Steps+1):
    p_eprime_s1[i] = exp(p_mac_ang_s1[i] * jay) * p_eqprime[i];
    //printf("eprime_s1 = %f + %fi\n", real(p_eprime_s1[i]), imag(p_eprime_s1[i]));

    // Calculate INorton_full
    p_INorton[i] = p_eprime_s1[i] / (p_dtr[i] * jay);
    //printf("INorton= %f + %fi\n", real(p_INorton[i]), imag(p_INorton[i]));
  }  
#else
  int i;
  for (i = 0; i < p_ngen; i++) {
    p_generators[i]->predictor(t_inc,flag);
  }
#endif
}

/**
 * Update values for vectors in each integration time step (Corrector)
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::DSFullBus::corrector_currentInjection(bool flag)
{
  if (p_ngen == 0) return;
  int i;
  for (i = 0; i < p_ngen; i++) {
    p_generators[i]->corrector_currentInjection(flag);
  }
}

/**
 * Update values for vectors in each integration time step (Corrector)
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::DSFullBus::corrector(double t_inc, bool flag)
{
  if (p_ngen == 0) return;
#if 0
  gridpack::ComplexType jay, curr, pelect;
  double sysFreq = 60.0;
  double pi = 4.0*atan(1.0);
  const double basrad = 2.0*pi*sysFreq;
  int i;
  // Evaluate updated values of machine parameters for integration
  jay = gridpack::ComplexType(0.0,1.0);
  p_INorton.clear();
  for (i = 0; i < p_ngen; i++) {
    // --------- CALL mac_em11(k,S_Steps) to calculate
    // terminal curr: curr = (eprime - volt) / GEN_dtr) -----------
    // curr:
    curr = p_eprime_s1[i];
    curr -= p_volt[i];
    //printf("\np_volt= %f + %fi\n", real(p_volt[i]), imag(p_volt[i]));
    curr = curr/(jay*p_dtr[i]);
    //printf("\ncurr = %f + %fi\n", real(curr), imag(curr));

    // CALL mac_em2(k, S_Steps+1):
    // pelect:
    imag(curr) = -imag(curr);
    pelect = real(p_eprime_s1[i] * curr);
    //printf("pelect = %f\n", pelect);
    // dmac_ang:
    p_dmac_ang_s1[i] = (p_mac_spd_s1[i] - 1.0) * basrad;
    //printf("dmac_ang_s0 = %f\n", p_dmac_ang_s0[i]);
    // dmac_spd:
    p_dmac_spd_s1[i] = (p_pmech[i] - pelect - p_d0[i] * (p_mac_spd_s1[i] - 1.0)) 
                        / (2.0 * p_h[i]);
    //printf("dmac_spd_s0 = %f\n", p_dmac_spd_s0[i]);

    p_mac_ang_s1[i] = p_mac_ang_s0[i] + (p_dmac_ang_s0[i] + p_dmac_ang_s1[i]) / 2.0 * t_inc;
    p_mac_spd_s1[i] = p_mac_spd_s0[i] + (p_dmac_spd_s0[i] + p_dmac_spd_s1[i]) / 2.0 * t_inc;
    //printf("mac_ang_s1 = %f\n", p_mac_ang_s1[i]);

    // CALL mac_em1(k, S_Steps+1):
    p_eprime_s1[i] = exp(p_mac_ang_s1[i] * jay) * p_eqprime[i];
    //printf("eprime_s1 = %f + %fi\n", real(p_eprime_s1[i]), imag(p_eprime_s1[i]));

    // Calculate INorton_full
    p_INorton[i] = p_eprime_s1[i] / (p_dtr[i] * jay);
    //printf("INorton= %f + %fi\n", real(p_INorton[i]), imag(p_INorton[i]));
    
    if (!flag) { 
      printf("%f %f\n", real(p_mac_ang_s0[i]), real(p_mac_spd_s0[i]));
    } else { 
      printf("%f %f\n", real(p_mac_ang_s1[i]), real(p_mac_spd_s1[i]));
    }
  }  
#else
  int i;
  for (i = 0; i < p_ngen; i++) {
    p_generators[i]->corrector(t_inc,flag);
  }
#endif
}

/**
 * Get roter angle of generators
 */
double gridpack::dynamic_simulation::DSFullBus::getAngle()
{
  if (p_ngen < 0) return 0.0;
  int i;
  for (i = 0; i < p_ngen; i++) {
    double angle = p_generators[i]->getAngle();
    return angle;
  }
  return 0.0;
}

/**
 * Set volt from volt_full
 */
void gridpack::dynamic_simulation::DSFullBus::setVolt(bool flag) 
{
  if (p_ngen > 0) {
    for (int i = 0; i < p_ngen; i++) {
#if 0
      p_volt[i] = p_volt_full;
      //printf("p_volt = %f + %fi\n", real(p_volt[i]), imag(p_volt[i]));
#else
      p_generators[i]->setVoltage(p_volt_full);
      if (flag) {
        //if (getOriginalIndex() == 6433) printf("bus id = %d: values[0] = %f, %f, mag = %f\n", getOriginalIndex(), p_volt_full, abs(p_volt_full)); // bpa 
        //if (getOriginalIndex() == 6) printf("bus id = %d: values[0] = %f, %f, mag = %f\n", getOriginalIndex(), p_volt_full, abs(p_volt_full)); // 179
      }
#endif
    }
  }
}

/**
 * Get values of YBus matrix. These can then be used in subsequent
 * calculations
 */
gridpack::ComplexType gridpack::dynamic_simulation::DSFullBus::getYBus(void)
{
  return YMBus::getYBus();
}

/**
 * Load values stored in DataCollection object into DSFullBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::dynamic_simulation::DSFullBus::load(
  const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  YMBus::load(data);

  p_sbase = 100.0;

  data->getValue(BUS_VOLTAGE_ANG, &p_angle);
  data->getValue(BUS_VOLTAGE_MAG, &p_voltage);
  //printf("p_voltage at bus %d: %f\n", getOriginalIndex(), p_voltage);

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
  } else if (p_type == 4) {
    p_isolated = true;
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
  bool has_ex, has_gov;
  DSFGeneratorFactory genFactory;
  p_generators.clear();
  int idx;
  data->getValue(BUS_NUMBER,&idx);
  if (data->getValue(GENERATOR_NUMBER, &p_ngen)) {
    std::string genid;
    //printf("p_ngen = %d\n", p_ngen);
    int icnt = 0;
    for (i=0; i<p_ngen; i++) { 
      int stat;
      data->getValue(GENERATOR_STAT, &stat, i);
      /* TBD: if (stat == 1 && GENERATOR_PG < 0) 
                modify Ybus
              */ 
      data->getValue(GENERATOR_PG, &pg, i);
      //printf("****p_pg[%d] = %f\n", i, p_pg[i]);
      /*if (stat == 1 && pg < 0) {
        p_pg_flag[i] = -1;
      } else {
        p_pg_flag[i] = 1;
      }*/
      std::string model;
      //if (data->getValue(GENERATOR_MODEL, &model, i) && stat == 1) {
      // TBD: if (data->getValue(GENERATOR_MODEL, &model, i) && stat == 1 && GENERATOR_PG >= 0) {
      if (data->getValue(GENERATOR_MODEL, &model, i) && stat == 1 && pg >= 0) {
        p_pg.push_back(pg);
        //std::cout << "generator: " << model << std::endl;
        DSFBaseGeneratorModel *generator
          = genFactory.createGeneratorModel(model);
        has_ex = false;
        has_gov = false;
        data->getValue(HAS_EXCITER, &has_ex, i);
        data->getValue(HAS_GOVERNOR, &has_gov, i);
        if (generator) {
          //boost::shared_ptr<DSFBaseGeneratorModel> tmp;
          //tmp.reset(generator);
          //boost::shared_ptr<DSFBaseGenerator> basegen;
          //basegen.reset(new DSFBaseGenerator);
          //basegen->setGeneratorModel(tmp);
          boost::shared_ptr<DSFBaseGeneratorModel> basegen;
          basegen.reset(generator);
          p_generators.push_back(basegen);
          data->getValue(GENERATOR_ID, &genid, i);
          p_genid.push_back(genid);
          if (has_ex) {
            if (data->getValue(EXCITER_MODEL, &model, i)) {
              //std::cout << "exciter: " << model << std::endl;
              DSFBaseExciterModel *exciter
                = genFactory.createExciterModel(model);
              boost::shared_ptr<DSFBaseExciterModel> ex;
              ex.reset(exciter);
              p_generators[icnt]->setExciter(ex);
            }
          }
          if (has_gov) {
            if (data->getValue(GOVERNOR_MODEL, &model, i)) {
              //std::cout << "governor: " << model << std::endl;
              DSFBaseGovernorModel *governor
                = genFactory.createGovernorModel(model);
              boost::shared_ptr<DSFBaseGovernorModel> gov;
              gov.reset(governor);
              p_generators[icnt]->setGovernor(gov);
            }
          }
        }
        p_generators[icnt]->load(data,i);
        if (has_gov) p_generators[icnt]->getGovernor()->load(data,i);
        if (has_ex) p_generators[icnt]->getExciter()->load(data,i);
        icnt++;
      } else if (stat == 1 && pg < 0) {
        // Evaluate correction to Y-bus
      }
    }
  }
  p_ngen = p_generators.size();
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::dynamic_simulation::DSFullBus::setMode(int mode)
{
  if (mode == YBUS || mode == YL || mode == PG || mode == jxd) {
    YMBus::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Return the value of the voltage magnitude on this bus
 * @return: voltage magnitude
 */
double gridpack::dynamic_simulation::DSFullBus::getVoltage(void)
{
  return 0.0;
}

/**
 * Return the value of the phase angle on this bus
 * @return: phase angle
 */
double gridpack::dynamic_simulation::DSFullBus::getPhase(void)
{
  return 0.0;
}

/**
 * Return whether or not a bus is isolated
 * @return true if bus is isolated
 */
bool gridpack::dynamic_simulation::DSFullBus::isIsolated(void) const
{
  return YMBus::isIsolated();
}

/**
 * Return the number of generators on this bus
 * @return number of generators on bus
 */
int gridpack::dynamic_simulation::DSFullBus::getNumGen(void)
{
  return p_ngen;
}

void gridpack::dynamic_simulation::DSFullBus::setIFunc(void)
{
}

void gridpack::dynamic_simulation::DSFullBus::setIJaco(void)
{
}

/**
 * Check to see if a fault event applies to this bus and set an internal
 * flag marking the bus as the "from" or "to" bus for the event
 * @param from_idx index of "from" bus for fault event
 * @param to_idx index of "to" bus for fault event
 */
void gridpack::dynamic_simulation::DSFullBus::setEvent(int from_idx, int to_idx,
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
void gridpack::dynamic_simulation::DSFullBus::clearEvent()
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
bool gridpack::dynamic_simulation::DSFullBus::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  if (p_ngen == 0) return false;
  int i;
  char buf[128];
  char *ptr = string;
  int idx = getOriginalIndex();
  if (signal == NULL) {
    return false;
  } else if (!strcmp(signal,"watch_header") ||
      !strcmp(signal,"watch")) {
    if (p_ngen == 0) return false;
    int i;
    char buf[128];
    char *ptr = string;
    int len = 0;
    bool ok;
    for (i=0; i<p_ngen; i++) {
      if (p_generators[i]->getWatch()) {
        ///printf("(DSFull::serialWrite) Got to 1\n");
        ok = p_generators[i]->serialWrite(buf,128,signal);
        ///printf("(DSFull::serialWrite) Got to 2\n");
        if (ok) {
          int slen = strlen(buf);
          if (len+slen < bufsize) sprintf(ptr,"%s",buf);
          len += slen;
          ptr += slen;
        }
      }
    }
    if (len > 0) return true;
  } else if (strlen(signal) > 0) {
    int i;
    char buf[128];
    int len = 0;
    bool ok = true;
    //printf("Writing for %d generators\n",p_ngen);
    for (i=0; i<p_ngen; i++) {
      p_generators[i]->serialWrite(buf,128,signal);
      int slen = strlen(buf);
      if (len+slen < bufsize) sprintf(ptr,"%s",buf);
      len += slen;
      ptr += slen;
    }
    if (len > 0) return true;
  }
  return false;
}

/**
 * Add constant impedance load admittance to diagonal elements of
 * Y-matrix
 */
void gridpack::dynamic_simulation::DSFullBus::addLoadAdmittance()
{
  p_ybusr = p_ybusr+p_pl/(p_voltage*p_voltage);
  p_ybusi = p_ybusi+(-p_ql)/(p_voltage*p_voltage);
  //printf("idx: %d %f %f\n", getOriginalIndex(), p_ybusr, p_ybusi);
}

/**
 * Set load on bus
 * @param pl real load
 * @param ql imaginary load
 */
void gridpack::dynamic_simulation::DSFullBus::setLoad(double pl, double ql)
{
  p_pl = pl;
  p_ql = ql;
}

/**
 * Get load on bus
 * @param pl real load
 * @param ql imaginary load
 */
void gridpack::dynamic_simulation::DSFullBus::getLoad(double *pl, double *ql)
{
  *pl = p_pl;
  *ql = p_ql;
}

#ifdef USE_FNCS
/**
 * Retrieve an opaque data item from component.
 * @param data item to retrieve from component
 * @param signal string to control behavior of routine
 * (currently ignored)
 * @return true if component is returning data item,
 * false otherwise
 */
bool gridpack::dynamic_simulation::DSFullBus::getDataItem(void *data, const char *signal)
{
  voltage_data *vdata = static_cast<voltage_data*>(data);
  vdata->busID = getOriginalIndex();
  vdata->voltage = gridpack::ComplexType(p_voltage*sin(p_angle),
      p_voltage*cos(p_angle));
}
#endif

/**
 * Set an internal parameter that specifies that the rotor speed and angle
 * for the generator corresponding to the string tag are to be printed to
 * output
 * @param tag 2-character identifier of generator
 * @param flag set to true to monitor generator
 */

void gridpack::dynamic_simulation::DSFullBus::setWatch(std::string tag, bool flag)
{
  int i;
  for (i=0; i<p_genid.size(); i++) {
    if (tag == p_genid[i]) {
      p_generators[i]->setWatch(flag);
      break;
    }
  }
}

/**
 *  Simple constructor
 */
gridpack::dynamic_simulation::DSFullBranch::DSFullBranch(void)
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
gridpack::dynamic_simulation::DSFullBranch::~DSFullBranch(void)
{
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::dynamic_simulation::DSFullBranch::matrixForwardSize(int *isize, int *jsize) const
{
  if (p_mode == YBUS || p_mode == YL || p_mode == PG || p_mode == onFY || p_mode == posFY
  || p_mode == jxd) { 
    return YMBranch::matrixForwardSize(isize,jsize);
  } else {
    return false;
  }
}
bool gridpack::dynamic_simulation::DSFullBranch::matrixReverseSize(int *isize, int *jsize) const
{
  if (p_mode == YBUS || p_mode == YL || p_mode == PG || p_mode == onFY || p_mode == posFY
  || p_mode == jxd) { 
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
bool gridpack::dynamic_simulation::DSFullBranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == YBUS || p_mode == YL || p_mode == PG || p_mode == jxd) {
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
}

bool gridpack::dynamic_simulation::DSFullBranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == YBUS || p_mode == YL || p_mode == PG || p_mode == jxd) {
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
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::dynamic_simulation::DSFullBranch::setYBus(void)
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
  gridpack::dynamic_simulation::DSFullBus *bus1 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get());
  gridpack::dynamic_simulation::DSFullBus *bus2 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get());
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
 * Load values stored in DataCollection object into DSFullBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::dynamic_simulation::DSFullBranch::load(
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
void gridpack::dynamic_simulation::DSFullBranch::setMode(int mode)
{
  if (mode == YBUS || mode == YL || mode == PG || mode == jxd) {
    YMBranch::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Return the complex admittance of the branch
 * @return: complex addmittance of branch
 */
gridpack::ComplexType gridpack::dynamic_simulation::DSFullBranch::getAdmittance(void)
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
gridpack::dynamic_simulation::DSFullBranch::getTransformer(gridpack::dynamic_simulation::DSFullBus *bus)
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
gridpack::dynamic_simulation::DSFullBranch::getShunt(gridpack::dynamic_simulation::DSFullBus *bus)
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

gridpack::ComplexType
gridpack::dynamic_simulation::DSFullBranch::getPosfy11YbusUpdateFactor(int sw2_2, int sw3_2)
{ 
  double retr, reti;
  int i;
  gridpack::dynamic_simulation::DSFullBus *bus1 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get());
  gridpack::dynamic_simulation::DSFullBus *bus2 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get());
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
  return gridpack::ComplexType(-999.0, -999.0);
}

gridpack::ComplexType 
gridpack::dynamic_simulation::DSFullBranch::getUpdateFactor()
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
void gridpack::dynamic_simulation::DSFullBranch::setEvent(const Event &event)
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
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
      (getBus1().get())->setEvent(idx1,idx2,this);
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
      (getBus2().get())->setEvent(idx1,idx2,this);
  }
}
