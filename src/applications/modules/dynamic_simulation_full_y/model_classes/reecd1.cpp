/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   reecd1.cpp
 * @author Shuangshuang Jin
 * @Created  January 04, 2023
 *
 * @brief
 * Renewable Energy Electrical Controller Model REECD1
 *
 */

#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

#include "base_exciter_model.hpp"
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "reecd1.hpp"

//bool gridpack::dynamic_simulation::Reecd1Model::getVoltageDip(double Vt)
bool gridpack::dynamic_simulation::Reecd1Model::getVoltageDip(double Vt_filter)
{
  // sj: Replace the following
  /*if(Vt < Vdip || Vt > Vup)
    return true;
  else return false;*/
  // sj: by below?
  if (Vt_filter < Vdip || Vt_filter > Vup)
      return true;
  else return false;
}

double gridpack::dynamic_simulation::Reecd1Model::computeS6()
{
    double s6;
    gridpack::ComplexType temp(rc, Xc);
    if (VCMPFLAG == 1) s6 = abs(Vt - temp * It);
    else s6 = Vt + Kc * Qelect;
    s6 = Vcmpflag_filter_blk.getoutput(Tr1);
    return s6;
    
}
/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Reecd1Model::Reecd1Model(void)
{
  omega_g = 1.0;
    Paux = 0;
    deltaPaux = 0;
    deltaQref = 0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Reecd1Model::~Reecd1Model(void) {}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::Reecd1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection> data, int idx) {
  data->getValue(CASE_SBASE, &p_sbase);
  data->getValue(BUS_NUMBER, &p_bus_num);
  data->getValue(GENERATOR_ID, &p_gen_id, idx);
  if (!data->getValue(GENERATOR_MBASE, &p_mbase, idx))
    p_mbase = 100.0; // Machine base

  if(!data->getValue(HAS_WIND_DRIVETRAIN,&p_has_drivetrain,idx))
    p_has_drivetrain = false;
  
  // Parameters
  // default values form
  // https://www.wecc.org/Reliability/WECC%20Wind%20Plant%20Dynamic%20Modeling%20Guidelines.pdf
  if (!data->getValue(GENERATOR_REECA_PFFLAG, &PFFLAG, idx))
    PFFLAG = 0;
  if (!data->getValue(GENERATOR_REECA_VFLAG, &VFLAG, idx))
    VFLAG = 0;
  if (!data->getValue(GENERATOR_REECA_QFLAG, &QFLAG, idx))
    QFLAG = 0;
  if (!data->getValue(GENERATOR_REECA_PFLAG, &PFLAG, idx))
    PFLAG = 0;
  if (!data->getValue(GENERATOR_REECA_PQFLAG, &PQFLAG, idx))
    PQFLAG = 0;
    if (!data->getValue(GENERATOR_REECD_VCMPFLAG, &VCMPFLAG, idx))
      VCMPFLAG = 0; // sj: set to 0?
    printf("flags:%d %d %d %d %d %d\n", PFFLAG, VFLAG, QFLAG, PFLAG, PQFLAG, VCMPFLAG);

  if (!data->getValue(GENERATOR_REECA_TRV, &Trv, idx))
    Trv = 0.0;
  if (!data->getValue(GENERATOR_REECA_DBD1, &dbd1, idx))
    dbd1 = -0.05;
  if (!data->getValue(GENERATOR_REECA_DBD2, &dbd2, idx))
    dbd2 = 0.05;
  if (!data->getValue(GENERATOR_REECA_VDIP, &Vdip, idx))
    Vdip = -99.0;
  if (!data->getValue(GENERATOR_REECA_VUP, &Vup, idx))
    Vup = 99.0;
  if (!data->getValue(GENERATOR_REECA_KQV, &Kqv, idx))
    Kqv = 0.0;
    
    printf("%f %f %f %f %f %f\n", Trv, dbd1, dbd2, Vdip, Vup, Kqv);

  if (!data->getValue(GENERATOR_REECA_IQH1, &Iqh1, idx))
    Iqh1 = 1.05;
  if (!data->getValue(GENERATOR_REECA_IQL1, &Iql1, idx))
    Iql1 = -1.05;
  if (!data->getValue(GENERATOR_REECA_VREF0, &Vref0, idx))
    Vref0 = 0.0;
  if (!data->getValue(GENERATOR_REECA_IQFRZ, &Iqfrz, idx))
    Iqfrz = 0.15;
  if (!data->getValue(GENERATOR_REECA_THLD, &Thld, idx))
    Thld = 0.0;
  if (!data->getValue(GENERATOR_REECA_THLD2, &Thld2, idx))
    Thld2 = 0.0;
  if (!data->getValue(GENERATOR_REECA_TP, &Tp, idx))
    Tp = 0.05;
  if (!data->getValue(GENERATOR_REECA_QMAX, &Qmax, idx))
    Qmax = 0.436;
  if (!data->getValue(GENERATOR_REECA_QMIN, &Qmin, idx))
    Qmin = -0.436;
  if (!data->getValue(GENERATOR_REECA_VMAX, &Vmax, idx))
    Vmax = 1.1;
  if (!data->getValue(GENERATOR_REECA_VMIN, &Vmin, idx))
    Vmin = 0.9;
  if (!data->getValue(GENERATOR_REECA_KQP, &Kqp, idx))
    Kqp = 0.0;
  if (!data->getValue(GENERATOR_REECA_KQI, &Kqi, idx))
    Kqi = 0.1;
  if (!data->getValue(GENERATOR_REECA_KVP, &Kvp, idx))
    Kvp = 0.0;
  if (!data->getValue(GENERATOR_REECA_KVI, &Kvi, idx))
    Kvi = 40.0;
  if (!data->getValue(GENERATOR_REECA_VBIAS, &Vbias, idx))
    Vbias = 0.0;

  if (!data->getValue(GENERATOR_REECA_IMAX, &Imax, idx))
    Imax = 1.82;
  if (!data->getValue(GENERATOR_REECA_TIQ, &Tiq, idx))
    Tiq = 0.02;
  if (!data->getValue(GENERATOR_REECA_DPMAX, &dPmax, idx))
    dPmax = 999.0;
  if (!data->getValue(GENERATOR_REECA_DPMIN, &dPmin, idx))
    dPmin = -999.0;
  if (!data->getValue(GENERATOR_REECA_PMAX, &Pmax, idx))
    Pmax = 1.0;
  if (!data->getValue(GENERATOR_REECA_PMIN, &Pmin, idx))
    Pmin = 0.0;
  if (!data->getValue(GENERATOR_REECA_TPORD, &Tpord, idx))
    Tpord = 0.02;
    printf("Imax=%f,Tiq=%f,dPmax=%f,dPmin=%f,Pmax=%f,Pmin=%f,Tpord=%f\n",
           Imax, Tiq, dPmax, dPmin, Pmax, Pmin, Tpord);

  if (!data->getValue(GENERATOR_REECA_VP1, &Vp1, idx))
    Vp1 = -1.0;
  if (!data->getValue(GENERATOR_REECA_IP1, &Ip1, idx))
    Ip1 = 1.1;
  if (!data->getValue(GENERATOR_REECA_VP2, &Vp2, idx))
    Vp2 = 0.0;
  if (!data->getValue(GENERATOR_REECA_IP2, &Ip2, idx))
    Ip2 = 1.1;
  if (!data->getValue(GENERATOR_REECA_VP3, &Vp3, idx))
    Vp3 = 1.0;
  if (!data->getValue(GENERATOR_REECA_IP3, &Ip3, idx))
    Ip3 = 1.1;
  if (!data->getValue(GENERATOR_REECA_VP4, &Vp4, idx))
    Vp4 = 2.0;
  if (!data->getValue(GENERATOR_REECA_IP4, &Ip4, idx))
    Ip4 = 1.1;
    
    printf("Vp1~4, Ip1~4:%f,%f,%f,%f,%f,%f,%f,%f\n", Vp1, Ip1, Vp2, Ip2, Vp3, Ip3, Vp4, Ip4);

    // sj: check initial values?
    if (!data->getValue(GENERATOR_REECD_VP5, &Vp5, idx))
      Vp5 = 3.0;
    if (!data->getValue(GENERATOR_REECD_IP5, &Ip5, idx))
      Ip5 = 1.1;
    if (!data->getValue(GENERATOR_REECD_VP6, &Vp6, idx))
      Vp6 = 4.0;
    if (!data->getValue(GENERATOR_REECD_IP6, &Ip6, idx))
      Ip6 = 1.1;
    if (!data->getValue(GENERATOR_REECD_VP7, &Vp7, idx))
      Vp7 = 5.0;
    if (!data->getValue(GENERATOR_REECD_IP7, &Ip7, idx))
      Ip7 = 1.1;
    if (!data->getValue(GENERATOR_REECD_VP8, &Vp8, idx))
      Vp8 = 6.0;
    if (!data->getValue(GENERATOR_REECD_IP8, &Ip8, idx))
      Ip8 = 1.1;
    if (!data->getValue(GENERATOR_REECD_VP9, &Vp9, idx))
      Vp9 = 7.0;
    if (!data->getValue(GENERATOR_REECD_IP9, &Ip9, idx))
      Ip9 = 1.1;
    if (!data->getValue(GENERATOR_REECD_VP10, &Vp10, idx))
      Vp10 = 8.0;
    if (!data->getValue(GENERATOR_REECD_IP10, &Ip10, idx))
      Ip10 = 1.1;
    
  if (!data->getValue(GENERATOR_REECA_VQ1, &Vq1, idx))
    Vq1 = -1.0;
  if (!data->getValue(GENERATOR_REECA_IQ1, &Iq1, idx))
    Iq1 = 1.1;
  if (!data->getValue(GENERATOR_REECA_VQ2, &Vq2, idx))
    Vq2 = 0.0;
  if (!data->getValue(GENERATOR_REECA_IQ2, &Iq2, idx))
    Iq2 = 1.1;
  if (!data->getValue(GENERATOR_REECA_VQ3, &Vq3, idx))
    Vq3 = 1.0;
  if (!data->getValue(GENERATOR_REECA_IQ3, &Iq3, idx))
    Iq3 = 1.1;
  if (!data->getValue(GENERATOR_REECA_VQ4, &Vq4, idx))
    Vq4 = 2.0;
  if (!data->getValue(GENERATOR_REECA_IQ4, &Iq4, idx))
    Iq4 = 1.1;
    
    printf("Vq1~4, Iq1~4:%f,%f,%f,%f,%f,%f,%f,%f\n", Vq1, Iq1, Vq2, Iq2, Vq3, Iq3, Vq4, Iq4);
    
    // sj: check initial values?
    if (!data->getValue(GENERATOR_REECD_VP5, &Vq5, idx))
      Vq5 = 3.0;
    if (!data->getValue(GENERATOR_REECD_IP5, &Iq5, idx))
      Iq5 = 1.1;
    if (!data->getValue(GENERATOR_REECD_VP6, &Vq6, idx))
      Vq6 = 4.0;
    if (!data->getValue(GENERATOR_REECD_IP6, &Iq6, idx))
      Iq6 = 1.1;
    if (!data->getValue(GENERATOR_REECD_VP7, &Vq7, idx))
      Vq7 = 5.0;
    if (!data->getValue(GENERATOR_REECD_IP7, &Iq7, idx))
      Iq7 = 1.1;
    if (!data->getValue(GENERATOR_REECD_VP8, &Vq8, idx))
      Vq8 = 6.0;
    if (!data->getValue(GENERATOR_REECD_IP8, &Iq8, idx))
      Iq8 = 1.1;
    if (!data->getValue(GENERATOR_REECD_VP9, &Vq9, idx))
      Vq9 = 7.0;
    if (!data->getValue(GENERATOR_REECD_IP9, &Iq9, idx))
      Iq9 = 1.1;
    if (!data->getValue(GENERATOR_REECD_VP10, &Vq10, idx))
      Vq10 = 8.0;
    if (!data->getValue(GENERATOR_REECD_IP10, &Iq10, idx))
      Iq10 = 1.1;
    
    // sj: check initial values?
    if (!data->getValue(GENERATOR_REECD_RC, &rc, idx))
      rc = 0.0;
    if (!data->getValue(GENERATOR_REECD_XC, &Xc, idx))
      Xc = 0.0;
    if (!data->getValue(GENERATOR_REECD_TR1, &Tr1, idx))
      Tr1 = 0.0;
    if (!data->getValue(GENERATOR_REECD_KC, &Kc, idx))
      Kc = 0.0;
    if (!data->getValue(GENERATOR_REECD_KE, &Ke, idx))
      Ke = 0.0;
    if (!data->getValue(GENERATOR_REECD_VBIKL, &Vblkl, idx))
      Vblkl = 0.0;
    if (!data->getValue(GENERATOR_REECD_VBIKH, &Vblkh, idx))
      Vblkh = 0.0;
    if (!data->getValue(GENERATOR_REECD_TBIK, &Tblk, idx))
      Tblk = 0.0;
    
    printf("Vp5 = %f, Ip5 = %f, Vq5 = %f, Iq5 = %f, Tblk = %f\n", Vp5, Ip5, Vq5, Iq5, Tblk);

  Ipmax = Imax;
  Ipmin = 0.0;
  Iqmax = Imax;
  Iqmin = -Iqmax;
  
  /* Set up blocks */
  // Vt_filter block
  Vt_filter_blk.setparams(1.0,Trv);
  // Voltage error deadband
  V_err_deadband.setparams(dbd1,dbd2);
  // Iqv_limit_blk
  Iqv_limit_blk.setparams(1.0,Iql1,Iqh1);

  // Iqcmd limiter blk
  Iqcmd_limit_blk.setparams(1.0,Iqmin,Iqmax);
  
  // Pe filter block
  Pe_filter_blk.setparams(1.0,Tp);
  // Q limiter block
  Qlim_blk.setparams(1.0,Qmin,Qmax);
  // Q PI control
  Q_PI_blk.setparams(Kqp,Kqi,Vmin,Vmax,-10000.0,10000.0);
  // Vlimiter block
  Vlim_blk.setparams(1.0,Vmin,Vmax);
  // Verr PI control
  Verr_PI_blk.setparams(Kvp,Kvi);
  // Iq lag block
  Iq_lag_blk.setparams(1.0,Tiq);
    
    Vcmpflag_filter_blk.setparams(1.0, Tr1);

  // Vt filter output used in division
  Vt_filter_lowcap_blk.setparams(1.0,0.01,2.0);

  // Initialize VDL1 (V-P) and VDL2 (V-Q)
  double uin[10], yin[10];
  uin[0] = Vp1; yin[0] = Ip1;
  uin[1] = Vp2; yin[1] = Ip2;
  uin[2] = Vp3; yin[2] = Ip3;
  uin[3] = Vp4; yin[3] = Ip4;
    uin[4] = Vp5; yin[4] = Ip5;
    uin[5] = Vp6; yin[5] = Ip6;
    uin[6] = Vp7; yin[6] = Ip7;
    uin[7] = Vp8; yin[7] = Ip8;
    uin[8] = Vp9; yin[8] = Ip9;
    uin[9] = Vp10; yin[9] = Ip10;

  VDL1.setparams(10,uin,yin);

  uin[0] = Vq1; yin[0] = Iq1;
  uin[1] = Vq2; yin[1] = Iq2;
  uin[2] = Vq3; yin[2] = Iq3;
  uin[3] = Vq4; yin[3] = Iq4;
    uin[4] = Vq5; yin[4] = Iq5;
    uin[5] = Vq6; yin[5] = Iq6;
    uin[6] = Vq7; yin[6] = Iq7;
    uin[7] = Vq8; yin[7] = Iq8;
    uin[8] = Vq9; yin[8] = Iq9;
    uin[9] = Vq10; yin[9] = Iq10;

  VDL2.setparams(10,uin,yin);

  // Pref limiter block
  Pref_limit_blk.setparams(1.0,dPmin,dPmax);

  // Pord block
  Pord_blk.setparams(1.0,Tpord,Pmin,Pmax,-1000.0,1000.0);

  // Ipcmd limiter block
  Ipcmd_limit_blk.setparams(1.0,Ipmin,Ipmax);

  Voltage_dip_prev = false;
  thld_timer = -1.0;
}

/**
 * Initialize generator model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 */
void gridpack::dynamic_simulation::Reecd1Model::init(double Vm, double Va,
                                                     double ts)
{
  double Pord_blk_in,Iq_Qflag;
  
  // Initialize Vt filter block
  Vt_filter = Vt_filter_blk.init_given_u(Vt);

  // Ipcmd related blocks initialization
  // Initial value of Ipcmd provided by generator controller
  //Pord = Vt_filter*Ipcmd;
    Pord = Vt_filter*Ipcmd-(Paux+deltaPaux);
  Pord_blk_in = Pord_blk.init_given_y(Pord);

  if(PFLAG == 0) Pref = Pord_blk_in;
  else {
    if(!p_has_drivetrain) {
      Pref = Pord_blk_in;
    } else {
      omega_g = 1.0;
      // *********
      // Need to do drive train initialization here to get omega_g
      // Using omega_g = 1.0 for now
      // *********
      Pref = Pord_blk_in/omega_g;
    }
  }

  // Iqcmd related blocks initialization
  // Initial value of Iqcmd provided by generator controller
  // Initialization of Iqinj path
  //Iqinj_sw = 0;

  // Is there a voltage dip?
  // Should not be during initialization
  //Voltage_dip = getVoltageDip(Vt);

  // Initialize Vref0 = Vt if Vref0 = 0 in the data file
  if(fabs(Vref0) < 1e-6) Vref0 = Vt;
  
  //if(Voltage_dip == 0) {
    //Iqinj = 0;
  //} else {
    V_err = V_err_deadband.getoutput(Vref0 - Vt_filter);
    Iqv = Kqv*V_err;
    Iqinj = Iqv_limit_blk.getoutput(Iqv);
  //}
  Iq_Qflag = Iqcmd - Iqinj;
  
  double Iq_lag_blk_in, Q_Pfflag;
  double V_err_PI_blk_in, Q_PI_blk_in;
  // Consider cases based on combination of VFLAG and QFLAG
  if(!QFLAG) {
    Iq_lag_blk_in = Iq_lag_blk.init_given_y(Iq_Qflag);
    Q_Pfflag = Vt_filter*Iq_lag_blk_in;
  } else {
    // QFLAG = 1
    V_err_PI_blk_in = Verr_PI_blk.init_given_y(Iq_Qflag);
      V_err_PI_blk_in += computeS6();
    if(!VFLAG) {
      Q_Pfflag = V_err_PI_blk_in + Vt_filter - Vbias;
    } else {
      Q_PI_blk_in = Q_PI_blk.init_given_y(V_err_PI_blk_in + Vt_filter);
      // ***** TO BE IMPLEMENTED ***
      // Need to get Qelec
    }
  }

  if(!PFFLAG) {
    Qref = Q_Pfflag;
  } else {
    // **** TO BE IMPLEMENTED
    // Need to get Pe
  }
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Reecd1Model::predictor(double t_inc,
                                                          bool flag)
{
  //Voltage_dip = getVoltageDip(Vt);
    Voltage_dip = getVoltageDip(Vt_filter);

  /*if(!Voltage_dip) {
    if(Voltage_dip_prev) {
      // Recovered from voltage dip
      // Check if thld == 0, in this case there is
      // no transition to state 2
      if(fabs(Thld) < 1e-6) Iqinj_sw = 0;
      else if(fabs(Thld) > 1e-6) {
    // Thld is positive, transition to state 2,
    // Set timer
    thld_timer = 0.0;
    Iqinj_sw = 2;
      }    else {
    // Thld is negative, transition to state 1
    thld_timer = 0.0;
    Iqinj_sw = 1;
      }
    } else {
      if(Iqinj_sw == 1) {
    thld_timer += t_inc;
    if(thld_timer > 1 - Thld) Iqinj_sw = 0;
      } else if(Iqinj_sw == 2) {
    thld_timer += t_inc;
    if(thld_timer > Thld) Iqinj_sw = 0;
      }
    }
  } else {
    if(Iqinj_sw == 0) {
      Iqinj_sw = 1;
      thld_timer = 0.0;
    }
  }*/

  //computeModel(Voltage_dip,Iqinj_sw,t_inc,PREDICTOR);
    computeModel(Voltage_dip,t_inc,PREDICTOR);
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Reecd1Model::corrector(double t_inc,
                                                          bool flag)
{
  Voltage_dip = getVoltageDip(Vt);

  //computeModel(Voltage_dip,Iqinj_sw,t_inc,CORRECTOR);
    computeModel(Voltage_dip,t_inc,CORRECTOR);

  Voltage_dip_prev = Voltage_dip; // Value of Voltage_dip at previous time-step

}


void gridpack::dynamic_simulation::Reecd1Model::computeModel(bool Voltage_dip,double t_inc,IntegrationStage int_flag)
{
  bool updateState = !Voltage_dip; // freezestate logic : updateState = 0 when no voltage dip, otherwise 1
  
  // Get terminal voltage measurement
  Vt_filter = Vt_filter_blk.getoutput(Vt,t_inc,int_flag,true);
  // Cap Vt_filter at 0.01, i.e., if Vt_filter < 0.01, Vt_filter = 0.01
  Vt_filter_lowcap_out = Vt_filter_lowcap_blk.getoutput(Vt_filter);

  // Current limit logic, calculate limits for Ipcmd and Iqcmd
  CurrentLimitLogic(PQFLAG,Vt_filter,Ipcmd,Iqcmd,&Ipmin,&Ipmax,&Iqmin,&Iqmax);

  // Pref limiter output
  Pref_limit_out = Pref_limit_blk.getoutput(Pref);

  double Pord_blk_in;

  // Pord block input based on PFLAG
  if(!PFLAG) {
    Pord_blk_in = Pref_limit_out;
  } else {
    Pord_blk_in = Pref_limit_out*omega_g;
  }

  // Pord block output, updateState = 0 (freezeState) if voltage_dip = 1
  Pord = Pord_blk.getoutput(Pord_blk_in,t_inc,int_flag,updateState);
    
    double temp = Pord + Paux + deltaPaux;

  // Ipcmd - output to generator controller
  //Ipcmd = Ipcmd_limit_blk.getoutput(Pord/Vt_filter_lowcap_out,Ipmin,Ipmax);
    Ipcmd = Ipcmd_limit_blk.getoutput(temp/Vt_filter_lowcap_out,Ipmin,Ipmax);

  // ****************
  // Iqcmd calculation
  // ****************
  double Qext;
  if(!PFFLAG) {
    //Qext = Qref;
      Qext = Qref + deltaQref;
  } else {
    // **** TO BE IMPLEMENTED
    // Need to get Pe
  }

  double Iq_lag_blk_in, Iq_Qflag;
  if(!QFLAG) {
    // Input to Iq_lag_blk
    Iq_lag_blk_in = Qext/Vt_filter_lowcap_out;
    Iq_Qflag = Iq_lag_blk.getoutput(Iq_lag_blk_in,t_inc,int_flag,updateState);
  } else {
    double Vlim_blk_in;
    if(!VFLAG) {
      Vlim_blk_in = Qext + Vbias;
    } else {
      // ******* Not implemented yet
      // Need to get Qelec
      // Calculate Vlim_blk_in
      // *************
    }
      Vlim_blk_in -= computeS6();
    double Verr_PI_blk_in;
    Verr_PI_blk_in = Vlim_blk.getoutput(Vlim_blk_in);

    Iq_Qflag = Verr_PI_blk.getoutput(Verr_PI_blk_in - Vt_filter,t_inc,Iqmin,Iqmax,-1000.0,1000.0,int_flag,updateState);
  }

  //if(Iqinj_sw == 0) {
    //Iqinj = 0;
  //} else if(Iqinj_sw == 1) {
    /* Get Iqinj */
    V_err = V_err_deadband.getoutput(Vref0 - Vt_filter);
    Iqv = Kqv*V_err;
    Iqinj = Iqv_limit_blk.getoutput(Iqv);
  //} else if(Iqinj_sw == 2) {
  //  Iqinj = Iqfrz;
  //}

  Iqcmd = Iq_Qflag + Iqinj;

  Iqcmd = Iqcmd_limit_blk.getoutput(Iqcmd,Iqmin,Iqmax);
}

/**
 * Implements current limiting logic updating limits for Ipcmd and Iqcmd
 * Logic from PowerWorld online documentation
 * https://www.powerworld.com/WebHelp/Content/TransientModels_HTML/Exciter%20REEC_A.htm
 **/
void gridpack::dynamic_simulation::Reecd1Model::CurrentLimitLogic(int PQFLAG,double Vt_filter, double Ipcmd, double Iqcmd,double *Ipmin_out, double *Ipmax_out, double *Iqmin_out, double *Iqmax_out)
{
  double Ipmax_temp,Iqmax_temp;
  double I_temp;

  *Ipmin_out = 0.0;
  
  // Iqmax look up from VDL1
  Iqmax_temp = VDL1.getoutput(Vt_filter);
  // Ipmax look up from VDL2
  Ipmax_temp = VDL2.getoutput(Vt_filter);

  if(!PQFLAG) { // Q priority
    if (Imax < Iqmax_temp) Iqmax_temp = Imax;
    // *************
    // Need to add timer logic here
    // *************
    I_temp = Imax*Imax - Iqcmd*Iqcmd;
    if(I_temp < 0) I_temp = 0;
    else I_temp = sqrt(I_temp);
    if(I_temp < Ipmax_temp) Ipmax_temp = I_temp;
  } else { // P priority
    if(Imax < Ipmax_temp) Ipmax_temp = Imax;
    I_temp = Imax*Imax - Ipcmd*Ipcmd;
    if(I_temp < 0) I_temp = 0.0;
    else I_temp = sqrt(I_temp);
    if(I_temp < Iqmax_temp) Iqmax_temp = I_temp;
  }

  *Ipmax_out =  Ipmax_temp;
  *Iqmin_out = -Iqmax_temp;
  *Iqmax_out =  Iqmax_temp;
}

/**
 * Set the value of the Vterminal
 * @return value of field current
 */
void gridpack::dynamic_simulation::Reecd1Model::setVterminal(double Vm) {
  Vt = Vm;
}

/**
 * Set the value of the omega_g
 *
 */
void gridpack::dynamic_simulation::Reecd1Model::setOmega(double omega)
{
  omega_g = omega;
}

void gridpack::dynamic_simulation::Reecd1Model::setIpcmdIqcmd(double Ipcmd_in,
                                                              double Iqcmd_in) {
  Ipcmd = Ipcmd_in;
  Iqcmd = Iqcmd_in;
}

void gridpack::dynamic_simulation::Reecd1Model::setPrefQext(double Pref_in,
                                                            double Qref_in) {
  Pref = Pref_in;
  Qref = Qref_in;
}

void gridpack::dynamic_simulation::Reecd1Model::setdeltaQref(double deltaQref_in) {
  deltaQref = deltaQref_in;
}

void gridpack::dynamic_simulation::Reecd1Model::setPaux(double Paux_in,
                                                              double deltaPaux_in) {
  Paux = Paux_in;
  deltaPaux = deltaPaux_in;
}

void gridpack::dynamic_simulation::Reecd1Model::setIt(double It_in) {
    It = It_in;
}

void gridpack::dynamic_simulation::Reecd1Model::setQelect(double Qelect_in) {
    Qelect = Qelect_in;
}

void gridpack::dynamic_simulation::Reecd1Model::setGeneratorPower(double Pg,
                                                            double Qg) {
  Pgen = Pg;
  Qgen = Qg;
}


double gridpack::dynamic_simulation::Reecd1Model::getPref() { return Pref; }

double gridpack::dynamic_simulation::Reecd1Model::getQext() { return Qref; }

double gridpack::dynamic_simulation::Reecd1Model::getIpcmd() { return Ipcmd; }

double gridpack::dynamic_simulation::Reecd1Model::getIqcmd() { return Iqcmd; }

  double gridpack::dynamic_simulation::Reecd1Model::getdeltaQref() { return deltaQref; }
  double gridpack::dynamic_simulation::Reecd1Model::getPaux() { return Paux; }
  double gridpack::dynamic_simulation::Reecd1Model::getdeltaPaux() { return deltaPaux; }


// Yuan added below 2020-6-23
/**
 * Set the exciter bus number
 * @return value of exciter bus number
 */
void gridpack::dynamic_simulation::Reecd1Model::setExtBusNum(int ExtBusNum) {
  p_bus_num = ExtBusNum;
}

/**
 * Set the exciter generator id
 * @return value of generator id
 */
void gridpack::dynamic_simulation::Reecd1Model::setExtGenId(
    std::string ExtGenId) {
  p_gen_id = ExtGenId;
}
// Yuan added above 2020-6-23

/**
 * Return Pord
 * @return Pord
 */
double gridpack::dynamic_simulation::Reecd1Model::getPord()
{
  return Pord;
}
