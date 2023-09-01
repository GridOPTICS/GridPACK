/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   repca1.cpp
 * @author Shrirang Abhyankar
 * @Updated November 20, 2022
 * 
 * @brief  
 *    Renewable plant controller model
 *    Used for both REPCA1 and REPCTA1
 */

#include <vector>
#include <iostream>
#include <stdio.h>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_plant_model.hpp"
#include "repca1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Repca1Model::Repca1Model(void)
{
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Repca1Model::~Repca1Model(void)
{
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * 
 */
void gridpack::dynamic_simulation::Repca1Model::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  data->getValue(CASE_SBASE,&p_sbase);
  data->getValue(BUS_NUMBER,&p_bus_num);
  data->getValue(GENERATOR_ID, &p_gen_id, idx);
  if(!data->getValue(GENERATOR_MBASE,&p_mbase,idx)) p_mbase = 100.0;
  
  if(!data->getValue(GENERATOR_REPCA_IREG,&ireg,idx)) ireg = p_bus_num;
  if(!data->getValue(GENERATOR_REPCA_RC,&Rc,idx)) Rc = 0.0;
  if(!data->getValue(GENERATOR_REPCA_XC,&Xc,idx)) Xc = 0.01;
  if(!data->getValue(GENERATOR_REPCA_BRCH_BUS_FROM,&fbus,idx)) fbus = 0;
  if(!data->getValue(GENERATOR_REPCA_BRCH_BUS_TO, &tbus, idx)) tbus = 0;
  if(!data->getValue(GENERATOR_REPCA_BRCH_CKT, &br_ckt,idx)) br_ckt = ' ';
  if(!data->getValue(GENERATOR_REPCA_VC_FLAG,&VCompFLAG,idx)) VCompFLAG = 0;
  if(!data->getValue(GENERATOR_REPCA_REF_FLAG,&RefFLAG,idx)) RefFLAG = 0;
  if(!data->getValue(GENERATOR_REPCA_F_FLAG,&FreqFLAG,idx)) FreqFLAG = 0;
  if (!data->getValue(GENERATOR_REPCA_KC,    &Kc     ,idx))   Kc =    0.02;
  if (!data->getValue(GENERATOR_REPCA_TFLTR, &Tfltr  ,idx))   Tfltr = 0.02;
  if (!data->getValue(GENERATOR_REPCA_DBD1,  &dbd1   , idx))  dbd1 =  -0.0;
  if (!data->getValue(GENERATOR_REPCA_DBD2,  &dbd2   , idx))  dbd2  = 0.0;
  if (!data->getValue(GENERATOR_REPCA_EMAX,  &Emax   , idx))  Emax  = 0.3;
  if (!data->getValue(GENERATOR_REPCA_EMIN,  &Emin   , idx))  Emin   = -0.3;
  if (!data->getValue(GENERATOR_REPCA_QMAX,  &Qmax   , idx))  Qmax  = 0.56;
  if (!data->getValue(GENERATOR_REPCA_QMIN,  &Qmin   , idx))  Qmin  = -0.56;
  if (!data->getValue(GENERATOR_REPCA_KP,    &Kp     , idx))  Kp   =  18.0;
  if (!data->getValue(GENERATOR_REPCA_KI,    &Ki     , idx))  Ki    = 5.0;
  if (!data->getValue(GENERATOR_REPCA_TFT,   &Tft    , idx))  Tft    = 0.0;
  if (!data->getValue(GENERATOR_REPCA_TFV,   &Tfv    , idx))  Tfv    = 0.05;
  if (!data->getValue(GENERATOR_REPCA_FDBD1, &fdbd1   , idx)) fdbd1  = 0.0;
  if (!data->getValue(GENERATOR_REPCA_FDBD2, &fdbd2   , idx)) fdbd2  = 0.0;
  if (!data->getValue(GENERATOR_REPCA_VFRZ,  &Vfrz     , idx)) Vfrz   = 0.5;
  if (!data->getValue(GENERATOR_REPCA_DDN,   &Ddn    , idx))  Ddn    = 20.0;
  if (!data->getValue(GENERATOR_REPCA_DUP,   &Dup    , idx))  Dup    =  -10.0;
  if (!data->getValue(GENERATOR_REPCA_TP,    &Tp     , idx))  Tp     = 0.05;
  if (!data->getValue(GENERATOR_REPCA_FEMAX, &femax   , idx)) femax  = 999.0;
  if (!data->getValue(GENERATOR_REPCA_FEMIN, &femin   , idx)) femin  = -999.0;
  if (!data->getValue(GENERATOR_REPCA_KPG,   &Kpg    , idx))  Kpg    = 0.1; 
  if (!data->getValue(GENERATOR_REPCA_KIG,   &Kig    , idx))  Kig    = 0.05;
  if (!data->getValue(GENERATOR_REPCA_PMAX,  &Pmax   , idx))  Pmax   = 1.5;
  if (!data->getValue(GENERATOR_REPCA_PMIN,  &Pmin   , idx))  Pmin   =  -1.5;
  if (!data->getValue(GENERATOR_REPCA_TG,    &Tg     , idx))  Tg     = 0.1;

  // Set constants
  Freq_ref = 1.0;
  
  // Set up model blocks
  // V filter block
  V_filter_blk.setparams(1.0,Tfltr);
  // Q branch filter block
  Qbranch_filter_blk.setparams(1.0,Tfltr);

  // VQerr deadband block
  VQerr_deadband.setparams(dbd1,dbd2);

  // VQerror limiter
  VQerr_limiter.setparams(1.0,Emin,Emax);

  // Qref PI controller
  Qref_PI_blk.setparams(Kp,Ki,Qmin,Qmax,-1000.0,1000.0);

  // Qref lead lag
  Qref_leadlag_blk.setparams(Tft,Tfv);

  // Frequency error deadband
  Freqerr_deadband.setparams(fdbd1,fdbd2);

  // Pbranch filter block
  Pbranch_filter_blk.setparams(1.0,Tp);

  // Frequency error limiter
  Freqerr_limiter.setparams(1.0,femin,femax);

  // Pref PI block
  Pref_PI_blk.setparams(Kpg,Kig,Pmin,Pmax,-1000.0,1000.0);

  // Pref filter block
  Pref_filter_blk.setparams(1.0,Tg);
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Repca1Model::init(double Vm, double Va, double ts)
{
  if(FreqFLAG) {
    double ferr;
    Pref_PI_blk_out = Pref_filter_blk.init_given_y(Pref);
    ferr = Pref_PI_blk.init_given_y(Pref_PI_blk_out);
    
    // ***********
    // Need to use Pbranch if given, using Pg
    // ***********

    Pbranch = Pg; // machine MVA base

    Pbranch_filter_blk_out = Pbranch_filter_blk.init_given_u(Pbranch); 
  }

  // Qref side now
  Qref_PI_blk_out = Qref_leadlag_blk.init_given_y(Qref);
  VQerr_limiter_out = Qref_PI_blk.init_given_y(Qref_PI_blk_out);

  if(!RefFLAG) {
    double temp;

    // ***********
    // Need to actual Qbranch if Q branch is given, using Qg
    // ***********

    Qbranch = Qg; // machine MVA base
    
    temp = Qbranch_filter_blk.init_given_u(Qbranch);
  } else {
    double V_VFLAG;
    if(!VCompFLAG) {

      Qbranch = Qg;
      V_VFLAG = Qbranch*Kc + Vt;
    } else {
      // Only considering Vt, line drop compensation
      // needs to be implemented in full implementation
      // This also means Rc = Xc = 0 in the file
      V_VFLAG = Vt;
    }
    V_filter_blk_out = V_filter_blk.init_given_u(V_VFLAG);
    Vref = V_filter_blk_out;
  }
}

void gridpack::dynamic_simulation::Repca1Model::computeModel(double t_inc,IntegrationStage int_flag, bool Vfreeze)
{
  bool updateState = !Vfreeze; // Do not update state (freeze) when Vfreeze is true
  
  // Pref part
  double ferr,dP;

  if(FreqFLAG) {
    ferr = Freq_ref - Freq;
    ferr = Freqerr_deadband.getoutput(ferr);
    dP = std::max(0.0,Ddn*ferr) + std::min(0.0,Dup*ferr);
    // ***********
    // Need to use Pbranch if given, using Pg
    // ***********
    Pbranch = Pg;
    
    Pbranch_filter_blk_out = Pbranch_filter_blk.getoutput(Pbranch,t_inc,int_flag,true);

    Freqerr_limiter_out = Freqerr_limiter.getoutput(Plant_ref - Pbranch_filter_blk_out + dP);

    Pref_PI_blk_out = Pref_PI_blk.getoutput(Freqerr_limiter_out,t_inc,int_flag,true);

    Pref = Pref_filter_blk.getoutput(Pref_PI_blk_out,t_inc,int_flag,true);
  }

  // Qref part
  double y_VCompFLAG=0.0; // value at VCompFLAG
  double y_RefFLAG = 0.0; // value at RefFLAG

  if(!VCompFLAG) {
    Qbranch = Qg;
    y_VCompFLAG = Qbranch*Kc + Vt;
  } else {
    // **************
    // Line drop compensation not implemented,
    // Only considering voltage control
    // *************
    y_VCompFLAG = Vt;
  }

  if(!RefFLAG) {
    // ***********
    // Need to actual Qbranch if Q branch is given, using Qg
    // ***********
    Qbranch = Qg; // machine MVA base

    Qbranch_filter_blk_out = Qbranch_filter_blk.getoutput(Qbranch,t_inc,int_flag,true);
    y_RefFLAG = Qref - Qbranch_filter_blk_out;
  } else {
    V_filter_blk_out = V_filter_blk.getoutput(y_VCompFLAG,t_inc,int_flag,true);
    y_RefFLAG = Vref - V_filter_blk_out;
  }

  // VQ error deadband block
  VQerr_deadband_out = VQerr_deadband.getoutput(y_RefFLAG);

  // VQ error limiter block
  VQerr_limiter_out = VQerr_limiter.getoutput(VQerr_deadband_out);

  // Qref PI control
  Qref_PI_blk_out = Qref_PI_blk.getoutput(VQerr_limiter_out,t_inc,int_flag,updateState);

  // Qref lead lag block
  Qref = Qref_leadlag_blk.getoutput(Qref_PI_blk_out,t_inc,int_flag,true);
  
}
/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Repca1Model::predictor(double t_inc, bool flag)
{
  if(Vt < Vfrz) {
    Vfreeze = true;
  }
  else Vfreeze = false;

  computeModel(t_inc,PREDICTOR,Vfreeze);
    
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Repca1Model::corrector(double t_inc, bool flag)
{

  computeModel(t_inc,CORRECTOR,Vfreeze);

}

void gridpack::dynamic_simulation::Repca1Model::setGenPQV(double P, double Q, double Vm)
{
  Pg = P;
  Qg = Q;
  Vt = Vm;
}

void gridpack::dynamic_simulation::Repca1Model::setBusFreq(double f)
{
  Freq  = f;
}

void gridpack::dynamic_simulation::Repca1Model::setPrefQext(double Pref_in, double Qref_in)
{
  Pref = Plant_ref = Pref_in;
  Qref = Qref_in;
}

double gridpack::dynamic_simulation::Repca1Model::getPref( )
{
  return Pref;
}

double gridpack::dynamic_simulation::Repca1Model::getQext( )
{
  return Qref;
}

void gridpack::dynamic_simulation::Repca1Model::setExtBusNum(int ExtBusNum)
{
  p_bus_num = ExtBusNum;
}	

void gridpack::dynamic_simulation::Repca1Model::setExtGenId(std::string ExtGenId)
{
  p_gen_id = ExtGenId;
}

/**
 * Set internal state parameter in plant controller
 * @param name character string corresponding to state variable
 * @param value new value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::Repca1Model::setState(std::string name,
    double value)
{
  bool ret = true;
  if(name == "S0") {
    V_filter_blk.setstate(value);
  } else if(name == "S1") {
    Qbranch_filter_blk.setstate(value);
  } else if(name == "S2") {
    Qref_PI_blk.setstate(value);
  } else if(name == "S3") {
    Qref_leadlag_blk.setstate(value);
  } else if(name == "S4") {
    Pbranch_filter_blk.setstate(value);
  } else if(name == "S5") {
    Pref_PI_blk.setstate(value);
  } else if(name == "S6") {
    Pref_filter_blk.setstate(value);
  } else ret = false;
  return ret;
}

/**
 * Get internal state parameter in plant controller
 * @param name character string corresponding to state variable
 * @param value current value for state parameter
 * @return false if no variable corresponding to name is found
 */
bool gridpack::dynamic_simulation::Repca1Model::getState(std::string name,
    double *value)
{
  bool ret = true;
  double state_value = -1E19; // Large negative value is set only when return flag is false, i.e., the name is not found

  if(name == "S0") {
    state_value = V_filter_blk.getstate();
  } else if(name == "S1") {
    state_value = Qbranch_filter_blk.getstate();
  } else if(name == "S2") {
    state_value = Qref_PI_blk.getstate();
  } else if(name == "S3") {
    state_value = Qref_leadlag_blk.getstate();
  } else if(name == "S4") {
    state_value = Pbranch_filter_blk.getstate();
  } else if(name == "S5") {
    state_value = Pref_PI_blk.getstate();
  } else if(name == "S6") {
    state_value = Pref_filter_blk.getstate();
  } else ret = false;

  *value = state_value;
  return ret;
}

/**
 * Write output from exciter to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if governor is contributing string to output, false otherwise
 */
bool gridpack::dynamic_simulation::Repca1Model::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  char *ptr = string;
  bool ret = false;

  if (!strcmp(signal,"watch_header")) {
    char buf[256];
    std::string tag;
    if(p_gen_id[0] != ' ') {
      tag = p_gen_id;
    } else {
      tag = p_gen_id[1];
    }
    sprintf(buf,", %d_%s_REPCA1_S0, %d_%s_REPCA1_S1, %d_%s_REPCA1_S2,%d_%s_REPCA1_S3, %d_%s_REPCA1_S4, %d_%s_REPCA1_S5,%d_%s_REPCA1_S6",p_bus_num,tag.c_str(),p_bus_num,tag.c_str(),p_bus_num,tag.c_str(),p_bus_num,tag.c_str(),p_bus_num,tag.c_str(),p_bus_num,tag.c_str(),p_bus_num,tag.c_str());
    if (strlen(buf) <= bufsize) {
      sprintf(string,"%s",buf);
      ret = true;
    } else {
      ret = false;
    }
  } else if (!strcmp(signal,"watch")) {
    double S0,S1,S2,S3,S4,S5,S6;
    S0 = V_filter_blk.getstate();
    S1 = Qbranch_filter_blk.getstate();
    S2 = Qref_PI_blk.getstate();
    S3 = Qref_leadlag_blk.getstate();
    S4 = Pbranch_filter_blk.getstate();
    S5 = Pref_PI_blk.getstate();
    S6 = Pref_filter_blk.getstate();

    char buf[256];
    sprintf(buf,",%f,%f,%f,%f,%f,%f,%f",
	    S0,S1,S2,S3,S4,S5,S6);
    if (strlen(buf) <= bufsize) {
      sprintf(string,"%s",buf);
      ret = true;
    } else {
      ret = false;
    }
  }

  return ret;
}
