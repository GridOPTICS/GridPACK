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

  // Qext PI controller
  Qext_PI_blk.setparams(Kp,Ki,Qmin,Qmax,-1000.0,1000.0);

  // Qext lead lag
  Qext_leadlag_blk.setparams(Tft,Tfv);

  // Frequency error deadband
  Freqerr_deadband.setparams(fdbd1,fdbd2);

  // Pbranch filter block
  Pbranch_filter_blk.setparams(1.0,Tp);

  // Frequency error limiter
  Freqerr_limiter.setparams(1.0,femin,femax);

  // Pext PI block
  Pext_PI_blk.setparams(Kpg,Kig,Pmin,Pmax,-1000.0,1000.0);

  // Pext filter block
  Pext_filter_blk.setparams(1.0,Tg);
}

/**
 * Initialize exciter model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 * @param ts time step 
 */
void gridpack::dynamic_simulation::Repca1Model::init(double Vm, double Va, double ts)
{
  // Assume Pext and Qext are 0
  Pext = Qext = 0.0;

  if(FreqFLAG) {
    double ferr;
    Pext_PI_blk_out = Pext_filter_blk.init_given_y(Pext);
    ferr = Pext_PI_blk.init_given_y(Pext_PI_blk_out);

    // ***********
    // Need to add Pbranch, assume = Pref for now
    // ***********
    Pbranch = Pref*p_mbase/p_sbase; // On system MVA base
  }

  // Qext side now
  Qext_PI_blk_out = Qext_leadlag_blk.init_given_y(Qext);
  VQerr_limiter_out = Qext_PI_blk.init_given_y(Qext_PI_blk_out);

  if(!RefFLAG) {
    double temp;
    Qbranch = Qref*p_mbase/p_sbase; // On system MVA base
    temp = Qbranch_filter_blk.init_given_u(Qref);
  } else {
    if(!VCompFLAG) {
      double V_VFLAG;
      Qbranch = Qref;
      V_VFLAG = Qbranch*Kc + Vt;

      V_filter_blk_out = V_filter_blk.init_given_u(V_VFLAG);
      Vref = V_filter_blk_out;
    } else {
      // **************
      // **** TO BE IMPLEMENTED
      // *************
    }
  }
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Repca1Model::predictor(double t_inc, bool flag)
{
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Repca1Model::corrector(double t_inc, bool flag)
{
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
  Pref = Pref_in;
  Qref = Qref_in;
}

double gridpack::dynamic_simulation::Repca1Model::getPref( )
{
  return Pref + Pext;
}

double gridpack::dynamic_simulation::Repca1Model::getQext( )
{
  return Qref + Qext;
}

void gridpack::dynamic_simulation::Repca1Model::setExtBusNum(int ExtBusNum)
{
  p_bus_num = ExtBusNum;
}	

void gridpack::dynamic_simulation::Repca1Model::setExtGenId(std::string ExtGenId)
{
  p_gen_id = ExtGenId;
}


