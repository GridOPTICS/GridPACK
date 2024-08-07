/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   acmotor.cpp
 * @author Shuangshuang Jin 
 * @Last modified:  Oct 21, 2016 
 * 
 * @brief  
 * 
 * 
 */

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstring>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "base_load_model.hpp"
#include "acmotor.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::AcmotorLoad::AcmotorLoad(void)
{ // with zero -> predictor step, without zero -> corrector
  bdebugprint = false;
  
  volt_measured0 = 0.0; // contractor and UV relays, voltage into the 1-ph A/C model is measured value with a delay of Tv, volt_measured = 1/(1+Tv*s) * vt
  freq_measured0 = 1.0; // for relay  models, frequency into the 1-ph A/C model is measured value with a delay of Tf, freq_measured = 1/(1+Tf*s) * freq
  temperatureA0 = 0.0; // for thermal protection, it is theta in thermal relay model
  temperatureB0 = 0.0; // for thermal protection, it is theta in thermal relay model
  volt_measured = 0.0; // contractor and UV relays
  freq_measured = 1.0; // for relay  models
  temperatureA = 0.0; // for thermal protection
  temperatureB = 0.0; // for thermal protection

  dv_dt0 = 0.0;
  dfreq_dt0 = 0.0;
  dThA_dt0 =0.0;
  dThB_dt0 =0.0;
  dv_dt = 0.0;
  dfreq_dt = 0.0;
  dThA_dt =0.0;
  dThB_dt =0.0;

  Tv = 0.02;
  Tf = 0.05;
  CompLF = 0.8;
  CompPF = 0.97;
  Vstall = 0.60;
  Rstall = 0.124;
  Xstall = 0.114;
  Tstall = 0.033;
  LFadj = 0.0;

  Kp1 = 0.0;
  Np1 = 1.0;
  Kq1 = 6.0;
  Nq1 = 2.0;
  Kp2 = 12.0;
  Np2 = 3.2;
  Kq2 = 11.0;
  Nq2 = 2.5;

  Vbrk = 0.86;
  Frst = 0.0;  
  Vrst = 0.9;
  Trst = 0.4;
  CmpKpf = 1.0;
  CmpKqf = -3.3;
  Vc1off = 0.45;
  Vc2off = 0.35;
  Vc1on = 0.50;
  Vc2on = 0.40;
  Tth = 10.0;
  Th1t = 1.3;
  Th2t = 4.3;
  Fuvr = 0.0; // fraction of AC motors equiped with UV relay
  Uvtr1 = 0.8;
  Ttr1 = 0.2;
  Uvtr2 = 0.9;
  Ttr2 = 5.0;
 
  volt  = 0.0; // terminal voltage
  freq = 1.0; // delta freq = freq - 1.0

  MVABase = 0.0;
  Pinit = 0.0; // unit in MW
  Qinit = 0.0; // unit in MW
  Pinit_pu = 0; // pu in Motor base
  Qinit_pu = 0; // pu in Motor base
  P0 = 0.0; // pu in Motor base
  Q0 = 0.0; // pu in Motor base
  PA = 0.0;
  QA = 0.0;
  PB = 0.0;
  QB = 0.0;
  Pmotor = 0.0;
  Qmotor = 0.0;
        
  Vstallbrk = 0.0;
        
  Gstall = 0.0;
  Bstall = 0.0;
  equivY =0.0; // motor base
  equivY_sysMVA =0.0; // system mva base
  
       
  INorton_sysMVA = 0.0; // system base
       
  Imotor_init = 0.0;
       
  statusA = 1; // 1- normal running, 0 stalling
  statusB = 1;
        
  stallTimer = 0.0; // timer for AC stalling, ac stalls if stallTimmer > Tstall
        
  restartTimer = 0.0;

  FthA = 1.0; // fraction of non-restarting part of load not tripped by thermal protection
  FthB = 1.0; // fraction of restarting part of load not tripped by thermal protection

  TthA =0.0; // temperature in non-restarting part of load
  TthB =0.0; // temperature in restarting part of load

  // thermal relay linear tripping equation
  thEqnA = 0.0;
  thEqnB = 0.0;

  // fraction coefficient of on-line AC motors for the UV relay and contractor
  Kuv = 1.0;
  Kcon = 1.0;
  fcon_trip = 0.0;
  isContractorActioned = 0;
       
  UVTimer1 = 0.0;
  UVTimer2 = 0.0;

  // time-varying equivalent admittance for network interfacing and calculating the current injection
  // it is calucated from the performance curve at the previous time step
  equivYpq_motorBase = 0.0;
        
  I_conv_factor_M2S = 0.0;
  
  Fonline = 1.0;
  acload_perc = 0.43; //ac motor load percentage to the total load, default 0.43;
  samebus_static_loadP = 0.0;
  samebus_static_loadQ = 0.0;
  samebus_static_equivY_sysMVA = gridpack::ComplexType(0.0, 0.0);
  
  

}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::AcmotorLoad::~AcmotorLoad(void)
{
}

/**
 * Load parameters from DataCollection object into load model
 * @param data collection of load parameters from input files
 * @param index of load on bus
 */
void gridpack::dynamic_simulation::AcmotorLoad::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx, double dloadP, double dloadQ, int ibCMPL)
{
  p_sbase = 100.0;

  //check with Qiuhua to see whether pl and ql are loaded correctly, is there a percentage?? renke??
  data->getValue(BUS_NUMBER,&p_bus_id);
  data->getValue(LOAD_ID,&p_loadid,idx);
  //if (!data->getValue(LOAD_PL, &p_pl, idx)) p_pl = 0.0;
  //if (!data->getValue(LOAD_QL, &p_ql, idx)) p_ql = 0.0;
  
  // set the information of dynamic load P, Q, id at the base class level too.
  p_pl = dloadP;
  p_ql = dloadQ;
  setDynLoadP(p_pl);
  setDynLoadQ(p_ql);
  setDynLoadID(p_loadid);
  
  if (bdebugprint) printf("AcmotorLoad::load: p_pl: %12.6f, p_ql: %12.6f \n", p_pl, p_ql);
  	
  if (!data->getValue(LOAD_TSTALL, &Tstall, idx))  		Tstall = 0.033;
  
  if ( ibCMPL==1 ){
	  if (!data->getValue(LOAD_TRST, &Trst 	, idx))  	Trst 	= 0.4;
  }else{
	  if (!data->getValue(LOAD_TRESTART, &Trst 	, idx))  	Trst 	= 0.4;
  }
  
  if (!data->getValue(LOAD_TV, &Tv 	, idx))  		Tv 	= 0.02;
  if (!data->getValue(LOAD_TF, &Tf 	, idx))  		Tf 	= 0.05;
  
  if ( ibCMPL==1 ){
	  if (!data->getValue(LOAD_LFM, &CompLF 	, idx))  	CompLF = 0.8;
  }else{
	  if (!data->getValue(LOAD_COMPLF, &CompLF, idx))  CompLF = 0.8;  
  }
  
  if (!data->getValue(LOAD_COMPPF, &CompPF, idx))  CompPF = 0.97;
  if (!data->getValue(LOAD_VSTALL, &Vstall, idx))  Vstall = 0.60;
  if (!data->getValue(LOAD_RSTALL, &Rstall, idx))  Rstall = 0.124;
  if (!data->getValue(LOAD_XSTALL, &Xstall, idx))  Xstall = 0.114;
  if (!data->getValue(LOAD_LFADJ, &LFadj , idx))  LFadj 	= 0.0;

  if (!data->getValue(LOAD_KP1, &Kp1, idx))  Kp1 = 0.0;
  if (!data->getValue(LOAD_NP1, &Np1, idx))  Np1 = 1.0;
  if (!data->getValue(LOAD_KQ1, &Kq1, idx))  Kq1 = 6.0;
  if (!data->getValue(LOAD_NQ1, &Nq1, idx))  Nq1 = 2.0;
  if (!data->getValue(LOAD_KP2, &Kp2, idx))  Kp2 = 12.0;
  if (!data->getValue(LOAD_NP2, &Np2, idx))  Np2 = 3.2;
  if (!data->getValue(LOAD_KQ2, &Kq2, idx))  Kq2 = 11.0;
  if (!data->getValue(LOAD_NQ2, &Nq2, idx))  Nq2 = 2.5;

  if (!data->getValue(LOAD_VBRK, &Vbrk 	, idx))  Vbrk 	= 0.86;
  if (!data->getValue(LOAD_FRST, &Frst 	, idx))  Frst 	= 0.0;  
  if (!data->getValue(LOAD_VRST, &Vrst 	, idx))  Vrst 	= 0.9;
  if (!data->getValue(LOAD_CMPKPF, &CmpKpf, idx))  CmpKpf = 1.0;
  if (!data->getValue(LOAD_CMPKQF, &CmpKqf, idx))  CmpKqf = -3.3;
  if (!data->getValue(LOAD_VC1OFF, &Vc1off, idx))  Vc1off = 0.45;
  if (!data->getValue(LOAD_VC2OFF, &Vc2off, idx))  Vc2off = 0.35;

  if (!data->getValue(LOAD_VC1ON, &Vc1on, idx))  Vc1on 	= 0.50;
  if (!data->getValue(LOAD_VC2ON, &Vc2on, idx))  Vc2on 	= 0.40;
  if (!data->getValue(LOAD_TTH, &Tth	, idx))  Tth	= 10.0;
  if (!data->getValue(LOAD_TH1T, &Th1t , idx))  Th1t 	= 1.3;  
  if (!data->getValue(LOAD_TH2T, &Th2t , idx))  Th2t 	= 4.3;  
  if (!data->getValue(LOAD_FUVR, &Fuvr , idx))  Fuvr 	= 0.0; // fraction of AC motors equiped with UV relay
  if (!data->getValue(LOAD_UVTR1, &Uvtr1, idx))  Uvtr1 	= 0.8;
  if (!data->getValue(LOAD_TTR1, &Ttr1 , idx))  Ttr1 	= 0.2;
  if (!data->getValue(LOAD_UVTR2, &Uvtr2, idx))  Uvtr2 	= 0.9;
  if (!data->getValue(LOAD_TTR2, &Ttr2 , idx))  Ttr2 	= 5.0;
 
    if (bdebugprint) printf("Tstall %f, Trst  %f, Tv %f, Tf %f, CompLF %f, CompPF %f, Vstall %f, Rstall %f, Xstall %f, LFadj %f \n", Tstall, Trst, Tv, Tf, CompLF, CompPF, Vstall, Rstall, Xstall, LFadj);
    if (bdebugprint) printf("Kp1 %f, Np1 %f, Kq1 %f, Nq1 %f, Kp2 %f, Np2 %f, Kq2 %f, Nq2 %f \n", Kp1, Np1, Kq1, Nq1, Kp2, Np2, Kq2, Nq2); 
    if (bdebugprint) printf ("Vbrk %f, Frst %f, Vrst %f, CmpKpf %f, CmpKqf %f, Vc1off %f, Vc2off %f \n", Vbrk, Frst, Vrst, CmpKpf, CmpKqf, Vc1off, Vc2off);
    if (bdebugprint) printf ("Vc1on  %f, Vc2on %f, Tth %f, Th1t %f, Th2t %f, Fuvr %f, Uvtr1 %f, Ttr1 %f, Uvtr2 %f, Ttr2 %f \n", Vc1on, Vc2on, Tth, Th1t, Th2t, Fuvr, Uvtr1, Ttr1, Uvtr2, Ttr2);
  
  
  // set the information of dynamic load P, Q, id at the base class level too.
  
  setDynLoadP(p_pl);
  setDynLoadQ(p_ql);
  setDynLoadID(p_loadid);
  
  if (bdebugprint) printf("AcmotorLoad::load: p_pl: %12.6f, p_ql: %12.6f \n", p_pl, p_ql);
  	
  if (!data->getValue(LOAD_TSTALL, &Tstall, idx))  		Tstall = 0.033;
  if (!data->getValue(LOAD_TRESTART, &Trst 	, idx))  	Trst 	= 0.4;
  if (!data->getValue(LOAD_TV, &Tv 	, idx))  		Tv 	= 0.02;
  if (!data->getValue(LOAD_TF, &Tf 	, idx))  		Tf 	= 0.05;
  if (!data->getValue(LOAD_COMPLF, &CompLF, idx))  CompLF = 0.8;
  if (!data->getValue(LOAD_COMPPF, &CompPF, idx))  CompPF = 0.97;
  if (!data->getValue(LOAD_VSTALL, &Vstall, idx))  Vstall = 0.60;
  if (!data->getValue(LOAD_RSTALL, &Rstall, idx))  Rstall = 0.124;
  if (!data->getValue(LOAD_XSTALL, &Xstall, idx))  Xstall = 0.114;
  if (!data->getValue(LOAD_LFADJ, &LFadj , idx))  LFadj 	= 0.0;

  if (!data->getValue(LOAD_KP1, &Kp1, idx))  Kp1 = 0.0;
  if (!data->getValue(LOAD_NP1, &Np1, idx))  Np1 = 1.0;
  if (!data->getValue(LOAD_KQ1, &Kq1, idx))  Kq1 = 6.0;
  if (!data->getValue(LOAD_NQ1, &Nq1, idx))  Nq1 = 2.0;
  if (!data->getValue(LOAD_KP2, &Kp2, idx))  Kp2 = 12.0;
  if (!data->getValue(LOAD_NP2, &Np2, idx))  Np2 = 3.2;
  if (!data->getValue(LOAD_KQ2, &Kq2, idx))  Kq2 = 11.0;
  if (!data->getValue(LOAD_NQ2, &Nq2, idx))  Nq2 = 2.5;

  if (!data->getValue(LOAD_VBRK, &Vbrk 	, idx))  Vbrk 	= 0.86;
  if (!data->getValue(LOAD_FRST, &Frst 	, idx))  Frst 	= 0.0;  
  if (!data->getValue(LOAD_VRST, &Vrst 	, idx))  Vrst 	= 0.9;
  if (!data->getValue(LOAD_CMPKPF, &CmpKpf, idx))  CmpKpf = 1.0;
  if (!data->getValue(LOAD_CMPKQF, &CmpKqf, idx))  CmpKqf = -3.3;
  if (!data->getValue(LOAD_VC1OFF, &Vc1off, idx))  Vc1off = 0.45;
  if (!data->getValue(LOAD_VC2OFF, &Vc2off, idx))  Vc2off = 0.35;

  if (!data->getValue(LOAD_VC1ON, &Vc1on, idx))  Vc1on 	= 0.50;
  if (!data->getValue(LOAD_VC2ON, &Vc2on, idx))  Vc2on 	= 0.40;
  if (!data->getValue(LOAD_TTH, &Tth	, idx))  Tth	= 10.0;
  if (!data->getValue(LOAD_TH1T, &Th1t , idx))  Th1t 	= 1.3;
  if (!data->getValue(LOAD_TH2T, &Th2t , idx))  Th2t 	= 4.3;
  if (!data->getValue(LOAD_FUVR, &Fuvr , idx))  Fuvr 	= 0.0; // fraction of AC motors equiped with UV relay
  if (!data->getValue(LOAD_UVTR1, &Uvtr1, idx))  Uvtr1 	= 0.8;
  if (!data->getValue(LOAD_TTR1, &Ttr1 , idx))  Ttr1 	= 0.2;
  if (!data->getValue(LOAD_UVTR2, &Uvtr2, idx))  Uvtr2 	= 0.9;
  if (!data->getValue(LOAD_TTR2, &Ttr2 , idx))  Ttr2 	= 5.0;
  if (!data->getValue(LOAD_AC_PERC, &acload_perc, idx)) acload_perc = 0.43;
  //printf("-------rkde,  ac motor acload_perc: %.2f\n", acload_perc);
  
  if (acload_perc< 0.01) acload_perc = 0.43; //if the LOAD_AC_PERC is not defined in the dyr file, as the normal PSS/E dyr file
  //printf("-------rkde after mod,  ac motor acload_perc: %.2f\n", acload_perc);
  
    if (bdebugprint) printf("Tstall %f, Trst  %f, Tv %f, Tf %f, CompLF %f, CompPF %f, Vstall %f, Rstall %f, Xstall %f, LFadj %f \n", Tstall, Trst, Tv, Tf, CompLF, CompPF, Vstall, Rstall, Xstall, LFadj);
    if (bdebugprint) printf("Kp1 %f, Np1 %f, Kq1 %f, Nq1 %f, Kp2 %f, Np2 %f, Kq2 %f, Nq2 %f \n", Kp1, Np1, Kq1, Nq1, Kp2, Np2, Kq2, Nq2); 
    if (bdebugprint) printf ("Vbrk %f, Frst %f, Vrst %f, CmpKpf %f, CmpKqf %f, Vc1off %f, Vc2off %f \n", Vbrk, Frst, Vrst, CmpKpf, CmpKqf, Vc1off, Vc2off);
    if (bdebugprint) printf ("Vc1on  %f, Vc2on %f, Tth %f, Th1t %f, Th2t %f, Fuvr %f, Uvtr1 %f, Ttr1 %f, Uvtr2 %f, Ttr2 %f \n", Vc1on, Vc2on, Tth, Th1t, Th2t, Fuvr, Uvtr1, Ttr1, Uvtr2, Ttr2);
  
}

/**
 * Update parameters in DataCollection object with current values from
 * load
 * @param data collection object for bus that hosts load
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::AcmotorLoad::updateData(
    boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  double pl = getDynLoadP();
  double ql = getDynLoadQ();
  if (!data->setValue(LOAD_PL_CURRENT, pl, idx)) {
    data->addValue(LOAD_PL_CURRENT, pl, idx);
  }
  if (!data->setValue(LOAD_QL_CURRENT, ql, idx)) {
    data->addValue(LOAD_QL_CURRENT, ql, idx);
  }

}

/**
 * Initialize load model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 */
void gridpack::dynamic_simulation::AcmotorLoad::init(double mag,
    double ang, double ts)
{
  //SJin: What are parameters systemMVABase, PintMW, and vt? From where to get them?
  //Fake declaration; 
  double systemMVABase, PintMW, vt;
  
  PintMW = p_pl;
  systemMVABase = 100.0;
  vt = mag;
  
  presentMag = mag;
  vt_complex = gridpack::ComplexType(mag*cos(ang), mag*sin(ang)); 

  if (CompLF ==0.0)
    CompLF = 1.0;

  // calculate the mva base
  Pinit = PintMW * acload_perc; // Yuan. please consider the ac motor/ static load percentage here
  MVABase = Pinit / CompLF;
  
  I_conv_factor_M2S = MVABase / systemMVABase;  // factor for Converting from motor to system base

  Gstall = Rstall / (Rstall * Rstall + Xstall * Xstall);
  Bstall = -Xstall / (Rstall * Rstall + Xstall * Xstall);

  equivY = gridpack::ComplexType(Gstall, Bstall); 

  equivY_sysMVA = equivY * MVABase / systemMVABase;  // Yuan, please modify the equivY_sysMVA based on the static load
  
  samebus_static_equivY_sysMVA = gridpack::ComplexType(0.0, 0.0);//Yuan, please modify this equivY to be the value of the corresponding 30% static load at the system MVA
  
  //printf("AcmotorLoad::init: equivY_sysMVA: %12.6f + j %12.6f \n", real(equivY_sysMVA), imag(equivY_sysMVA));

  // initial P and Q in motor MVA base

  Pinit_pu = Pinit / MVABase;
  Qinit_pu = Pinit_pu * tan(acos(CompPF));
  
  //printf("AcmotorLoad::init: Pinit_pu: %12.6f, Qinit_pu: %12.6f \n", Pinit_pu, Qinit_pu);

  gridpack::ComplexType tmp(Pinit_pu, -Qinit_pu);
  equivYpq_motorBase = tmp / vt / vt; 

  // initial temparitureA/B
  Imotor_init = abs(tmp) / vt;

  temperatureA =  Imotor_init * Imotor_init * Rstall;
  temperatureB =  temperatureA;

  volt = vt;
            
  volt_measured = vt;
            
  // Apply adjustments
  Vstall = Vstall * (1.0 + LFadj * (CompLF - 1.0));
  Vbrk = Vbrk * (1.0 + LFadj * (CompLF -1.0));
  
  //printf("AcmotorLoad::init: Vstall: %12.6f, Vbrk: %12.6f \n", Vstall, Vbrk);

  // calculate P0 and Q0 for the algebraic P/Q curve

  P0 = Pinit_pu - Kp1 * pow(volt - Vbrk, Np1);
  Q0 =  Qinit_pu - Kq1 * pow(volt- Vbrk, Nq1);
  
  //printf("AcmotorLoad::init: P0: %12.6f, Q0: %12.6f \n", P0, Q0);

  double v = 0.4;
  while (v <= Vbrk) {
    double pst = Gstall * v * v; //renke??
    double p_comp = P0 + Kp1 * pow(v - Vbrk, Np1); // original
	//double p_comp = P0 + Kp2 * pow(v - Vbrk, Np2);  // Yuan changed on 20200709

    if (p_comp < pst) {
      Vstallbrk = v;
      break;
    }
    v += 0.001;
  }

  if (Vstallbrk > Vstall)
    Vstallbrk = Vstall;

  // initialize the thermal tripping equation

  if (Th1t >0.0 && Th2t > Th1t) {
    thEqnA = -1.0/(Th2t - Th1t);
    thEqnB = Th2t/(Th2t - Th1t); // for the thermal relay model, if Th1t<theta<Th2t, fth = 1.0/(Th2t - Th1t) * (Th2t - theta) = thEqnA * theta + thEqnB
  }

  // check parameter Vc2off and Vc1off, Vc2off should be no larger than Vc1ff
  if (Vc2off > 0.0 && Vc2off > Vc1off) 
    Vc2off  = 0.0; // fix bad data

  // Vc2on should be no larger than Vc1on
  if (Vc2on > 0.0 && Vc2on > Vc1on)
    Vc2on  = Vc1on - 0.1; // fix bad data

  volt_measured0 = volt_measured; // contracttor and UV relays
  freq_measured0 = freq_measured; // for relay  models
  temperatureA0 = temperatureA; // for thermal protection
  temperatureB0 = temperatureB; // for thermal protection
  
  if (bdebugprint) printf("AcmotorLoad::Init(), I_conv_factor_M2S: %12.6f, MVABase: %12.6f  \n", I_conv_factor_M2S, MVABase);
  if (bdebugprint) printf("AcmotorLoad::Init(), equivY_sysMVA: %12.6f + j %12.6f  \n", real(equivY_sysMVA), imag(equivY_sysMVA) );
  if (bdebugprint) printf("AcmotorLoad::Init(), equivYpq_motorBase: %12.6f + j %12.6f  \n", real(equivYpq_motorBase), imag(equivYpq_motorBase) );
  if (bdebugprint) printf("AcmotorLoad::Init(), Pinit_pu: %12.6f, Qinit_pu: %12.6f, Imotor_init: %12.6f,  temperatureA: %12.6f\n", Pinit_pu, Qinit_pu, Imotor_init, temperatureA);
  if (bdebugprint) printf("AcmotorLoad::Init(), volt: %12.6f, volt_measured: %12.6f, P0: %12.6f, Q0: %12.6f  \n", volt, volt_measured, P0, Q0);
  if (bdebugprint) printf("AcmotorLoad::Init(), Vstall: %12.6f, Vbrk: %12.6f, Vstallbrk: %12.6f  \n", Vstall, Vbrk, Vstallbrk);
  if (bdebugprint) printf("AcmotorLoad::Init(), thEqnA: %12.6f, thEqnB: %12.6f, Vc2off: %12.6f, Vc2on: %12.6f  \n", thEqnA, thEqnB, Vc2off, Vc2on);
  
  setDynLoadP(Pinit);
  setDynLoadQ(getInitReactivePower()); // original, 20200714
  
  Fonline = 1.0;
  //acload_perc = 0.7;
  
}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType gridpack::dynamic_simulation::AcmotorLoad::INorton()
{
  //SJIN: Matlab getINorton?
  return p_INorton;
}


/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType gridpack::dynamic_simulation::AcmotorLoad::NortonImpedence()
{
  //SJIN: Matlab getNortonImpedance?
  //return ??; // refer to gensal.cpp
  double systemMVABase, PintMW, vt;
  
  PintMW = p_pl;
  systemMVABase = 100.0;

  if (CompLF ==0.0)
    CompLF = 1.0;

  // calculate the mva base
  Pinit = PintMW * acload_perc; //SJin: what is PintMW
  MVABase = Pinit / CompLF;
  I_conv_factor_M2S = MVABase / systemMVABase;  // factor for Converting from motor to system base

  Gstall = Rstall / (Rstall * Rstall + Xstall * Xstall);
  Bstall = -Xstall / (Rstall * Rstall + Xstall * Xstall);

  equivY = gridpack::ComplexType(Gstall, Bstall); 

  equivY_sysMVA = equivY * MVABase / systemMVABase; 
  
  return equivY_sysMVA;
  
}

/**
 * Predict part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::AcmotorLoad::predictor_currentInjection(bool flag)
{
  if (!flag) {
    volt_measured0 = volt_measured;
    freq_measured0 = freq_measured;
    temperatureA0 = temperatureA;
    temperatureB0 = temperatureB;
  }
  gridpack::ComplexType Imotor_motorBase = equivYpq_motorBase * vt_complex;
  
  if (bdebugprint) printf("AcmotorLoad::predictor_currentInjection, equivY_sysMVA: %12.6f +j %12.6f\n", real(equivY_sysMVA), imag(equivY_sysMVA));
  if (bdebugprint) printf("AcmotorLoad::predictor_currentInjection, vt_complex: %12.6f +j %12.6f\n", real(vt_complex), imag(vt_complex));
  if (bdebugprint) printf("AcmotorLoad::predictor_currentInjection, Imotor_motorBase: %12.6f +j %12.6f\n", real(Imotor_motorBase), imag(Imotor_motorBase));
  if (bdebugprint) printf("AcmotorLoad::predictor_currentInjection, I_conv_factor_M2S: %12.6f \n", I_conv_factor_M2S);
  
  //INorton_sysMVA = equivY_sysMVA * vt_complex - Imotor_motorBase * I_conv_factor_M2S; //original one without considering loadshedding of the same bus static load
  INorton_sysMVA = equivY_sysMVA * vt_complex - Imotor_motorBase * I_conv_factor_M2S + samebus_static_equivY_sysMVA*vt_complex*(1.0-Fonline); // original
  p_INorton = INorton_sysMVA; // SJin: Correct? 
} 

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::AcmotorLoad::predictor(
    double t_inc, bool flag)
{
  //SJin: What are parameters vt, freq, and dt? From where to get them?
  //Fake declaration; 
  double vt, freq;
  double dt = t_inc;
  vt = presentMag;
  freq = presentFreq;
         
  // step-1, we need to update the state variable
  if (!flag) {
    volt_measured0 = volt_measured;
    freq_measured0 = freq_measured;
    temperatureA0 = temperatureA;
    temperatureB0 = temperatureB;
  }
         
  // step -2 dx/dt
  dv_dt0 = (vt - volt_measured0) / Tv; // only used for relay and protection, not in dynamic calculation, in dv_dt0, v represents volt_measured. This is the dyn eqns of voltage sensor.
  dfreq_dt0 = (freq - freq_measured0) / Tf; // only used for relay and protection, not in dynamic calculation. This is the dyn eqns of frequency sensor.
  // the thermal proection temperature is only measured when the AC is stalled
  //double i = sqrt(-1);
  gridpack::ComplexType ImotorA_pu(0.0, 0.0);
  if (statusA == 0)
    ImotorA_pu = vt * equivY ;
  else {
    gridpack::ComplexType tmp(PA, -QA);
    ImotorA_pu = tmp / vt;
  }
             
  dThA_dt0 = (abs(ImotorA_pu) * abs(ImotorA_pu) * Rstall - temperatureA0) / Tth;
            
  gridpack::ComplexType ImotorB_pu(0.0, 0.0);
  if (statusB == 0)
    ImotorB_pu = vt * equivY ;
  else {
    gridpack::ComplexType tmp(PB, -QB);
    ImotorA_pu = tmp / vt;
  }
               
  dThB_dt0 = (abs(ImotorB_pu) * abs(ImotorB_pu) * Rstall - temperatureB0) / Tth;
         
  // step-3 integration 
  volt_measured = volt_measured0 + dv_dt0 * dt;
  freq_measured = freq_measured0 + dfreq_dt0 * dt;
  temperatureA = temperatureA0 +  dThA_dt0 * dt;
  temperatureB = temperatureB0 +  dThB_dt0 * dt;
}

/**
 * Corrector part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::AcmotorLoad::corrector_currentInjection(bool flag)
{
  gridpack::ComplexType Imotor_motorBase = equivYpq_motorBase * vt_complex;
  
  //INorton_sysMVA = equivY_sysMVA * vt_complex - Imotor_motorBase * I_conv_factor_M2S; //original one without considering loadshedding of the same bus static load
  INorton_sysMVA = equivY_sysMVA * vt_complex - Imotor_motorBase * I_conv_factor_M2S + samebus_static_equivY_sysMVA*vt_complex*(1.0-Fonline); // original
  p_INorton = INorton_sysMVA; // SJin: Correct? 
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::AcmotorLoad::corrector(
    double t_inc, bool flag)
{
  //SJin: What are parameters vt, freq, and dt? From where to get them?
  //Fake declaration; 
  double vt, freq;
  double dt = t_inc;
  vt = presentMag;
  freq = presentFreq;
        
  // step -1 dx'/dt
         
  dv_dt = (vt - volt_measured) / Tv ; // only used for relay and protection, not in dynamic calculation
  dfreq_dt = (freq - freq_measured) / Tf; // only used for relay and protection, not in dynamic calculation
  // the thermal proection temperature is only measured when the AC is stalled
  //double i = sqrt(-1);
  gridpack::ComplexType ImotorA_pu(0.0, 0.0);
  if (statusA == 0)
    ImotorA_pu = vt * equivY ;
  else {
    gridpack::ComplexType tmp(PA, -QA);
    ImotorA_pu = tmp / vt;
  }
             
  dThA_dt = (abs(ImotorA_pu) * abs(ImotorA_pu) * Rstall - temperatureA) / Tth;
            
  gridpack::ComplexType ImotorB_pu(0.0, 0.0);
  if (statusB == 0)
    ImotorB_pu = vt * equivY ;
  else {
    gridpack::ComplexType tmp(PB, -QB);
    ImotorA_pu = tmp / vt;
  }
               
  dThB_dt = (abs(ImotorB_pu) * abs(ImotorB_pu) * Rstall -temperatureB) / Tth;
               
  // step-2 integration 
  volt_measured = volt_measured0 + 0.5 * (dv_dt0 + dv_dt) * dt;
  freq_measured = freq_measured0 + 0.5 * (dfreq_dt0 + dfreq_dt) * dt;
  temperatureA = temperatureA0 +  0.5 * (dThA_dt0 + dThA_dt) * dt;
  temperatureB = temperatureB0 +  0.5 * (dThB_dt0 + dThB_dt) * dt;
}

/**
 * post process for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::AcmotorLoad::dynamicload_post_process(
    double t_inc, bool flag)
{
  double vt, freq;
  double dt = t_inc;
  vt = presentMag;
  freq = presentFreq;
  
  // UV Relay
  // Two levels of undervoltage load shedding can be represented: If the voltage drops
  // below uvtr1 for ttr1 seconds, the fraction “fuvtr” of the load is tripped; If the voltage drops below uvtr2
  // for ttr2 seconds, the fraction “fuvtr” of the load is tripped

  if (volt_measured < Uvtr1 && Fuvr > 0.0){
	UVTimer1 = UVTimer1 + dt;
  }else {
	 UVTimer1 = 0.0; 
  }
  
  if ( volt_measured <  Uvtr2 &&  Fuvr > 0.0){
	 UVTimer2 =  UVTimer2 + dt;
  }else{
	 UVTimer2 = 0.0;  
  }
          
  if ( UVTimer1 >  Ttr1){
	 Kuv = 1.0 -  Fuvr;  
  }
      
  if ( UVTimer2 >  Ttr2){
	 Kuv = 1.0 -  Fuvr;  
  }
  
  // contractor
  // Contactor – If the voltage drops to below Vc2off, all of the load is tripped; if the voltage is between
  // Vc1off and Vc2off, a linear fraction of the load is tripped. If the voltage later recovers to above Vc1on, all
  // of the motor is reconnected; if the voltage recovers to between Vc2on and Vc1on, a linear fraction of the
  // load is reconnected.
  
  if ( volt_measured <  Vc2off){
     Kcon = 0.0;  // fraction not tripped by contactor
     fcon_trip = 1.0;
  }
  else if ( volt_measured >=  Vc2off &&  
		 volt_measured <  Vc1off){
	 Kcon =   ( volt_measured -  Vc2off)/( Vc1off -  Vc2off);
     fcon_trip = 1.0 -   Kcon;
  }  
  
  // !!! this part has issue, should recover the parts being tripped
  // AC reconnects after voltge recovery 
  if ( volt_measured >=   Vc1on){
     Kcon = 1.0;
  } else if  ( volt_measured <  Vc1on &&   volt_measured >=  Vc2on) {  
    double Frecv =  ( volt_measured -  Vc2on)/( Vc1on -  Vc2on);
     Kcon = 1.0 -  fcon_trip*(1.0- Frecv);  // Kcon = (1 - fcon_trip) + fcon_trip * Frecv, fraction not tripped by contactor
  } 
  // if (volt_measured > Vc1off and volt_measured < Vc2on) then No change on Kcon
  
  // update timer for AC thermal protection
  if  ( statusA == 1) {
     if (vt <  Vstall){
          stallTimer  =   stallTimer + dt;
	 }
     else {// reset timer
          stallTimer = 0.0;
     }
  }
  
  if  ( statusB == 0) {
      if (vt >  Vrst){
           restartTimer =  restartTimer + dt;
	  }
      else {// reset timer
           restartTimer =0.0;
      }
  }
  
  // update the status of the motor. transition from running to
  // stalling is the same for the equivalent Motor A and B
  if( stallTimer >  Tstall &&  statusA == 1){
        statusA = 0;
        statusB = 0;
  }
  
  
  // considering the AC restarting
  if ( Frst > 0.0 &&  statusB == 0 &&  restartTimer >  Trst){
        statusB = 1;
  } 
  
  // check whether AC motor will be trip next step, and what will be
  // the remaining fraction; 
         
  if ( thEqnA < 0.0){
         if ( temperatureA >  Th1t) {
             // it only trips under the stalled condition???
             if ( statusA == 0) {
                 FthA =  temperatureA* thEqnA +  thEqnB;
			 
				if ( FthA <0.0) {
                      FthA = 0.0;
				}
             }   
         } 
        
        if ( temperatureB >  Th1t) {
             // it only trips under the stalled condition???
             if ( statusB == 0) {
                  FthB =  temperatureB* thEqnA +  thEqnB;

                 if ( FthB <0.0) {
                      FthB = 0.0;
                 }
             }           
		}   
  } // end of if ( thEqnA < 0.0)
   
  // calculate the AC motor power 
  if (bdebugprint) printf("AcmotorLoad::dynamicload_post_process, P0: %12.6f, Q0: %12.6f \n", P0, Q0);
   if ( statusA == 1 ) {// MotorA running
       
           if (vt >=  Vbrk) {

                PA = ( P0 +  Kp1*pow(vt -  Vbrk, Np1) )*(1.0 +  CmpKpf*(freq - 1.0));  // delta_f = freq - 1 or 1 - freq?
                QA = ( Q0 +  Kq1*pow(vt -  Vbrk, Nq1) )*(1.0 +  CmpKqf*(freq - 1.0));
		   }
           else if (vt <  Vbrk && vt >  Vstallbrk) {

                PA = ( P0 +  Kp2*pow(Vbrk- vt, Np2) )*(1.0 +  CmpKpf*(freq - 1.0));
                QA = ( Q0 +  Kq2*pow(Vbrk -vt, Nq2) )*(1.0 +  CmpKqf*(freq - 1.0));
		   }
           else  {

                PA =  Gstall*vt*vt;
                QA = - Bstall*vt*vt; // motor oriented--draw power from system as positive 
           }
   } else {// MotorA stalled
        
            PA =  Gstall*vt*vt;
            QA = - Bstall*vt*vt;   // motor oriented--draw power from system as positive 
        
   } // end of if ( statusA == 1 ) {// MotorA running
    
    if ( Frst > 0.0  ) { // power is calcuated only when Frst > 0.0
        if ( statusB == 1) { // MotorB running

               if (vt >  Vbrk) {

                    PB = ( P0 +  Kp1*pow(vt -  Vbrk, Np1) )*(1.0 +  CmpKpf*(freq - 1.0));
                    QB = ( Q0 +  Kq1*pow(vt -  Vbrk, Nq1) )*(1.0 +  CmpKqf*(freq - 1.0));
			   }
               else if (vt <  Vbrk && vt >  Vstallbrk) {

                    PB = ( P0 +  Kp2*pow(Vbrk- vt, Np2) )*(1.0 +  CmpKpf*(freq - 1.0));
                    QB = ( Q0 +  Kq2*pow(Vbrk- vt, Nq2) )*(1.0 +  CmpKqf*(freq - 1.0));
			   }
               else  {

                    PB =  Gstall*vt*vt;
                    QB = - Bstall*vt*vt; // motor oriented--draw power from system as positive 
               }
		} else {// MotorA stalled

              PB =  Gstall*vt*vt;
              QB = - Bstall*vt*vt; // motor oriented--draw power from system as positive 

		}
    }
  
    // MOTOR A part -- non-restartable
    // MOTOR B part -- restartable
     Pmotor =  PA*(1.0- Frst)* FthA +  PB*  Frst* FthB;
     Qmotor =  QA*(1.0- Frst)* FthA +  QB*  Frst* FthB;
  
    // consider the UR relay and contractor
	 //Fonline = Kuv* Kcon* Fonline;
     Pmotor =  Kuv*Kcon*Fonline* Pmotor;
     Qmotor =  Kuv*Kcon*Fonline* Qmotor;
	 
	 gridpack::ComplexType tmpcplx(Pmotor, -Qmotor);
	 equivYpq_motorBase = tmpcplx/vt/vt;
     
    //equivYpq_motorBase = ( Pmotor - i* Qmotor)/vt/vt;

    if (bdebugprint) printf("dynamic load output step: ,%8d, %2s, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %8d, %8d, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f,\n",
          p_bus_id, p_loadid.c_str(), volt_measured, freq_measured, temperatureA, temperatureB, presentMag, presentFreq,
                  statusA, statusB, Pmotor, Qmotor, FthA, FthB, Fonline, Kuv, Kcon);

}

/**
 * Set voltage on each load
 */
void gridpack::dynamic_simulation::AcmotorLoad::setVoltage(
    gridpack::ComplexType voltage)
{
  if (bdebugprint) printf("AcmotorLoad::setVoltage, %12.6f + j %12.6f \n", real(voltage), imag(voltage));	
  presentMag = abs(voltage);
  presentAng = atan2(imag(voltage), real(voltage));  
  vt_complex = voltage;
  
}

/**
 * Set terminal voltage frequency on each load
 */
void gridpack::dynamic_simulation::AcmotorLoad::setFreq(double dFreq)
{
  presentFreq = dFreq;
}

/**
 * get intialized reactive power of the dynamic load model
 */
double gridpack::dynamic_simulation::AcmotorLoad::getInitReactivePower(void)
{
	return Qinit_pu*MVABase;
}

/**
 * get the variable Fonline 
 */
double gridpack::dynamic_simulation::AcmotorLoad::getFonline(void)
{
	return Fonline;
}

/** 
 * Get the value of the field voltage parameter
 * @return value of field voltage
 */
/*double gridpack::dynamic_simulation::AcmotorLoad::getFieldVoltage()
{
  return Efd;
}*/

/**
 * Write out load state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
bool gridpack::dynamic_simulation::AcmotorLoad::serialWrite(
    char* string, const int bufsize, const char *signal)
{
    string[0] = '\0';
    if (!strcmp(signal,"standard")) {
    //sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f\n",
    //    p_bus_id,p_ckt.c_str(),real(p_mac_ang_s1),real(p_mac_spd_s1),real(p_mech),
    //    real(p_pelect));
    sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f\n",
          p_bus_id, p_loadid.c_str(), volt_measured, freq_measured, temperatureA, temperatureB);
    return true;
  } else if (!strcmp(signal,"init_debug")) {
    sprintf(string," %8d  %2s Something\n",p_bus_id,p_loadid.c_str());
    return true;
  } else if (!strcmp(signal,"debug_initial")) {
    sprintf(string,"Bus: %d PLOAD: %f QLOAD: %f\n",p_bus_id,p_pl,p_ql);
    return true;
  } else if (!strcmp(signal,"load_watch")) {
    if (getWatch()) {
      char buf[128];
//    sprintf(buf,", %f, %f",real(p_mac_ang_s1),real(p_mac_spd_s1));
      sprintf(string,",%8d, %2s, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %8d, %8d, %12.6f, %12.6f, %12.6f, %12.6f, \n",
          p_bus_id, p_loadid.c_str(), volt_measured, freq_measured, temperatureA, temperatureB, presentMag, presentFreq, 
		  statusA, statusB, Pmotor, Qmotor, FthA, FthB);
      return true;
/*      if (strlen(buf) <= bufsize) {
        sprintf(string,"%s",buf);
        return true; 
      } else {
        printf ("size watch problem at gensal serialWrite() \n");
        return false;
      }
    } else {
        return false;
*/
    }
  } else if (!strcmp(signal,"debug_initial")) {
    return false;
  }	
		
  return false;
}


/**
 * return true if load change is enabled
 * @param percentageFactor: the fraction (percentage) of load that is changed. Negative: load reduction, Positive: load increase
 */
bool gridpack::dynamic_simulation::AcmotorLoad::changeLoad(double percentageFactor)
{
	if (percentageFactor < -1.0) {
		if (bdebugprint) printf("percentageFactor < -1.0, this change will not be applied.  \n");
		return false;
	}
	
	if (Fonline == 0.0){
		if (bdebugprint) printf("-------!!!! gridpack warning: the dynamic load at bus %d with ID %s has 0 percent of load and could not be shedding the percentange of %f !!!! \n", 
		                                              p_bus_id, p_loadid.c_str(), percentageFactor);
		//percentageFactor = 0.0;
		//Fonline = Fonline + percentageFactor;
		if (bdebugprint) printf("----renke debug load shed, AcmotorLoad::changeLoad, Fonline: %f \n", Fonline);
		return true;
	}
	
	if ( (Fonline + percentageFactor) < 0.0 ){
		if (bdebugprint) printf("-------!!!! gridpack warning: the dynamic load at bus %d with ID %s has %f percent of load and could not be shedding the percentange of %f !!!! \n", 
		                                              p_bus_id, p_loadid.c_str(), Fonline, percentageFactor);
		Fonline = 0.0;
		if (bdebugprint) printf("----renke debug load shed, AcmotorLoad::changeLoad, Fonline: %f \n", Fonline);
		return true;
	}
	
	Fonline = Fonline + percentageFactor;
	// Yuan added below 20200709
	if (Fonline < 0.0) Fonline=0.0;
	// Yuan added above 20200709
	
	if (bdebugprint) printf("----renke debug load shed, AcmotorLoad::changeLoad, the dynamic load at bus %d with ID %s, percentageFactor: %f, remaining Fonline: %f \n", 
																				p_bus_id, p_loadid.c_str(), percentageFactor, Fonline);
	return true;
}
/**
 * Set same bus static load p and q for load shedding action usage
 */
void gridpack::dynamic_simulation::AcmotorLoad::setSameBusStaticLoadPQ(double static_pl, double static_ql, double mag)
{
  samebus_static_loadP = static_pl;
  samebus_static_loadQ = static_ql;
  
  double samebus_static_load_yr = samebus_static_loadP/(mag*mag);
  double samebus_static_load_yi = (-samebus_static_loadQ)/(mag*mag);
  samebus_static_equivY_sysMVA = samebus_static_equivY_sysMVA + gridpack::ComplexType(samebus_static_load_yr, samebus_static_load_yi);
  if (bdebugprint) printf("AcmotorLoad::setSameBusStaticLoadPQ, samebus_static_equivY_sysMVA: %12.6f +j %12.6f\n", real(samebus_static_equivY_sysMVA), imag(samebus_static_equivY_sysMVA));
}
