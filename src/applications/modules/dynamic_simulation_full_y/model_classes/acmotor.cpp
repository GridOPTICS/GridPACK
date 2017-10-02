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
{
  volt_measured0 = 0.0; // contracttor and UV relays
  freq_measured0 = 1.0; // for relay  models
  temperatureA0 = 0.0; // for thermal protection
  temperatureB0 = 0.0; // for thermal protection
  volt_measured = 0.0; // contracttor and UV relays
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
  
  printf("AcmotorLoad::load: p_pl: %12.6f, p_ql: %12.6f \n", p_pl, p_ql);
  	
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
 
 printf("Tstall %f, Trst  %f, Tv %f, Tf %f, CompLF %f, CompPF %f, Vstall %f, Rstall %f, Xstall %f, LFadj %f \n", Tstall, Trst, Tv, Tf, CompLF, CompPF, Vstall, Rstall, Xstall, LFadj);
    printf("Kp1 %f, Np1 %f, Kq1 %f, Nq1 %f, Kp2 %f, Np2 %f, Kq2 %f, Nq2 %f \n", Kp1, Np1, Kq1, Nq1, Kp2, Np2, Kq2, Nq2); 
    printf ("Vbrk %f, Frst %f, Vrst %f, CmpKpf %f, CmpKqf %f, Vc1off %f, Vc2off %f \n", Vbrk, Frst, Vrst, CmpKpf, CmpKqf, Vc1off, Vc2off);
    printf ("Vc1on  %f, Vc2on %f, Tth %f, Th1t %f, Th2t %f, Fuvr %f, Uvtr1 %f, Ttr1 %f, Uvtr2 %f, Ttr2 %f \n", Vc1on, Vc2on, Tth, Th1t, Th2t, Fuvr, Uvtr1, Ttr1, Uvtr2, Ttr2);
  
  
  // set the information of dynamic load P, Q, id at the base class level too.
  
  setDynLoadP(p_pl);
  setDynLoadQ(p_ql);
  setDynLoadID(p_loadid);
  
  printf("AcmotorLoad::load: p_pl: %12.6f, p_ql: %12.6f \n", p_pl, p_ql);
  	
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
 
 printf("Tstall %f, Trst  %f, Tv %f, Tf %f, CompLF %f, CompPF %f, Vstall %f, Rstall %f, Xstall %f, LFadj %f \n", Tstall, Trst, Tv, Tf, CompLF, CompPF, Vstall, Rstall, Xstall, LFadj);
    printf("Kp1 %f, Np1 %f, Kq1 %f, Nq1 %f, Kp2 %f, Np2 %f, Kq2 %f, Nq2 %f \n", Kp1, Np1, Kq1, Nq1, Kp2, Np2, Kq2, Nq2); 
    printf ("Vbrk %f, Frst %f, Vrst %f, CmpKpf %f, CmpKqf %f, Vc1off %f, Vc2off %f \n", Vbrk, Frst, Vrst, CmpKpf, CmpKqf, Vc1off, Vc2off);
    printf ("Vc1on  %f, Vc2on %f, Tth %f, Th1t %f, Th2t %f, Fuvr %f, Uvtr1 %f, Ttr1 %f, Uvtr2 %f, Ttr2 %f \n", Vc1on, Vc2on, Tth, Th1t, Th2t, Fuvr, Uvtr1, Ttr1, Uvtr2, Ttr2);
  
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

  if (CompLF ==0.0)
    CompLF = 1.0;

  // calculate the mva base
  Pinit = PintMW; //SJin: what is PintMW
  MVABase = Pinit / CompLF;
  I_conv_factor_M2S = MVABase / systemMVABase;  // factor for Converting from motor to system base

  Gstall = Rstall / (Rstall * Rstall + Xstall * Xstall);
  Bstall = -Xstall / (Rstall * Rstall + Xstall * Xstall);

  equivY = gridpack::ComplexType(Gstall, Bstall); 

  equivY_sysMVA = equivY * MVABase / systemMVABase;
  
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
    double p_comp = P0 + Kp1 * pow(v - Vbrk, Np1);

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
    thEqnB = Th2t/(Th2t - Th1t);
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
  
  printf("AcmotorLoad::Init(), I_conv_factor_M2S: %12.6f, MVABase: %12.6f  \n", I_conv_factor_M2S, MVABase);
  printf("AcmotorLoad::Init(), equivY_sysMVA: %12.6f + j %12.6f  \n", real(equivY_sysMVA), imag(equivY_sysMVA) );
  printf("AcmotorLoad::Init(), equivYpq_motorBase: %12.6f + j %12.6f  \n", real(equivYpq_motorBase), imag(equivYpq_motorBase) );
  printf("AcmotorLoad::Init(), Pinit_pu: %12.6f, Qinit_pu: %12.6f, Imotor_init: %12.6f,  temperatureA: %12.6f\n", Pinit_pu, Qinit_pu, Imotor_init, temperatureA);
  printf("AcmotorLoad::Init(), volt: %12.6f, volt_measured: %12.6f, P0: %12.6f, Q0: %12.6f  \n", volt, volt_measured, P0, Q0);
  printf("AcmotorLoad::Init(), Vstall: %12.6f, Vbrk: %12.6f, Vstallbrk: %12.6f  \n", Vstall, Vbrk, Vstallbrk);
  printf("AcmotorLoad::Init(), thEqnA: %12.6f, thEqnB: %12.6f, Vc2off: %12.6f, Vc2on: %12.6f  \n", thEqnA, thEqnB, Vc2off, Vc2on);
  
  setDynLoadQ(getInitReactivePower());
 
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
  Pinit = PintMW; //SJin: what is PintMW
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
  
  printf("AcmotorLoad::predictor_currentInjection, equivY_sysMVA: %12.6f +j %12.6f\n", real(equivY_sysMVA), imag(equivY_sysMVA));
  printf("AcmotorLoad::predictor_currentInjection, vt_complex: %12.6f +j %12.6f\n", real(vt_complex), imag(vt_complex));
  printf("AcmotorLoad::predictor_currentInjection, Imotor_motorBase: %12.6f +j %12.6f\n", real(Imotor_motorBase), imag(Imotor_motorBase));
  printf("AcmotorLoad::predictor_currentInjection, I_conv_factor_M2S: %12.6f \n", I_conv_factor_M2S);
  
  INorton_sysMVA = equivY_sysMVA * vt_complex - Imotor_motorBase * I_conv_factor_M2S;
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
  dv_dt0 = (vt - volt_measured0) / Tv; // only used for relay and protection, not in dynamic calculation
  dfreq_dt0 = (freq - freq_measured0) / Tf; // only used for relay and protection, not in dynamic calculation
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
  INorton_sysMVA = equivY_sysMVA * vt_complex - Imotor_motorBase * I_conv_factor_M2S;
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
  // Vc1off and Vc2off, a linear fraction of the load is tripped. If the voltage later recovers to above Vc2on, all
  // of the motor is reconnected; if the voltage recovers to between Vc2on and Vc1on, a linear fraction of the
  // load is reconnected.
  
  if ( volt_measured <  Vc2off){
     Kcon = 0.0;
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
     Kcon = 1.0 -  fcon_trip*(1.0- Frecv);
  } 
  
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
  printf("AcmotorLoad::dynamicload_post_process, P0: %12.6f, Q0: %12.6f \n", P0, Q0);
   if ( statusA == 1 ) {// MotorA running
       
           if (vt >=  Vbrk) {

                PA = ( P0 +  Kp1*pow(vt -  Vbrk, Np1) )*(1.0 +  CmpKpf*(freq - 1.0));
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
     Pmotor =  Kuv* Kcon* Pmotor;
     Qmotor =  Kuv* Kcon* Qmotor;
	 
	 gridpack::ComplexType tmpcplx(Pmotor, -Qmotor);
	 equivYpq_motorBase = tmpcplx/vt/vt;
     
    //equivYpq_motorBase = ( Pmotor - i* Qmotor)/vt/vt;

    printf("dynamic load output step: ,%8d, %2s, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %8d, %8d, %12.6f, %12.6f, %12.6f, %12.6f, \n",
          p_bus_id, p_loadid.c_str(), volt_measured, freq_measured, temperatureA, temperatureB, presentMag, presentFreq,
                  statusA, statusB, Pmotor, Qmotor, FthA, FthB);

}

/**
 * Set voltage on each load
 */
void gridpack::dynamic_simulation::AcmotorLoad::setVoltage(
    gridpack::ComplexType voltage)
{
  printf("AcmotorLoad::setVoltage, %12.6f + j %12.6f \n", real(voltage), imag(voltage));	
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
