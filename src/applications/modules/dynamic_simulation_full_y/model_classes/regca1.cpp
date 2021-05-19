/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   regca1.cpp
 * @author Renke Huang
 * @Last modified:   May. 16, 2021
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
#include "base_generator_model.hpp"
#include "regca1.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Regca1Generator::Regca1Generator(void)
{
    x1ip_0 = 0.0;
	x2iq_0 = 0.0;
	x3vmeas_0 = 0.0;
	
    x1ip_1 = 0.0;
	x2iq_1 = 0.0;
	x3vmeas_1 = 0.0;
	
    dx1ip_0 = 0.0;
	dx2iq_0 = 0.0;
	dx3vmeas_0 = 0.0;
	
    dx1ip_1 = 0.0;
	dx2iq_1 = 0.0;
	dx3vmeas_1 = 0.0;
	
	presentMag = 1.0;
	presentAng = 0.0;
	
	busfreq = 1.0;
	
	p_tripped = false;
	bmodel_debug = false;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Regca1Generator::~Regca1Generator(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::Regca1Generator::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  p_sbase = 100.0;

  data->getValue(BUS_NUMBER,&p_bus_id);
  data->getValue(GENERATOR_ID,&p_ckt,idx);
  if (!data->getValue(GENERATOR_PG, &p_pg,idx)) p_pg = 0.0;
  if (!data->getValue(GENERATOR_QG, &p_qg,idx)) p_qg = 0.0;
  if (!data->getValue(GENERATOR_STAT, &p_status,idx)) p_status = 0;
  
  if (bmodel_debug)
  {
     printf("Regca1Generator::load p_pg = %f, p_qg = %f\n", p_pg, p_qg);
  }
  p_pg *= p_sbase;
  p_qg *= p_sbase;

  if (!data->getValue(GENERATOR_MBASE, &MVABase, idx)) MVABase = 100.0; // MVABase
  if (!data->getValue(GENERATOR_REGCA_LVPLSW , &lvplsw, idx)) lvplsw = 0; 
  if (!data->getValue(GENERATOR_REGCA_TG  ,    &tg, idx))     tg = 0.02; 
  if (!data->getValue(GENERATOR_REGCA_RRPWR  , &rrpwr, idx))  rrpwr = 10.0; 
  if (!data->getValue(GENERATOR_REGCA_BRKPT  , &brkpt, idx))  brkpt = 0.9; 
  if (!data->getValue(GENERATOR_REGCA_ZEROX  , &zerox, idx))  zerox = 0.4; 
  if (!data->getValue(GENERATOR_REGCA_LVPL1  , &lvpl1, idx))  lvpl1 = 1.22; 
  if (!data->getValue(GENERATOR_REGCA_VOLIM  , &volim, idx))  volim = 1.2; 
  if (!data->getValue(GENERATOR_REGCA_LVPNT1 , &lvpnt1, idx)) lvpnt1 = 0.8; 
  if (!data->getValue(GENERATOR_REGCA_LVPNT0 , &lvpnt0, idx)) lvpnt0 = 0.4; 
  if (!data->getValue(GENERATOR_REGCA_LOLIM  , &lolim, idx))  lolim = -1.3; 
  if (!data->getValue(GENERATOR_REGCA_TFLTR  , &tfltr, idx))  tfltr = 0.02; 
  if (!data->getValue(GENERATOR_REGCA_KHV  ,   &khv, idx))    khv = 0.0; 
  if (!data->getValue(GENERATOR_REGCA_LQRMAX , &iqrmax, idx)) iqrmax = 999.0; 
  if (!data->getValue(GENERATOR_REGCA_LQRMIN , &iqrmin, idx)) iqrmin = -999.0; 
  if (!data->getValue(GENERATOR_REGCA_ACCEL  , &accel, idx))  accel = 0.7; 
  
  double tmp = sqrt(p_pg*p_pg +p_qg*p_qg);
  if ( tmp > MVABase) {
       //MVABase = tmp*1.3;
      //printf("-----------generator at bus %d  has P: %f, Q: %f, S: %f, MVABASE: %f  \n", p_bus_id, p_pg, p_qg, tmp, MVABase);
	  MVABase = tmp*1.2;
  }

  if (bmodel_debug){
	printf("\n--------Regca1Generator parameters: MVABase = %12.6f, tg = %12.6f, rrpwr = %12.6f, brkpt = %12.6f, zerox = %12.6f, lvpl1 = %12.6f, volim = %12.6f, lvpnt1 = %12.6f, lvpnt0 = %12.6f, lolim = %12.6f, tfltr = %12.6f, khv = %12.6f, iqrmax = %12.6f, iqrmin = %12.6f, accel = %12.6f  \n", MVABase, tg, rrpwr, brkpt, zerox, lvpl1, volim, lvpnt1, lvpnt0, lolim, tfltr, khv, iqrmax, iqrmin, accel);
  }
  
}

/**
 * Initialize generator model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 */
void gridpack::dynamic_simulation::Regca1Generator::init(double mag,
    double ang, double ts)
{
  //printf("Step0 gen%d mag = %f\n", p_bus_id, mag);
  Vterm = mag;
  presentMag = mag;
  Theta = ang;
  presentAng = ang;

  double P = p_pg / MVABase;
  double Q = p_qg / MVABase;
  genP = P;
  genQ = Q;
  busfreq = 1.0;
  
  if (bmodel_debug){
	printf("Regca1Generator p_pg = %f, p_qg = %f, MVABase = %f\n", p_pg, p_qg, MVABase);
	printf("Regca1Generator Vterm = %f, Theta = %f, P = %f, Q = %f\n", Vterm, Theta, P, Q);
  }
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  
  //ip = (P * Vrterm + Q * Viterm) / (Vterm * Vterm);
  //iq = (P * Viterm - Q * Vrterm) / (Vterm * Vterm);
  ip = genP/Vterm;
  iq = genQ/Vterm;

  //initialize x2iq block
  iqcmd = delayblocklimit_x2iq.init( iq, tg, iqrmax, iqrmin);

  //initialize x3vmeas block
  double lvpl;
  lvpl = delayblock_x3vmeas.init(Vterm, tfltr);
  
  //initialize x1ip block
  ipcmd = delayblocklimit_x1ip.init(ip, tg, rrpwr, -9999.0);
  
  //initialize the REEC and REPC model here---------------
  if (p_hasExciter){
	p_exciter = getExciter();
	p_exciter->setVterminal(Vterm); 
	p_exciter->setIpcmdIqcmd(ipcmd, iqcmd); 


	p_exciter->init(mag, ang, ts);
	
	pref = p_exciter->getPref( );
	qext = p_exciter->getQext( );

  }
  
  if (p_hasPlantController){
	p_plant = getPlantController();
	p_plant->setGenPQV(genP, genQ, Vterm);
	p_plant->setPrefQext(pref, qext);
	p_plant->setExtBusNum(p_bus_id); // Yuan added on 2020-6-23

	p_plant->init(mag, ang, ts);	

  }
  
    x1ip_0 = delayblocklimit_x1ip.x0;
	x2iq_0 = delayblocklimit_x2iq.x0;
	x3vmeas_0 = delayblock_x3vmeas.x0;
	
    x1ip_1 = x1ip_0;
	x2iq_1 = x2iq_0;
	x3vmeas_1 = x3vmeas_0;
  
  if (bmodel_debug){
    printf("Regca1Generator Ip = %f, Iq = %f, Ipcmd = %f, Iqcmd = %f, pref = %f, qext = %f\n", ip, iq, ipcmd, iqcmd, pref, qext);
	printf("Regca1Generator x1ip_0 = %f, x2iq_0 = %f, x3vmeas_0 = %f, \n", x1ip_0, x2iq_0, x3vmeas_0);
  }
  
  p_Norton_Ya = NortonImpedence();
  if (bmodel_debug){
	printf("------renke debug in Regca1Generator::init, p_Norton_Ya = %f, %f \n", real(p_Norton_Ya), imag(p_Norton_Ya));
  }

}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType gridpack::dynamic_simulation::Regca1Generator::INorton()
{
  return p_INorton;
}


/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType gridpack::dynamic_simulation::Regca1Generator::NortonImpedence()
{
  /*double ra = Ra * p_sbase / MVABase;
  double xd = Xdpp * p_sbase / MVABase;
  B = -xd / (ra + xd);
  G = ra / (ra + xd);
  gridpack::ComplexType Y_a(B, G);
  return Y_a;*/
  //double ra = 0.0; //Ra * p_sbase / MVABase;
  //double xd = XL * p_sbase / MVABase;
  //printf("Ra = %f, Xdpp = %f\n", Ra, Xdpp);
  //printf("ra = %f, xd = %f, p_sbase = %f, MVABase = %f\n", ra, xd, p_sbase, MVABase);
  double B = 0.0; //-xd / (ra * ra + xd * xd);
  double G = 0.0; //ra / (ra * ra + xd * xd);
  //printf("B = %f, G = %f\n", B, G);
  gridpack::ComplexType Y_a(G, B);
  return Y_a;
}


/**
 * Predict part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Regca1Generator::predictor_currentInjection(bool flag)
{
  if (!flag) {
    x1ip_0 = x1ip_1;
    x2iq_0 = x2iq_1;
    x3vmeas_0 = x3vmeas_1;

  }  
  
  Vterm = presentMag;
  //printf("Regca1Generator predictor_currentInjection: %d %f\n", p_bus_id, Vterm);
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  
  //note, need to convert current to system mva base!!!!!!!!!!
  double Ireal = ip * MVABase / p_sbase;
  double Iimag = -iq * MVABase / p_sbase;  /// important, here Iimag has a minus sign!!!!!!
  
  // need to rotate the current by Theta to go to the system common reference frame!!!!!!!!!!!!!!
  IrNorton = Ireal*cos(Theta) - Iimag*sin(Theta);
  IiNorton = Ireal*sin(Theta) + Iimag*cos(Theta);
  
  gridpack::ComplexType vt_complex_tmp = gridpack::ComplexType(Vrterm, Viterm); 
  
  if (getGenStatus()){
	  if (p_tripped){
		  p_INorton = p_Norton_Ya*vt_complex_tmp;
	  }else{
		  p_INorton = gridpack::ComplexType(IrNorton, IiNorton);	
	  }		
  }else {
	  p_INorton = gridpack::ComplexType(0.0, 0.0);
  }
  //printf("Regca1Generator::predictor_currentInjection: presentMag = %f, presentAng = %f \n", presentMag, presentAng);
  if (bmodel_debug){
	printf("Regca1Generator::predictor_currentInjection: presentMag = %f, presentAng = %f, ip = %f, iq = %f, Ireal = %f, Iimag = %f \n", presentMag, presentAng, ip, iq, Ireal, Iimag);
	printf("Regca1Generator::predictor_currentInjuction: p_INorton = %f, %f \n", real(p_INorton), imag(p_INorton));
  }
} 

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Regca1Generator::predictor(
    double t_inc, bool flag)
{
  ///printf("\n***** GEN %d Predicator:\n", p_bus_id);
  if (getGenStatus()){

	if (!flag) {
		x1ip_0 = x1ip_1;
		x2iq_0 = x2iq_1;
		x3vmeas_0 = x3vmeas_1;
	}    	
	
	//compute generator P and Q first 
	Vterm = presentMag;
	Theta = presentAng;
	double Vrterm = Vterm * cos(Theta);
	double Viterm = Vterm * sin(Theta);
	
	//genP = Vrterm*ip + Viterm*iq;
	//genQ = Viterm*ip - Vrterm*iq;
	
	genP = Vterm*ip;
	genQ = Vterm*iq;
	
	if (bmodel_debug){
		printf("------renke debug in Regca1Generator::predictor, Vterm = %f, Theta = %f, Vrterm= %f, Viterm= %f, Ip= %f, Iq= %f\n", Vterm, Theta, Vrterm, Viterm, ip, iq);
		printf("------renke debug in Regca1Generator::predictor, genP = %f, genQ = %f, busfreq = %f \n", genP, genQ, busfreq);
		printf("------renke debug in Regca1Generator::predictor, presentMag, presentAng = %12.6f, %12.6f \n", presentMag, presentAng);
    }
	
	
	//----------Set Vt, Ip and Iq, genP, gen Q, busfreq for the REECA1 model first
	
	//get pref and qext from the plant controller REPCA1 MODEL
	if (p_hasPlantController){
		p_plant = getPlantController();
		p_plant->setGenPQV(genP, genQ, Vterm);
		p_plant->setBusFreq(busfreq);
		
		p_plant->predictor(t_inc, flag);
		
		pref = p_plant->getPref();
		qext = p_plant->getQext();
	}
	
   //get ipcmd and iqcmd from the exctier REECA1 MODEL
    if (p_hasExciter){
		p_exciter = getExciter();
		p_exciter->setPrefQext(pref, qext);
		p_exciter->setVterminal(Vterm);
		p_exciter->predictor(t_inc, flag);
		
		ipcmd = p_exciter->getIpcmd();
		iqcmd = p_exciter->getIqcmd();		
	}

	///----------------------------------------------------
	 iq= delayblocklimit_x2iq.predictor (iqcmd, t_inc, flag);
	 
	 double vlimitout;
	 vlimitout = delayblock_x3vmeas.predictor (Vterm, t_inc, flag);
	 
	 //!!!!!!!!!!!!!!!! linear block for vlimitout
	 //delayblocklimit_x1ip.Max = vlimitout;
	 ip= delayblocklimit_x1ip.predictor (ipcmd, t_inc, flag);
	
	//------------!!!!!!!!need to update ip and iq here to connect to current injection function
	

	if (bmodel_debug){	
		printf("---Regca1Generator predictor: pbusid: %d, dx: %12.6f, %12.6f, %12.6f\n", p_bus_id, delayblocklimit_x1ip.dx0, delayblocklimit_x2iq.dx0, delayblock_x3vmeas.dx0);
	}
	
	if (bmodel_debug){
		printf("------renke debug in Regca1Generator::predictor, presentMag, presentAng = %12.6f, %12.6f \n", presentMag, presentAng);
		printf("------Regca1Generator test 1 predictor output x results: %8d, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f  \n",
          p_bus_id,  delayblocklimit_x1ip.x0, delayblocklimit_x2iq.x0, delayblock_x3vmeas.x0,  Vterm, Theta, genP, genQ);
		printf("------Regca1Generator test 2 predictor other: %8d, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f,\n",
          p_bus_id,  ipcmd,  ip, iqcmd, iq, Vterm);
	}
	
	if (p_tripped){
		delayblocklimit_x2iq.x0 = 0.0;
		delayblocklimit_x1ip.x0 = 0.0;
		delayblock_x3vmeas.x0 = 0.0;
		
		delayblocklimit_x2iq.x1 = 0.0;
		delayblocklimit_x1ip.x1 = 0.0;
		delayblock_x3vmeas.x1 = 0.0;
	
	}
	
  }else {
    delayblocklimit_x2iq.x0 = 0.0;
	delayblocklimit_x1ip.x0 = 0.0;
	delayblock_x3vmeas.x0 = 0.0;
		
	delayblocklimit_x2iq.x1 = 0.0;
	delayblocklimit_x1ip.x1 = 0.0;
	delayblock_x3vmeas.x1 = 0.0;
  }//pair with getGenStatus()

}

/**
 * Corrector part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Regca1Generator::corrector_currentInjection(bool flag)
{
  Vterm = presentMag;
  //printf("Regca1Generator predictor_currentInjection: %d %f\n", p_bus_id, Vterm);
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  
  //note, need to convert current to system mva base!!!!!!!!!!
  double Ireal = ip * MVABase / p_sbase;
  double Iimag = -iq * MVABase / p_sbase;  /// important, here Iimag has a minus sign!!!!!!
  
  // need to rotate the current by Theta to go to the system common reference frame!!!!!!!!!!!!!!
  IrNorton = Ireal*cos(Theta) - Iimag*sin(Theta);
  IiNorton = Ireal*sin(Theta) + Iimag*cos(Theta);
  
  gridpack::ComplexType vt_complex_tmp = gridpack::ComplexType(Vrterm, Viterm);   
  
  if (getGenStatus()){
	  if (p_tripped){
		  p_INorton = p_Norton_Ya*vt_complex_tmp;
	  }else{
		  p_INorton = gridpack::ComplexType(IrNorton, IiNorton);	
	  }		
  }else {
	  p_INorton = gridpack::ComplexType(0.0, 0.0);
  }
  //printf("Regca1Generator::predictor_currentInjection: presentMag = %f, presentAng = %f \n", presentMag, presentAng);
  if (bmodel_debug){
	printf("Regca1Generator::correct_currentInjuction: p_INorton = %f, %f \n", real(p_INorton), imag(p_INorton));
  }
  
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Regca1Generator::corrector(
    double t_inc, bool flag)
{
  ///printf("\n***** GEN %d Corrector:\n", p_bus_id);
  if (getGenStatus()){
	
	//compute generator P and Q first 
	Vterm = presentMag;
	Theta = presentAng;
	double Vrterm = Vterm * cos(Theta);
	double Viterm = Vterm * sin(Theta);
	
	//genP = Vrterm*ip + Viterm*iq;
	//genQ = Viterm*ip - Vrterm*iq;
	
	genP = Vterm*ip;
	genQ = Vterm*iq;
	
	if (bmodel_debug){
		printf("------renke debug in Regca1Generator::corrector, Vterm = %f, Theta = %f, Vrterm= %f, Viterm= %f, Ip= %f, Iq= %f\n", Vterm, Theta, Vrterm, Viterm, ip, iq);
		printf("------renke debug in Regca1Generator::corrector, genP = %f, genQ = %f \n", genP, genQ);
		printf("------renke debug in Regca1Generator::corrector, presentMag, presentAng = %12.6f, %12.6f \n", presentMag, presentAng);
    }
	
	
	//----------Set Vt, Ip and Iq, genP, gen Q, busfreq for the REECA1 model first
	
	//get pref and qext from the plant controller REPCA1 MODEL
	if (p_hasPlantController){
		p_plant = getPlantController();
		p_plant->setGenPQV(genP, genQ, Vterm);
		p_plant->setBusFreq(busfreq);
		
		p_plant->corrector(t_inc, flag);
		
		pref = p_plant->getPref();
		qext = p_plant->getQext();
	}
	
   //get ipcmd and iqcmd from the exctier REECA1 MODEL
    if (p_hasExciter){
		p_exciter = getExciter();
		p_exciter->setPrefQext(pref, qext);
		p_exciter->setVterminal(Vterm);
		
		p_exciter->corrector(t_inc, flag);
		
		ipcmd = p_exciter->getIpcmd();
		iqcmd = p_exciter->getIqcmd();		
	}
	
	
	///----------------------------------------------------
	 iq= delayblocklimit_x2iq.corrector (iqcmd, t_inc, flag);
	 
	 double vlimitout;
	 vlimitout = delayblock_x3vmeas.corrector (Vterm, t_inc, flag);
	 
	 //!!!!!!!!!!!!!!!! linear block for vlimitout
	 //delayblocklimit_x1ip.Max = vlimitout;
	 ip= delayblocklimit_x1ip.corrector (ipcmd, t_inc, flag);
	
	//------------!!!!!!!!need to update ip and iq here to connect to current injection function
	

	if (bmodel_debug){	
		printf("---Regca1Generator corrector: pbusid: %d, dx: %12.6f, %12.6f, %12.6f\n", p_bus_id, delayblocklimit_x1ip.dx1, delayblocklimit_x2iq.dx1, delayblock_x3vmeas.dx1);
	}
	
	if (bmodel_debug){
		printf("------renke debug in Regca1Generator::corrector, presentMag, presentAng = %12.6f, %12.6f \n", presentMag, presentAng);
		printf("------------corrector output results: %8d, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f  \n",
          p_bus_id,  delayblocklimit_x1ip.x1, delayblocklimit_x2iq.x1, delayblock_x3vmeas.x1,  Vterm, Theta, genP, genQ);
	}
	
	if (p_tripped){
		delayblocklimit_x2iq.x0 = 0.0;
		delayblocklimit_x1ip.x0 = 0.0;
		delayblock_x3vmeas.x0 = 0.0;
		
		delayblocklimit_x2iq.x1 = 0.0;
		delayblocklimit_x1ip.x1 = 0.0;
		delayblock_x3vmeas.x1 = 0.0;
	
	}
	
  }else {
    delayblocklimit_x2iq.x0 = 0.0;
	delayblocklimit_x1ip.x0 = 0.0;
	delayblock_x3vmeas.x0 = 0.0;
		
	delayblocklimit_x2iq.x1 = 0.0;
	delayblocklimit_x1ip.x1 = 0.0;
	delayblock_x3vmeas.x1 = 0.0;
  }//pair with getGenStatus()
  
}

bool gridpack::dynamic_simulation::Regca1Generator::tripGenerator()
{
	p_tripped = true;
	
	return true;
}

/**
* return true if modify the generator parameters successfully
* input controlTyp: 0: GFI mp adjust; 1: GFI mq adjust; 2: GFI Pset adjust; 3: GFI Qset adjust; others: invalid 
* input newParValScaletoOrg:  GFI new parameter scale factor to the very initial parameter value at the begining of dynamic simulation
* 
*/
bool gridpack::dynamic_simulation::Regca1Generator::applyGeneratorParAdjustment(int controlType, double newParValScaletoOrg){
	
	return true;
	
}

/**
 * Set voltage on each generator
 */
void gridpack::dynamic_simulation::Regca1Generator::setVoltage(
    gridpack::ComplexType voltage)
{
  presentMag = abs(voltage);
  presentAng = atan2(imag(voltage), real(voltage));  
}

/**
 * Set frequency on each generator, frequency is perunit
 */
void gridpack::dynamic_simulation::Regca1Generator::setFreq(double dFreq)
{
   busfreq = dFreq;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
bool gridpack::dynamic_simulation::Regca1Generator::serialWrite(
    char* string, const int bufsize, const char *signal)
{
  if (!strcmp(signal,"standard")) {
    //sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f\n",
    //    p_bus_id,p_ckt.c_str(),real(p_mac_ang_s1),real(p_mac_spd_s1),real(p_mech),
    //    real(p_pelect));
    sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    \n",
          p_bus_id, p_ckt.c_str(), delayblocklimit_x1ip.x1, delayblocklimit_x2iq.x1, delayblock_x3vmeas.x1);
    return true;
  } else if (!strcmp(signal,"init_debug")) {
    sprintf(string," %8d  %2s Something\n",p_bus_id,p_ckt.c_str());
    return true;
  } else if (!strcmp(signal,"debug_initial")) {
    sprintf(string,"Bus: %d PG: %f QG: %f\n",p_bus_id,p_pg,p_qg);
    return true;
  } else if (!strcmp(signal,"watch")) {
    if (getWatch()) {
      char buf[256];
//    sprintf(buf,", %f, %f",real(p_mac_ang_s1),real(p_mac_spd_s1));
      sprintf(string,",%8d, %2s, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f , %12.6f ",
          p_bus_id, p_ckt.c_str(), delayblocklimit_x1ip.x1, delayblocklimit_x2iq.x1, delayblock_x3vmeas.x1, presentMag, presentAng, genP, genQ, busfreq);
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
}

/**
 * return a vector containing any generator values that are being
 * watched
 * @param vals vector of watched values
 */
void gridpack::dynamic_simulation::Regca1Generator::getWatchValues(
    std::vector<double> &vals)
{
  vals.clear();
  vals.push_back(delayblocklimit_x1ip.x1);
  vals.push_back(delayblocklimit_x2iq.x1);
  vals.push_back(genP*MVABase/p_sbase); //output at system mva base
  vals.push_back(genQ*MVABase/p_sbase); //output at system mva base
}