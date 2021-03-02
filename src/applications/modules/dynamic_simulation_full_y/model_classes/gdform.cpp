/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   gdform.hpp
 * @author Renke Huang
 * @Last modified:   Jan. 16, 2021
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
#include "gdform.hpp"

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::GridFormingGenerator::GridFormingGenerator(void)
{
    x1E_0 = 0.0;
	x2d_0 = 0.0;
	x3_0 = 0.0;
	x4_0 = 0.0;
	xvterm_0= 0.0;
    xp_0= 0.0; 
	xq_0= 0.0;
	
	x1E_1 = 0.0;
	x2d_1 = 0.0;
	x3_1 = 0.0;
	x4_1 = 0.0;
	xvterm_1= 0.0;
    xp_1= 0.0; 
	xq_1= 0.0;
	
	dx1E_0 = 0.0;
	dx2d_0 = 0.0;
	dx3_0 = 0.0;
	dx4_0 = 0.0;
	dxvterm_0= 0.0;
    dxp_0= 0.0; 
	dxq_0= 0.0;
	
	dx1E_1 = 0.0;
	dx2d_1 = 0.0;
	dx3_1 = 0.0;
	dx4_1 = 0.0;
	dxvterm_1= 0.0;
    dxp_1= 0.0; 
	dxq_1= 0.00;
	
	presentMag = 1.0;
	presentAng = 0.0;
	fset = 60.0;
	Ra = 0.0;
	omega = 1.0;
	delta_omega_lim = 999.0;
	Poutctrl = 0.0;
	Ts = 0.01;
	
	p_tripped = false;
	bmodel_debug = false;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::GridFormingGenerator::~GridFormingGenerator(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::GridFormingGenerator::load(
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
     printf("GridFormingGenerator::load p_pg = %f, p_qg = %f\n", p_pg, p_qg);
  }
  p_pg *= p_sbase;
  p_qg *= p_sbase;

  if (!data->getValue(GENERATOR_MBASE, &MVABase, idx)) MVABase = 1000.0; // MVABase
  if (!data->getValue(GENERATOR_XL  , &XL, idx)) XL = 0.075; 
  if (!data->getValue(GENERATOR_VSET, &Ts, idx)) Ts = 0.01666; // D
  if (!data->getValue(GENERATOR_MQ, &mq, idx)) mq=0.05; // 
  if (!data->getValue(GENERATOR_KPV, &kpv, idx)) kpv=0.0; // 
  if (!data->getValue(GENERATOR_KIV, &kiv, idx)) kiv=5.86; // 
  if (!data->getValue(GENERATOR_EMAX, &Emax, idx)) Emax = 2.0; // 
  if (!data->getValue(GENERATOR_EMIN, &Emin, idx)) Emin = -2.0; // 
  if (!data->getValue(GENERATOR_MP, &mp, idx)) mp=3.77; // 
  if (!data->getValue(GENERATOR_KPPMAX, &kppmax, idx)) kppmax=3.0; // 
  if (!data->getValue(GENERATOR_KIPMAX, &kipmax, idx)) kipmax=30.0; // 
  if (!data->getValue(GENERATOR_PSET, &Pset, idx)) Pset=0.44; // 
  if (!data->getValue(GENERATOR_PMAX, &Pmax, idx)) Pmax=1.0; // 
  if (!data->getValue(GENERATOR_PMIN, &Pmin, idx)) Pmin=0.0; // 
  
  mp_org = mp;
  mq_org = mq; 
  
  double tmp = sqrt(p_pg*p_pg +p_qg*p_qg);
  if ( tmp > MVABase) {
       //MVABase = tmp*1.3;
      //printf("-----------generator at bus %d  has P: %f, Q: %f, S: %f, MVABASE: %f  \n", p_bus_id, p_pg, p_qg, tmp, MVABase);
	  MVABase = tmp*1.2;
  }

  if (bmodel_debug){
	printf("\n--------gdform parameters: MVABase = %12.6f, XL = %12.6f, Ts = %12.6f, mq = %12.6f, kpv = %12.6f, kiv = %12.6f, Emax = %12.6f, Emin = %12.6f, mp = %12.6f, kppmax = %12.6f, kipmax = %12.6f, Pset = %12.6f, Pmax = %12.6f, Pmin = %12.6f  \n", MVABase, XL, Ts, mq, kpv, kiv, Emax, Emin, mp, kppmax, kipmax, Pset, Pmax, Pmin);
  }
  
}

/**
 * Initialize generator model before calculation
 * @param mag voltage magnitude
 * @param ang voltage angle
 */
void gridpack::dynamic_simulation::GridFormingGenerator::init(double mag,
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
  if (bmodel_debug){
	printf("p_pg = %f, p_qg = %f, MVABase = %f\n", p_pg, p_qg, MVABase);
	printf("Vterm = %f, Theta = %f, P = %f, Q = %f\n", Vterm, Theta, P, Q);
  }
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  Ir = (P * Vrterm + Q * Viterm) / (Vterm * Vterm);
  Ii = (P * Viterm - Q * Vrterm) / (Vterm * Vterm);
  
  //compute E and delta
  double Etermr = Vrterm - Ii*XL;
  double Etermi = Viterm + Ir*XL;
  
  x1E_0 = sqrt(Etermr*Etermr+Etermi*Etermi);
  x1E_1 = x1E_0;
  E_term = x1E_0;
  x2d_0 = atan2(Etermi, Etermr);  // from -pi to pi, rads
  x2d_1 = x2d_0;
  xvterm_0 = Vterm;
  xp_0 = genP;
  xq_0 = genQ;
  xvterm_1 = Vterm;
  xp_1 = genP;
  xq_1 = genQ;
  
  if (bmodel_debug){
    printf("Ir = %f, Ii = %f\n", Ir, Ii);
	printf("Etermr = %f, Etermi = %f\n", Etermr, Etermi);
	printf("x1E_0 = %f, x2d_0 = %f\n", x1E_0, x2d_0);
  }
  
  Vset = Vterm + Q*mq;
  Pset = P;
  
  mp_org = mp;
  mq_org = mq;
  Vset_org = Vset;
  Pset_org = Pset;
  
  if (bmodel_debug){
	 printf("Vset = %f, Pset = %f\n", Vset, Pset);
  }
	
  p_Norton_Ya = NortonImpedence();
  if (bmodel_debug){
	printf("------renke debug in GridFormingGenerator::init, p_Norton_Ya = %f, %f \n", real(p_Norton_Ya), imag(p_Norton_Ya));
  }

}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType gridpack::dynamic_simulation::GridFormingGenerator::INorton()
{
  return p_INorton;
}


/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType gridpack::dynamic_simulation::GridFormingGenerator::NortonImpedence()
{
  /*double ra = Ra * p_sbase / MVABase;
  double xd = Xdpp * p_sbase / MVABase;
  B = -xd / (ra + xd);
  G = ra / (ra + xd);
  gridpack::ComplexType Y_a(B, G);
  return Y_a;*/
  double ra = 0.0; //Ra * p_sbase / MVABase;
  double xd = XL * p_sbase / MVABase;
  //printf("Ra = %f, Xdpp = %f\n", Ra, Xdpp);
  //printf("ra = %f, xd = %f, p_sbase = %f, MVABase = %f\n", ra, xd, p_sbase, MVABase);
  B = -xd / (ra * ra + xd * xd);
  G = ra / (ra * ra + xd * xd);
  //printf("B = %f, G = %f\n", B, G);
  gridpack::ComplexType Y_a(G, B);
  return Y_a;
}


/**
 * Predict part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GridFormingGenerator::predictor_currentInjection(bool flag)
{
  if (!flag) {
    x1E_0 = x1E_1;
    x2d_0 = x2d_1;
    x3_0 = x3_1;
    x4_0 = x4_1;
	xvterm_0 = xvterm_1;
	xp_0=xp_1;
	xq_0=xq_1;

  }  
  // Calculate INorton_full
  // Admittance
  B = -XL / (Ra * Ra + XL * XL);
  G = Ra / (Ra * Ra + XL * XL);
  
  double Etermr = E_term * cos(x2d_0);
  double Etermi = E_term * sin(x2d_0);
  double Irnorton_gen = Etermr * G - Etermi * B;
  double Iinorton_gen = Etermr * B + Etermi * G;
  
  double IrNorton = Irnorton_gen * MVABase / p_sbase; 
  double IiNorton = Iinorton_gen * MVABase / p_sbase;
    
  Vterm = presentMag;
  //printf("Gensal predictor_currentInjection: %d %f\n", p_bus_id, Vterm);
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  
  Ir = (Etermr - Vrterm) * G - (Etermi - Viterm) * B;
  Ii = (Etermr - Vrterm) * B + (Etermi - Viterm) * G;
  //printf("Xdpp = %f, Ra = %f, B = %f, G = %f\n", Xdpp, Ra, B, G);
  
  gridpack::ComplexType vt_complex_tmp = gridpack::ComplexType(Vrterm, Viterm); 
  //printf("x5Psiqpp_0 = %f, x3Eqp_0 = %f, Xl = %f, Xdp = %f, Psidpp = %f\n", x5Psiqpp_0, x3Eqp_0, Xl, Xdp, Psidpp);
  
  if (getGenStatus()){
	  if (p_tripped){
		  p_INorton = p_Norton_Ya*vt_complex_tmp;
	  }else{
		  p_INorton = gridpack::ComplexType(IrNorton, IiNorton);	
	  }		
  }else {
	  p_INorton = gridpack::ComplexType(0.0, 0.0);
  }
  //printf("gensal::predictor_currentInjection: presentMag = %f, presentAng = %f \n", presentMag, presentAng);
  if (bmodel_debug){
	printf("gridforming::predictor_currentInjuction: p_INorton = %f, %f \n", real(p_INorton), imag(p_INorton));
  }
} 

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GridFormingGenerator::predictor(
    double t_inc, bool flag)
{
  ///printf("\n***** GEN %d Predicator:\n", p_bus_id);
  if (getGenStatus()){

	if (!flag) {
	x1E_0 = x1E_1;
	x2d_0 = x2d_1;
	x3_0 = x3_1;
	x4_0 = x4_1;
	xvterm_0 = xvterm_1;
	xp_0=xp_1;
	xq_0=xq_1;
	}  
	
	//compute generator P and Q first 
	Vterm = presentMag;
	Theta = presentAng;
	double Vrterm = Vterm * cos(Theta);
	double Viterm = Vterm * sin(Theta);
	double Etermr = E_term * cos(x2d_0);
    double Etermi = E_term * sin(x2d_0);
	
	B = -XL / (Ra * Ra + XL * XL);
    G = Ra / (Ra * Ra + XL * XL);
	
	Ir = (Etermr - Vrterm) * G - (Etermi - Viterm) * B;
	Ii = (Etermr - Vrterm) * B + (Etermi - Viterm) * G;
	
	genP = Vrterm*Ir + Viterm*Ii;
	genQ = Viterm*Ir - Vrterm*Ii;
	
	if (bmodel_debug){
		printf("------renke debug in GridFormingGenerator::predictor, Vterm = %f, Theta = %f, E_term= %f, x2d_0= %f, Vrterm= %f, Viterm= %f, Etermr= %f, Etermi= %f, B= %f, G= %f, Ir= %f, Ii= %f\n", Vterm, Theta, E_term, x2d_0, Vrterm, Viterm, Etermr, Etermi, B, G, Ir, Ii);
		printf("------renke debug in GridFormingGenerator::predictor, genP = %f, genQ = %f \n", genP, genQ);
		printf("------renke debug in GridFormingGenerator::predictor, presentMag, presentAng = %12.6f, %12.6f \n", presentMag, presentAng);
    }
	
	//--------------add delay function here for genP, genQ, and Vterm----------------
	if (Ts< 4*t_inc){
		dxvterm_0 = 0.0;
		dxp_0 = 0.0;
		dxq_0 = 0.0;
		Vterm_delay = Vterm;
		genP_delay = genP;
		genQ_delay = genQ;
		xvterm_0 = Vterm;
		xp_0 = genP;
		xq_0 = genQ;
	}else{
		
		dxvterm_0 = (Vterm-xvterm_0)/Ts;
		dxp_0 = (genP-xp_0)/Ts;
		dxq_0 = (genQ-xq_0)/Ts;
		
		Vterm_delay = xvterm_0;
		genP_delay = xp_0;
		genQ_delay = xq_0;
	}
	
	double tmpin = -mq*genQ_delay - Vterm_delay + Vset;
	double tmpmax = Emax - kpv*tmpin;
	double tmpmin = Emin - kpv*tmpin;
	
	if (x1E_0 > tmpmax) x1E_0 = tmpmax;
	if (x1E_0 < tmpmin) x1E_0 = tmpmin;
	
	dx1E_0 = kiv*tmpin;
	
	if( dx1E_0>0.0 && x1E_0>=tmpmax ) dx1E_0 = 0.0;
	if( dx1E_0<0.0 && x1E_0<=tmpmin ) dx1E_0 = 0.0;
	
	Poutctrl = 0.0;
	double out3 = 0.0;
	double out4 = 0.0;
	//compute s3
	if ( genP_delay>Pmax ){
		tmpin = Pmax - genP_delay;
		tmpmax = 0.0 - kppmax*tmpin;
		tmpmin = -delta_omega_lim - kppmax*tmpin;
		
		if (x3_0 > tmpmax) x3_0 = tmpmax;
		if (x3_0 < tmpmin) x3_0 = tmpmin;
		
		dx3_0 = kipmax*tmpin;
		
		if( dx3_0>0.0 && x3_0>=tmpmax ) dx3_0 = 0.0;
		if( dx3_0<0.0 && x3_0<=tmpmin ) dx3_0 = 0.0;
		out3 = kipmax * tmpin + x3_0;
		if (out3 > 0.0) out3 = 0.0;
		if (out3 < -delta_omega_lim) out3 = -delta_omega_lim;
	}else{
		 x3_0 = 0.0;
		dx3_0 = 0.0;
	}
	
    //compute s4
	if ( genP_delay<Pmin ){
		tmpin = Pmin - genP_delay;
		tmpmax = delta_omega_lim - kppmax*tmpin;
		tmpmin = 0.0 - kppmax*tmpin;
		
		if (x4_0 > tmpmax) x4_0 = tmpmax;
		if (x4_0 < tmpmin) x4_0 = tmpmin;
		
		dx4_0 = kipmax*tmpin;
		
		if( dx4_0>0.0 && x4_0>=tmpmax ) dx4_0 = 0.0;
		if( dx4_0<0.0 && x4_0<=tmpmin ) dx4_0 = 0.0;
		
		out4 = kipmax * tmpin + x4_0;
		if (out4 > delta_omega_lim) out4 = delta_omega_lim;
		if (out4 < 0.0) out4 = 0.0;
	}else{
		 x4_0 = 0.0;
		dx4_0 = 0.0;
	}
	
	Poutctrl = out3 + out4;
	
	if (bmodel_debug){
		printf("---gdform predictor x3x4 debug, pbusid, %d,  x3_0, %12.6f, x4_0, %12.6f, dx3_0, %12.6f, dx4_0, %12.6f, out3, %12.6f, out4, %12.6f, Poutctrl, %12.6f, \n", p_bus_id, x3_0, x4_0, dx3_0, dx4_0, out3, out4, Poutctrl);
	}
	// compute s2
	double pi = 4.0*atan(1.0);
	tmpin = (Pset-genP_delay)*mp + Poutctrl;
	tmpin = tmpin + 2*pi*fset - 2*pi*60.0;
	dx2d_0 = tmpin/1.0;
	omega = tmpin/(2*pi*60.0)+1.0;
	
	x1E_1 = x1E_0 + dx1E_0 * t_inc;
	x2d_1 = x2d_0 + dx2d_0 * t_inc;
	x3_1 = x3_0 + dx3_0 * t_inc;
	x4_1 = x4_0 + dx4_0 * t_inc;
	xvterm_1 = xvterm_0 + dxvterm_0 * t_inc;
	xp_1 = xp_0 + dxp_0 * t_inc;
	xq_1 = xq_0 + dxq_0 * t_inc;
	
	E_term = tmpin*kpv + x1E_1;
	
	if (E_term > Emax) E_term = Emax;
	if (E_term < Emin) E_term = Emin;

	if (bmodel_debug){	
		printf("---gdform predictor: pbusid: %d, dx: %12.6f, %12.6f, %12.6f, %12.6f \n", p_bus_id, dx1E_0, dx2d_0, dx3_0, dx4_0);
	}
	///printf("gensal x: %f\t%f\t%f\t%f\t%f\n", x1d_1, x2w_1, x3Eqp_1, x4Psidp_1, x5Psiqpp_1);
	//printf ("gensal predictor presentMag=%f\n", presentMag);

	
	if (p_tripped){
	x1E_0 = 0.0;
	x2d_0 = 0.0;
	x3_0 = 0.0;
	x4_0 = 0.0;
	xvterm_0 = 0.0;
	xp_0 = 0.0;
	xq_0 = 0.0;
	
	x1E_1 = 0.0;
	x2d_1 = 0.0;
	x3_1 = 0.0;
	x4_1 = 0.0;
	xvterm_1 = 0.0;
	xp_1 = 0.0;
	xq_1 = 0.0;
	genP = 0.0;
	genQ = 0.0;
	omega= 1.0;
	}
	
  }else {
    x1E_0 = 0.0;
	x2d_0 = 0.0;
	x3_0 = 0.0;
	x4_0 = 0.0;
	xvterm_0 = 0.0;
	xp_0 = 0.0;
	xq_0 = 0.0;
	
	x1E_1 = 0.0;
	x2d_1 = 0.0;
	x3_1 = 0.0;
	x4_1 = 0.0;
	xvterm_1 = 0.0;
	xp_1 = 0.0;
	xq_1 = 0.0;
	genP = 0.0;
	genQ = 0.0;
	omega= 1.0;
  }//pair with getGenStatus()

}

/**
 * Corrector part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GridFormingGenerator::corrector_currentInjection(bool flag)
{
  // Calculate INorton_full
  // Admittance
  // Calculate INorton_full
  // Admittance
  B = -XL / (Ra * Ra + XL * XL);
  G = Ra / (Ra * Ra + XL * XL);
  
  double Etermr = E_term * cos(x2d_1);
  double Etermi = E_term * sin(x2d_1);
  double Irnorton_gen = Etermr * G - Etermi * B;
  double Iinorton_gen = Etermr * B + Etermi * G;
  
  double IrNorton = Irnorton_gen * MVABase / p_sbase; 
  double IiNorton = Iinorton_gen * MVABase / p_sbase;
    
  Vterm = presentMag;
  //printf("Gensal predictor_currentInjection: %d %f\n", p_bus_id, Vterm);
  Theta = presentAng;
  double Vrterm = Vterm * cos(Theta);
  double Viterm = Vterm * sin(Theta);
  
  Ir = (Etermr - Vrterm) * G - (Etermi - Viterm) * B;
  Ii = (Etermr - Vrterm) * B + (Etermi - Viterm) * G;
  //printf("Xdpp = %f, Ra = %f, B = %f, G = %f\n", Xdpp, Ra, B, G);
  
  gridpack::ComplexType vt_complex_tmp = gridpack::ComplexType(Vrterm, Viterm); 
  //printf("x5Psiqpp_0 = %f, x3Eqp_0 = %f, Xl = %f, Xdp = %f, Psidpp = %f\n", x5Psiqpp_0, x3Eqp_0, Xl, Xdp, Psidpp);
  
  if (getGenStatus()){
	  if (p_tripped){
		  p_INorton = p_Norton_Ya*vt_complex_tmp;
	  }else{
		  p_INorton = gridpack::ComplexType(IrNorton, IiNorton);	
	  }		
  }else {
	  p_INorton = gridpack::ComplexType(0.0, 0.0);
  }
  
  if (bmodel_debug){
	printf("gridforming::corrector_currentInjuction: p_INorton = %d, %f, %f \n", p_bus_id, real(p_INorton), imag(p_INorton));
  }
  
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::GridFormingGenerator::corrector(
    double t_inc, bool flag)
{
  ///printf("\n***** GEN %d Corrector:\n", p_bus_id);
  if (getGenStatus()){
	  
	  //compute generator P and Q first 
	Vterm = presentMag;
	Theta = presentAng;
	double Vrterm = Vterm * cos(Theta);
	double Viterm = Vterm * sin(Theta);
	double Etermr = E_term * cos(x2d_1);
    double Etermi = E_term * sin(x2d_1);
	
	B = -XL / (Ra * Ra + XL * XL);
    G = Ra / (Ra * Ra + XL * XL);
	
	Ir = (Etermr - Vrterm) * G - (Etermi - Viterm) * B;
	Ii = (Etermr - Vrterm) * B + (Etermi - Viterm) * G;
	
	genP = Vrterm*Ir + Viterm*Ii;
	genQ = Viterm*Ir - Vrterm*Ii;
	
	if (bmodel_debug){
		printf("------renke debug in GridFormingGenerator::corrector, genP = %f, genQ = %f, Vterm = %f, x1E_1 = %f\n", genP, genQ, Vterm, x1E_1);
    }
	
	//--------------add delay function here for genP, genQ, and Vterm----------------
	if (Ts< 4*t_inc){
		dxvterm_1 = 0.0;
		dxp_1 = 0.0;
		dxq_1 = 0.0;
		Vterm_delay = Vterm;
		genP_delay = genP;
		genQ_delay = genQ;
		xvterm_1 = Vterm;
		xp_1 = genP;
		xq_1 = genQ;
	}else{
		
		dxvterm_1 = (Vterm-xvterm_1)/Ts;
		dxp_1 = (genP-xp_1)/Ts;
		dxq_1 = (genQ-xq_1)/Ts;
		
		Vterm_delay = xvterm_1;
		genP_delay = xp_1;
		genQ_delay = xq_1;
	}
	
	double tmpin = -mq*genQ_delay - Vterm_delay + Vset;
	double tmpmax = Emax - kpv*tmpin;
	double tmpmin = Emin - kpv*tmpin;
	
	if (x1E_1 > tmpmax) x1E_1 = tmpmax;
	if (x1E_1 < tmpmin) x1E_1 = tmpmin;
	
	if (bmodel_debug){
		printf("------renke debug in GridFormingGenerator::corrector, tmpmax = %f, tmpmin = %f, tmpin = %f, x1E_1 = %f \n", tmpmax, tmpmin, tmpin, x1E_1);
    }
	
	dx1E_1 = kiv*tmpin;
	if (bmodel_debug){
		printf("------renke debug in GridFormingGenerator::corrector, dx1E_1 = %12.6f \n", dx1E_1);
    }
	
	if( dx1E_1>0.0 && x1E_1>=tmpmax ) dx1E_1 = 0.0;
	if( dx1E_1<0.0 && x1E_1<=tmpmin ) dx1E_1 = 0.0;
	
	Poutctrl = 0.0;
	double out3 = 0.0;
	double out4 = 0.0;
	
	//compute s3
	if (genP_delay > Pmax){
		tmpin = Pmax - genP_delay;
		tmpmax = 0.0 - kppmax*tmpin;
		tmpmin = -delta_omega_lim - kppmax*tmpin;
		
		if (x3_1 > tmpmax) x3_1 = tmpmax;
		if (x3_1 < tmpmin) x3_1 = tmpmin;
		
		dx3_1 = kipmax*tmpin;
		
		if( dx3_1>0.0 && x3_1>=tmpmax ) dx3_1 = 0.0;
		if( dx3_1<0.0 && x3_1<=tmpmin ) dx3_1 = 0.0;
		double out3 = kipmax * tmpin + x3_1;
		if (out3 > 0.0) out3 = 0.0;
		if (out3 < -delta_omega_lim) out3 = -delta_omega_lim;
	}else{
		 x3_1 = 0.0;
		dx3_1 = 0.0;
	}
	
    //compute s4
	if (genP_delay < Pmin){
		tmpin = Pmin - genP_delay;
		tmpmax = delta_omega_lim - kppmax*tmpin;
		tmpmin = 0.0 - kppmax*tmpin;
		
		if (x4_1 > tmpmax) x4_1 = tmpmax;
		if (x4_1 < tmpmin) x4_1 = tmpmin;
		
		dx4_1 = kipmax*tmpin;
		
		if( dx4_1>0.0 && x4_1>=tmpmax ) dx4_1 = 0.0;
		if( dx4_1<0.0 && x4_1<=tmpmin ) dx4_1 = 0.0;
		
		double out4 = kipmax * tmpin + x4_1;
		if (out4 > delta_omega_lim) out4 = delta_omega_lim;
		if (out4 < 0.0) out4 = 0.0;
	}else{
		 x4_1 = 0.0;
		dx4_1 = 0.0;
	}
	
	Poutctrl = out3 + out4;
	
	if (bmodel_debug){
		printf("---gdform corrector x3x4 debug, pbusid, %d,  x3_1, %12.6f, x4_1, %12.6f, dx3_1, %12.6f, dx4_1, %12.6f, out3, %12.6f, out4, %12.6f, Poutctrl, %12.6f, \n", p_bus_id, x3_1, x4_1, dx3_1, dx4_1, out3, out4, Poutctrl);
	}
	
	// compute s2
	double pi = 4.0*atan(1.0);
	tmpin = (Pset-genP_delay)*mp + Poutctrl;
	tmpin = tmpin + 2*pi*fset - 2*pi*60.0;
	dx2d_1 = tmpin/1.0;
	omega = tmpin/(2*pi*60.0)+1.0;
	
	
	if (bmodel_debug){	
		printf("---gdform corrector: pbusid: %d, dx: %12.6f, %12.6f, %12.6f, %12.6f \n", p_bus_id, dx1E_1, dx2d_1, dx3_1, dx4_1);
	}
  
  
	x1E_1 = x1E_0 + (dx1E_0+dx1E_1)/2.0 * t_inc;
	x2d_1 = x2d_0 + (dx2d_0+dx2d_1)/2.0 * t_inc;
	x3_1 = x3_0 + (dx3_0+dx3_1)/2.0 * t_inc;
	x4_1 = x4_0 + (dx4_0+dx4_1)/2.0 * t_inc;
	xvterm_1 = xvterm_0 + (dxvterm_0+dxvterm_1)/2.0 * t_inc;
	xp_1 = xp_0 + (dxp_0+dxp_1)/2.0 * t_inc;
	xq_1 = xq_0 + (dxq_0+dxq_1)/2.0 * t_inc;
	
	E_term = tmpin*kpv + x1E_1;
	
	if (E_term > Emax) E_term = Emax;
	if (E_term < Emin) E_term = Emin;
	
	if (bmodel_debug){
		printf("------renke debug in GridFormingGenerator::corrector, presentMag, presentAng = %12.6f, %12.6f \n", presentMag, presentAng);
		printf("------------corrector output results: %8d, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f  \n",
          p_bus_id, x1E_1, x2d_1, x3_1, x4_1, E_term, omega, Vterm, Theta, genP, genQ);
	}

 	if (p_tripped){
	x1E_0 = 0.0;
	x2d_0 = 0.0;
	x3_0 = 0.0;
	x4_0 = 0.0;
	xvterm_0 = 0.0;
	xp_0 = 0.0;
	xq_0 = 0.0;
	x1E_1 = 0.0;
	x2d_1 = 0.0;
	x3_1 = 0.0;
	x4_1 = 0.0;
	xvterm_1 = 0.0;
	xp_1 = 0.0;
	xq_1 = 0.0;
	genP = 0.0;
	genQ = 0.0;
	omega= 1.0;
	}
	
  }else {
    x1E_0 = 0.0;
	x2d_0 = 0.0;
	x3_0 = 0.0;
	x4_0 = 0.0;
	xvterm_0 = 0.0;
	xp_0 = 0.0;
	xq_0 = 0.0;
	x1E_1 = 0.0;
	x2d_1 = 0.0;
	x3_1 = 0.0;
	x4_1 = 0.0;
	xvterm_1 = 0.0;
	xp_1 = 0.0;
	xq_1 = 0.0;
	genP = 0.0;
	genQ = 0.0;
	omega= 1.0;
  }//pair with getGenStatus()
  
}

bool gridpack::dynamic_simulation::GridFormingGenerator::tripGenerator()
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
bool gridpack::dynamic_simulation::GridFormingGenerator::applyGeneratorParAdjustment(int controlType, double newParValScaletoOrg){
	
	if ( controlType == 0){
		mp = mp_org * newParValScaletoOrg;
		//printf ("GFI at bus %d, apply newParValScaletoOrg %f to mp, org mp value %f, new mp value %f \n", p_bus_id, newParValScaletoOrg, mp_org, mp);
		return true;
	}else if(controlType == 1){
		mq = mq_org * newParValScaletoOrg;
		//printf ("GFI at bus %d, apply newParValScaletoOrg %f to mq, org mq value %f, new mq value %f \n", p_bus_id, newParValScaletoOrg, mq_org, mq);
		return true;
	}else if(controlType == 2) {
		Pset = Pset_org  * newParValScaletoOrg;
		//printf ("GFI at bus %d, apply newParValScaletoOrg %f to Pset, org Pset value %f, new Pset value %f \n", p_bus_id, newParValScaletoOrg, Pset_org, Pset);
		return true;
	}else if(controlType == 3) {
		Vset = Vset_org  * newParValScaletoOrg;
		//printf ("GFI at bus %d, apply newParValScaletoOrg %f to Vset, org Vset value %f, new Vset value %f \n", p_bus_id, newParValScaletoOrg, Vset_org, Vset);
		return true;
	}else{
		
		//printf ("GFI controlType %d not defined correctly:  \n", controlType);
		return false;
		
	}
	
}

/**
 * Set voltage on each generator
 */
void gridpack::dynamic_simulation::GridFormingGenerator::setVoltage(
    gridpack::ComplexType voltage)
{
  presentMag = abs(voltage);
  presentAng = atan2(imag(voltage), real(voltage));  
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
bool gridpack::dynamic_simulation::GridFormingGenerator::serialWrite(
    char* string, const int bufsize, const char *signal)
{
  if (!strcmp(signal,"standard")) {
    //sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f\n",
    //    p_bus_id,p_ckt.c_str(),real(p_mac_ang_s1),real(p_mac_spd_s1),real(p_mech),
    //    real(p_pelect));
    sprintf(string,"      %8d            %2s    %12.6f    %12.6f    %12.6f    %12.6f	%12.6f\n",
          p_bus_id, p_ckt.c_str(), x1E_1, x2d_1, x3_1, x4_1, omega);
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
      sprintf(string,",%8d, %2s, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f ",
          p_bus_id, p_ckt.c_str(), x3_1, x4_1, E_term, omega, presentMag, presentAng, genP, genQ);
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
void gridpack::dynamic_simulation::GridFormingGenerator::getWatchValues(
    std::vector<double> &vals)
{
  vals.clear();
  vals.push_back(x2d_1);
  vals.push_back(omega);
  vals.push_back(genP*MVABase/p_sbase); //output at system mva base
  vals.push_back(genQ*MVABase/p_sbase); //output at system mva base
}