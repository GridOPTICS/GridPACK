/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -----------------------------------------------------------
/**
 * @file   regcb1.cpp
 * @author Shrirang Abhyankar
 * @Added: January 1, 2023
 * 
 * @brief  WECC generic generator/converter REGCB1 model
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
#include "regcb1.hpp"

/**
 *   TRANSFORM -
 *   Does a reference frame transformation, i.e.
 *   Converter values from one reference frame to another
 *   through a rotation angle theta
 *   Input
 *     xin - x coord (real part) of input reference frame
 *     yin - y coord (imag part) of input reference frame
 * theta - rotation angle
 *   
 *   Output
 *     xout - x coord (real part) of output reference frame
 *     yout - y coord (imag part) of output reference frame
 *
 *    The reference frame transformation is given by
      xout = xin*cos(theta) - yin*sin(theta)
      yout = xin*sin(theta) + yin*cos(theta)
   
 **/
static void transform(double xin, double yin, double theta,double *xout, double *yout)
{
  *xout =  xin*cos(theta) - yin*sin(theta);
  *yout =  xin*sin(theta) + yin*cos(theta);
}

/**
 *  Basic constructor
 */
gridpack::dynamic_simulation::Regcb1Generator::Regcb1Generator(void)
{
  domega_g = 0.0;
  Thetapitch = 0.0;
  omega_ref = 1.0;
}

/**
 *  Basic destructor
 */
gridpack::dynamic_simulation::Regcb1Generator::~Regcb1Generator(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 */
void gridpack::dynamic_simulation::Regcb1Generator::load(
    boost::shared_ptr<gridpack::component::DataCollection>
    data, int idx)
{
  data->getValue(BUS_NUMBER,&p_bus_num);
  data->getValue(CASE_SBASE,&p_sbase);
  data->getValue(GENERATOR_ID,&p_gen_id,idx);
  
  if (!data->getValue(GENERATOR_ZSOURCE, &Zsource, idx)) Zsource = gridpack::ComplexType(0.01,0.02);
  if (!data->getValue(GENERATOR_PG, &p_pg,idx)) p_pg = 0.0;
  if (!data->getValue(GENERATOR_QG, &p_qg,idx)) p_qg = 0.0;
  if (!data->getValue(GENERATOR_STAT, &p_status,idx)) p_status = 0;
  if (!data->getValue(GENERATOR_MBASE, &p_mbase, idx))  p_mbase = 100.0; // MBase
  if(fabs(p_mbase) < 1e-6) p_mbase = p_pg; // Set mbase = p_pg if it is not given in file

  if (!data->getValue(GENERATOR_REGC_RATEFLAG , &RateFlag, idx)) RateFlag = 0;
  if (!data->getValue(GENERATOR_REGC_PQFLAG , &PQFlag, idx)) PQFlag = 0; 
  if (!data->getValue(GENERATOR_REGC_TG  ,    &Tg, idx))     Tg = 0.02; 
  if (!data->getValue(GENERATOR_REGC_RRPWR  , &Rrpwr, idx))  Rrpwr = 10.0; 
  if (!data->getValue(GENERATOR_REGC_TFLTR  , &Tfltr, idx))  Tfltr = 0.02;
  if (!data->getValue(GENERATOR_REGC_TE  , &Te, idx))  Te = 0.02; 
  if (!data->getValue(GENERATOR_REGC_IQRMAX , &Iqrmax, idx)) Iqrmax = 999.0; 
  if (!data->getValue(GENERATOR_REGC_IQRMIN , &Iqrmin, idx)) Iqrmin = -999.0;
  if (!data->getValue(GENERATOR_REGC_IMAX , &Imax, idx)) Imax = 1.6; 

  // Set up blocks

  // transfer function blocks
  Ip_blk.setparams(1.0,Tg);
  
  Iq_blk.setparams(-1.0,Tg);
  
  if(p_qg > 0.0) {
    // Upper limit active with Qg > 0
    Iq_blk.setdxlimits(-1000.0,Iqrmax);
  } else {
    // Lower limit active when Qg < 0
    Iq_blk.setdxlimits(Iqrmin,1000.0);
  }
  // Lower limit of 0.01 applied on output of Vt_filter block
  Vt_filter_blk.setparams(1.0,Tfltr,-1000.0,1000.0,0.01,2.0);

  Ed_blk.setparams(1.0,Te);
  Eq_blk.setparams(1.0,Te);
}

/**
 * Initialize generator model before calculation
 * @param Vm voltage magnitude
 * @param Va voltage angle
 */
void gridpack::dynamic_simulation::Regcb1Generator::init(double Vm,
    double Va, double ts)
{
  double Vq;

  theta = Va; // save to local variable for later use

  // Change of base from system MVA base to Machine base
  p_pg *= p_sbase/p_mbase;
  p_qg *= p_sbase/p_mbase;

  // Real and imaginary components of voltage
  VR = Vm*cos(Va);
  VI = Vm*sin(Va);

  // Transformation such that Vm aligns with VR
  // Rotation by -Va
  transform(VR,VI,-Va,&Vt,&Vq);

  Ip =  p_pg/Vt;
  Iq = -p_qg/Vt;

  // Internal voltage source
  gridpack::ComplexType E,I,V;
  V = gridpack::ComplexType(Vt,0.0);
  I = gridpack::ComplexType(Ip,Iq);
  E = V + Zsource*I;

  Ed = E.real();
  Eq = E.imag();
  
  // Initialize blocks
  // Assume no limits are hit
  if(RateFlag == 0) {
    Ipcmd = Ip_blk.init_given_y(Ip);
  } else {
    Ipcmd = Ip_blk.init_given_y(Ip*Vt)/Vt;
  }
  
  Iqcmd = Iq_blk.init_given_y(Iq);

  Vt_filter = Vt_filter_blk.init_given_u(Vt);

  Ed = Ed_blk.init_given_y(Ed);
  Eq = Eq_blk.init_given_y(Eq);
  
  double Qext, Pord;
  double domega_t = 0.0;

  // Initialize electrical controller model
  if (p_hasExciter){
    p_exciter = getExciter();
    p_exciter->setVterminal(Vm);
    p_exciter->setGeneratorPower(p_pg,p_qg);
    p_exciter->setIpcmdIqcmd(Ipcmd, Iqcmd); 

    p_exciter->init(Vm, Va, ts);
	
    Pref = p_exciter->getPref( );
    Qext = p_exciter->getQext( );
    Pord = p_exciter->getPord( );
  }

  // Initialize plant controller model
  if (p_hasPlantController){
    p_plant = getPlantController();
    p_plant->setGenPQV(p_pg, p_qg, Vm);
    p_plant->setPrefQext(Pref, Qext);
    p_plant->setExtBusNum(p_bus_num);
    
    p_plant->init(Vm, Va, ts);
  }

  // Initialize torque controller model
  if(p_hasTorqueController) {
    p_torquecontroller = getTorqueController();
    p_torquecontroller->setPref0(Pref);
    p_torquecontroller->setPelec(Pref);

    p_torquecontroller->init(Vm, Va, ts);
  }

  // Initialize drive train model
  if(p_hasDriveTrainModel) {
    if(p_hasTorqueController) {
      omega_ref = p_torquecontroller->getOmegaref();
    }
    p_drivetrainmodel = getDriveTrainModel();
    
    p_drivetrainmodel->setOmegaref(omega_ref);
    p_drivetrainmodel->setTelec(p_pg/omega_ref);

    p_drivetrainmodel->init(Vm, Va, ts);
  }

  double Theta0 = 0.0;
  // Initialize aerodynamic model
  if(p_hasAeroDynamicModel) {
    p_aerodynamicmodel = getAeroDynamicModel();

    double Tm;
    if(p_hasDriveTrainModel) {
      domega_t = p_drivetrainmodel->getTurbineSpeedDeviation();
      Tm = p_drivetrainmodel->getTmech();
    }
    p_aerodynamicmodel->setTurbineSpeedDeviation(domega_t);
    p_aerodynamicmodel->setPmech(Tm*(1+domega_t));

    p_aerodynamicmodel->init(Vm, Va, ts);

    Taero = p_aerodynamicmodel->getTaero();
    Theta0 = p_aerodynamicmodel->getTheta();
  }

  // Initialize pitch controller model
  if(p_hasPitchController) {
    p_pitchcontroller = getPitchController();
    p_pitchcontroller->setPord(p_pg);
    p_pitchcontroller->setPord0(p_pg);

    if(p_hasDriveTrainModel) {
      domega_t = p_drivetrainmodel->getTurbineSpeedDeviation();
    }
    if(p_hasTorqueController) {
      omega_ref = p_torquecontroller->getOmegaref();
    }
    p_pitchcontroller->setTurbineSpeedDeviation(domega_t);
    p_pitchcontroller->setOmegaref(omega_ref);
    p_pitchcontroller->setTheta(Theta0);

    p_pitchcontroller->init(Vm, Va, ts);
  }
}

/**
 * Return contribution to Norton current
 * @return contribution to Norton vector
 */
gridpack::ComplexType gridpack::dynamic_simulation::Regcb1Generator::INorton()
{
  gridpack::ComplexType E,V,Iout;

  E = gridpack::ComplexType(Ed,Eq);
  V = gridpack::ComplexType(Vt,0.0);

  Iout = (E - V)/Zsource;
  
  Ipout = Iout.real();
  Iqout = Iout.imag();
  
  transform(Ipout,Iqout,theta,&Irout,&Iiout);

  // Scaled to system MVAbase
  Irout *= p_mbase/p_sbase;
  Iiout *= p_mbase/p_sbase;

  p_INorton = gridpack::ComplexType(Irout,Iiout);
    
  return p_INorton;
}


/**
 * Return Norton impedence
 * @return value of Norton impedence
 */
gridpack::ComplexType gridpack::dynamic_simulation::Regcb1Generator::NortonImpedence()
{
  double B = 0.0;
  double G = 0.0;
  gridpack::ComplexType Y_a(G, B);
  return Y_a;
}


/**
 * Predict part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Regcb1Generator::predictor_currentInjection(bool flag)
{
  p_INorton = INorton();
}

void gridpack::dynamic_simulation::Regcb1Generator::computeModel(double t_inc, IntegrationStage int_flag,bool flag)
{
  double Pg,Qg;
  double Qext;
  bool Vdip; // Voltage dip flag for electrical controller
  
  double domega_t = 0.0;
  double Pord = 0.0,Preal;
  double Pref_plant;
  double Ip_blk_in;

  // Pg, Qg on machine MVAbase
  Pg = Vt*Ipout;
  Qg = -Vt*Iqout;

  if(p_hasDriveTrainModel) {
    p_drivetrainmodel = getDriveTrainModel();
    domega_g = p_drivetrainmodel->getGeneratorSpeedDeviation();
    domega_t = p_drivetrainmodel->getTurbineSpeedDeviation();
  }

  if (p_hasPlantController){
    p_plant = getPlantController();
    p_plant->setGenPQV(Pg, Qg, Vt);
    p_plant->setBusFreq(busfreq);

    if(int_flag == PREDICTOR) {
      p_plant->predictor(t_inc, flag);
    } else {
      p_plant->corrector(t_inc,flag);
    }
		
    Pref = Pref_plant = p_plant->getPref();
    Qext = p_plant->getQext();
  }

  if(p_hasTorqueController) {
    p_torquecontroller = getTorqueController();
    p_torquecontroller->setPref0(Pref);
    p_torquecontroller->setPelec(Pg);
    p_torquecontroller->setGeneratorSpeedDeviation(domega_g);

    Vdip = getExciter()->getVoltageDip(Vt);

    p_torquecontroller->setVdip(Vdip);

    if(int_flag == PREDICTOR) {
      p_torquecontroller->predictor(t_inc,flag);
    } else {
      p_torquecontroller->corrector(t_inc,flag);
    }
    
    Pref = p_torquecontroller->getPref();
    omega_ref = p_torquecontroller->getOmegaref();
  }

  //get ipcmd and iqcmd from the exctier REECA1 MODEL
  if (p_hasExciter){
    p_exciter = getExciter();
    p_exciter->setPrefQext(Pref, Qext);
    p_exciter->setVterminal(Vt);

    p_exciter->setGeneratorPower(Pg,Qg);

    if(int_flag == PREDICTOR) {
      p_exciter->predictor(t_inc, flag);
    } else {
      p_exciter->corrector(t_inc,flag);
    }
		
    Ipcmd = p_exciter->getIpcmd();
    Iqcmd = p_exciter->getIqcmd();

    Pord = p_exciter->getPord();
  }

  if(p_hasPitchController) {
    p_pitchcontroller = getPitchController();
    p_pitchcontroller->setPord(Pord);
    p_pitchcontroller->setPord0(Pref_plant);
    p_pitchcontroller->setTurbineSpeedDeviation(domega_g); // using domega_g instead of domega_g gives better response (why?)
    p_pitchcontroller->setOmegaref(omega_ref);

    if(int_flag == PREDICTOR) {
      p_pitchcontroller->predictor(t_inc,flag);
    } else {
      p_pitchcontroller->corrector(t_inc,flag);
    }

    Thetapitch = p_pitchcontroller->getTheta();
  }

  if(p_hasAeroDynamicModel) {
    p_aerodynamicmodel = getAeroDynamicModel();
    p_aerodynamicmodel->setTheta(Thetapitch);
    p_aerodynamicmodel->setTurbineSpeedDeviation(domega_t);

    if(int_flag == PREDICTOR) {
      p_aerodynamicmodel->predictor(t_inc,flag);
    } else {
      p_aerodynamicmodel->corrector(t_inc,flag);
    }

    Taero = p_aerodynamicmodel->getTaero();
  }

  if(p_hasDriveTrainModel) {
    p_drivetrainmodel = getDriveTrainModel();

    p_drivetrainmodel->setTelec(Pg/(1+domega_g));
    p_drivetrainmodel->setTmech(Taero);

    if(int_flag == PREDICTOR) {
      p_drivetrainmodel->predictor(t_inc,flag);
    } else {
      p_drivetrainmodel->corrector(t_inc,flag);
    }      
  }

  Vt_filter = Vt_filter_blk.getoutput(Vt,t_inc,int_flag,true);

  if(RateFlag) {
    Ip_blk_in = Ipcmd*Vt_filter;
    if(Preal > 0) {
      Ip_blk.setdxlimits(-1000.0,Rrpwr);
    } else {
      Ip_blk.setdxlimits(-Rrpwr,1000.0);
    } 
    Preal = Ip_blk.getoutput(Ip_blk_in,t_inc,int_flag,true);
    Ip = Preal/Vt_filter;
  } else {
    if(Ipcmd > 0.0) {
      Ip_blk.setdxlimits(-1000.0,Rrpwr);
    } else {
      Ip_blk.setdxlimits(-Rrpwr,1000.0);
    }
    Ip = Ip_blk.getoutput(Ipcmd,t_inc,int_flag,true);
  }
				   
  Iq = Iq_blk.getoutput(Iqcmd,t_inc,int_flag,true);

  gridpack::ComplexType V = gridpack::ComplexType(Vt,0.0);
  gridpack::ComplexType Ipq = gridpack::ComplexType(Ip,Iq);

  gridpack::ComplexType E;

  E = V + Zsource*Ipq;

  Ed = Ed_blk.getoutput(E.real(),t_inc,int_flag,true);
  Eq = Eq_blk.getoutput(E.imag(),t_inc,int_flag,true);
  
}

/**
 * Predict new state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Regcb1Generator::predictor(
    double t_inc, bool flag)
{
  computeModel(t_inc,PREDICTOR,flag);
}

/**
 * Corrector part calculate current injections
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Regcb1Generator::corrector_currentInjection(bool flag)
{
  p_INorton = INorton();
}

/**
 * Correct state variables for time step
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::Regcb1Generator::corrector(
    double t_inc, bool flag)
{
  computeModel(t_inc,CORRECTOR,flag);
}

bool gridpack::dynamic_simulation::Regcb1Generator::tripGenerator()
{
	return false;
}

/**
* return true if modify the generator parameters successfully
* input controlTyp: 0: GFI mp adjust; 1: GFI mq adjust; 2: GFI Pset adjust; 3: GFI Qset adjust; others: invalid 
* input newParValScaletoOrg:  GFI new parameter scale factor to the very initial parameter value at the begining of dynamic simulation
* 
*/
bool gridpack::dynamic_simulation::Regcb1Generator::applyGeneratorParAdjustment(int controlType, double newParValScaletoOrg){
	
  return false;
	
}

/**
 * Set voltage on each generator
 */
void gridpack::dynamic_simulation::Regcb1Generator::setVoltage(
    gridpack::ComplexType voltage)
{
  Vt    = abs(voltage);
  theta = atan2(imag(voltage), real(voltage));
  VR    = real(voltage);
  VI    = imag(voltage);
}

/**
 * Set frequency on each generator, frequency is perunit
 */
void gridpack::dynamic_simulation::Regcb1Generator::setFreq(double dFreq)
{
   busfreq = dFreq;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
bool gridpack::dynamic_simulation::Regcb1Generator::serialWrite(
    char* string, const int bufsize, const char *signal)
{
  bool ret = false;
  string[0] = '\0';
  if (!strcmp(signal,"watch")) {
    if(getWatch()) {
      double Pg,Qg;

      // Pg, Qg on system MVAbase
      Pg = Vt*Ipout*p_mbase/p_sbase; 
      Qg = -Vt*Iqout*p_mbase/p_sbase;

      sprintf(string,",%12.6f,%12.6f, %12.6f, %12.6f, %12.6f ",
          Vt,Pg, Qg, busfreq,Pref);
      if(p_hasPitchController) {
        sprintf(string + strlen(string),",%12.6f,%12.6f ", Thetapitch,omega_ref);
      }
      if(p_hasDriveTrainModel) {
        sprintf(string + strlen(string),",%12.6f ", 1+domega_g);
      }
      if(p_hasAeroDynamicModel) {
        sprintf(string + strlen(string),",%12.6f ", Taero);
      }

      ret = true; 
    }
  } else if(!strcmp(signal,"watch_header")) {
    if(getWatch()) {
      char buf[128];
      std::string tag;
      if(p_gen_id[0] != ' ') {
        tag = p_gen_id;
      } else {
        tag = p_gen_id[1];
      }
      sprintf(buf,", %d_%s_V,%d_%s_Pg, %d_%s_Qg, %d_%s_freq, %d_%s_Pref",
          p_bus_num,tag.c_str(),p_bus_num,tag.c_str(),
          p_bus_num,tag.c_str(),p_bus_num,tag.c_str(), p_bus_num, tag.c_str());
      if(p_hasPitchController) {
        sprintf(buf + strlen(buf),", %d_%s_Tpitch,%d_%s_Omegaref",
            p_bus_num, tag.c_str(),p_bus_num, tag.c_str());
      }
      if(p_hasDriveTrainModel) {
        sprintf(buf + strlen(buf),", %d_%s_Omega_g",p_bus_num, tag.c_str());
      }
      if(p_hasAeroDynamicModel) {
        sprintf(buf + strlen(buf),", %d_%s_Taero",p_bus_num, tag.c_str());
      }


      if (strlen(buf) <= bufsize) {
        sprintf(string,"%s",buf);
        ret = true;
      } else {
        ret = false;
      }
    }
  }
  return ret;
}

/**
 * return a vector containing any generator values that are being
 * watched
 * @param vals vector of watched values
 */
void gridpack::dynamic_simulation::Regcb1Generator::getWatchValues(
    std::vector<double> &vals)
{
  vals.clear();
  vals.push_back(Ipout);
  vals.push_back(Iqout);

  double Pg,Qg;

  // Pg, Qg on system MVAbase
  Pg = Vt*Ipout*p_mbase/p_sbase;
  Qg = -Vt*Iqout*p_mbase/p_sbase;
  
  if (p_generatorObservationPowerSystemBase){
    vals.push_back(Pg);  //output at system mva base
    vals.push_back(Qg);  //output at system mva base
  }else{
    vals.push_back(Pg);  //output at generator mva base
    vals.push_back(Qg);  //output at generator mva base
  }
}

/**
 * Set internal state parameter in generator
 * @param name character string corresponding to state variable
 * @param value new value for state parameter
 * @return false if no variable corresponding to name is found
 */   
bool gridpack::dynamic_simulation::Regcb1Generator::setState(std::string name,
    double value)
{
  return false;
}

/**
 * Get internal state parameter in generator
 * @param name character string corresponding to state variable
 * @param value current value for state parameter
 * @return false if no variable corresponding to name is found
 */   
bool gridpack::dynamic_simulation::Regcb1Generator::getState(std::string name,
    double *value)
{
  return false;
}
