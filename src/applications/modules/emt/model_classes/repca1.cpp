/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   repca1.cpp
 *  
 * @brief REPCA1 model implementation 
 *
 *
 */

#include <repca1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Repca1::Repca1(void)
{
  nxplant = 0;
}

Repca1::~Repca1(void)
{
}

void Repca1::getnvar(int *nvar)
{
  *nvar = nxplant;
}

void Repca1::getPrefQext(double *Prefout, double *Qextout)
{
  *Prefout = Pref;
  *Qextout = Qref;
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to BaseExciterModel
 */
void Repca1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTPlantControllerModel::load(data,idx); // load parameters in base exciter model
  
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
 * @param [output] values - array where initialized exciter variables should be set
 */
void Repca1::init(gridpack::RealType* xin) 
{
  gridpack::RealType *x = xin+offsetb; // exciter array starts from this location
  double Vt = sqrt(VD*VD + VQ*VQ);

  // Get Initial Pref and Qref
  getElectricalController()->getInitialPrefQext(&Pref, &Qref);
  Plant_ref = Pref;

  // Get machine power
  getGenerator()->getInitialPower(&Pg, &Qg);

  // Convert to machine mvabase
  Pg *= sbase/mbase;
  Qg *= sbase/mbase;

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

/**
 * Write output from exciters to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Repca1::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out exciter state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Repca1::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto exciters
 * @param values array containing exciter state variables
*/
void Repca1::setValues(gridpack::RealType *val)
{
  gridpack::RealType *values = val+offsetb; // exciter array starts from this location

  if(p_mode == XVECTOBUS) {
  } else if(p_mode == XDOTVECTOBUS) {
  }
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Repca1::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // exciter array starts from this location
}

/**
   Prestep function
*/
void Repca1::preStep(double time ,double timestep)
{
  double Pgen,Qgen, Vd, Vq, Vt;
  double vabc[3],vdq0[3],theta;
  
  if(integrationtype != EXPLICIT) return;

  vabc[0] = p_va;
  vabc[1] = p_vb;
  vabc[2] = p_vc;

  theta = getGenerator()->getAngle();
  
  abc2dq0(vabc,time,theta,vdq0);

  Vd = vdq0[0];
  Vq = vdq0[1];

  Vt = sqrt(Vd*Vd + Vq*Vq);

  if(Vt < Vfrz) {
    Vfreeze = true;
  }
  else Vfreeze = false;

  bool updateState = !Vfreeze; // Do not update state (freeze) when Vfreeze is true
  
  // Pref part
  double ferr,dP;

  Freq = getGenerator()->getFreq();

  getGenerator()->getPower(time,&Pg,&Qg);
  // Convert power to mbase
  Pg *= sbase/mbase;
  Qg *= sbase/mbase;
  
  if(FreqFLAG) {
    ferr = Freq_ref - Freq;
    ferr = Freqerr_deadband.getoutput(ferr);
    dP = std::max(0.0,Ddn*ferr) + std::min(0.0,Dup*ferr);
    // ***********
    // Need to use Pbranch if given, using Pg
    // ***********
    Pbranch = Pg;
    
    Pbranch_filter_blk_out = Pbranch_filter_blk.getoutput(Pbranch,timestep,true);

    Freqerr_limiter_out = Freqerr_limiter.getoutput(Plant_ref - Pbranch_filter_blk_out + dP);

    Pref_PI_blk_out = Pref_PI_blk.getoutput(Freqerr_limiter_out,timestep,true);

    Pref = Pref_filter_blk.getoutput(Pref_PI_blk_out,timestep,true);
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

    Qbranch_filter_blk_out = Qbranch_filter_blk.getoutput(Qbranch,timestep,true);
    y_RefFLAG = Qref - Qbranch_filter_blk_out;
  } else {
    V_filter_blk_out = V_filter_blk.getoutput(y_VCompFLAG,timestep,true);
    y_RefFLAG = Vref - V_filter_blk_out;
  }

  // VQ error deadband block
  VQerr_deadband_out = VQerr_deadband.getoutput(y_RefFLAG);

  // VQ error limiter block
  VQerr_limiter_out = VQerr_limiter.getoutput(VQerr_deadband_out);

  // Qref PI control
  Qref_PI_blk_out = Qref_PI_blk.getoutput(VQerr_limiter_out,timestep,updateState);

  // Qref lead lag block
  Qref = Qref_leadlag_blk.getoutput(Qref_PI_blk_out,timestep,true);

}

/**
   Poststep function
*/
void Repca1::postStep(double time)
{
}

/**
 * Get number of matrix values contributed by exciter
 * @return number of matrix values
 */
int Repca1::matrixNumValues()
{
  return 0;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Repca1::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;
  *nvals = ctr;		     
}

/**
 * Update the event function values
 */
void Repca1::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues)
{
  int offset    = getLocalOffset();
} 

/**
 * Event handler
 */
void Repca1::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  int offset    = getLocalOffset();
}

/**
 * Set event
 */
void Repca1::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  gridpack::math::RealDAESolver::EventPtr e(new Repca1Event(this));

  eman->add(e);
}

void Repca1Event::p_update(const double& t,gridpack::RealType *state)
{
  p_plant->eventFunction(t,state,p_current);
}

void Repca1Event::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  p_plant->eventHandlerFunction(triggered,t,state);
}
