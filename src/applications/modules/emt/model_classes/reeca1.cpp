/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   reeca1.cpp
 *  
 * @brief REECA1 model implementation 
 *
 *
 */

#include <reeca1.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

Reeca1::Reeca1(void)
{
  omega_g = 1.0;
  nxexc = 0;
}

Reeca1::~Reeca1(void)
{
}

bool Reeca1::getVoltageDip(double Vt)
{
  if(Vt < Vdip || Vt > Vup)
    return true;
  else return false;
}

void Reeca1::getnvar(int *nvar)
{
  *nvar = nxexc;
}

/**
 * Load parameters from DataCollection object into exciter model
 * @param data collection of exciter parameters from input files
 * @param index of exciter on bus
 * TODO: might want to move this functionality to BaseExciterModel
 */
void Reeca1::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTExcModel::load(data,idx); // load parameters in base exciter model
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

  // Vt filter output used in division
  Vt_filter_lowcap_blk.setparams(1.0,0.01,2.0);

  // Initialize VDL1 (V-P) and VDL2 (V-Q)
  double uin[4], yin[4];
  uin[0] = Vp1; yin[0] = Ip1;
  uin[1] = Vp2; yin[1] = Ip2;
  uin[2] = Vp3; yin[2] = Ip3;
  uin[3] = Vp4; yin[3] = Ip4;

  VDL1.setparams(4,uin,yin);

  uin[0] = Vq1; yin[0] = Iq1;
  uin[1] = Vq2; yin[1] = Iq2;
  uin[2] = Vq3; yin[2] = Iq3;
  uin[3] = Vq4; yin[3] = Iq4;

  VDL2.setparams(4,uin,yin);

  // Pref limiter block
  Pref_limit_blk.setparams(1.0,-1000.0,1000.0,dPmin,dPmax);

  // Pord block
  Pord_blk.setparams(1.0,Tpord,Pmin,Pmax,-1000.0,1000.0);

  // Ipcmd limiter block
  Ipcmd_limit_blk.setparams(1.0,Ipmin,Ipmax);

  Voltage_dip_prev = false;
  thld_timer = -1.0;
}

/**
 * Initialize exciter model before calculation
 * @param [output] values - array where initialized exciter variables should be set
 */
void Reeca1::init(gridpack::RealType* xin) 
{
  gridpack::RealType *x = xin+offsetb; // exciter array starts from this location
  double Vt = sqrt(VD*VD + VQ*VQ);

  double Pord_blk_in,Iq_Qflag;
  
  // Initialize Vt filter block
  Vt_filter = Vt_filter_blk.init_given_u(Vt);

  // Ipcmd related blocks initialization
  // Initial value of Ipcmd provided by generator controller
  Pord = Vt_filter*Ipcmd;
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
  Iqinj_sw = 0;

  // Is there a voltage dip?
  // Should not be during initialization
  Voltage_dip = getVoltageDip(Vt);

  // Initialize Vref0 = Vt if Vref0 = 0 in the data file
  if(fabs(Vref0) < 1e-6) Vref0 = Vt;
  
  if(Voltage_dip == 0) {
    Iqinj = 0;
  } else {
    V_err = V_err_deadband.getoutput(Vref0 - Vt_filter);
    Iqv = Kqv*V_err;
    Iqinj = Iqv_limit_blk.getoutput(Iqv);
  }
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
    if(!VFLAG) {
      Q_Pfflag = V_err_PI_blk_in + Vt_filter - Vbias;
    } else {
      Q_PI_blk_in = Q_PI_blk.init_given_y(V_err_PI_blk_in + Vt_filter);
      Q_Pfflag = Q_PI_blk_in + Qgen;
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
 * Write output from exciters to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Reeca1::serialWrite(char *string, const int bufsize,const char *signal)
{
  return false;
}

/**
 * Write out exciter state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Reeca1::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto exciters
 * @param values array containing exciter state variables
*/
void Reeca1::setValues(gridpack::RealType *val)
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
void Reeca1::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // exciter array starts from this location
}

/**
   Prestep function
*/
void Reeca1::preStep(double time ,double timestep)
{

}

/**
   Poststep function
*/
void Reeca1::postStep(double time)
{
}

/**
 * Get number of matrix values contributed by exciter
 * @return number of matrix values
 */
int Reeca1::matrixNumValues()
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
void Reeca1::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr = 0;
  *nvals = ctr;		     
}

/**
 * Update the event function values
 */
void Reeca1::eventFunction(const double&t,gridpack::RealType *state,std::vector<gridpack::RealType >& evalues)
{
  int offset    = getLocalOffset();
} 

/**
 * Event handler
 */
void Reeca1::eventHandlerFunction(const bool *triggered, const double& t, gridpack::RealType *state)
{
  int offset    = getLocalOffset();
}

/**
 * Set event
 */
void Reeca1::setEvent(gridpack::math::RealDAESolver::EventManagerPtr eman)
{
  gridpack::math::RealDAESolver::EventPtr e(new Reeca1Event(this));

  eman->add(e);
}

void Reeca1Event::p_update(const double& t,gridpack::RealType *state)
{
  p_exc->eventFunction(t,state,p_current);
}

void Reeca1Event::p_handle(const bool *triggered, const double& t, gridpack::RealType *state)
{
  p_exc->eventHandlerFunction(triggered,t,state);
}
