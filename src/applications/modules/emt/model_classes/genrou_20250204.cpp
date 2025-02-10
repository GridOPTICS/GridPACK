#include <genrou_20250204.hpp>
// #include <genrou.hpp>
#include <gridpack/include/gridpack.hpp>
#include <constants.hpp>

/**
   Sauer and Pai - Power System Dynamics and Stability book
   Eqs. 3.148 - 3.159 - Page 42 with the following modifications
   - Stator flux differential equations ignored
   - Additional algebraic equations for machine three-phase currents
   - Speed sensitivity for mechanical torque Tm
*/
Genrou::Genrou(void)
{
  nxgen   = 12; // Number of variables for this model
  flux_speed_sensitivity = 0;
}

void Genrou::getnvar(int *nvar)
{
  if(integrationtype == EXPLICIT) nxgen = 6;
  *nvar = nxgen;
}


Genrou::~Genrou(void)
{
}

/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseGeneratorModel
 */
void Genrou::load(const boost::shared_ptr<gridpack::component::DataCollection> data, int idx)
{
  BaseEMTGenModel::load(data,idx); // load parameters in base generator model

  gridpack::ComplexType Zsource;
  
  // load parameters for the model type
  data->getValue(BUS_NUMBER, &bid);
    // load parameters for the model type
  if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H, &H, idx)) H = 0.0; // H
  if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0, &D, idx)) D = 0.0; // D
  data->getValue(GENERATOR_ZSOURCE,&Zsource,idx);
  Ra = real(Zsource);

  if (!data->getValue(GENERATOR_RESISTANCE, &Ra, idx)) Ra=0.0; // Ra

  if (!data->getValue(GENERATOR_XD, &Xd, idx)) Xd=0.0; // Xd
  if (!data->getValue(GENERATOR_XQ, &Xq, idx)) Xq=0.0; // Xq
  if (!data->getValue(GENERATOR_XDP, &Xdp, idx)) Xdp=0.0; // Xdp
  if (!data->getValue(GENERATOR_XDPP, &Xdpp, idx)) Xdpp=0.0; // Xdpp
  if (!data->getValue(GENERATOR_XL, &Xl, idx)) Xl=0.0; // Xl
  if (!data->getValue(GENERATOR_TDOP, &Tdop, idx)) Tdop=0.0; // Tdop
  if (!data->getValue(GENERATOR_TDOPP, &Tdopp, idx)) Tdopp=0.0; // Tdopp
  if (!data->getValue(GENERATOR_TQOPP, &Tqopp, idx)) Tqopp=0.0; // Tqopp
  if (!data->getValue(GENERATOR_S1, &S10, idx)) S10=0.067; // S10 TBD: check parser
  if (!data->getValue(GENERATOR_S12, &S12, idx)) S12=0.579; // S12 TBD: check parser
  if (!data->getValue(GENERATOR_XQP, &Xqp, idx)) Xqp=0.0; // Xqp
  if (!data->getValue(GENERATOR_XDPP, &Xqpp, idx)) Xqpp=Xdpp; // Xqpp 
  if (!data->getValue(GENERATOR_TQOP, &Tqop, idx)) Tqop=0.0; // Tqop

  L = Xdpp/OMEGA_S;
  
  // parameter conversion into EMT format
  L_ad = Xd - Xl;
  L_aq = Xq - Xl;
  L_l = Xl;
  L_d = L_ad + L_l;
  L_q = L_aq + L_l;
  L_fd = (Xdp - Xl) * L_ad / (L_ad - (Xdp - Xl));
  L_1q = (Xqp - Xl) * L_aq / (L_aq - (Xqp - Xl));

  double z = Xdpp - Xl;
  double y = L_ad * L_fd / (L_ad + L_fd);
  L_1d = y * z / (y - z);
  // ec_R1d = (y + ec_L1d) / Td0pp; // orig
  R_1d = (y + L_1d) / Tdopp / w_base; // Tdopp here is in second, the numerator is in pu


  z = Xdpp - Xl;
  y = L_aq * L_1q / (L_aq + L_1q);
  L_2q = y * z / (y - z);

  R_2q =  (y + L_2q) / Tqopp / w_base; 
  R_fd = (L_ad + L_fd) / Tdop / w_base; 
  R_1q = (L_aq + L_1q) / Tqop / w_base; 

  R_a = Ra;

  L_f1d = L_ad;

  L_ffd = L_ad * L_ad / (Xd - Xdp);
  L_11d =  L_ad * L_ad / (Xd - Xdpp);
  L_11q = L_aq * L_aq / (Xq - Xqp);
  L_22q = L_aq * L_aq / (Xq - Xdpp);

  w_r = 1.0; // rotor speed in pu

}

/**
 * Initialize generator model before calculation
 * @param [output] xin - array where initialized generator variables should be set
 */
void Genrou::init(gridpack::RealType* xin)
{
  double IGD,IGQ; // Machine currents in cartesian coordinates
  double Pg, Qg;  // Generator real and reactive power
  double dw0=0.0;  // Initial machine speed deviation
  gridpack::RealType *x = xin+offsetb; // generator array starts from this location

  Pg = pg/mbase;
  Qg = qg/mbase;

  VD = p_Vm0*cos(p_Va0);
  VQ = p_Vm0*sin(p_Va0);

  gridpack::ComplexType V = gridpack::ComplexType(VD,VQ);
  gridpack::ComplexType S = gridpack::ComplexType(Pg,Qg);
  gridpack::ComplexType I;
  gridpack::ComplexType Z = gridpack::ComplexType(Ra,Xdpp);
  gridpack::ComplexType Z1 = gridpack::ComplexType(Ra,Xq); // used in delta calculation
  gridpack::ComplexType E,E1;

  // Machine current
  I = conj(S/V);
  double Im = abs(I);
  double Ia = arg(I);
  
  double ia,ib,ic;
  ia = Im*sin(OMEGA_S*p_time + Ia);
  ib = Im*sin(OMEGA_S*p_time + Ia - TWOPI_OVER_THREE);
  ic = Im*sin(OMEGA_S*p_time + Ia + TWOPI_OVER_THREE);

  iabc[0] = ia;
  iabc[1] = ib;
  iabc[2] = ic;

  vabc[0] = p_va;
  vabc[1] = p_vb;
  vabc[2] = p_vc;
  
  // Calculation of machine angle delta and coordinate transformation
  // angle theta
  E1 = V + I*Z1;
  delta = arg(E1);
  // theta is behind delta by 90 degrees
  double theta = delta - PI/2.0;

  // Generator internal voltage on network reference frame
  E = V + I*Z;
  double Em = abs(E);
  double Eang = arg(E);

  // Network to machine reference frame transformation
  abc2dq0(vabc,p_time,theta,vdq0);
  abc2dq0(iabc,p_time,theta,idq0);
  
  double Vd, Vq, V0, Id, Iq, I0;

  Vd = vdq0[0];
  Vq = vdq0[1];
  V0 = vdq0[2];

  Id = idq0[0];
  Iq = idq0[1];
  I0 = idq0[2];

  psi_q = (-R_a * Id - Vd) / w_r; // initially w_r==1.0
  psi_d = (R_a * Iq + Vq) / w_r;

  i_fd = ( psi_d + (L_ad + L_l)*Id ) / L_ad;
  e_fd = R_fd * i_fd;

  psi_fd = L_ffd*i_fd - L_ad*Id;
  psi_1d = L_f1d*i_fd - L_ad*Id;
  psi_1q = L_11q*i_1q - L_aq*Iq;
  psi_2q = - L_aq*Iq;

  TE = k_te * (psi_d * Iq - psi_q * Id);
  TM = TE;
  delt_w_r = 0.0;

  // psid = Ra*Iq + Vq;
  // psiq = -Ra*Id - Vd;
  // psi0 = 0.0;

  // dw = dw0;

  // psi1d = psid + Xl*Id;
  // Eqp   = psi1d + (Xdp - Xl)*Id;

  // psi2q = psiq + Xl*Iq;
  // Edp = -psi2q - (Xqp - Xl)*Iq;

  // TM = psid*Iq - psiq*Id;

  // double param, LadIfd;
  // double dpsi1ddt;

  // param = (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl));
  
  // dpsi1ddt = psi1d + (Xdp - Xl)*Id - Eqp;

  // LadIfd = -Eqp - (Xd - Xdp)*(Id - param*dpsi1ddt);

  // Efd = -LadIfd;

  if(integrationtype != EXPLICIT) {
    // Initialized state variables
    x[0] = psid;
    x[1] = psiq;
    x[2] = psi0;
    x[3] = Eqp;
    x[4] = psi1d;
    x[5] = Edp;
    x[6] = psi2q;
    x[7] = delta;
    x[8] = dw;
    x[9] = iabc[0]*mbase/sbase;
    x[10] = iabc[1]*mbase/sbase;
    x[11] = iabc[2]*mbase/sbase;
  } else {
    x[0] = psid;
    x[1] = psiq;
    x[2] = psi0;
    x[3] = iabc[0]*mbase/sbase;
    x[4] = iabc[1]*mbase/sbase;
    x[5] = iabc[2]*mbase/sbase;

  }
}

/**
 * Write output from generators to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Genrou::serialWrite(char *string, const int bufsize,const char *signal)
{
  if(!strcmp(signal,"header")) {
    /* Print output header */
    sprintf(string,", %d_%s_V,%d_%s_Pg,%d_%s_delta, %d_%s_dw",busnum,id.c_str(),busnum,id.c_str(),busnum,id.c_str(),busnum,id.c_str());
    return true;
  } else if(!strcmp(signal,"monitor")) {
    /* Print output */
    double Vm, Pgen,dspd;
    if(p_online) {
      Vm = sqrt(vdq0[0]*vdq0[0] + vdq0[1]*vdq0[1]);
      dspd = dw;
    } else {
      Vm = p_Vm0;
      dspd = 0.0;
    }
    Pgen = p_online*(vdq0[0]*idq0[0] + vdq0[1]*idq0[1])*mbase/sbase;
    sprintf(string,", %6.5f,%6.5f,%6.5f, %6.5f",Vm,Pgen,delta,dspd);
    return true;
  }
  return false;
}

/**
 * Write out generator state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Genrou::write(const char* signal, char* string)
{
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto generators
 * @param values array containing generator state variables
*/
void Genrou::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values+offsetb; // generator array starts from this location

  if(integrationtype == EXPLICIT) {
    if(p_mode == XVECTOBUS) {
      psid  = x[0];
      psiq  = x[1];
      psi0  = x[2];
      iabc[0] = x[3];
      iabc[1] = x[4];
      iabc[2] = x[5];
    } else if(p_mode == XDOTVECTOBUS) {
      dpsid  = x[0];
      dpsiq  = x[1];
      dpsi0  = x[2];
      diabc[0] = x[3];
      diabc[1] = x[4];
      diabc[2] = x[5];
    }
  } else {
    if(p_mode == XVECTOBUS) {
      psid  = x[0];
      psiq  = x[1];
      psi0  = x[2];
      Eqp   = x[3];
      psi1d   = x[4];
      Edp     = x[5];
      psi2q = x[6];
      delta = x[7];
      dw = x[8];
      iabc[0] = x[9];
      iabc[1] = x[10];
      iabc[2] = x[11];
    } else if(p_mode == XDOTVECTOBUS) {
      dpsid  = x[0];
      dpsiq  = x[1];
      dpsi0  = x[2];
      dEqp   = x[3];
      dpsi1d   = x[4];
      dEdp     = x[5];
      dpsi2q = x[6];
      ddelta = x[7];
      ddw = x[8];
      diabc[0] = x[9];
      diabc[1] = x[10];
      diabc[2] = x[11];
    }
  }
}

/**
 * Return the values of the generator vector block
 * @param values: pointer to vector values
 * @return: false if generator does not contribute
 *        vector element
 */
void Genrou::vectorGetValues(gridpack::RealType *values)
{
  gridpack::RealType *f = values+offsetb; // generator array starts from this location

  if(p_mode == RESIDUAL_EVAL) {
    double theta = delta - PI/2.0;

    double tempd1,tempd2,tempq1,tempq2;
    tempd1 = (Xdpp - Xl)/(Xdp - Xl);
    tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
    tempq1 = (Xdpp - Xl)/(Xqp - Xl);
    tempq2 = (Xqp - Xdpp)/(Xqp - Xl);

    vabc[0] = p_va;
    vabc[1] = p_vb;
    vabc[2] = p_vc;

    // Network to machine reference frame transformation
    abc2dq0(vabc,p_time,theta,vdq0);
  
    double Vd, Vq, V0, Id, Iq, I0;

    Vd = vdq0[0];
    Vq = vdq0[1];
    V0 = vdq0[2];

    Id = (psid - tempd1*Eqp - tempd2*psi1d)/-Xdpp;
    Iq = (psiq + tempq1*Edp - tempq2*psi2q)/-Xdpp;
    I0 = psi0/-Xl;

    idq0[0] = Id;
    idq0[1] = Iq;
    idq0[2] = I0;

    if(hasExciter()) {
      Efd = getExciter()->getFieldVoltage();
    }
		    
    double dpsi1ddt;
    double param1 = (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl));

    dpsi1ddt = -psi1d + Eqp - (Xdp - Xl)*Id;

    double dpsi2qdt;
    double param2 = (Xqp - Xdpp)/((Xqp - Xl)*(Xqp - Xl));

    dpsi2qdt = -psi2q - Edp - (Xqp - Xl)*Iq;

    if(hasGovernor()) {
      TM = getGovernor()->getMechanicalPower();
    }

    double igen[3];
    dq02abc(idq0,p_time,theta,igen);

    if(integrationtype != EXPLICIT) {
      f[0] = OMEGA_S*(Ra*Id + (1 + flux_speed_sensitivity*dw)*psiq + Vd) - dpsid;
      f[1] = OMEGA_S*(Ra*Iq - (1 + flux_speed_sensitivity*dw)*psid + Vq) - dpsiq;
      f[2] = OMEGA_S*(Ra*I0 + V0) - dpsi0;
      f[3] = (-Eqp - (Xd - Xdp)*(Id - param1*-dpsi1ddt) + Efd)/Tdop - dEqp;
      f[4] = dpsi1ddt/Tdopp - dpsi1d;
      f[5] = (-Edp + (Xq - Xqp)*(Iq - param2*-dpsi2qdt))/Tqop - dEdp;
      f[6] = dpsi2qdt/Tqopp - dpsi2q;
      f[7] = dw*OMEGA_S - ddelta;
      f[8] = 1 / (2 *H) * ((TM - D*dw)/(1+dw) - (psid*Iq - psiq*Id)) - ddw;
      f[9] = igen[0]*mbase/sbase - iabc[0];
      f[10] = igen[1]*mbase/sbase - iabc[1];
      f[11] = igen[2]*mbase/sbase - iabc[2];
    } else {
      if(p_online) {
	f[0] = OMEGA_S*(Ra*Id + (1 + flux_speed_sensitivity*dw)*psiq + Vd) - dpsid;
	f[1] = OMEGA_S*(Ra*Iq - (1 + flux_speed_sensitivity*dw)*psid + Vq) - dpsiq;
	f[2] = OMEGA_S*(Ra*I0 + V0) - dpsi0;
	f[3] = igen[0]*mbase/sbase - iabc[0];
	f[4] = igen[1]*mbase/sbase - iabc[1];
	f[5] = igen[2]*mbase/sbase - iabc[2];
      } else {
	f[0] = psid;
	f[1] = psiq;
	f[2] = psi0;
	f[3] = iabc[0];
	f[4] = iabc[1];
	f[5] = iabc[2];
      }
    }
  }
}

  /**
   * Return the generator current injection 
   * @param [output] ia - phase a current
   * @param [output] ib - phase b current
   * @param [output] ic - phase c current
   */
void Genrou::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = p_online*iabc[0];
  *ib = p_online*iabc[1];
  *ic = p_online*iabc[2];
}

/**
 * Return the global location for the generator current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Genrou::getCurrentGlobalLocation(int *i_gloc)
{
  if(integrationtype == IMPLICIT) {
    *i_gloc = p_gloc + 9;
  } else {
    *i_gloc = p_gloc + 3;
  }
}


/**
 * Get number of matrix values contributed by generator
 * @return number of matrix values

 Non-zero pattern of the Jacobian (x denotes non-zero value)

          psid    psiq    psi0    Eqp    psi1d    Edp    psi2q    delta    omega    ia    ib    ic    va    vb    vc    Efd    PM
 eq. 0 |   x       x               x       x                       x         x                        x      x     x 
 eq. 1 |   x       x                               x       x       x         x                        x      x     x
 eq. 2 |                   x                                       x                                  x      x     x 
 eq. 3 |   x                       x       x                                                                             x                 
 eq. 4 |   x                       x       x         
 eq. 5 |           x                               x       x
 eq. 6 |           x                               x       x
 eq. 7 |                                                           x         x
 eq. 8 |   x       x               x       x       x       x                 x                                                x
 eq. 9 |   x       x       x       x       x       x       x       x                x
 eq.10 |   x       x       x       x       x       x       x       x                      x
 eq.11 |   x       x       x       x       x       x       x       x                            x
 
 Number of non-zero values = 9 + 9 + 5 + 4 + 3 + 3 + 3 + 2 + 7 + 9 + 9 + 9 = 69
 */
int Genrou::matrixNumValues()
{
  int numVals = 0;
  if(integrationtype == IMPLICIT) {
    numVals = 69 + 10;
    if(hasExciter()) numVals += 1;
    if(hasGovernor()) numVals += 1;
  } else if(integrationtype == EXPLICIT) {
    numVals = 15 + 12;
  }
  return numVals;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Genrou::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  *nvals = 0;
  int ctr = 0;
  if(integrationtype == EXPLICIT) {
    // Set up some indices
    int psid_idx = p_gloc;
    int psiq_idx = p_gloc+1;
    int psi0_idx = p_gloc+2;
    int ia_idx = p_gloc+3;
    int ib_idx = p_gloc+4;
    int ic_idx = p_gloc+5;

    if(!p_online) {
      rows[ctr]   = psid_idx; cols[ctr] = psid_idx;
      rows[ctr+1] = psiq_idx; cols[ctr+1] = psiq_idx;
      rows[ctr+2] = psi0_idx; cols[ctr+2] = psi0_idx;
      rows[ctr+3] = ia_idx;  cols[ctr+3] = ia_idx;
      rows[ctr+4] = ib_idx;  cols[ctr+4] = ib_idx;
      rows[ctr+5] = ic_idx;  cols[ctr+5] = ic_idx;

      values[ctr] = 1.0;
      values[ctr+1] = 1.0;
      values[ctr+2] = 1.0;
      values[ctr+3] = 1.0;
      values[ctr+4] = 1.0;
      values[ctr+5] = 1.0;

      ctr += 6;

      *nvals = ctr;

      return;
    }
      
    double dId_dpsid;
    double dIq_dpsiq;
    double dI0_dpsi0;

    dId_dpsid  = -1/Xdpp;
    dIq_dpsiq  = -1/Xdpp;
    dI0_dpsi0  = -1/Xl;
    
    double vabc[3];
    vabc[0] = p_va;
    vabc[1] = p_vb;
    vabc[2] = p_vc;

    double Tdq0[3][3];

    double theta = delta - PI/2.0;
    getTdq0(p_time,theta,Tdq0);

        // Derivatives w.r.t dpisd_dt
    // psid
    rows[ctr] = psid_idx; cols[ctr] = psid_idx;
    values[ctr] = OMEGA_S*Ra*dId_dpsid - shift;
  
    rows[ctr+1] = psid_idx; cols[ctr+1] = psiq_idx;
    values[ctr+1] = OMEGA_S*(1 + flux_speed_sensitivity*dw);
    
    ctr += 2;
    
    rows[ctr]   = psid_idx;
    cols[ctr]   = p_glocvoltage;
    values[ctr] = OMEGA_S*Tdq0[0][0];

    rows[ctr+1]   = psid_idx;
    cols[ctr+1]   = p_glocvoltage+1;
    values[ctr+1] = OMEGA_S*Tdq0[0][1];

    rows[ctr+2]   = psid_idx;
    cols[ctr+2]   = p_glocvoltage+2;
    values[ctr+2] = OMEGA_S*Tdq0[0][2];

    ctr += 3;
    
    assert(ctr == 5);
    
    // Derivatives w.r.t dpisq_dt

    rows[ctr] = psiq_idx; cols[ctr] = psid_idx;
    values[ctr] = OMEGA_S*-(1 + flux_speed_sensitivity*dw);
    
    rows[ctr+1] = psiq_idx; cols[ctr+1] = psiq_idx;
    values[ctr+1] = OMEGA_S*Ra*dIq_dpsiq - shift;
    
    ctr += 2;
    
    rows[ctr]   = psiq_idx;
    cols[ctr]   = p_glocvoltage;
    values[ctr] = OMEGA_S*Tdq0[1][0];
    
    rows[ctr+1]   = psiq_idx;
    cols[ctr+1]   = p_glocvoltage+1;
    values[ctr+1] = OMEGA_S*Tdq0[1][1];
    
    rows[ctr+2]   = psiq_idx;
    cols[ctr+2]   = p_glocvoltage+2;
    values[ctr+2] = OMEGA_S*Tdq0[1][2];
    
    ctr += 3;
    
    assert(ctr == 10);
    
    // Derivatives w.r.t dpsi0_dt
    
    rows[ctr] = psi0_idx;  cols[ctr] = psi0_idx;
    values[ctr] = OMEGA_S*Ra*dI0_dpsi0 - shift;
    
    ctr += 1;
    
    rows[ctr]   = psi0_idx;
    cols[ctr]   = p_glocvoltage;
    values[ctr] = OMEGA_S*Tdq0[2][0];
    
    rows[ctr+1]   = psi0_idx;
    cols[ctr+1]   = p_glocvoltage+1;
    values[ctr+1] = OMEGA_S*Tdq0[2][1];
    
    rows[ctr+2]   = psi0_idx;
    cols[ctr+2]   = p_glocvoltage+2;
    values[ctr+2] = OMEGA_S*Tdq0[2][2];
    
    ctr += 3;
    
    assert(ctr == 14);

    rows[ctr]   = ia_idx; cols[ctr]   = psid_idx;
    rows[ctr+1] = ia_idx; cols[ctr+1] = psiq_idx;
    rows[ctr+2] = ia_idx; cols[ctr+2] = psi0_idx;
    rows[ctr+3] = ia_idx; cols[ctr+3] = ia_idx;

    double Tdq0inv[3][3];
    getTdq0inv(p_time,theta,Tdq0inv);

    double scal = mbase/sbase;

    values[ctr]   = scal*(Tdq0inv[0][0]*dId_dpsid);
    values[ctr+1] = scal*(Tdq0inv[0][1]*dIq_dpsiq);
    values[ctr+2] = scal*(Tdq0inv[0][2]*dI0_dpsi0);
    values[ctr+3] = -1.0;

    ctr += 4;

    rows[ctr]   = ib_idx; cols[ctr]   = psid_idx;
    rows[ctr+1] = ib_idx; cols[ctr+1] = psiq_idx;
    rows[ctr+2] = ib_idx; cols[ctr+2] = psi0_idx;
    rows[ctr+3] = ib_idx; cols[ctr+3] = ib_idx;

    values[ctr]   = scal*(Tdq0inv[1][0]*dId_dpsid);
    values[ctr+1] = scal*(Tdq0inv[1][1]*dIq_dpsiq);
    values[ctr+2] = scal*(Tdq0inv[1][2]*dI0_dpsi0);
    values[ctr+3] = -1.0;

    ctr += 4;

    rows[ctr]   = ic_idx; cols[ctr]   = psid_idx;
    rows[ctr+1] = ic_idx; cols[ctr+1] = psiq_idx;
    rows[ctr+2] = ic_idx; cols[ctr+2] = psi0_idx;
    rows[ctr+3] = ic_idx; cols[ctr+3] = ic_idx;

    values[ctr]   = scal*(Tdq0inv[2][0]*dId_dpsid);
    values[ctr+1] = scal*(Tdq0inv[2][1]*dIq_dpsiq);
    values[ctr+2] = scal*(Tdq0inv[2][2]*dI0_dpsi0);
    values[ctr+3] = -1.0;

    ctr += 4;
    
  } else {
    // Set up some indices
    int psid_idx = p_gloc;
    int psiq_idx = p_gloc+1;
    int psi0_idx = p_gloc+2;
    int Eqp_idx  = p_gloc+3;
    int psi1d_idx = p_gloc+4;
    int Edp_idx  = p_gloc+5;
    int psi2q_idx = p_gloc+6;
    int delta_idx = p_gloc+7;
    int dw_idx  = p_gloc+8;
    int ia_idx  = p_gloc+9;
    int ib_idx = p_gloc+10;
    int ic_idx = p_gloc+11;

    double tempd1,tempd2,tempq1,tempq2;
    tempd1 = (Xdpp - Xl)/(Xdp - Xl);
    tempd2 = (Xdp - Xdpp)/(Xdp - Xl);
    tempq1 = (Xdpp - Xl)/(Xqp - Xl);
    tempq2 = (Xqp - Xdpp)/(Xqp - Xl);

    double dId_dpsid, dId_dEqp, dId_dpsi1d;
    double dIq_dpsiq, dIq_dEdp, dIq_dpsi2q;
    double dI0_dpsi0;

    dId_dpsid  = -1/Xdpp;
    dId_dEqp   =  tempd1/Xdpp;
    dId_dpsi1d =  tempd2/Xdpp;
    
    dIq_dpsiq  = -1/Xdpp;
    dIq_dEdp   = -tempq1/Xdpp;
    dIq_dpsi2q =  tempq2/Xdpp;
    
    dI0_dpsi0  = -1/Xl;
    
    double vabc[3];
    vabc[0] = p_va;
    vabc[1] = p_vb;
    vabc[2] = p_vc;

    double dvdq0_ddelta[3];
    double Tdq0[3][3],dTdq0_ddelta[3][3];

    double theta = delta - PI/2.0;
    getTdq0(p_time,theta,Tdq0);
    getdTdq0dtheta(p_time,theta,dTdq0_ddelta);

    matvecmult3x3(dTdq0_ddelta,vabc,dvdq0_ddelta);

    // Derivatives w.r.t dpisd_dt
    // psid
    rows[ctr] = psid_idx; cols[ctr] = psid_idx;
    values[ctr] = OMEGA_S*Ra*dId_dpsid - shift;
  
    rows[ctr+1] = psid_idx; cols[ctr+1] = psiq_idx;
    values[ctr+1] = OMEGA_S*(1 + flux_speed_sensitivity*dw);
    
    ctr += 2;
    
    rows[ctr] = psid_idx; cols[ctr] = Eqp_idx;
    values[ctr] = OMEGA_S*Ra*dId_dEqp;

    rows[ctr+1] = psid_idx; cols[ctr+1] = psi1d_idx;
    values[ctr+1] = OMEGA_S*Ra*dId_dpsi1d;
    
    ctr += 2;
    
    rows[ctr] = psid_idx; cols[ctr] = delta_idx;
    values[ctr] = OMEGA_S*dvdq0_ddelta[0];

    rows[ctr+1] = psid_idx; cols[ctr+1] = dw_idx;
    values[ctr+1] = OMEGA_S*flux_speed_sensitivity*psiq;
  
    ctr += 2;
  
    rows[ctr]   = psid_idx;
    cols[ctr]   = p_glocvoltage;
    values[ctr] = OMEGA_S*Tdq0[0][0];

    rows[ctr+1]   = psid_idx;
    cols[ctr+1]   = p_glocvoltage+1;
    values[ctr+1] = OMEGA_S*Tdq0[0][1];

    rows[ctr+2]   = psid_idx;
    cols[ctr+2]   = p_glocvoltage+2;
    values[ctr+2] = OMEGA_S*Tdq0[0][2];

    ctr += 3;
    
    assert(ctr == 9);
    
    // Derivatives w.r.t dpisq_dt

    rows[ctr] = psiq_idx; cols[ctr] = psid_idx;
    values[ctr] = OMEGA_S*-(1 + flux_speed_sensitivity*dw);
    
    rows[ctr+1] = psiq_idx; cols[ctr+1] = psiq_idx;
    values[ctr+1] = OMEGA_S*Ra*dIq_dpsiq - shift;
    
    ctr += 2;
    
    rows[ctr] = psiq_idx; cols[ctr] = Edp_idx;
    values[ctr] = OMEGA_S*Ra*dIq_dEdp;
    
    rows[ctr+1] = psiq_idx; cols[ctr+1] = psi2q_idx;
    values[ctr+1] = OMEGA_S*Ra*dIq_dpsi2q;
    
    ctr += 2;
    
    rows[ctr] = psiq_idx; cols[ctr] = delta_idx;
    values[ctr] = OMEGA_S*dvdq0_ddelta[1];
    
    rows[ctr+1] = psiq_idx; cols[ctr+1] = dw_idx;
    values[ctr+1] = OMEGA_S*-psid*flux_speed_sensitivity;
    
    ctr += 2;
    
    rows[ctr]   = psiq_idx;
    cols[ctr]   = p_glocvoltage;
    values[ctr] = OMEGA_S*Tdq0[1][0];
    
    rows[ctr+1]   = psiq_idx;
    cols[ctr+1]   = p_glocvoltage+1;
    values[ctr+1] = OMEGA_S*Tdq0[1][1];
    
    rows[ctr+2]   = psiq_idx;
    cols[ctr+2]   = p_glocvoltage+2;
    values[ctr+2] = OMEGA_S*Tdq0[1][2];
    
    ctr += 3;
    
    assert(ctr == 18);
    
    // Derivatives w.r.t dpsi0_dt
    
    rows[ctr] = psi0_idx;  cols[ctr] = psi0_idx;
    values[ctr] = OMEGA_S*Ra*dI0_dpsi0 - shift;
    
    ctr += 1;
    
    rows[ctr] = psi0_idx; cols[ctr] = delta_idx;
    values[ctr] = OMEGA_S*dvdq0_ddelta[2];
    
    ctr += 1;
    
    rows[ctr]   = psi0_idx;
    cols[ctr]   = p_glocvoltage;
    values[ctr] = OMEGA_S*Tdq0[2][0];
    
    rows[ctr+1]   = psi0_idx;
    cols[ctr+1]   = p_glocvoltage+1;
    values[ctr+1] = OMEGA_S*Tdq0[2][1];
    
    rows[ctr+2]   = psi0_idx;
    cols[ctr+2]   = p_glocvoltage+2;
    values[ctr+2] = OMEGA_S*Tdq0[2][2];
    
    ctr += 3;
    
    assert(ctr == 23);
    
    // Derivatives w.r.t. dEqp_dt
    double dpsi1ddt_dpsi1d, dpsi1ddt_dEqp, dpsi1ddt_dpsid;
    double param1 = (Xdp - Xdpp)/((Xdp - Xl)*(Xdp - Xl));
    
    dpsi1ddt_dpsid = -(Xdp - Xl)*dId_dpsid;
    dpsi1ddt_dEqp  = 1.0 -(Xdp - Xl)*dId_dEqp;
    dpsi1ddt_dpsi1d = -1.0 -(Xdp - Xl)*dId_dpsi1d;
    
    rows[ctr] = Eqp_idx;  cols[ctr] = psid_idx;
    values[ctr] = (-(Xd - Xdp)*(dId_dpsid - param1*-dpsi1ddt_dpsid))/Tdop;
    
    rows[ctr+1] = Eqp_idx; cols[ctr+1] = Eqp_idx;
    values[ctr+1] = (-1.0 -(Xd - Xdp)*(dId_dEqp - param1*(-dpsi1ddt_dEqp)))/Tdop -shift;
    
    rows[ctr+2] = Eqp_idx; cols[ctr+2] = psi1d_idx;
    values[ctr+2] = (-(Xd - Xdp)*(dId_dpsi1d -param1*-dpsi1ddt_dpsi1d))/Tdop;
    
    ctr += 3;
    if(hasExciter()) {
      int Efd_idx;
      double Efd;
      Efd = getExciter()->getFieldVoltage(&Efd_idx);
      if(Efd_idx >= 0) { // Exciter is using implicit integration
	rows[ctr] = Eqp_idx; cols[ctr] = Efd_idx;
	values[ctr] = 1.0/Tdop;
	ctr += 1;
      }
    }
    
    // Derivative of dpsi1d_dt
    rows[ctr] = psi1d_idx;  cols[ctr] = psid_idx;
    values[ctr] = (dpsi1ddt_dpsid)/Tdopp;
    
    rows[ctr+1] = psi1d_idx; cols[ctr+1] = Eqp_idx;
    values[ctr+1] = (dpsi1ddt_dEqp)/Tdopp;
    
    rows[ctr+2] = psi1d_idx; cols[ctr+2] = psi1d_idx;
    values[ctr+2] = (dpsi1ddt_dpsi1d)/Tdopp - shift;
    
    ctr += 3;
    
    // Derivatives w.r.t. dEdp_dt
    double dpsi2qdt_dpsi2q, dpsi2qdt_dEdp, dpsi2qdt_dpsiq;
    double param2 = (Xqp - Xdpp)/((Xqp - Xl)*(Xqp - Xl));
    
    dpsi2qdt_dpsiq = -(Xqp - Xl)*dIq_dpsiq;
    dpsi2qdt_dEdp  = -1.0 -(Xqp - Xl)*dIq_dEdp;
    dpsi2qdt_dpsi2q = -1.0 -(Xqp - Xl)*dIq_dpsi2q;
    
    rows[ctr] = Edp_idx;  cols[ctr] = psiq_idx;
    values[ctr] = ((Xq - Xqp)*(dIq_dpsiq - param2*-dpsi2qdt_dpsiq))/Tqop;
    
    rows[ctr+1] = Edp_idx; cols[ctr+1] = Edp_idx;
    values[ctr+1] = (-1.0 +(Xq - Xqp)*(dIq_dEdp - param2*(-dpsi2qdt_dEdp)))/Tqop - shift;
    
    rows[ctr+2] = Edp_idx; cols[ctr+2] = psi2q_idx;
    values[ctr+2] = ((Xq - Xqp)*(dIq_dpsi2q -param2*-dpsi2qdt_dpsi2q))/Tqop;
    
    ctr += 3;
    
    // Derivative of dpsi2q_dt
    rows[ctr] = psi2q_idx;  cols[ctr] = psiq_idx;
    values[ctr] = (dpsi2qdt_dpsiq)/Tqopp;
    
    rows[ctr+1] = psi2q_idx; cols[ctr+1] = Edp_idx;
    values[ctr+1] = (dpsi2qdt_dEdp)/Tqopp;
    
    rows[ctr+2] = psi2q_idx; cols[ctr+2] = psi2q_idx;
    values[ctr+2] = (dpsi2qdt_dpsi2q)/Tqopp - shift;
    
    ctr += 3;
    
    // Derivative of ddelta_dt
    rows[ctr] = delta_idx; cols[ctr] = delta_idx;
    values[ctr] = -shift;
    
    rows[ctr+1] = delta_idx; cols[ctr+1] = dw_idx;
    values[ctr+1] = OMEGA_S;
    
    ctr += 2;
    
    // derivative of ddw_dt
    double Minv = 1 / (2*H);
    double Id,Iq,I0;
    
    Id = (psid - tempd1*Eqp - tempd2*psi1d)/-Xdpp;
    Iq = (psiq + tempq1*Edp - tempq2*psi2q)/-Xdpp;
    I0 = psi0/-Xl;
    
    rows[ctr] = dw_idx; cols[ctr] = psid_idx;
    values[ctr] = Minv*-(Iq - psiq*dId_dpsid);
    
    rows[ctr+1] = dw_idx; cols[ctr+1] = psiq_idx;
    values[ctr+1] = Minv*-(psid*dIq_dpsiq - Id);
    
    rows[ctr+2] = dw_idx; cols[ctr+2] = Eqp_idx;
    values[ctr+2] = Minv*-(-psiq*dId_dEqp);
    
    rows[ctr+3] = dw_idx; cols[ctr+3] = psi1d_idx;
    values[ctr+3] = Minv*-(-psiq*dId_dpsi1d);
    
    rows[ctr+4] = dw_idx; cols[ctr+4] = Edp_idx;
    values[ctr+4] = Minv*-(psid*dIq_dEdp);
    
    rows[ctr+5] = dw_idx; cols[ctr+5] = psi2q_idx;
    values[ctr+5] = Minv*-(psid*dIq_dpsi2q);
    
    int TM_idx;
    if(hasGovernor()) {
      TM = getGovernor()->getMechanicalPower(&TM_idx);
    }
    
    rows[ctr+6] = dw_idx; cols[ctr+6] = dw_idx;
    values[ctr+6] = Minv*(-D*(1/(1+dw) - (TM - D*dw)/((1+dw)*(1+dw)))) -shift;
    
    ctr += 7;
    
    if(hasGovernor()) {
      // Partial derivatives w.r.t Governor
      if(TM_idx >= 0) {
	// >=0 indiciates governor is using implicit method
	// so we also need to set the derivative
	rows[ctr] = dw_idx; cols[ctr] = TM_idx;
	values[ctr] = Minv*(1/(1+dw));
	ctr += 1;
      }
    }
    
    // derivative of currents iabc
    
    rows[ctr]   = ia_idx; cols[ctr]   = psid_idx;
    rows[ctr+1] = ia_idx; cols[ctr+1] = psiq_idx;
    rows[ctr+2] = ia_idx; cols[ctr+2] = psi0_idx;
    rows[ctr+3] = ia_idx; cols[ctr+3] = Eqp_idx;
    rows[ctr+4] = ia_idx; cols[ctr+4] = psi1d_idx;
    rows[ctr+5] = ia_idx; cols[ctr+5] = Edp_idx;
    rows[ctr+6] = ia_idx; cols[ctr+6] = psi2q_idx;
    rows[ctr+7] = ia_idx; cols[ctr+7] = delta_idx;
    rows[ctr+8] = ia_idx; cols[ctr+8] = ia_idx;
    
    double Tdq0inv[3][3],dTdq0inv_ddelta[3][3];
    double scal = mbase/sbase;
    
    getTdq0inv(p_time,theta,Tdq0inv);
    getdTdq0invdtheta(p_time,theta,dTdq0inv_ddelta);
    
    values[ctr]   = scal*(Tdq0inv[0][0]*dId_dpsid);
    values[ctr+1] = scal*(Tdq0inv[0][1]*dIq_dpsiq);
    values[ctr+2] = scal*(Tdq0inv[0][2]*dI0_dpsi0);
    
    values[ctr+3] = scal*(Tdq0inv[0][0]*dId_dEqp);
    values[ctr+4] = scal*(Tdq0inv[0][0]*dId_dpsi1d);
    
    values[ctr+5] = scal*(Tdq0inv[0][1]*dIq_dEdp);
    values[ctr+6] = scal*(Tdq0inv[0][1]*dIq_dpsi2q);
    
    values[ctr+7] = scal*(dTdq0inv_ddelta[0][0]*Id + dTdq0inv_ddelta[0][1]*Iq + dTdq0inv_ddelta[0][2]*I0);
    
    values[ctr+8] = -1.0;
    
    ctr += 9;
    
    rows[ctr]   = ib_idx; cols[ctr]   = psid_idx;
    rows[ctr+1] = ib_idx; cols[ctr+1] = psiq_idx;
    rows[ctr+2] = ib_idx; cols[ctr+2] = psi0_idx;
    rows[ctr+3] = ib_idx; cols[ctr+3] = Eqp_idx;
    rows[ctr+4] = ib_idx; cols[ctr+4] = psi1d_idx;
    rows[ctr+5] = ib_idx; cols[ctr+5] = Edp_idx;
    rows[ctr+6] = ib_idx; cols[ctr+6] = psi2q_idx;
    rows[ctr+7] = ib_idx; cols[ctr+7] = delta_idx;
    rows[ctr+8] = ib_idx; cols[ctr+8] = ib_idx;
    
    values[ctr]   = scal*(Tdq0inv[1][0]*dId_dpsid);
    values[ctr+1] = scal*(Tdq0inv[1][1]*dIq_dpsiq);
    values[ctr+2] = scal*(Tdq0inv[1][2]*dI0_dpsi0);
    
    values[ctr+3] = scal*(Tdq0inv[1][0]*dId_dEqp);
    values[ctr+4] = scal*(Tdq0inv[1][0]*dId_dpsi1d);
    
    values[ctr+5] = scal*(Tdq0inv[1][1]*dIq_dEdp);
    values[ctr+6] = scal*(Tdq0inv[1][1]*dIq_dpsi2q);
    
    values[ctr+7] = scal*(dTdq0inv_ddelta[1][0]*Id + dTdq0inv_ddelta[1][1]*Iq + dTdq0inv_ddelta[1][2]*I0);
    
    values[ctr+8] = -1.0;
    
    ctr += 9;
    
    rows[ctr]   = ic_idx; cols[ctr]   = psid_idx;
    rows[ctr+1] = ic_idx; cols[ctr+1] = psiq_idx;
    rows[ctr+2] = ic_idx; cols[ctr+2] = psi0_idx;
    rows[ctr+3] = ic_idx; cols[ctr+3] = Eqp_idx;
    rows[ctr+4] = ic_idx; cols[ctr+4] = psi1d_idx;
    rows[ctr+5] = ic_idx; cols[ctr+5] = Edp_idx;
    rows[ctr+6] = ic_idx; cols[ctr+6] = psi2q_idx;
    rows[ctr+7] = ic_idx; cols[ctr+7] = delta_idx;
    rows[ctr+8] = ic_idx; cols[ctr+8] = ic_idx;
    
    values[ctr]   = scal*(Tdq0inv[2][0]*dId_dpsid);
    values[ctr+1] = scal*(Tdq0inv[2][1]*dIq_dpsiq);
    values[ctr+2] = scal*(Tdq0inv[2][2]*dI0_dpsi0);
    
    values[ctr+3] = scal*(Tdq0inv[2][0]*dId_dEqp);
    values[ctr+4] = scal*(Tdq0inv[2][0]*dId_dpsi1d);
    
    values[ctr+5] = scal*(Tdq0inv[2][1]*dIq_dEdp);
    values[ctr+6] = scal*(Tdq0inv[2][1]*dIq_dpsi2q);
    
    values[ctr+7] = scal*(dTdq0inv_ddelta[2][0]*Id + dTdq0inv_ddelta[2][1]*Iq + dTdq0inv_ddelta[2][2]*I0);
    
    values[ctr+8] = -1.0;
    
    ctr += 9;
  }
  *nvals = ctr;
}


/**
 * Set Jacobian values
 * @param values a 2-d array of Jacobian block for the bus
 */
bool Genrou::setJacobian(gridpack::RealType **values)
{
  return true;
}

/**
 * Returns the initial field voltage (Efd(t0))
 * @param [out] Efd0 - Initial field voltage
 */
double Genrou::getInitialFieldVoltage()
{
  return Efd;
}

/**
   Prestep function
*/
void Genrou::preStep(double time ,double timestep)
{
  if(integrationtype != EXPLICIT) {
    return;
  }

  double x[9],f[9]; // Yuan: why initializing 9-dim array instead of 5 or 7 depending on the model's order?

  x[0] = psi_d;
  x[1] = psi_q;
  x[2] = psi_fd;
  x[3] = psi_1d;
  x[4] = psi_1q;
  x[5] = psi_2q;
  x[6] = delt_w_r;
  x[7] = delta;

  double theta = delta - PI/2.0;
  
  vabc[0] = p_va;
  vabc[1] = p_vb;
  vabc[2] = p_vc;

  // Network to machine reference frame transformation
  abc2dq0(vabc,time,theta,vdq0);
  
  double Vd, Vq, V0, Id, Iq, I0;

  Vd = vdq0[0];
  Vq = vdq0[1];
  V0 = vdq0[2];
  
  // Id = (psid - tempd1*Eqp - tempd2*psi1d)/-Xdpp;
  // Iq = (psiq + tempq1*Edp - tempq2*psi2q)/-Xdpp;
  // I0 = psi0/-Xl;
  
  // idq0[0] = Id;
  // idq0[1] = Iq;
  // idq0[2] = I0;

  if(hasExciter()) {
    Efd = getExciter()->getFieldVoltage();
  }
		    
       
  // Step-2: calculate predictor dx/dt       
  // d_psi_d = w_base * (Vd + psi_q*w_r - (R_a*(L_f1d^2*psi_d + 
  //           L_11d*L_ad*psi_fd - L_11d*L_ffd*psi_d - 
  //               L_ad*L_f1d*psi_1d - L_ad*L_f1d*psi_fd + L_ad*L_ffd*psi_1d)) / 
  //               (L_11d*L_ad^2 + L_ad*L_f1d^2 - 2*L_ad^2*L_f1d + L_ad^2*L_ffd + 
  //               L_f1d^2*L_l - L_11d*L_ad*L_ffd - L_11d*L_ffd*L_l) );
  
  f[0] = w_base * (Vd + psi_q*w_r - (R_a*(L_f1d*L_f1d *psi_d + 
    L_11d*L_ad*psi_fd - L_11d*L_ffd*psi_d - 
        L_ad*L_f1d*psi_1d - L_ad*L_f1d*psi_fd + L_ad*L_ffd*psi_1d)) / 
        (L_11d*L_ad*L_ad + L_ad*L_f1d*L_f1d - 2*L_ad*L_ad*L_f1d + L_ad*L_ad*L_ffd + 
        L_f1d*L_f1d*L_l - L_11d*L_ad*L_ffd - L_11d*L_ffd*L_l) );

  // d_psi_q = w_base * ( Vq - psi_d*w_r - (R_a*(L_aq^2*psi_q - 
  //               L_aq^2*psi_1q - L_aq^2*psi_2q - L_11q*L_22q*psi_q + 
  //               L_11q*L_aq*psi_2q + L_22q*L_aq*psi_1q)) / 
  //               (L_11q*L_aq^2 + L_22q*L_aq^2 + L_aq^2*L_l - L_aq^3 - 
  //               L_11q*L_22q*L_aq - L_11q*L_22q*L_l) );

  f[1] = w_base * ( Vq - psi_d*w_r - (R_a*(L_aq*L_aq*psi_q - 
    L_aq*L_aq*psi_1q - L_aq*L_aq*psi_2q - L_11q*L_22q*psi_q + 
    L_11q*L_aq*psi_2q + L_22q*L_aq*psi_1q)) / 
    (L_11q*L_aq*L_aq + L_22q*L_aq*L_aq + L_aq*L_aq*L_l - L_aq*L_aq*L_aq - 
    L_11q*L_22q*L_aq - L_11q*L_22q*L_l) );

  // d_psi_fd = w_base * ( e_fd + (R_fd*(L_ad^2*psi_1d - L_ad^2*psi_fd - 
  //               L_11d*L_ad*psi_d + L_11d*L_ad*psi_fd + L_11d*L_l*psi_fd + 
  //               L_ad*L_f1d*psi_d - L_ad*L_f1d*psi_1d - L_f1d*L_l*psi_1d)) / 
  //               (L_11d*L_ad^2 + L_ad*L_f1d^2 - 2*L_ad^2*L_f1d + L_ad^2*L_ffd + 
  //               L_f1d^2*L_l - L_11d*L_ad*L_ffd - L_11d*L_ffd*L_l) );

  f[2] = w_base * ( e_fd + (R_fd*(L_ad*L_ad*psi_1d - L_ad*L_ad*psi_fd - 
    L_11d*L_ad*psi_d + L_11d*L_ad*psi_fd + L_11d*L_l*psi_fd + 
    L_ad*L_f1d*psi_d - L_ad*L_f1d*psi_1d - L_f1d*L_l*psi_1d)) / 
    (L_11d*L_ad*L_ad + L_ad*L_f1d*L_f1d - 2*L_ad*L_ad*L_f1d + L_ad*L_ad*L_ffd + 
    L_f1d*L_f1d*L_l - L_11d*L_ad*L_ffd - L_11d*L_ffd*L_l) );

  // d_psi_1d = w_base * ( -(R_1d*(L_ad^2*psi_1d - L_ad^2*psi_fd - 
  //               L_ad*L_f1d*psi_d + L_ad*L_ffd*psi_d + 
  //               L_ad*L_f1d*psi_fd - L_ad*L_ffd*psi_1d + 
  //               L_f1d*L_l*psi_fd - L_ffd*L_l*psi_1d)) / 
  //               (L_11d*L_ad^2 + L_ad*L_f1d^2 - 2*L_ad^2*L_f1d + 
  //               L_ad^2*L_ffd + L_f1d^2*L_l - L_11d*L_ad*L_ffd - 
  //               L_11d*L_ffd*L_l) );
  
  f[3] = w_base * ( -(R_1d*(L_ad*L_ad*psi_1d - L_ad*L_ad*psi_fd - 
    L_ad*L_f1d*psi_d + L_ad*L_ffd*psi_d + 
    L_ad*L_f1d*psi_fd - L_ad*L_ffd*psi_1d + 
    L_f1d*L_l*psi_fd - L_ffd*L_l*psi_1d)) / 
    (L_11d*L_ad*L_ad + L_ad*L_f1d*L_f1d - 2*L_ad*L_ad*L_f1d + 
    L_ad*L_ad*L_ffd + L_f1d*L_f1d*L_l - L_11d*L_ad*L_ffd - 
    L_11d*L_ffd*L_l) );

  // d_psi_1q = w_base * ( (R_1q*(L_aq^2*psi_q - L_aq^2*psi_1q - 
  //               L_22q*L_aq*psi_q + L_22q*L_aq*psi_1q + 
  //               L_22q*L_l*psi_1q - L_aq*L_l*psi_2q)) / 
  //               (L_11q*L_aq^2 + L_22q*L_aq^2 + L_aq^2*L_l - 
  //               L_aq^3 - L_11q*L_22q*L_aq - L_11q*L_22q*L_l) );

  f[4] = w_base * ( (R_1q*(L_aq*L_aq*psi_q - L_aq*L_aq*psi_1q - 
    L_22q*L_aq*psi_q + L_22q*L_aq*psi_1q + 
    L_22q*L_l*psi_1q - L_aq*L_l*psi_2q)) / 
    (L_11q*L_aq*L_aq + L_22q*L_aq*L_aq + L_aq*L_aq*L_l - 
    L_aq*L_aq*L_aq - L_11q*L_22q*L_aq - L_11q*L_22q*L_l) );

  // d_psi_2q = w_base * ( (R_2q*(L_aq^2*psi_q - L_aq^2*psi_2q - 
  //               L_11q*L_aq*psi_q + L_11q*L_aq*psi_2q + 
  //               L_11q*L_l*psi_2q - L_aq*L_l*psi_1q)) / 
  //               (L_11q*L_aq^2 + L_22q*L_aq^2 + L_aq^2*L_l - 
  //               L_aq^3 - L_11q*L_22q*L_aq - L_11q*L_22q*L_l) );

  f[5] = w_base * ( (R_2q*(L_aq*L_aq*psi_q - L_aq*L_aq*psi_2q - 
    L_11q*L_aq*psi_q + L_11q*L_aq*psi_2q + 
    L_11q*L_l*psi_2q - L_aq*L_l*psi_1q)) / 
    (L_11q*L_aq*L_aq + L_22q*L_aq*L_aq + L_aq*L_aq*L_l - 
    L_aq*L_aq*L_aq - L_11q*L_22q*L_aq - L_11q*L_22q*L_l) );

  f[6] = (T_m - T_e - K_D * delt_w_r) / (2*H);

  if(hasGovernor()) {
    TM = getGovernor()->getMechanicalPower();
  }
  
  f[5] = 1 / (2 *H) * ((TM - D*dw)/(1+dw) - (psid*Iq - psiq*Id)); 

  for(int i=0; i < 6; i++) {
    x[i] += timestep*f[i]; // Forward Euler update
  }

  Eqp   = x[0];
  psi1d = x[1];
  Edp   = x[2];
  psi2q = x[3];
  delta = x[4];
  dw    = x[5];
  
}

/**
   Post step function
*/
void Genrou::postStep(double time)
{
  if(integrationtype != EXPLICIT) {
    return;
  }

}


