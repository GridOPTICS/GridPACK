#include <transformer.hpp>
#include <gridpack/include/gridpack.hpp>

Transformer::Transformer(void)
{
  nxbranch  = 6;
}

Transformer::~Transformer(void)
{
}

/**
 * Return the number of variables
 * @param [output] nvar - number of variables
 */
void Transformer::getnvar(int *nvar)
{
  *nvar = nxbranch;
}


/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseEMTGeneratorModel
 */
void Transformer::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int i)
{
  // Get line
  BaseEMTBranchModel::load(data,i);

  data->getValue(BRANCH_R,&R,i);
  data->getValue(BRANCH_X,&X,i);
  data->getValue(BRANCH_TAP,&tap,i);

  if(fabs(R) > 1e-6) {
    p_hasResistance = true;
  }

  if(fabs(X) > 1e-6) {
    p_hasInductance = true;
  }

  R1 = R; L1 = X/OMEGA_S;
  R0 = 3*R1; L0 = 3*L1;

  double Rs = R1;
  double Rm = 0.0;
  p_R[0][0] = p_R[1][1] = p_R[2][2] = Rs;
  p_R[0][1] = p_R[1][0] = Rm;
  p_R[0][2] = p_R[2][0] = Rm;
  p_R[1][2] = p_R[2][1] = Rm;
  
  double Ls = L1;
  double Lm = 0.0;
  p_L[0][0] = p_L[1][1] = p_L[2][2] = Ls;
  p_L[0][1] = p_L[1][0] = Lm;
  p_L[0][2] = p_L[2][0] = Lm;
  p_L[1][2] = p_L[2][1] = Lm;
}

/**
 * Set up branch model before calculation
 */
void Transformer::setup()
{

}

/**
 * Initialize branch model before calculation
 * @param [output] values - array where initialized branch variables should be set
 */
void Transformer::init(gridpack::RealType *values)
{
  double *x = values + offsetb;
  
  double Vmf,Vaf,Vmt,Vat;
  double VDf,VQf,VDt,VQt;
  fbus->getInitialVoltage(&Vmf,&Vaf);
  tbus->getInitialVoltage(&Vmt,&Vat);

  VDf = Vmf*cos(Vaf);
  VQf = Vmf*sin(Vaf);
  VDt = Vmt*cos(Vat);
  VQt = Vmt*sin(Vat);

  gridpack::ComplexType Vf = gridpack::ComplexType(VDf,VQf);
  gridpack::ComplexType Vt = gridpack::ComplexType(VDt,VQt);

  gridpack::ComplexType Zline = gridpack::ComplexType(R,X);

  gridpack::ComplexType Iline = (Vf/tap - Vt)/Zline;

  double Ilinem,Ilinea;
  
  Ilinem = abs(Iline);
  Ilinea = arg(Iline);

  x[0] = Ilinem*sin(Ilinea)/tap;
  x[1] = Ilinem*sin(Ilinea - 2.0*PI/3.0)/tap;
  x[2] = Ilinem*sin(Ilinea + 2.0*PI/3.0)/tap;

  x[3] = Ilinem*sin(Ilinea);
  x[4] = Ilinem*sin(Ilinea - 2.0*PI/3.0);
  x[5] = Ilinem*sin(Ilinea + 2.0*PI/3.0);
}

/**
 * Write output from branchs to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Transformer::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

/**
 * Write out branch state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Transformer::write(const char* signal, char* string)
{
}

/**
 * Return the branch from bus current injection 
 * @param [output] ia - phase a current
 * @param [output] ib - phase b current
 * @param [output] ic - phase c current
 */
void Transformer::getFromBusCurrent(double *ia, double *ib, double *ic)
{
  *ia = ibr_from[0];
  *ib = ibr_from[1];
  *ic = ibr_from[2];
}

/**
 * Return the branch to bus current injection 
 * @param [output] ia - phase a current
 * @param [output] ib - phase b current
 * @param [output] ic - phase c current
 */
void Transformer::getToBusCurrent(double *ia, double *ib, double *ic)
{
  *ia = ibr_to[0];
  *ib = ibr_to[1];
  *ic = ibr_to[2];
}


/**
 * Return the location for the current in the local branch array
 * @param [output] i_loc - location for the first current variable in the local branch array
 */
void Transformer::getCurrentLocalLocation(int *i_loc)
{
  *i_loc = offsetb;
}

/**
 * Return the global location for the branch current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Transformer::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc;
}

/**
 * Get number of matrix values contributed by branch
 * @return number of matrix values
 */
int Transformer::matrixNumValues()
{
  int numvals = 0;
  if(p_hasInductance) {
    // fval = (vf - R*i - vt) - L*di_dt
    // dfval_dvf = I => 3 entries
    // dfval_di  = -R - sL => 9 entries
    // dfval_dvt = -I => 3 entries
    // 6 entries for to bus current equation
    // Total 21 entries
    numvals += 21;
    
  } else {
    // fval = vf - R*i - vt
    // dfval_dvf = I => 3 entries
    // dfval_di = -R => 9 entries
    // dfval_dvt = -I => 3 entries
    // six entries for to bus current equation
    // Total 21 entries
    numvals += 21;
  }

  return numvals;
}

/**
 * Return values from Jacobian matrix
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void Transformer::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
{
  int ctr=0;
  int j,k;
  int vf_gloc,vt_gloc,i_gloc;

  getCurrentGlobalLocation(&i_gloc);
    
  if(p_hasInductance) {
    for(j=0; j < 3; j++) {
      for(k=0; k < 3; k++) {
	rows[ctr]   = i_gloc + j;
	cols[ctr]   = i_gloc + k;
	values[ctr] = -p_R[j][k] - shift*p_L[j][k];
	ctr++;
      }
    }
  } else {
    for(j=0; j < 3; j++) {
      for(k=0; k < 3; k++) {
	rows[ctr]   = i_gloc + j;
	cols[ctr]   = i_gloc + k;
	values[ctr] = -p_R[j][k];
	ctr++;
      }
    }
  }

  fbus->getVoltageGlobalLocation(&vf_gloc);
  tbus->getVoltageGlobalLocation(&vt_gloc);


  // Partial derivatives w.r..t voltages
  for(j=0; j < 3; j++) {
    rows[ctr] = i_gloc + j;
    cols[ctr] = vf_gloc + j;
    values[ctr] = 1.0/tap;

    rows[ctr+1] = i_gloc + j;
    cols[ctr+1] = vt_gloc + j;
    values[ctr+1] = -1.0;

    ctr += 2;
  }

  // Partials for to bus current equations
  for(j = 0; j < 3; j++) {
    // partial derivative w.r.t. from bus current
    rows[ctr] = i_gloc + 3 + j;
    cols[ctr] = i_gloc + j;
    values[ctr] = -1.0;

    ctr++;

    // partial derivative w.r.t to bus current
    rows[ctr] = i_gloc + 3 + j;
    cols[ctr] = i_gloc + 3 + j;
    values[ctr] = 1.0/tap;

    ctr++;
  }
  
  *nvals = ctr;

}

/**
 * Return vector values from the branch model 
 * @param values - array of returned values
 *
 * Note: This function is used to return the entries in vector,
 * for e.g., the entries in the residual vector from the branch
 * object
   */
void Transformer::vectorGetValues(gridpack::RealType *values)
{
  double *f = values + offsetb;

  if(p_mode == RESIDUAL_EVAL) {
    double vf[3],vt[3];

    fbus->getVoltages(&vf[0],&vf[1],&vf[2]);
    tbus->getVoltages(&vt[0],&vt[1],&vt[2]);

    double vf_minus_vt[3];
    
    vf_minus_vt[0] = vf[0]/tap - vt[0];
    vf_minus_vt[1] = vf[1]/tap - vt[1];
    vf_minus_vt[2] = vf[2]/tap - vt[2];
    
    if(p_hasInductance) {
      double Ribr[3],Ldidt[3];
      // vf - R*ibr_from - vt - L*didt = 0
      matvecmult3x3(p_R, ibr_from,Ribr); // fval1 = R*ibr
      matvecmult3x3(p_L,dibrf_dt,Ldidt); // fval2 = L*didt
      f[0] = vf_minus_vt[0] - Ribr[0] - Ldidt[0];
      f[1] = vf_minus_vt[1] - Ribr[1] - Ldidt[1];
      f[2] = vf_minus_vt[2] - Ribr[2] - Ldidt[2];

      // To bus current equation
      // To bus current is same as from bus current
      f[3] = ibr_to[0]/tap - ibr_from[0];
      f[4] = ibr_to[1]/tap - ibr_from[1];
      f[5] = ibr_to[2]/tap - ibr_from[2];
    } else {
      double Ribr[3];

      matvecmult3x3(p_R,ibr_from,Ribr);
      
      f[0] = -Ribr[0] + vf_minus_vt[0];
      f[1] = -Ribr[1] + vf_minus_vt[1];
      f[2] = -Ribr[2] + vf_minus_vt[2];

      // To bus current equation
      f[3] = ibr_to[0]/tap - ibr_from[0];
      f[4] = ibr_to[1]/tap - ibr_from[1];
      f[5] = ibr_to[2]/tap - ibr_from[2];
    }
  }
}

/**
 * Pass solution vector values to the branch object
 * @param values - array of returned values
 *
 * Note: This function is used to pass the entries in vector
 * to the branch object,
 * for e.g., the state vector values for this branch
 */
void Transformer::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values + offsetb;

  if(p_mode == XVECTOBUS) {
    ibr_from[0] = x[0];
    ibr_from[1] = x[1];
    ibr_from[2] = x[2];
    ibr_to[0]   = x[3];
    ibr_to[1]   = x[4];
    ibr_to[2]   = x[5];
  } else if(p_mode == XDOTVECTOBUS) {
    dibrf_dt[0] = x[0];
    dibrf_dt[1] = x[1];
    dibrf_dt[2] = x[2];
    dibrt_dt[0] = x[3];
    dibrt_dt[1] = x[4];
    dibrt_dt[2] = x[5];
  }
}

/**
   Prestep function
*/
void Transformer::preStep(double time, double timestep)
{
}

/**
   Poststep function
*/
void Transformer::postStep(double time)
{
}
