#include <lumpedline.hpp>
#include <gridpack/include/gridpack.hpp>

Lumpedline::Lumpedline(void)
{
  nxbranch  = 3;
}

Lumpedline::~Lumpedline(void)
{
}

/**
 * Return the number of variables
 * @param [output] nvar - number of variables
 */
void Lumpedline::getnvar(int *nvar)
{
  *nvar = nxbranch;
}


/**
 * Load parameters from DataCollection object into generator model
 * @param data collection of generator parameters from input files
 * @param index of generator on bus
 * TODO: might want to move this functionality to BaseEMTGeneratorModel
 */
void Lumpedline::load(const boost::shared_ptr<gridpack::component::DataCollection>
        data, int i)
{
  // Get line
  BaseEMTBranchModel::load(data,i);

  data->getValue(BRANCH_R,&R,i);
  data->getValue(BRANCH_X,&X,i);
  data->getValue(BRANCH_B,&Bc,i);

  if(abs(R) > 1e-6) {
    p_hasResistance = true;
  }

  if(abs(X) > 1e-6) {
    p_hasInductance = true;
  }

  R1 = R; L1 = X/OMEGA_S; C1 = Bc/OMEGA_S;
  R0 = 3*R1; L0 = 3*L1; C0 = 3*C1;

  double Rs = (2*R1 + R0)/3.0;
  double Rm = (R0 - R1)/3.0;
  p_R[0][0] = p_R[1][1] = p_R[2][2] = Rs;
  p_R[0][1] = p_R[1][0] = Rm;
  p_R[0][2] = p_R[2][0] = Rm;
  p_R[1][2] = p_R[2][1] = Rm;
  
  double Ls = (2*L1 + L0)/3.0;
  double Lm = (L0 - L1)/3.0;
  p_L[0][0] = p_L[1][1] = p_L[2][2] = Ls;
  p_L[0][1] = p_L[1][0] = Lm;
  p_L[0][2] = p_L[2][0] = Lm;
  p_L[1][2] = p_L[2][1] = Lm;
  
  double Cp = (2*C1 + C0)/3.0;
  double Cg = (C0 - C1)/3.0;
  p_C[0][0] = p_C[1][1] = p_C[2][2] = Cp;
  p_C[0][1] = p_C[1][0] = Cg;
  p_C[0][2] = p_C[2][0] = Cg;
  p_C[1][2] = p_C[2][1] = Cg;
}

/**
 * Initialize branch model before calculation
 * @param [output] values - array where initialized branch variables should be set
 */
void Lumpedline::init(gridpack::RealType *values)
{
  double *x = values + offsetb;
  
  fbus->addLumpedLineCshunt(p_C,0.5);
  tbus->addLumpedLineCshunt(p_C,0.5);

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

  gridpack::ComplexType Iline = (Vf - Vt)/Zline;

  double Ilinem,Ilinea;
  
  Ilinem = abs(Iline);
  Ilinea = atan2(imag(Iline),real(Iline));

  x[0] = Ilinem*sin(Ilinea);
  x[1] = Ilinem*sin(Ilinea - 2.0*PI/3.0);
  x[2] = Ilinem*sin(Ilinea + 2.0*PI/3.0);

}

/**
 * Write output from branchs to a string.
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool Lumpedline::serialWrite(char *string, const int bufsize,
			       const char *signal)
{
  return false;
}

/**
 * Write out branch state
 * @param signal character string used to determine behavior
 * @param string buffer that contains output
 */
void Lumpedline::write(const char* signal, char* string)
{
}

/**
 * Return the branch current injection 
 * @param [output] ia - phase a current
 * @param [output] ib - phase b current
 * @param [output] ic - phase c current
 */
void Lumpedline::getCurrent(double *ia, double *ib, double *ic)
{
  *ia = ibr[0];
  *ib = ibr[1];
  *ic = ibr[2];
}

/**
 * Return the location for the current in the local branch array
 * @param [output] i_loc - location for the first current variable in the local branch array
 */
void Lumpedline::getCurrentLocalLocation(int *i_loc)
{
  *i_loc = offsetb;
}

/**
 * Return the global location for the branch current injection 
 * @param [output] i_gloc - global location for the first current variable
 */
void Lumpedline::getCurrentGlobalLocation(int *i_gloc)
{
  *i_gloc = p_gloc;
}

/**
 * Get number of matrix values contributed by branch
 * @return number of matrix values
 */
int Lumpedline::matrixNumValues()
{
  int numvals = 0;
  if(p_hasInductance) {
    // fval = (vf - R*i - vt) - L*di_dt
    // dfval_dvf = I => 3 entries
    // dfval_di  = -R - sL => 9 entries
    // dfval_dvt = -I => 3 entries
    // Total 15 entries
    numvals += 15;
  } else {
    // fval = vf - R*i - vt
    // dfval_dvf = I => 3 entries
    // dfval_di = -R => 9 entries
    // dfval_dvt = -I => 3 entries
    // Total 15 entries
    numvals += 15;
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
void Lumpedline::matrixGetValues(int *nvals, gridpack::RealType *values, int *rows, int *cols)
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
    values[ctr] = 1.0;

    rows[ctr+1] = i_gloc + j;
    cols[ctr+1] = vt_gloc + j;
    values[ctr+1] = -1.0;

    ctr += 2;
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
void Lumpedline::vectorGetValues(gridpack::RealType *values)
{
  double *f = values + offsetb;

  if(p_mode == RESIDUAL_EVAL) {
    double vf[3],vt[3];

    fbus->getVoltages(&vf[0],&vf[1],&vf[2]);
    tbus->getVoltages(&vt[0],&vt[1],&vt[2]);

    double vf_minus_vt[3];
    
    vf_minus_vt[0] = vf[0] - vt[0];
    vf_minus_vt[1] = vf[1] - vt[1];
    vf_minus_vt[2] = vf[2] - vt[2];
    
    if(p_hasInductance) {
      double Ribr[3],Ldidt[3];
      // vf - R*ibr - vt - L*didt = 0
      matvecmult3x3(p_R, ibr,Ribr); // fval1 = R*ibr
      matvecmult3x3(p_L,dibr_dt,Ldidt); // fval2 = L*didt
      f[0] = vf_minus_vt[0] - Ribr[0] - Ldidt[0];
      f[1] = vf_minus_vt[1] - Ribr[1] - Ldidt[1];
      f[2] = vf_minus_vt[2] - Ribr[2] - Ldidt[2];
    } else {
      double Ribr[3];
      
      matvecmult3x3(p_R,ibr,Ribr);
      
      f[0] = -Ribr[0] + vf_minus_vt[0];
      f[1] = -Ribr[1] + vf_minus_vt[1];
      f[2] = -Ribr[2] + vf_minus_vt[2];
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
void Lumpedline::setValues(gridpack::RealType *values)
{
  gridpack::RealType *x = values + offsetb;

  if(p_mode == XVECTOBUS) {
    ibr[0] = x[0];
    ibr[1] = x[1];
    ibr[2] = x[2];
  } else if(p_mode == XDOTVECTOBUS) {
    dibr_dt[0] = x[0];
    dibr_dt[1] = x[1];
    dibr_dt[2] = x[2];
  }
}

/**
   Prestep function
*/
void Lumpedline::preStep(double time, double timestep)
{
}

/**
   Poststep function
*/
void Lumpedline::postStep(double time)
{
}
