/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   emtnetwork.cpp
 * @brief  EMT network implementation
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>
#include <emtnetwork.hpp>
#include <gridpack/include/gridpack.hpp>
#include <gridpack/utilities/complex.hpp>
#include <constants.hpp>
#include <gridpack/math/dae_solver.hpp>

/**
 *  Simple constructor
 */
EmtBus::EmtBus(void)
{
  p_gl = p_bl = 0.0;
  p_pl = p_ql = 0.0;
  p_ngen = p_nactivegen = 0;
  p_isolated = false;
  p_mode = NONE;
  p_vptr = NULL;
  p_nvar = 3;
  p_neqsgen = NULL;
  p_neqsexc= NULL;
  p_neqsgov= NULL;
  p_gen     = NULL;
  p_num_vals = -1;
  p_hasCapacitiveShunt = false;
  p_hasResistiveShunt = false;
  p_hasInductiveShunt = false;
  p_nphases = 3;

  p_Gshunt[0][0] = p_Gshunt[0][1] = p_Gshunt[0][2] = 0.0;
  p_Gshunt[1][0] = p_Gshunt[1][1] = p_Gshunt[1][2] = 0.0;
  p_Gshunt[2][0] = p_Gshunt[2][1] = p_Gshunt[2][2] = 0.0;

  p_Lshunt[0][0] = p_Lshunt[0][1] = p_Lshunt[0][2] = 0.0;
  p_Lshunt[1][0] = p_Lshunt[1][1] = p_Lshunt[1][2] = 0.0;
  p_Lshunt[2][0] = p_Lshunt[2][1] = p_Lshunt[2][2] = 0.0;

  p_Cshunt[0][0] = p_Cshunt[0][1] = p_Cshunt[0][2] = 0.0;
  p_Cshunt[1][0] = p_Cshunt[1][1] = p_Cshunt[1][2] = 0.0;
  p_Cshunt[2][0] = p_Cshunt[2][1] = p_Cshunt[2][2] = 0.0;


}

/**
 *  Simple destructor
 */
EmtBus::~EmtBus(void)
{
  if(p_ngen) {
    for(int i=0; i < p_ngen; i++) {
      if(p_gen[i]) delete(p_gen[i]);
    }
    free(p_neqsgen);
    free(p_neqsexc);
    free(p_neqsgov);
    free(p_gen);
  }
}

/**
  Set the shift value provided by TS
*/
void EmtBus::setTSshift(double shift)
{
  p_TSshift = shift;
}

/**
  *  Check if the bus is isolated. Returns true if the bus is isolated

*/
bool EmtBus::isIsolated(void) const
{
  bool isol = false;
  if(p_isolated) isol = true;
  return isol;
}

/**
 * Set local offset for the bus (starting location of variables for this bus in the state vector)
   and propogate that to the models on this bus. Used for events 
*/
void EmtBus::setLocalOffset(int offset)
{
  int i;

  p_offset = offset;

  for(i=0; i < p_ngen; i++) {
    if(!p_gen[i]->getGenStatus()) continue;

    p_gen[i]->setBusLocalOffset(offset);
      
    if(p_gen[i]->hasExciter()) {
      p_gen[i]->getExciter()->setBusLocalOffset(offset);
    }
    
    if(p_gen[i]->hasGovernor()) {
      p_gen[i]->getGovernor()->setBusLocalOffset(offset);
    }
  }
}

/**
 * Reset limiter flags after event has occured. Only called when the network
 * is resolved
*/
void EmtBus::resetEventFlags()
{
  int i;

  for(i=0; i < p_ngen; i++) {
    if(!p_gen[i]->getGenStatus()) continue;

    p_gen[i]->resetEventFlags();
    
    if(p_gen[i]->hasExciter()) {
      p_gen[i]->getExciter()->resetEventFlags();
    }
    
    if(p_gen[i]->hasGovernor()) {
      p_gen[i]->getGovernor()->resetEventFlags();
    }
  }
}

/**
 * Get voltages va, vb, vc from the exchange buffer
 * @param double va - phase a voltage
 * @param double vb - phase b voltage
 * @param double vc - phase c voltage
 */
void EmtBus::getVoltages(double *va,double *vb,double *vc) const
{
  *va = *p_vptr;
  *vb = *(p_vptr+1);
  *vc = *(p_vptr+2);
}

/**
 * Add events
 * @eman - EventManager pointer
 */
void EmtBus::setEvent(gridpack::math::DAESolver::EventManagerPtr eman)
{
  int i;
  bool has_ex=false,has_gov=false;

  for(i=0; i < p_ngen; i++) {
    if(!p_gen[i]->getGenStatus()) continue;

    has_ex = p_gen[i]->hasExciter();
    if(has_ex) {
      p_gen[i]->getExciter()->setEvent(eman);
    }
    has_gov = static_cast<bool>(p_gen[i]->getGovernor());
    if(has_gov) {
      p_gen[i]->getGovernor()->setEvent(eman);
    }
  }
}

/**
 * Load values stored in DataCollection object into EmtBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 *
 * Note: It is assumed that the network data file read has data for "solved"
 * power flow so we can do the initialization when loading the data
 */
void EmtBus::load(const
         boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  double Vm, Va;  // Voltage magnitude and angle from the case file
  double pi=PI;
  int    i=0;
  double Pg, Qg,mbase,Rs,Xdp,H,D; // temp. variables to hold machine data
  int    gstatus;
  double delta,dw=0.0;
  double Pm,Ep;
  double IGD,IGQ; /* D and Q components of generator currents */

  p_sbase = DEFAULT_MVABASE; // Default MVA base

  // Get MVAbase
  data->getValue(CASE_SBASE,&p_sbase);

  // Get Voltage magnitude and angle
  data->getValue(BUS_VOLTAGE_ANG,&Va); // This is in degress
  data->getValue(BUS_VOLTAGE_MAG,&Vm);

  /* Convert from degrees to radians */
  Va *= pi/180.0;

  // Assume balanced conditions
  // Three phase voltages
  double p_va = Vm*sin(Va);
  double p_vb = Vm*sin(Va - 2*PI/3.0);
  double p_vc = Vm*sin(Va + 2*PI/3.0);

  if(p_vptr) {
    *p_vptr     = p_va;
    *(p_vptr+1) = p_vb;
    *(p_vptr+2) = p_vc;
  }

  /* Save the voltage magnitude and angle so that we can use it later for building loads */
  p_Vm0 = Vm;
  p_Va0 = Va;

  // Get the bus type
  data->getValue(BUS_TYPE,&p_bustype);
  if(p_bustype == 4) p_isolated = true;

  if(p_isolated) return;

  // Read shunts
  data->getValue(BUS_SHUNT_GL,&p_gl);
  data->getValue(BUS_SHUNT_BL,&p_bl);
  p_gl /= p_sbase;
  p_bl /= p_sbase;

  if(abs(p_gl) > 1e-6) {
    p_hasResistiveShunt = true;
    p_Gshunt[0][0] = p_gl;
    p_Gshunt[1][1] = p_gl;
    p_Gshunt[2][2] = p_gl;
  }

  if(abs(p_bl) > 1e-6) {
    if(p_bl > 0) {
      p_hasCapacitiveShunt = true;
      p_Cshunt[0][0] = p_bl/OMEGA_S;
      p_Cshunt[1][1] = p_bl/OMEGA_S;
      p_Cshunt[2][2] = p_bl/OMEGA_S;
    } else {
      p_hasInductiveShunt = true;
      p_Lshunt[0][0] = -p_bl/OMEGA_S;
      p_Lshunt[1][1] = -p_bl/OMEGA_S;
      p_Lshunt[2][2] = -p_bl/OMEGA_S;
    }
  }
      
  gridpack::utility::StringUtils util;
  bool has_ex;
  bool has_gv;

  // Read Generators 
  // Get number of generators incident on this bus
  data->getValue(GENERATOR_NUMBER, &p_ngen);
  if(p_ngen) {
    p_gen = (BaseGenModel**)malloc(p_ngen*sizeof(BaseGenModel*));
    p_neqsgen = (int*)malloc(p_ngen*sizeof(int));
    p_neqsexc= (int*)malloc(p_ngen*sizeof(int));
    p_neqsgov= (int*)malloc(p_ngen*sizeof(int));

    for(i=0; i < p_ngen; i++) {
      data->getValue(GENERATOR_STAT,&gstatus,i); // Generator status
      if(!gstatus) {
        p_gen[i] = new BaseGenModel;
        p_gen[i]->setGenStatus(0);
        continue;
      }

      p_gen[i] = NULL;
      p_neqsgen[i] = 0;
      p_neqsexc[i] = 0;
      p_neqsgov[i] = 0;
      std::string model;
      data->getValue(GENERATOR_MODEL,&model,i);

      std::string type = util.trimQuotes(model);
      util.toUpper(type);
 

      // Read generator data stored in data collection objects
      p_gen[i]->load(data,i); // load data
      has_ex = p_gen[i]->hasExciter();
      has_gv = p_gen[i]->hasGovernor();
      if (has_ex) p_gen[i]->getExciter()->load(data,i); // load exciter data
      if (has_gv) p_gen[i]->getGovernor()->load(data,i); // load governor model

      // Set number of equations for this generator
      p_gen[i]->vectorSize(&p_neqsgen[i]);
      // Set the offset for the first variable in the bus variable array 
      p_gen[i]->setBusOffset(p_nvar);
      if (has_ex) {
        p_gen[i]->getExciter()->vectorSize(&p_neqsexc[i]);
        p_gen[i]->getExciter()->setBusOffset(p_nvar+p_neqsgen[i]);
      }
      if (has_gv) {
        p_gen[i]->getGovernor()->vectorSize(&p_neqsgov[i]);
        p_gen[i]->getGovernor()->setBusOffset(p_nvar+p_neqsgen[i]+p_neqsexc[i]);
      }

      /* Update number of variables for this bus */
      p_nvar += p_neqsgen[i]+p_neqsexc[i]+p_neqsgov[i];
    }
  }

  data->getValue(LOAD_NUMBER, &p_nload);
  if(p_nload > 1) printf("Warning: More than one load detected on bus\n");
  
  for(i=0; i < p_nload; i++) {
    // Get active and reactive power load
    data->getValue(LOAD_PL,&p_pl,i);
    data->getValue(LOAD_QL,&p_ql,i);
    p_pl /= p_sbase;
    p_ql /= p_sbase;
  }

}

/**
 * Return number of rows (dependent variables) that bus contributes
 * to matrix
 * @return number of dependent variables (equations)
 */
int EmtBus::matrixNumRows()
{
  return p_nvar;
}

/**
 * Return number of columns (independent variables) that bus contributes
 * to matrix
 * @return number of independent variables
 */
int EmtBus::matrixNumCols()
{
  return p_nvar;
}

/** 
 * Set global row index
 * @param irow local row index
 * @param idx global row index
 */
void EmtBus::matrixSetRowIndex(int irow, int idx)
{
  if (p_rowidx.size() == 0) {
    p_rowidx.resize(p_nvar);
    int i;
    for (i=0; i<p_nvar; i++) p_rowidx[i] = -1;
  }
  p_rowidx[irow] = idx;
}

/** 
 * Set global column index
 * @param icol local column index
 * @param idx global column index
 */
void EmtBus::matrixSetColIndex(int icol, int idx)
{
  if (p_colidx.size() == 0) {
    p_colidx.resize(p_nvar);
    int i;
    for (i=0; i<p_nvar; i++) p_colidx[i] = -1;
  }
  p_colidx[icol] = idx;
}

/**
 * Return global row index given local row index
 * @param irow local row index
 * @return global row index
 */
int EmtBus::matrixGetRowIndex(int irow)
{
  return p_rowidx[irow];
}

/**
 * Return global column index given local column index
 * @param icol local column index
 * @return global column index
 */
int EmtBus::matrixGetColIndex(int icol)
{
  return p_colidx[icol];
}

/**
 * Total number of matrix elements returned by bus
 * @return number of matrix elements
 */
int EmtBus::matrixNumValues()
{
  if (p_isolated) {
    p_num_vals = 3;
  } else {
    
  }

  return p_num_vals;
}

/**
 * Return values from a matrix block
 * @param nvals: number of values to be inserted
 * @param values: pointer to matrix block values
 * @param rows: pointer to matrix block rows
 * @param cols: pointer to matrix block cols
 */
void EmtBus::matrixGetValues(int *nvals, gridpack::ComplexType *values,
    int *rows, int *cols)
{

}

/**
 * Set the model to control what matrices and vectors and built when using the
 * mapper
 * @param mode: enumerated constant for different modes
 */
void EmtBus::setMode(int mode)
{
  p_mode = (EMTMode)mode;
}

int EmtBus::getXCBufSize(void)
{
  return 3*sizeof(double);
}

void EmtBus::setXCBuf(void *buf)
{
  p_vptr = static_cast<double*>(buf);

}

/**
 * Set value of global index for corresponding local index
 * @param ielem local index for element
 * @param idx global index of element
 */
void EmtBus::vectorSetElementIndex(int ielem, int idx)
{
}

/**
 * Return a set of element indices that map the local indices to
 * global indices
 * @param idx array of global indices
 */
void EmtBus::vectorGetElementIndices(int *idx)
{

}

/**
 * Return number elements contributed by this bus
 * @return number of elements
 */
int EmtBus::vectorNumElements() const
{
  return p_nvar;
}

/**
 * Return the elements and their global indices in the vector
 * @param values array of element values
 * @param idx array of element indices
 */
void EmtBus::vectorGetElementValues(gridpack::ComplexType *values, int *idx)
{
}

/**
 * Set network elements based on values in vector
 * @param array containing vector values
 */
void EmtBus::vectorSetElementValues(gridpack::ComplexType *values)
{
}

/**
 * Write output from buses to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool EmtBus::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  return false;
}

/**
 *  Simple constructor
 */
EmtBranch::EmtBranch(void)
{
  p_nparlines = 0;
  p_nvar = 3;
  p_num_vals = 0;
  p_mode = NONE;
}

/**
 *  Simple destructor
 */
EmtBranch::~EmtBranch(void)
{
}

/**
 * Load values stored in DataCollection object into EmtBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 *
 * Note: GridPACK stores the data for all parallel branches in the same
 *       DataCollection object. 
 */
void EmtBranch::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{

}

/**
 * Set the model to control what matrices and vectors and built when using the
 * mapper
 * @param mode: enumerated constant for different modes
 */
void EmtBranch::setMode(int mode)
{
  p_mode = mode;
}

/**
 * Return number of rows (dependent variables) that branch contributes
 * to matrix
 * @return number of dependent variables (equations)
 */
int EmtBranch::matrixNumRows()
{
  return p_nvar;
}

/**
 * Return number of columns (independent variables) that branch contributes
 * to matrix
 * @return number of independent variables
 */
int EmtBranch::matrixNumCols()
{
  return p_nvar;
}

/** 
 * Set global row index
 * @param irow local row index
 * @param global row index
 */
void EmtBranch::matrixSetRowIndex(int irow, int idx)
{

}

/** 
 * Set global column index
 * @param icol local column index
 * @param global column index
 */
void EmtBranch::matrixSetColIndex(int icol, int idx)
{
  if (p_colidx.size() == 0) {
    p_colidx.resize(p_nvar);
    int i;
    for (i=0; i<p_nvar; i++) p_colidx[i] = -1;
  }
  p_colidx[icol] = idx;
}

/**
 * Return global row index given local row index
 * @param irow local row index
 * @return global row index
 */
int EmtBranch::matrixGetRowIndex(int irow)
{
  return p_rowidx[irow];
}

/**
 * Return global column index given local column index
 * @param icol local column index
 * @return global column index
 */
int EmtBranch::matrixGetColIndex(int icol)
{
  return p_colidx[icol];
}

/**
 * Total number of matrix elements returned by branch
 * @return number of matrix elements
 */
int EmtBranch::matrixNumValues()
{
  return p_num_vals;
}

/**
 * Return list of matrix values and their locations generated by the branch
 * @params nvals number of values inserted
 * @param values list of matrix values
 * @param rows list of row indices
 * @param cols list of column indices
 */
void EmtBranch::matrixGetValues(int *nvals,gridpack::ComplexType *values,
    int *rows, int *cols)
{

}

/**
 * Write output from branches to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if branch is contributing string to output, false otherwise
 */
bool EmtBranch::serialWrite(char *string, const int
    bufsize,  const char *signal)
{
  return false;
}
