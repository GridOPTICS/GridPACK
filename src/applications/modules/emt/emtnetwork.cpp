/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   emtnetwork.cpp
 * @brief  EMT network implementation
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>
#include <emtnetwork.hpp>
#include <gridpack/include/gridpack.hpp>
#include <gridpack/utilities/complex.hpp>
#include <constants.hpp>
#include <gridpack/math/dae_solver.hpp>
#include <model_classes/gencls.hpp>
#include <model_classes/constantimpedance.hpp>
//#include <model_classes/lumpedline.hpp>


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
  p_nvar = 0;
  p_nvarbus = 0;
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

  if(p_nload) {
    for(int i=0; i < p_nload; i++) {
      if(p_load[i]) delete(p_load[i]);
    }
    free(p_neqsload);
  }
  if(p_nvar) {
    delete(p_vecidx);
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
  Set the shift value provided by TS
*/
void EmtBus::setTime(double time)
{
  p_time = time;
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
    if(!p_gen[i]->getStatus()) continue;

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
    if(!p_gen[i]->getStatus()) continue;

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
    if(!p_gen[i]->getStatus()) continue;

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
   Sets up the bus component

   - Set up arrays and calculates number of variables
**/
void EmtBus::setup()
{
  int i,j;

  p_nvar = 0; // Initialize number of variables
  p_nvarbus = 3;
  /* Check if there are any capactive shunts */
  for(i=0; i < 3; i++) {
    if(p_hasCapacitiveShunt) break;
    for(j=0; j < 3; j++) {
      if(abs(p_Cshunt[i][j]) > 1e-6) {
	p_hasCapacitiveShunt = true;
	break;
      }
    }
  }

  if(p_hasCapacitiveShunt) {
    // Invert the capacitive shunt
    inverse3x3(p_Cshunt,p_Cshuntinv);
  }
  
  for(i=0; i < 3; i++) {
    if(p_hasInductiveShunt) break;
    for(j=0; j < 3; j++) {
      if(abs(p_Lshunt[i][j]) > 1e-6) {
	p_nvarbus += 3; // Bus also has inductive shunt
	p_hasInductiveShunt = true;
	break;
      }
    }
  }

  p_nvar += p_nvarbus;

  // Set up generators
  if(p_ngen) {
    p_neqsgen = (int*)malloc(p_ngen*sizeof(int));
    p_neqsexc= (int*)malloc(p_ngen*sizeof(int));
    p_neqsgov= (int*)malloc(p_ngen*sizeof(int));
  }

  for(i=0; i < p_ngen; i++) {
    if(!p_gen[i]->getStatus()) {
      p_neqsgen[i] = 0;
      p_neqsexc[i] = 0;
      p_neqsgov[i] = 0;
    } else {

      // Set number of equations for this generator
      p_gen[i]->getnvar(&p_neqsgen[i]);
      // Set the offset for the first variable in the bus variable array 
      p_gen[i]->setBusOffset(p_nvar);

      bool has_ex = p_gen[i]->hasExciter();
      bool has_gv = p_gen[i]->hasGovernor();

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

  if(p_nload) {
    p_neqsload = (int*)malloc(p_nload*sizeof(int));
  }

  for(i=0; i < p_nload; i++) {
    if(!p_load[i]->getStatus()) {
      p_neqsload[i] = 0;
    } else {
      // Get number of equations for this load
      p_load[i]->getnvar(&p_neqsload[i]);
      // Set the offset for the first variable in the bus variable array 
      p_load[i]->setBusOffset(p_nvar);
      
      p_nvar += p_neqsload[i];
    }
  }
  if(p_nvar) {
    p_vecidx = new int[p_nvar];
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
  int    gstatus,lstatus;
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
  double va = Vm*sin(Va);
  double vb = Vm*sin(Va - 2*PI/3.0);
  double vc = Vm*sin(Va + 2*PI/3.0);

  if(p_vptr) {
    *p_vptr     = va;
    *(p_vptr+1) = vb;
    *(p_vptr+2) = vc;
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
    p_gen = (BaseEMTGenModel**)malloc(p_ngen*sizeof(BaseEMTGenModel*));

    for(i=0; i < p_ngen; i++) {
      data->getValue(GENERATOR_STAT,&gstatus,i); // Generator status
      if(!gstatus) {
        p_gen[i] = new BaseEMTGenModel;
        p_gen[i]->setStatus(0);
        continue;
      }

      p_gen[i] = NULL;
      std::string model;
      data->getValue(GENERATOR_MODEL,&model,i);

      std::string type = util.trimQuotes(model);
      util.toUpper(type);

      if(type == "GENCLS") {
        Gencls *clgen;
        clgen = new Gencls;
        p_gen[i] = clgen;
      }
      
      // Read generator data stored in data collection objects
      p_gen[i]->load(data,i); // load data
      has_ex = p_gen[i]->hasExciter();
      has_gv = p_gen[i]->hasGovernor();
      if (has_ex) p_gen[i]->getExciter()->load(data,i); // load exciter data
      if (has_gv) p_gen[i]->getGovernor()->load(data,i); // load governor model
    }
  }

  data->getValue(LOAD_NUMBER, &p_nload);
  if(p_nload) {
    p_load = (BaseEMTLoadModel**)malloc(p_nload*sizeof(BaseEMTLoadModel*));
  }
  
  for(i=0; i < p_nload; i++) {
    data->getValue(LOAD_STATUS,&lstatus,i); // Generator status
    if(!lstatus) {
      p_load[i] = new BaseEMTLoadModel;
      p_load[i]->setStatus(0);
      continue;
    }

    std::string lmodel;
    data->getValue(LOAD_MODEL,&lmodel,i);

    if(lmodel.empty()) {
      // No load model given, so use a R-L impedance model
      Constantimpedance *ciload;
      ciload = new Constantimpedance;
      p_load[i] = ciload;
    } else {
      std::string type = util.trimQuotes(lmodel);
      util.toUpper(type);
    }

    p_load[i]->load(data,i);

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

void EmtBus::matrixGetValues(gridpack::math::Matrix &matrix)
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
  p_vecidx[ielem] = idx;
}

/**
 * Return a set of element indices that map the local indices to
 * global indices
 * @param idx array of global indices
 */
void EmtBus::vectorGetElementIndices(int *idx)
{
  int i;
  for (i=0; i< p_nvar; i++) {
    idx[i] = p_vecidx[i];
  }
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
  bool has_ex=false, has_gv=false;
  int i;
  gridpack::ComplexType *x,*f;
  for(i=0; i < p_nvar; i++) {
    idx[i] = p_vecidx[i];
  }
  if(p_mode == INIT_X) { /* Initialization of values */
    int i;
    double va = p_Vm0*sin(p_Va0);
    double vb = p_Vm0*sin(p_Va0 - 2*PI/3.0);
    double vc = p_Vm0*sin(p_Va0 + 2*PI/3.0);

    x = values;
    
    if(p_isolated) {
      x[0] = va;
      x[1] = vb;
      x[2] = vc;

      return;
    }

    x[0] = va;
    x[1] = vb;
    x[2] = vc;

    if(p_hasInductiveShunt) {
      /* Calculate the current in the inductive shunt */
      double VR,VI;
      VR = p_Vm0*cos(p_Va0);
      VI = p_Vm0*sin(p_Va0);
      gridpack::ComplexType V = gridpack::ComplexType(VR,VI);
      gridpack::ComplexType S = gridpack::ComplexType(p_gl,p_bl);
      gridpack::ComplexType Ishunt = S/V;
      double Ishuntmag,Ishuntang;
      Ishuntmag = abs(Ishunt);
      Ishuntang = atan2(imag(Ishunt),real(Ishunt));

      x[3] = Ishuntmag*sin(Ishuntang);
      x[4] = Ishuntmag*sin(Ishuntang-2*PI/3.0);
      x[5] = Ishuntmag*sin(Ishuntang+2*PI/3.0);
    }

    for(i=0; i < p_ngen; i++) {
      if(!p_gen[i]->getStatus()) continue;

      p_gen[i]->setVoltage(va,vb,vc);
      p_gen[i]->setInitialVoltage(p_Vm0,p_Va0);
      p_gen[i]->init(x);
    }

    for(i=0; i < p_nload; i++) {
      if(!p_load[i]->getStatus()) continue;

      p_load[i]->setVoltage(va,vb,vc);
      p_load[i]->setInitialVoltage(p_Vm0,p_Va0);
      p_load[i]->init(x);
    }
  } else if(p_mode == RESIDUAL_EVAL) {
    int i,j;
    double *v = p_vptr;
    double i_gen[3]; // Total generator current
    double i_geni[3]; // Current from generator i
    double i_load[3]; // Total load current
    double i_loadi[3]; // Current from load i
    double i_br[3];   // Total current from branches
    double i_bri[3]; // Current from branch i
    double i_mis[3]; // i_gen - i_br - i_load

    f = values;
    if(p_isolated) {
      f[0] = v[0] - p_Vm0*sin(OMEGA_S*p_time + p_Va0);
      f[1] = v[1] - p_Vm0*sin(OMEGA_S*p_time + p_Va0 - TWOPI_OVER_THREE);
      f[2] = v[2] - p_Vm0*sin(OMEGA_S*p_time + p_Va0 + TWOPI_OVER_THREE);

      return;
    }

    if(p_hasInductiveShunt) {
    }

    std::vector<boost::shared_ptr<BaseComponent> > branches;
    //Get the edges connected to this bus
    gridpack::component::BaseBusComponent::getNeighborBranches(branches);
    int nconnbranch = branches.size();

    i_br[0] = i_br[1] = i_br[2] = 0.0;
    for(i=0; i < nconnbranch; i++) {
      EmtBranch *branch = dynamic_cast<EmtBranch*>(branches[i].get());
      
      int nparlines = branch->getNumParallelLines();
      
      for(j=0; j < nparlines; j++) {
	if(!branch->getStatus(j)) continue;

	i_bri[0] = i_bri[1] = i_bri[2] = 0.0;
	branch->getCurrent(j,&i_bri[0],&i_bri[1],&i_bri[2]);

	i_br[0] += i_bri[0];
	i_br[1] += i_bri[1];
	i_br[2] += i_bri[2];

      }
    }
    
    i_gen[0] = i_gen[1] = i_gen[2] = 0.0;
    for(i=0; i < p_ngen; i++) {
      if(!p_gen[i]->getStatus()) continue;

      p_gen[i]->setVoltage(v[0],v[1],v[2]);
      p_gen[i]->setTime(p_time);

      /* Get generator current */
      i_geni[0] = i_geni[1] = i_geni[2] = 0.0;
      p_gen[i]->getCurrent(&i_geni[0],&i_geni[1],&i_geni[2]);

      i_gen[0] += i_geni[0];
      i_gen[1] += i_geni[1];
      i_gen[2] += i_geni[2];

      /* Generator residual evaluation */
      p_gen[i]->setMode(p_mode);
      p_gen[i]->vectorGetValues(f);
    }

    i_load[0] = i_load[1] = i_load[2] = 0.0;
    for(i=0; i < p_nload; i++) {
      if(!p_load[i]->getStatus()) continue;

      p_load[i]->setVoltage(v[0],v[1],v[2]);
      p_load[i]->setTime(p_time);

      /* Get load current */
      i_loadi[0] = i_loadi[1] = i_loadi[2] = 0.0;
      p_load[i]->getCurrent(&i_loadi[0],&i_loadi[1],&i_loadi[2]);

      i_load[0] += i_loadi[0];
      i_load[1] += i_loadi[1];
      i_load[2] += i_loadi[2];

      /* Load residual evaluation */
      p_load[i]->setMode(p_mode);
      p_load[i]->vectorGetValues(f);
    }

    i_mis[0] = i_gen[0] - i_br[0] - i_load[0];
    i_mis[1] = i_gen[1] - i_br[1] - i_load[1];
    i_mis[2] = i_gen[2] - i_br[2] - i_load[2];

    if(p_hasCapacitiveShunt) {
      double fval[3];
      matvecmult3x3(p_Cshuntinv,i_mis,fval);
      f[0] = fval[0] - p_dvdt[0];
      f[1] = fval[1] - p_dvdt[1];
      f[2] = fval[2] - p_dvdt[2];
    } else {
      f[0] = i_mis[0];
      f[1] = i_mis[1];
      f[2] = i_mis[2];
    }
    
  } 
}

/**
 * Set network elements based on values in vector
 * @param array containing vector values
 */
void EmtBus::vectorSetElementValues(gridpack::ComplexType *values)
{
  int i;
  double va,vb,vc;
  
  if(p_mode == XVECTOBUS) {
    *p_vptr     = va = real(values[0]);
    *(p_vptr+1) = vb = real(values[1]);
    *(p_vptr+2) = vc = real(values[2]);

    if(p_isolated) return;

    for(i=0; i < p_ngen; i++) {
      if(!p_gen[i]->getStatus()) continue;

      p_gen[i]->setMode(p_mode);
      p_gen[i]->setVoltage(va,vb,vc);
      p_gen[i]->setValues(values);
    }

    for(i=0; i < p_nload; i++) {
      if(!p_load[i]->getStatus()) continue;

      p_load[i]->setMode(p_mode);
      p_load[i]->setVoltage(va,vb,vc);
      p_load[i]->setValues(values);
    }

    
  } else if(p_mode == XDOTVECTOBUS) {
    p_dvdt[0] = real(values[0]);
    p_dvdt[1] = real(values[1]);
    p_dvdt[2] = real(values[2]);

    if(p_isolated) return;

    for(i=0; i < p_ngen; i++) {
      if(!p_gen[i]->getStatus()) continue;
      
      p_gen[i]->setMode(p_mode);
      p_gen[i]->setValues(values);
    }
  }
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
  p_iptr = NULL;
  p_nvar = 3;
  p_num_vals = 0;
  p_mode = NONE;
}

/**
 *  Simple destructor
 */
EmtBranch::~EmtBranch(void)
{
  if(p_nvar) {
    delete [] p_ibr;
    delete [] p_didt;
    delete [] p_vecidx;
  }
}

void EmtBranch::setup()
{
  int i;
  p_nvar = 0;

  for(i = 0; i < p_nparlines; i++) {
    p_localoffset.push_back(p_nvar);
    if(!p_status[i]) continue;

    p_nvar += 3;
  }

  if(p_nvar) {
    p_vecidx = new int[p_nvar];
    p_ibr  = new double[3*p_nparlines];
    p_didt = new double[3*p_nparlines];
  }
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
  double pi=PI;
  int    status,i;
  std::string cktid;
  double R,X, Bc;

  data->getValue(BRANCH_NUM_ELEMENTS,&p_nparlines);

  try {
    if(p_nparlines > 1) {
      throw(p_nparlines);
    }
  } catch(int nlines) {
    std::cout << "Parallel lines not allowed\n";
  }

  for(i=0; i < p_nparlines; i++) {
    //    LumpedLineModel *lumpedline;
    //    lumpedpi = new LumpedLineModel;
    //    p_impl[i] = lumpedline;

    // Get line parameters
    data->getValue(BRANCH_STATUS,&status,i);
    data->getValue(BRANCH_CKT,&cktid,i);
    // Positive sequence values
    data->getValue(BRANCH_R,&R,i);
    data->getValue(BRANCH_X,&X,i);
    data->getValue(BRANCH_B,&Bc,i);

    p_status.push_back(status);
    p_cktid.push_back(cktid);
    p_lineR.push_back(R);
    p_lineX.push_back(X);
    
    // Assuming zero sequence values as 3 times of +ve sequence
    // we get the following phase matrices
    double R1,L1,C1;
    double R0,L0,C0;

    if(abs(R) > 1e-6) {
      p_hasResistance = true;
    }

    if(abs(X) > 1e-6) {
      p_hasInductance = true;
    }
      
    R1 = R; L1 = X/OMEGA_S; C1 = Bc/OMEGA_S;
    R0 = 3*R1; L0 = 3*L1; C0 = 3*C1;

    p_R[0][0] = p_R[1][1] = p_R[2][2] = R1 + R0;
    p_R[0][1] = p_R[1][0] = R1;
    p_R[0][2] = p_R[2][0] = R1;
    p_R[1][2] = p_R[2][1] = R1;

    p_L[0][0] = p_L[1][1] = p_L[2][2] = L1 + L0;
    p_L[0][1] = p_L[1][0] = L1;
    p_L[0][2] = p_L[2][0] = L1;
    p_L[1][2] = p_L[2][1] = L1;

    if(p_hasInductance) {
      // L^-1
      inverse3x3(p_L,p_Linv);

      // -L^-1*R
      scaledmatmatmult3x3(p_Linv,p_R,p_minusLinvR,-1.0);
    }
    
    p_C[0][0] = p_C[1][1] = p_C[2][2] = C1 + C0;
    p_C[0][1] = p_C[1][0] = C1;
    p_C[0][2] = p_C[2][0] = C1;
    p_C[1][2] = p_C[2][1] = C1;

    /* Add Capacitive shunt to from and to buses */
    EmtBus *busf = dynamic_cast<EmtBus*>((getBus1()).get());
    EmtBus *bust = dynamic_cast<EmtBus*>((getBus2()).get());
    busf->addLumpedLineCshunt(p_C,0.5);
    bust->addLumpedLineCshunt(p_C,0.5);

  }
  
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
  Set the shift value provided by TS
*/
void EmtBranch::setTSshift(double shift)
{
  p_TSshift = shift;
}

/**
  Set the shift value provided by TS
*/
void EmtBranch::setTime(double time)
{
  p_time = time;
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

void EmtBranch::matrixGetValues(gridpack::math::Matrix &matrix)
{

}

int EmtBranch::getXCBufSize(void)
{
  return p_nvar*sizeof(double);
}

void EmtBranch::setXCBuf(void *buf)
{
  p_iptr = static_cast<double*>(buf);

}

/**
 * Set value of global index for corresponding local index
 * @param ielem local index for element
 * @param idx global index of element
 */
void EmtBranch::vectorSetElementIndex(int ielem, int idx)
{
  p_vecidx[ielem] = idx;
}

/**
 * Return a set of element indices that map the local indices to
 * global indices
 * @param idx array of global indices
 */
void EmtBranch::vectorGetElementIndices(int *idx)
{
  int i;
  for (i=0; i<p_nvar; i++) {
    idx[i] = p_vecidx[i];
    //    printf("Branch %d -- %d: p_vecidx[%d] = %d\n",this->getBus1OriginalIndex(),this->getBus2OriginalIndex(),i,p_vecidx[i]);
  }
}

/**
 * Return number elements contributed by this bus
 * @return number of elements
 */
int EmtBranch::vectorNumElements() const
{
  return p_nvar;
}

/**
 * Return the elements and their global indices in the vector
 * @param values array of element values
 * @param idx array of element indices
 */
void EmtBranch::vectorGetElementValues(gridpack::ComplexType *values, int *idx)
{
  int i;
  for(i=0; i < p_nvar; i++) {
    idx[i] = p_vecidx[i];
  }
  if(p_mode == INIT_X) { /* Initialization of values */
    gridpack::ComplexType *x = values;
    for(i=0; i < p_nparlines; i++) {
      double *ibr = p_ibr + 3*i;
      ibr[0] = ibr[1] = ibr[2] = 0.0;
      if(!p_status[i]) continue;
      
      EmtBus *busf = dynamic_cast<EmtBus*>((getBus1()).get());
      EmtBus *bust = dynamic_cast<EmtBus*>((getBus2()).get());

      double Vmf,Vaf,Vmt,Vat;
      double VDf,VQf,VDt,VQt;
      busf->getInitialVoltage(&Vmf,&Vaf);
      bust->getInitialVoltage(&Vmt,&Vat);

      VDf = Vmf*cos(Vaf);
      VQf = Vmf*sin(Vaf);
      VDt = Vmt*cos(Vat);
      VQt = Vmt*sin(Vat);

      gridpack::ComplexType Vf = gridpack::ComplexType(VDf,VQf);
      gridpack::ComplexType Vt = gridpack::ComplexType(VDt,VQt);

      gridpack::ComplexType Zline = gridpack::ComplexType(p_lineR[i],p_lineX[i]);

      gridpack::ComplexType Iline = (Vf - Vt)/Zline;

      double Ilinem,Ilinea;

      Ilinem = abs(Iline);
      Ilinea = atan2(imag(Iline),real(Iline));

      x += p_localoffset[i];

      x[0] = ibr[0] = Ilinem*sin(Ilinea);
      x[1] = ibr[1] = Ilinem*sin(Ilinea - 2.0*PI/3.0);
      x[2] = ibr[2] = Ilinem*sin(Ilinea + 2.0*PI/3.0);
    }
  } else if(p_mode == RESIDUAL_EVAL) {
    gridpack::ComplexType *f = values;
    for(i=0; i < p_nparlines; i++) {
      double *ibr = p_ibr + 3*i;
      if(!p_status[i]) continue;
      
      EmtBus *busf = dynamic_cast<EmtBus*>((getBus1()).get());
      EmtBus *bust = dynamic_cast<EmtBus*>((getBus2()).get());

      double vf[3],vt[3];

      busf->getVoltages(&vf[0],&vf[1],&vf[2]);
      bust->getVoltages(&vt[0],&vt[1],&vt[2]);

      double vf_minus_vt[3];
      
      vf_minus_vt[0] = vf[0] - vt[0];
      vf_minus_vt[1] = vf[1] - vt[1];
      vf_minus_vt[2] = vf[2] - vt[2];
      
      if(p_hasInductance) {
	double fval1[3],fval2[3];

	matvecmult3x3(p_Linv,vf_minus_vt,fval1); // fval1 = L^-1*(vf - vt)
	matvecmult3x3(p_minusLinvR,ibr,fval2); // fval2 = -L^-1*R*i_br
	f[0] = fval1[0] + fval2[0] - p_didt[0];
	f[1] = fval1[1] + fval2[1] - p_didt[1];
	f[2] = fval1[2] + fval2[2] - p_didt[2];
      } else {
	double fval[3];

	scaledmatvecmult3x3(p_R,ibr,fval,-1.0);

	f[0] = fval[0] + vf_minus_vt[0];
	f[1] = fval[1] + vf_minus_vt[1];
	f[2] = fval[2] + vf_minus_vt[2];
      }
      
      f += p_localoffset[i];
    }
  }    
}

/**
 * Set network elements based on values in vector
 * @param array containing vector values
 */
void EmtBranch::vectorSetElementValues(gridpack::ComplexType *values)
{
  int i,l=0;
  double *ibr    = p_ibr;
  double *dibrdt = p_didt;
  gridpack::ComplexType *x = values;

  if(p_mode == XVECTOBUS) {
    for(i=0; i < p_nparlines; i++) {
      if(p_status[i]) {
	*p_iptr     = ibr[0] = real(x[0]);
	*(p_iptr+1) = ibr[1] = real(x[1]);
	*(p_iptr+2) = ibr[2] = real(x[2]);

	p_iptr += 3;
	x += 3;
      }
      ibr += 3;
    }
  } else if(p_mode == XDOTVECTOBUS) {
    for(i=0; i < p_nparlines; i++) {
      if(p_status[i]) {
	dibrdt[0] = real(x[0]);
	dibrdt[1] = real(x[1]);
	dibrdt[2] = real(x[2]);

	x += 3;
      }
      dibrdt += 3;
    }
  }
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
