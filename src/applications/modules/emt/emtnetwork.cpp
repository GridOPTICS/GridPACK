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
#include <model_classes/classical_gen_model.hpp>
#include <model_classes/genrou.hpp>
#include <model_classes/gensal.hpp>
#include <model_classes/exdc1.hpp>
#include <model_classes/esst1a.hpp>
#include <model_classes/wsieg1.hpp>
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
  p_VDQptr = NULL;
  p_nvar = 2;
  p_nrows = 2;
  p_ncols = 2;
  p_neqsgen = NULL;
  p_neqsexc= NULL;
  p_neqsgov= NULL;
  p_gen     = NULL;
  p_num_vals = -1;
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
    free(p_matvalues);
  }
}

/**
   Get number of equations (variables) for this bus
*/
void EmtBus::getNvar(int *nvar) const
{
  *nvar = p_nvar;
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
 * Get voltages in the rectangular form VD, VQ
 * @param double VD - real part of complex voltage at this bus
 * @param double VQ - imaginary part of complex voltage at this bus
 */
void EmtBus::getVoltagesRectangular(double *VD,double *VQ) const
{
  *VD = *p_VDQptr;
  *VQ = *(p_VDQptr+1);
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

  /* Convert from polar to rectangular */
  Va *= pi/180.0;

  double p_VD = Vm*cos(Va);
  double p_VQ = Vm*sin(Va);
  if(p_VDQptr) {
    *p_VDQptr = p_VD;
    *(p_VDQptr+1) = p_VQ;
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

      if(type == "GENCLS") {
        ClassicalGen *clgen;
        clgen = new ClassicalGen;
        p_gen[i] = clgen;
      } else if (type == "GENROU") {
        GenrouGen *grgen;
        grgen = new GenrouGen;
        p_gen[i] = grgen;
      } else if (type == "GENSAL") {
        GensalGen *gsgen;
        gsgen = new GensalGen;
        p_gen[i] = gsgen;
      }

      has_ex = false;
      data->getValue(HAS_EXCITER, &has_ex, i);
      if (has_ex) {
        has_ex = true;
        if (data->getValue(EXCITER_MODEL, &model, i)) {
          type = util.trimQuotes(model);
          if (type == "EXDC1") { 
            Exdc1Exc *e1Exc;
            e1Exc = new Exdc1Exc;
            e1Exc->setGenerator(p_gen[i]);
            boost::shared_ptr<BaseExcModel> ex;
            ex.reset(e1Exc);
            p_gen[i]->setExciter(ex);
          } else if (type == "ESST1A") { 
            Esst1aExc *e1Exc;
            e1Exc = new Esst1aExc;
            e1Exc->setGenerator(p_gen[i]);
            boost::shared_ptr<BaseExcModel> ex;
            ex.reset(e1Exc);
            p_gen[i]->setExciter(ex);
          }
        }
      } 

      has_gv = false;
      data->getValue(HAS_GOVERNOR, &has_gv, i);
      if (has_gv) {
        has_gv = true;
        if (data->getValue(GOVERNOR_MODEL, &model, i)) {
          type = util.trimQuotes(model);
          if (type == "WSIEG1") { 
            Wsieg1Gov *w1Gov;
            w1Gov = new Wsieg1Gov;
            w1Gov->setGenerator(p_gen[i]);
            boost::shared_ptr<BaseGovModel> gv;
            gv.reset(w1Gov);
            p_gen[i]->setGovernor(gv);
          } 
        }
      } 

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
      p_nrows += p_neqsgen[i]+p_neqsexc[i]+p_neqsgov[i];
      p_ncols += p_neqsgen[i]+p_neqsexc[i]+p_neqsgov[i];
      //      printf("p_nvar = %d (ngen:%d,nexc:%d,ngov:%d)\n",p_nvar,p_neqsgen[i],p_neqsexc[i],p_neqsgov[i]);
    }
  }
  p_matvalues = (gridpack::ComplexType**)malloc(p_nvar*sizeof(gridpack::ComplexType*));

  // Get active and reactive power load
  data->getValue(LOAD_PL,&p_pl);
  data->getValue(LOAD_QL,&p_ql);
  // Convert to p.u.
  p_pl /= p_sbase;
  p_ql /= p_sbase;

}

/**
 * Return size of matrix block on the diagonal contributed by component
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix
 *        element
 */
bool EmtBus::matrixDiagSize(int *isize, int *jsize) const
{
  *isize = p_nvar;
  *jsize = *isize;

  return true;
}

/**
 * Return number of rows (dependent variables) that bus contributes
 * to matrix
 * @return number of dependent variables (equations)
 */
int EmtBus::matrixNumRows()
{
  return p_nrows;
}

/**
 * Return number of columns (independent variables) that bus contributes
 * to matrix
 * @return number of independent variables
 */
int EmtBus::matrixNumCols()
{
  return p_ncols;
}

/** 
 * Set global row index
 * @param irow local row index
 * @param idx global row index
 */
void EmtBus::matrixSetRowIndex(int irow, int idx)
{
  if (p_rowidx.size() == 0) {
    p_rowidx.resize(p_nrows);
    int i;
    for (i=0; i<p_nrows; i++) p_rowidx[i] = -1;
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
    p_colidx.resize(p_ncols);
    int i;
    for (i=0; i<p_ncols; i++) p_colidx[i] = -1;
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
//  if (p_num_vals == -1) {
    if (p_isolated) {
      p_num_vals = 2;
    } else {
      std::map<std::pair<int,int>,gridpack::ComplexType> value_map;
      matrixGetValueMap(value_map);
      p_num_vals = value_map.size();
    }
//  }
  return p_num_vals;
}


/**
 * Get value map for generalized matrix interface
 * @param value_map standard map structure that contains all non-zero
 *                  matrix elements
 */
void EmtBus::matrixGetValueMap(std::map<std::pair<int,int>,
    gridpack::ComplexType> &value_map)
{
  int i,j;
  int nvar = p_nvar;
  int VD_col_start = 0; // Starting col index for VD
  int VQ_col_start = nvar; // Starting col index for VQ
  int VD_idx=0; // Location of VD in the solution vector for this bus
  int VQ_idx=1; // Location of VQ in the solution vector for this bus
  double p_VD,p_VQ;

  // Make sure value_map is empty
  value_map.clear();

  getVoltagesRectangular(&p_VD,&p_VQ);

#define MAP_PAIR(_i, _j, _val)                                     \
{                                                                  \
  std::map<std::pair<int,int>,gridpack::ComplexType>::iterator it; \
  it = value_map.find(std::pair<int,int>(_i,_j));                  \
  if (it != value_map.end()) {                                     \
    it->second += _val;                                            \
  } else {                                                         \
    value_map.insert(                                              \
      std::pair<std::pair<int,int>,gridpack::ComplexType>(         \
      std::pair<int,int>(_i,_j),_val));                            \
  }                                                                \
}

  if(p_isolated) {
    MAP_PAIR(0,0,1.0);
    MAP_PAIR(1,1,1.0);
    return;
  }

  std::vector<boost::shared_ptr<BaseComponent> > branches;
  // Get the edges connected to this bus
  gridpack::component::BaseBusComponent::getNeighborBranches(branches);
  // Number of branches connected to this
  int nconnbranch = branches.size();
  for(i=0; i < nconnbranch; i++) {
    EmtBranch *branch = dynamic_cast<EmtBranch*>(branches[i].get());
    EmtBus *busf = dynamic_cast<EmtBus*>((branch->getBus1()).get());
    EmtBus *bust = dynamic_cast<EmtBus*>((branch->getBus2()).get());

    if(this == busf) {
      /* This bus is the from bus of branch[i] */
      double Gff=0.0,Bff=0.0;
      branch->getForwardSelfAdmittance(&Gff,&Bff);
      MAP_PAIR(0,0,-Bff);
      MAP_PAIR(0,1,-Gff);
      MAP_PAIR(1,0,-Gff);
      MAP_PAIR(1,1, Bff);
    } else {
      double Gtt=0.0,Btt=0.0;
      /* This bus is the to bus of branch[i] */
      branch->getReverseSelfAdmittance(&Gtt,&Btt);
      MAP_PAIR(0,0,-Btt);
      MAP_PAIR(0,1,-Gtt);
      MAP_PAIR(1,0,-Gtt);
      MAP_PAIR(1,1, Btt);
    }
  }
  // Partials of contributions from shunt elements
  MAP_PAIR(0,0,-p_bl);
  MAP_PAIR(0,1,-p_gl);
  MAP_PAIR(1,0,-p_gl);
  MAP_PAIR(1,1, p_bl);

  // Partials of contributions from load
  // Assuming constant impedance load
  double yp,yq;
  yp = p_pl/(p_Vm0*p_Vm0);
  yq = p_ql/(p_Vm0*p_Vm0);

  MAP_PAIR(0,0, yq);
  MAP_PAIR(0,1,-yp);
  MAP_PAIR(1,0,-yp);
  MAP_PAIR(1,1,-yq);

  // Partials of generator equations and contributions to the network<->generator 
  EMTMode dsmode;
  for(i=0; i < p_ngen; i++) {
    if(!p_gen[i]->getGenStatus()) continue;

    p_gen[i]->setMode(p_mode);
    p_gen[i]->setVoltage(p_VDQptr[0],p_VDQptr[1]);
    p_gen[i]->setTSshift(p_TSshift);

    p_gen[i]->setJacobian(value_map);

    if(p_gen[i]->hasExciter()) {
      p_gen[i]->getExciter()->setMode(p_mode);
      p_gen[i]->getExciter()->setTSshift(p_TSshift);
      p_gen[i]->getExciter()->setJacobian(value_map);
    }

    if(p_gen[i]->hasGovernor()) {
      p_gen[i]->getGovernor()->setMode(p_mode);
      p_gen[i]->getGovernor()->setTSshift(p_TSshift);
      p_gen[i]->getGovernor()->setJacobian(value_map);
    }
  } 
#undef MAP_PAIR
}

/**
 * Return list of matrix values and their locations generated by the bus
 * @param values list of matrix values
 * @param rows list of row indices
 * @param cols list of column indices
 */
void EmtBus::matrixGetValues(gridpack::ComplexType *values,
    int *rows, int *cols)
{
  int i,j;
  // Get value map for all non-zero indices
  std::map<std::pair<int,int>,gridpack::ComplexType> value_map;
  matrixGetValueMap(value_map);

  int nvals = value_map.size();
  std::map<std::pair<int,int>,gridpack::ComplexType>::iterator it;
  it = value_map.begin();
  int ncnt = 0;
  while (it != value_map.end()) {
    rows[ncnt] = p_rowidx[it->first.first];
    cols[ncnt] = p_colidx[it->first.second];
    values[ncnt] = it->second;
    ncnt++;
    it++;
  }
}

/**
 * Return the values of for a diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute
 *        matrix element
 */
bool EmtBus::matrixDiagValues(gridpack::ComplexType *values)
{
  int i,j;
  int nvar = p_nvar;
  int VD_col_start = 0; // Starting col index for VD
  int VQ_col_start = nvar; // Starting col index for VQ
  int VD_idx=0; // Location of VD in the solution vector for this bus
  int VQ_idx=1; // Location of VQ in the solution vector for this bus
  double p_VD,p_VQ;
  // matvalues is a 2-d representation of the 1-d values array
  // Note that values is column-major, so each array of matvalues
  // points a column of the matrix. Hence, we need to store the values
  // in the transpose of the locations.
  gridpack::ComplexType **matvalues = p_matvalues;

  getVoltagesRectangular(&p_VD,&p_VQ);

  // Zero out values first in case they haven't been zeroed out.
  for(i=0; i < nvar; i++) {
    matvalues[i] = values + i*p_nvar;
    for(j=0; j < nvar; j++) {
      values[nvar*i + j] = 0.0;
    }
  }

  if(p_isolated) {
    matvalues[0][0] = matvalues[1][1] = 1.0;
    return true;
  }

  std::vector<boost::shared_ptr<BaseComponent> > branches;
  // Get the edges connected to this bus
  gridpack::component::BaseBusComponent::getNeighborBranches(branches);
  // Number of branches connected to this
  int nconnbranch = branches.size();
  for(i=0; i < nconnbranch; i++) {
    EmtBranch *branch = dynamic_cast<EmtBranch*>(branches[i].get());
    EmtBus *busf = dynamic_cast<EmtBus*>((branch->getBus1()).get());
    EmtBus *bust = dynamic_cast<EmtBus*>((branch->getBus2()).get());

    if(this == busf) {
      /* This bus is the from bus of branch[i] */
      double Gff=0.0,Bff=0.0;
      branch->getForwardSelfAdmittance(&Gff,&Bff);
      matvalues[0][0] += -Bff;
      matvalues[0][1] += -Gff;
      matvalues[1][0] += -Gff;
      matvalues[1][1] +=  Bff;
    } else {
      double Gtt=0.0,Btt=0.0;
      /* This bus is the to bus of branch[i] */
      branch->getReverseSelfAdmittance(&Gtt,&Btt);
      matvalues[0][0] += -Btt;
      matvalues[0][1] += -Gtt;
      matvalues[1][0] += -Gtt;
      matvalues[1][1] +=  Btt;
    }
  }
  // Partials of contributions from shunt elements
  matvalues[0][0] += -p_bl;
  matvalues[0][1] += -p_gl;
  matvalues[1][0] += -p_gl;
  matvalues[1][1] +=  p_bl;

  // Partials of contributions from load
  // Assuming constant impedance load
  double yp,yq;
  yp = p_pl/(p_Vm0*p_Vm0);
  yq = p_ql/(p_Vm0*p_Vm0);

  matvalues[0][0] +=  yq;
  matvalues[0][1] += -yp;
  matvalues[1][0] += -yp;
  matvalues[1][1] += -yq;

  // Partials of generator equations and contributions to the network<->generator 
  EMTMode dsmode;
  for(i=0; i < p_ngen; i++) {
    if(!p_gen[i]->getGenStatus()) continue;

    p_gen[i]->setMode(p_mode);
    p_gen[i]->setVoltage(p_VDQptr[0],p_VDQptr[1]);
    p_gen[i]->setTSshift(p_TSshift);

    p_gen[i]->setJacobian(matvalues);

    if(p_gen[i]->hasExciter()) {
      p_gen[i]->getExciter()->setMode(p_mode);
      p_gen[i]->getExciter()->setTSshift(p_TSshift);
      p_gen[i]->getExciter()->setJacobian(matvalues);
    }

    if(p_gen[i]->hasGovernor()) {
      p_gen[i]->getGovernor()->setMode(p_mode);
      p_gen[i]->getGovernor()->setTSshift(p_TSshift);
      p_gen[i]->getGovernor()->setJacobian(matvalues);
    }
  } 

  return true;
}

/**
 * Return size of vector block contributed by component
 * @param isize: number of vector elements
 * @return: false if network component does not contribute
 *        vector element
 */
bool EmtBus::vectorSize(int *isize) const
{
  *isize = p_nvar;
  return true;

}

/**
 * Return the values of the vector block
 * @param values: pointer to vector values
 * @return: false if network component does not contribute
 *        vector element
 */
bool EmtBus::vectorValues(gridpack::ComplexType *values)
{
  gridpack::ComplexType *fbus = values;  // First two locations are for the bus current balance equations.
  gridpack::ComplexType *fgen = values + 2;
  bool has_ex=false,has_gv=false;
  if(p_mode == INIT_X) { /* Values go in X */
    int i;
    if(p_isolated) {
      fbus[0] = p_Vm0*cos(p_Va0); // VD
      fbus[1] = p_Vm0*sin(p_Va0); // VQ
      return true;
    }
    else {
      fbus[0]   = p_Vm0*cos(p_Va0); // VD
      fbus[1] = p_Vm0*sin(p_Va0); // VQ

      p_VDQptr[0] = real(values[0]);
      p_VDQptr[1] = real(values[1]);

      for(i=0; i < p_ngen; i++) {
        if(!p_gen[i]->getGenStatus()) continue;

        p_gen[i]->setVoltage(p_VDQptr[0],p_VDQptr[1]);
        p_gen[i]->init(fgen);

        has_ex = p_gen[i]->hasExciter();
        if (has_ex) {
          fgen += p_neqsgen[i];
          p_gen[i]->getExciter()->setVoltage(p_VDQptr[0],p_VDQptr[1]);
          p_gen[i]->getExciter()->init(fgen);
        }	

        has_gv = p_gen[i]->hasGovernor();
        if (has_gv) {
          if (has_ex) fgen += p_neqsexc[i];
          else fgen += p_neqsgen[i];
          p_gen[i]->getGovernor()->init(fgen);
        }	
      }
    }	  
  } else if(p_mode == RESIDUAL_EVAL || p_mode == FAULT_EVAL) { /* Values go in F */
    int i;
    int VD_idx=0; // Location of VD in the solution vector for this bus
    int VQ_idx=1; // Location of VQ in the solution vector for this bus
    double p_VD,p_VQ;
    getVoltagesRectangular(&p_VD,&p_VQ);

    if(p_isolated) {
      fbus[VD_idx] = p_VD - p_Vm0*cos(p_Va0); // f(VD)
      fbus[VQ_idx] = p_VQ - p_Vm0*sin(p_Va0); // f(VQ)
      return true;
    }

    fbus[VD_idx] = fbus[VQ_idx] = 0.0;
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    // Get the edges connected to this bus
    gridpack::component::BaseBusComponent::getNeighborBranches(branches);
    // Number of branches connected to this
    int nconnbranch = branches.size();
    double IbrD=0.0,IbrQ=0.0;
    int thisbusnum = this->getOriginalIndex();

    for(i=0; i < nconnbranch; i++) {
      EmtBranch *branch = dynamic_cast<EmtBranch*>(branches[i].get());
      EmtBus *busf = dynamic_cast<EmtBus*>((branch->getBus1()).get());
      EmtBus *bust = dynamic_cast<EmtBus*>((branch->getBus2()).get());

      int busfnum = branch->getBus1OriginalIndex();
      int bustnum = branch->getBus2OriginalIndex();

      if(busfnum == thisbusnum) { 	  /* This bus is the from bus of branch[i] */
        double VDf,VQf,VDt,VQt;
        VDf = p_VD; VQf = p_VQ;
        /* Get to bus voltages */
        bust->getVoltagesRectangular(&VDt,&VQt);

        double Gff=0.0,Bff=0.0,Gft=0.0,Bft=0.0;
        branch->getForwardSelfAdmittance(&Gff,&Bff);
        branch->getForwardTransferAdmittance(&Gft,&Bft);

        IbrD += Gff*VDf - Bff*VQf + Gft*VDt - Bft*VQt;
        IbrQ += Bff*VDf + Gff*VQf + Bft*VDt + Gft*VQt;
      } else { 	/* This bus is the to bus of branch[i] */
        double VDf,VQf,VDt,VQt;
        VDt = p_VD; VQt = p_VQ;
        /* Get from bus voltages */
        busf->getVoltagesRectangular(&VDf,&VQf);
        double Gtt=0.0,Btt=0.0,Gtf=0.0,Btf=0.0;
        branch->getReverseSelfAdmittance(&Gtt,&Btt);
        branch->getReverseTransferAdmittance(&Gtf,&Btf);

        IbrD += Gtt*VDt - Btt*VQt + Gtf*VDf - Btf*VQf;
        IbrQ += Btt*VDt + Gtt*VQt + Btf*VDf + Gtf*VQf;
      }
    }

    // Current injection from shunts
    double IshuntD,IshuntQ;
    IshuntD = p_gl*p_VD - p_bl*p_VQ;
    IshuntQ = p_bl*p_VD + p_gl*p_VQ;

    // Load current injections
    // Assuming constant impedance load
    double yp,yq;
    double IloadD,IloadQ;
    yp = p_pl/(p_Vm0*p_Vm0);
    yq = p_ql/(p_Vm0*p_Vm0);

    IloadD = yp*p_VD  + yq*p_VQ;
    IloadQ = -yq*p_VD + yp*p_VQ;

    // Generator equations and current contributions from generators
    int ctr=0;
    double IgenD=0.0,IgenQ=0.0;
    for(i=0; i < p_ngen; i++) {
      double IGD=0.0, IGQ=0.0;
      if(!p_gen[i]->getGenStatus()) continue;

      p_gen[i]->setMode(p_mode);
      p_gen[i]->setVoltage(p_VDQptr[0],p_VDQptr[1]);
      p_gen[i]->vectorValues(fgen);
      fgen += p_neqsgen[i];
      p_gen[i]->getCurrent(&IGD,&IGQ);

      has_ex = p_gen[i]->hasExciter();
      if (has_ex) {
        p_gen[i]->getExciter()->setMode(p_mode);
        p_gen[i]->getExciter()->setVoltage(p_VDQptr[0],p_VDQptr[1]);
        p_gen[i]->getExciter()->vectorValues(fgen);
        fgen += p_neqsexc[i];
      }	

      has_gv = p_gen[i]->hasGovernor();
      if (has_gv) {
        p_gen[i]->getGovernor()->setMode(p_mode);
        p_gen[i]->getGovernor()->vectorValues(fgen);
        fgen += p_neqsgov[i];
      }	

      IgenD += IGD;
      IgenQ += IGQ;
    }

    fbus[VD_idx] = IgenQ - IshuntQ - IbrQ - IloadQ;
    fbus[VQ_idx] = IgenD - IshuntD - IbrD - IloadD;

    //    printf("Bus\n");
    //    printf("VD = %f, VQ = %f, fbus[VD_idx] = %f, fbus[VQ_idx] = %f\n",p_VDQptr[0],p_VDQptr[1],real(fbus[VD_idx]),real(fbus[VQ_idx]));
  }
  return true;
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
  return 2*sizeof(double);
}

void EmtBus::setXCBuf(void *buf)
{
  p_VDQptr = static_cast<double*>(buf);

}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto buses 
 * @param values array containing voltage magnitude and angle
 */
void EmtBus::setValues(gridpack::ComplexType *values)
{
  int i,ctr=2;
  gridpack::ComplexType *genvals;
  EMTMode mode;
  bool has_ex=false,has_gv=false;

  if(p_isolated) return;

  mode = p_mode;

  if(p_mode == XVECTOBUS) { // Push values from X vector back on the bus
    *p_VDQptr = real(values[0]);
    *(p_VDQptr+1) = real(values[1]);
  }

  genvals = values + 2;
  for(i=0; i < p_ngen; i++) {
    if(!p_gen[i]->getGenStatus()) continue;

    p_gen[i]->setMode(mode);
    p_gen[i]->setVoltage(p_VDQptr[0],p_VDQptr[1]);
    p_gen[i]->setValues(genvals);
    genvals += p_neqsgen[i];

    has_ex = p_gen[i]->hasExciter();
    if (has_ex) {
      p_gen[i]->getExciter()->setMode(p_mode);
      p_gen[i]->getExciter()->setValues(genvals);
      genvals += p_neqsexc[i];
    }	

    has_gv = p_gen[i]->hasGovernor();
    if (has_gv) {
      p_gen[i]->getGovernor()->setMode(p_mode);
      p_gen[i]->getGovernor()->setValues(genvals);
      genvals += p_neqsgov[i];
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
  p_nrows = 0;
  p_ncols = 0;
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
  double tap, shift;
  double R,X,Bc,G,B,Zm,tap2,tapr,tapi;
  double Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  int    status,i;
  double pi=PI;
  std::string cktid;

  // Get the number of parallel lines
  data->getValue(BRANCH_NUM_ELEMENTS,&p_nparlines);

  p_cktid.reserve(p_nparlines);
  p_Gff.reserve(p_nparlines);
  p_Bff.reserve(p_nparlines);
  p_Gft.reserve(p_nparlines);
  p_Bft.reserve(p_nparlines);
  p_Gtf.reserve(p_nparlines);
  p_Btf.reserve(p_nparlines);
  p_Gtt.reserve(p_nparlines);
  p_Btt.reserve(p_nparlines);
  p_status.reserve(p_nparlines);

  for(i=0; i < p_nparlines; i++) {
    // Get line parameters
    data->getValue(BRANCH_STATUS,&status,i);
    data->getValue(BRANCH_R,&R,i);
    data->getValue(BRANCH_X,&X,i);
    data->getValue(BRANCH_B,&Bc,i);
    data->getValue(BRANCH_TAP,&tap,i);
    data->getValue(BRANCH_SHIFT,&shift,i);
    data->getValue(BRANCH_CKT,&cktid,i);

    p_cktid.push_back(cktid);
    p_status.push_back(status);
    // We calculate the two-port network elements for the branch irrespecitive
    // of whether the branch is on or off. This is done so that it can be
    // reused later if branch switching is done

    Zm = R*R + X*X;
    G  = R/Zm;
    B  = -X/Zm;

    if(tap == 0.0) tap = 1.0;
    shift *= pi/180.0;
    tap2 = tap*tap;
    tapr = tap*cos(shift);
    tapi = tap*sin(shift);

    Gff = G/tap2;
    Bff = (B+Bc/2.0)/tap2;

    Gft = -(G*tapr - B*tapi)/tap2;
    Bft = -(B*tapr + G*tapi)/tap2;

    Gtf = -(G*tapr + B*tapi)/tap2;
    Btf = -(B*tapr - G*tapi)/tap2;

    Gtt = G;
    Btt = B+Bc/2.0;

    p_Gff.push_back(Gff);
    p_Bff.push_back(Bff);
    p_Gft.push_back(Gft);
    p_Bft.push_back(Bft);
    p_Gtf.push_back(Gtf);
    p_Btf.push_back(Btf);
    p_Gtt.push_back(Gtt);
    p_Btt.push_back(Btt);

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
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool EmtBranch::matrixForwardSize(int *isize, int *jsize) const
{
  // Get from and to buses
  EmtBus *busf = dynamic_cast<EmtBus*>(getBus1().get());
  EmtBus *bust = dynamic_cast<EmtBus*>(getBus2().get());
  // If either the from or to bus is isolated then there is no contribution to the matrix as there
  // is no flow on the branch
  if(busf->isIsolated() || bust->isIsolated()) return false;

  int nvarf,nvart;
  busf->getNvar(&nvarf);
  bust->getNvar(&nvart);
  *isize = nvarf;
  *jsize = nvart;

  return true;
}

/**
 * Return number of rows (dependent variables) that branch contributes
 * to matrix
 * @return number of dependent variables (equations)
 */
int EmtBranch::matrixNumRows()
{
  return p_nrows;
}

/**
 * Return number of columns (independent variables) that branch contributes
 * to matrix
 * @return number of independent variables
 */
int EmtBranch::matrixNumCols()
{
  return p_ncols;
}

/** 
 * Set global row index
 * @param irow local row index
 * @param global row index
 */
void EmtBranch::matrixSetRowIndex(int irow, int idx)
{
  if (p_rowidx.size() == 0) {
    p_rowidx.resize(p_nrows);
    int i;
    for (i=0; i<p_nrows; i++) p_rowidx[i] = -1;
  }
  p_rowidx[irow] = idx;
}

/** 
 * Set global column index
 * @param icol local column index
 * @param global column index
 */
void EmtBranch::matrixSetColIndex(int icol, int idx)
{
  if (p_colidx.size() == 0) {
    p_colidx.resize(p_ncols);
    int i;
    for (i=0; i<p_ncols; i++) p_colidx[i] = -1;
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
  if (p_num_vals == 0) {
  }
  return p_num_vals;
}

/**
 * Return list of matrix values and their locations generated by the branch
 * @param values list of matrix values
 * @param rows list of row indices
 * @param cols list of column indices
 */
void EmtBranch::matrixGetValues(gridpack::ComplexType *values,
    int *rows, int *cols)
{
}

bool EmtBranch::matrixReverseSize(int *isize, int *jsize) const
{
  // Get from and to buses
  EmtBus *busf = dynamic_cast<EmtBus*>(getBus1().get());
  EmtBus *bust = dynamic_cast<EmtBus*>(getBus2().get());
  // If either the from or to bus is isolated then there is no contribution to the matrix as there
  // is no flow on the branch
  if(busf->isIsolated() || bust->isIsolated()) return false;

  int nvarf,nvart;
  busf->getNvar(&nvarf);
  bust->getNvar(&nvart);
  *isize = nvart;
  *jsize = nvarf;

  return true;
}

/**
 * Return the values of the forward/reverse matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool EmtBranch::matrixForwardValues(gridpack::ComplexType *values)
{
  int i;
  // Get from and to buses
  EmtBus *busf = dynamic_cast<EmtBus*>(getBus1().get());
  EmtBus *bust = dynamic_cast<EmtBus*>(getBus2().get());

  if(busf->isGhost()) return false;

  // If either the from or to bus is isolated then there is no contribution to the matrix as there
  // is no flow on the branch
  if(busf->isIsolated() || bust->isIsolated()) return false;

  int nvarf,nvart;
  busf->getNvar(&nvarf);
  bust->getNvar(&nvart);

  values[0] = values[1] = values[2] = values[3] = 0.0;
  if(p_mode == INIT_X || p_mode == RESIDUAL_EVAL || p_mode == FAULT_EVAL) {
    for(i=0; i < p_nparlines;i++) {
      if(p_status[i]) {
        values[0]       += -p_Bft[i];
        values[nvarf]   += -p_Gft[i];
        values[1]       += -p_Gft[i];
        values[nvarf+1] +=  p_Bft[i];
      }
    }
  }  
  return true;
}

bool EmtBranch::matrixReverseValues(gridpack::ComplexType *values)
{
  int i;
  // Get from and to buses
  EmtBus *busf = dynamic_cast<EmtBus*>(getBus1().get());
  EmtBus *bust = dynamic_cast<EmtBus*>(getBus2().get());
  // If either the from or to bus is isolated then there is no contribution to the matrix as there
  // is no flow on the branch
  if(busf->isIsolated() || bust->isIsolated()) return false;

  if(bust->isGhost()) return false;

  int nvarf,nvart;
  busf->getNvar(&nvarf);
  bust->getNvar(&nvart);


  if(p_mode == INIT_X || p_mode == RESIDUAL_EVAL || p_mode == FAULT_EVAL) {
    for(i=0; i < p_nparlines;i++) {
      if(p_status[i]) {
        values[0]       += -p_Btf[i];
        values[nvart]   += -p_Gtf[i];
        values[1]       += -p_Gtf[i];
        values[nvart+1] +=  p_Btf[i];
      }
    }
  }  
  return true;
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

bool EmtBranch::getForwardSelfAdmittance(double *Gff,double *Bff)
{
  int i;
  *Gff = *Bff = 0.0;
  for(i=0; i < p_nparlines; i++) {
    if(p_status[i]) {
      *Gff += p_Gff[i];
      *Bff += p_Bff[i];
    }
  }
  return true;
}

bool EmtBranch::getReverseSelfAdmittance(double *Gtt,double *Btt)
{
  int i;
  *Gtt = *Btt = 0.0;
  for(i=0; i < p_nparlines; i++) {
    if(p_status[i]) {
      *Gtt += p_Gtt[i];
      *Btt += p_Btt[i];
    }
  }
  return true;
}

bool EmtBranch::getForwardTransferAdmittance(double *Gft,double *Bft)
{
  int i;
  *Gft = *Bft = 0.0;
  for(i=0; i < p_nparlines; i++) {
    if(p_status[i]) {
      *Gft += p_Gft[i];
      *Bft += p_Bft[i];
    }
  }
  return true;
}

bool EmtBranch::getReverseTransferAdmittance(double *Gtf,double *Btf)
{
  int i;
  *Gtf = *Btf = 0.0;
  for(i=0; i < p_nparlines; i++) {
    if(p_status[i]) {
      *Gtf += p_Gtf[i];
      *Btf += p_Btf[i];
    }
  }
  return true;
}
