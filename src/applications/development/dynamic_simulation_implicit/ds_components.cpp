/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_components.cpp
 * @author Shrirang Abhyankar
 * @date   2016-07-14 14:29:57 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "ds_components.hpp"
#include "gridpack/parser/dictionary.hpp"

/**
 *  Simple constructor
 */
gridpack::dsimplicit::DSBus::DSBus(void)
{
  p_gl = p_bl = 0.0;
  p_pl = p_ql = 0.0;
  p_ngen = p_nactivegen = 0;
  p_isolated = false;
  p_mode = INIT_X;
  p_ws = 2*22.0/7.0*60.0;
}

/**
 *  Simple destructor
 */
gridpack::dsimplicit::DSBus::~DSBus(void)
{
}

/**
   Get number of equations (variables) for this bus
*/
void gridpack::dsimplicit::DSBus::getNvar(int *nvar) const
{
  *nvar = 2 + 2*p_nactivegen;
}

/**
  Set the shift value provided by TS
*/
void gridpack::dsimplicit::DSBus::setTSshift(double shift)
{
  p_TSshift = shift;
}

/**
  *  Check if the bus is isolated. Returns true if the bus is isolated

*/
bool gridpack::dsimplicit::DSBus::isIsolated(void) const
{
  bool isol = false;
  if(p_isolated) isol = true;
  return isol;
}

/**
 * Get voltages in the rectangular form VD, VQ
 * @param double VD - real part of complex voltage at this bus
 * @param double VQ - imaginary part of complex voltage at this bus
 */
void gridpack::dsimplicit::DSBus::getVoltagesRectangular(double *VD,double *VQ) const
{
  *VD = p_VD;
  *VQ = p_VQ;
}


/**
 * Load values stored in DataCollection object into DSBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 *
 * Note: It is assumed that the network data file read has data for "solved"
 * power flow so we can do the initialization when loading the data
 */
void gridpack::dsimplicit::DSBus::load(const
         boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  double MVAbase; // Base MVA 
  double Vm, Va;  // Voltage magnitude and angle from the case file
  double pi=22.0/7.0;
  int    i=0;
  double Pg, Qg,mbase,Rs,Xdp,H,D; // temp. variables to hold machine data
  int    gstatus;
  double delta,dw=0.0;
  double Pm,Ep;
  double IGD,IGQ; /* D and Q components of generator currents */
  int    bustype; // Type of bus

  MVAbase = 100.0;

  // Get MVAbase
  data->getValue(CASE_SBASE,&MVAbase);

  // Get Voltage magnitude and angle
  data->getValue(BUS_VOLTAGE_ANG,&Va); // This is in degress
  data->getValue(BUS_VOLTAGE_MAG,&Vm);

  // 
  /* Each bus has variables VD and VQ,
     In addition if active generators are incident on the bus
     then each generator has additional variables delta and dw
     Note that all generators are modeled as GENCLS
  */
  /* Convert from polar to rectangular */
  Va *= pi/180.0;

  p_VD = Vm*cos(Va);
  p_VQ = Vm*sin(Va);

  /* Save the voltage magnitude and angle so that we can use it later for building loads */
  p_Vm0 = Vm;
  p_Va0 = Va;

  // Get the bus type
  data->getValue(BUS_TYPE,&bustype);
  if(bustype == 4) p_isolated = true;

  if(p_isolated) return;

  // Get active and reactive power load
  data->getValue(LOAD_PL,&p_pl);
  data->getValue(LOAD_QL,&p_ql);
  // Convert to p.u.
  p_pl /= MVAbase;
  p_ql /= MVAbase;

  // Read shunts
  data->getValue(BUS_SHUNT_GL,&p_gl);
  data->getValue(BUS_SHUNT_BL,&p_bl);
  p_gl /= MVAbase;
  p_bl /= MVAbase;

  // Get number of generators incident on this bus
  data->getValue(GENERATOR_NUMBER, &p_ngen);
  if(p_ngen) {
    // Allocate containers
    p_gstatus.reserve(p_ngen);
    p_pg.reserve(p_ngen);   
    p_qg.reserve(p_ngen);
    p_mbase.reserve(p_ngen);
    p_Rs.reserve(p_ngen);
    p_Xdp.reserve(p_ngen);
    p_H.reserve(p_ngen);
    p_D.reserve(p_ngen);
    p_Pm.reserve(p_ngen);
    p_Ep.reserve(p_ngen);
    p_delta.reserve(p_ngen);
    p_dw.reserve(p_ngen);
    p_deltadot.reserve(p_ngen);
    p_dwdot.reserve(p_ngen);
    
    // Read generator data stored in data collection objects
    for(i=0; i < p_ngen; i++) {
      // Generator real and reactive power
      data->getValue(GENERATOR_PG,&Pg,i);
      data->getValue(GENERATOR_QG,&Qg,i);
      Pg /= MVAbase;
      Qg /= MVAbase;

      // Generator parameters
      data->getValue(GENERATOR_STAT,&gstatus,i);
      data->getValue(GENERATOR_MBASE,&mbase,i);
      data->getValue(GENERATOR_RESISTANCE,&Rs,i);
      data->getValue(GENERATOR_TRANSIENT_REACTANCE,&Xdp,i);
      data->getValue(GENERATOR_INERTIA_CONSTANT_H,&H,i);
      data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&D,i);

      // Convert generator parameters from machine base to MVA base
      H *= mbase/MVAbase;
      D *= mbase/MVAbase;
      Xdp /= mbase/MVAbase;

      p_pg.push_back(Pg);
      p_qg.push_back(Qg);
      p_gstatus.push_back(gstatus);
      p_mbase.push_back(mbase);
      p_Rs.push_back(Rs);
      p_Xdp.push_back(Xdp);
      p_H.push_back(H);
      p_D.push_back(D);

      // Initialize the generator dynamic variables delta and dw and controller constants
      if(gstatus) {
	p_nactivegen++; // Increment the number of active generators
	IGD = (p_VD*Pg + p_VQ*Qg)/(Vm*Vm);
	IGQ = (p_VQ*Pg - p_VD*Qg)/(Vm*Vm);

	delta = atan2(p_VQ + Xdp*IGD,p_VD-Xdp*IGQ);
	p_delta.push_back(delta);
	p_dw.push_back(dw);
	
	Ep = sqrt(pow((p_VD - Xdp*IGQ),2) + pow((p_VQ + Xdp*IGD),2));
	p_Ep.push_back(Ep);
	p_Pm.push_back(Pg);
      } else {
	p_delta.push_back(0.0);
	p_dw.push_back(0.0);
	p_Pm.push_back(0.0);
	p_Ep.push_back(0.0);
      }
      p_deltadot.push_back(0.0);
      p_dwdot.push_back(0.0);
    }
  }
}

/**
 * Return size of matrix block on the diagonal contributed by component
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix
 *        element
 */
bool gridpack::dsimplicit::DSBus::matrixDiagSize(int *isize, int *jsize) const
{
  *isize = 2+ 2*p_nactivegen;
  *jsize = *isize;

  return true;
}

/**
 * Return the values of for a diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute
 *        matrix element
 */
bool gridpack::dsimplicit::DSBus::matrixDiagValues(ComplexType *values)
{
  int i,j;
  int nvar = 2 + 2*p_nactivegen;
  int VD_col_start = 0; // Starting col index for VD
  int VQ_col_start = nvar; // Starting col index for VQ
  int delta_col_start, dw_col_start;
  int VD_idx=0; // Location of VD in the solution vector for this bus
  int VQ_idx=1; // Location of VQ in the solution vector for this bus
  int delta_idx; // Location of delta for geni in the soluton vector for this bus
  int dw_idx;   // Location of dw for geni in the solution vector for this bus
  
  // Zero out values first in case they haven't been zeroed out.
  for(i=0; i < nvar; i++) {
    for(j=0; j < nvar; j++) {
      values[nvar*i + j] = 0.0;
    }
  }
   
 if(p_isolated) {
   values[VD_col_start+VD_idx] = values[VQ_col_start+VQ_idx] = 1.0;
   values[VD_col_start+VQ_idx] = values[VQ_col_start+VD_idx] = 0.0;
   return true;
 }
 
 std::vector<boost::shared_ptr<BaseComponent> > branches;
 // Get the edges connected to this bus
 gridpack::component::BaseBusComponent::getNeighborBranches(branches);
 // Number of branches connected to this
 int nconnbranch = branches.size();
 for(i=0; i < nconnbranch; i++) {
   gridpack::dsimplicit::DSBranch *branch = dynamic_cast<gridpack::dsimplicit::DSBranch*>(branches[i].get());
   gridpack::dsimplicit::DSBus *busf = dynamic_cast<gridpack::dsimplicit::DSBus*>((branch->getBus1()).get());
   gridpack::dsimplicit::DSBus *bust = dynamic_cast<gridpack::dsimplicit::DSBus*>((branch->getBus2()).get());
   
   if(this == busf) {
     /* This bus is the from bus of branch[i] */
     double Gff=0.0,Bff=0.0;
     branch->getForwardSelfAdmittance(&Gff,&Bff);
     values[VD_col_start + VD_idx] += -Bff;
     values[VQ_col_start + VD_idx] += -Gff;
     values[VD_col_start + VQ_idx] += -Gff;
     values[VQ_col_start + VQ_idx] +=  Bff;
   } else {
     double Gtt=0.0,Btt=0.0;
     /* This bus is the to bus of branch[i] */
     branch->getReverseSelfAdmittance(&Gtt,&Btt);
     values[VD_col_start + VD_idx] += -Btt;
     values[VQ_col_start + VD_idx] += -Gtt;
     values[VD_col_start + VQ_idx] += -Gtt;
     values[VQ_col_start + VQ_idx] +=  Btt;
   }
 }
 // Partials of contributions from shunt elements
 values[VD_col_start + VD_idx] += -p_bl;
 values[VQ_col_start + VD_idx] += -p_gl;
 values[VD_col_start + VQ_idx] += -p_gl;
 values[VQ_col_start + VQ_idx] +=  p_bl;
 
 // Partials of contributions from load
 // Assuming constant impedance load
 double yp,yq;
 yp = p_pl/(p_Vm0*p_Vm0);
 yq = p_ql/(p_Vm0*p_Vm0);
 values[VD_col_start + VD_idx] +=  yq;
 values[VQ_col_start + VD_idx] += -yp;
 values[VD_col_start + VQ_idx] += -yp;
 values[VQ_col_start + VQ_idx] += -yq;
 
 // Partials of generator equations and contributions to the network<->generator
 int ctr=0;
 for(i=0; i < p_ngen; i++) {
   if(p_gstatus[i]) {
     delta_col_start = (2+ctr)*nvar;
     dw_col_start    = (2+ctr+1)*nvar;
     delta_idx       = 2 + 2*ctr;
     dw_idx          = 2 + 2*ctr + 1;
     
     if(p_mode == FAULT_EVAL) {
       values[delta_col_start+delta_idx] = 1.0;
       values[dw_col_start+dw_idx]       = 1.0;

       // Partials of generator current injections into the network w.r.t VD, VQ
       values[VD_col_start + VD_idx] +=  1/p_Xdp[i];
       values[VQ_col_start + VQ_idx] +=  -1/p_Xdp[i];
     } else {
     // Partials of generator equations w.r.t generator variables
       values[delta_col_start+delta_idx] = -p_TSshift; 
       values[dw_col_start+delta_idx]    = 1.0/p_ws;
       values[delta_col_start+dw_idx]    = (-p_VD*p_Ep[i]*cos(p_delta[i])/p_Xdp[i] - p_VQ*p_Ep[i]*sin(p_delta[i])/p_Xdp[i])/(2*p_H[i]);
       values[dw_col_start+dw_idx]       = -p_TSshift - p_D[i]/(2*p_H[i]);
       
       // Partials of generator equations w.r.t VD, VQ
       values[VD_col_start+dw_idx] = (-p_Ep[i]*sin(p_delta[i])/p_Xdp[i])/(2*p_H[i]);
       values[VQ_col_start+dw_idx] = (p_Ep[i]*cos(p_delta[i])/p_Xdp[i])/(2*p_H[i]);
       
       // Partials of generator current injections into the network w.r.t generator variables
       values[delta_col_start + VD_idx] = p_Ep[i]*sin(p_delta[i])/p_Xdp[i];
       values[delta_col_start + VQ_idx] = p_Ep[i]*cos(p_delta[i])/p_Xdp[i];
       
       // Partials of generator current injections into the network w.r.t VD, VQ
       values[VD_col_start + VD_idx] +=  1/p_Xdp[i];
       values[VQ_col_start + VQ_idx] +=  -1/p_Xdp[i];
     }
     ctr += 2;
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
bool gridpack::dsimplicit::DSBus::vectorSize(int *isize) const
{
  // [V_D, V_Q] + if(p_nactivegen) => [delta_i,dw] for each gen
  *isize = 2 + 2*p_nactivegen;
  return true;
}

/**
 * Return the values of the vector block
 * @param values: pointer to vector values
 * @return: false if network component does not contribute
 *        vector element
 */
bool gridpack::dsimplicit::DSBus::vectorValues(ComplexType *values)
{
  if(p_mode == INIT_X) { /* Values go in X */
    int i,ctr=0;
    if(p_isolated) {
      values[0] = p_Vm0*cos(p_Va0); // VD
      values[1] = p_Vm0*sin(p_Va0); // VQ
      return true;
    } else {
      values[ctr] = p_Vm0*cos(p_Va0); // VD
      values[ctr+1] = p_Vm0*sin(p_Va0); // VQ
      ctr += 2;
      for(i=0; i < p_ngen; i++) {
	if(p_gstatus[i]) {
	  values[ctr] = p_delta[i];
	  values[ctr+1] = p_dw[i];
	  ctr += 2;
	}
      }
    } 
  } else if(p_mode == RESIDUAL_EVAL || p_mode == FAULT_EVAL) { /* Values go in F */
    int i;
    int VD_idx=0; // Location of VD in the solution vector for this bus
    int VQ_idx=1; // Location of VQ in the solution vector for this bus
    int delta_idx; // Location of delta for geni in the soluton vector for this bus
    int dw_idx;   // Location of dw for geni in the solution vector for this bus
    
    if(p_isolated) {
      values[VD_idx] = p_VD - p_Vm0*cos(p_Va0); // f(VD)
      values[VQ_idx] = p_VQ - p_Vm0*sin(p_Va0); // f(VQ)
      return true;
    }
      
    values[VD_idx] = values[VQ_idx] = 0.0;
    std::vector<boost::shared_ptr<BaseComponent> > branches;
    // Get the edges connected to this bus
    gridpack::component::BaseBusComponent::getNeighborBranches(branches);
    // Number of branches connected to this
    int nconnbranch = branches.size();
    //    printf("Nconnbranch = %d\n",nconnbranch);
    double IbrD=0.0,IbrQ=0.0;
    for(i=0; i < nconnbranch; i++) {
      gridpack::dsimplicit::DSBranch *branch = dynamic_cast<gridpack::dsimplicit::DSBranch*>(branches[i].get());
      gridpack::dsimplicit::DSBus *busf = dynamic_cast<gridpack::dsimplicit::DSBus*>((branch->getBus1()).get());
      gridpack::dsimplicit::DSBus *bust = dynamic_cast<gridpack::dsimplicit::DSBus*>((branch->getBus2()).get());
      
      if(this == busf) { 	  /* This bus is the from bus of branch[i] */
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
      if(p_gstatus[i]) {
	delta_idx       = 2 + 2*ctr;
	dw_idx          = 2 + 2*ctr + 1;
	
	if(p_mode == FAULT_EVAL) {
	  values[delta_idx] = values[dw_idx] = 0.0;
	} else {
	  // Generator equations
	  values[delta_idx] = p_dw[i]/p_ws - p_deltadot[i];
	  values[dw_idx]    = (p_Pm[i] - p_VD*p_Ep[i]*sin(p_delta[i])/p_Xdp[i] + p_VQ*p_Ep[i]*cos(p_delta[i])/p_Xdp[i] - p_D[i]*p_dw[i])/(2*p_H[i]) - p_dwdot[i];
	}

	// Generator current injections in the network
	IgenD += (-p_VQ + p_Ep[i]*sin(p_delta[i]))/p_Xdp[i];
	IgenQ += (p_VD - p_Ep[i]*cos(p_delta[i]))/p_Xdp[i];

	ctr += 2;
      }
    }

    values[VD_idx] = IgenQ - IshuntQ - IbrQ - IloadQ;
    values[VQ_idx] = IgenD - IshuntD - IbrD - IloadD;

    //    printf("Bus %d:misQ = %4.3f  IgenQ = %4.3f  IshuntQ = %4.3f  IbrQ = %4.3f  IloadQ = %4.3f\n",getOriginalIndex(),real(values[VD_idx]),IgenQ,IshuntQ,IbrQ,IloadQ);
    //    printf("Bus %d:misD = %4.3f  IgenD = %4.3f  IshuntD = %4.3f  IbrD = %4.3f  IloadD = %4.3f\n",getOriginalIndex(),real(values[VQ_idx]),IgenD,IshuntD,IbrD,IloadD);

  }
  return true;
}

/**
 * Set the model to control what matrices and vectors and built when using the
 * mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::dsimplicit::DSBus::setMode(int mode)
{
  p_mode = mode;
}

int gridpack::dsimplicit::DSBus::getXCBufSize(void)
{
  return (2+2*p_nactivegen)*sizeof(double);
}

void gridpack::dsimplicit::DSBus::setXCBuf(void *buf)
{
  double* ptr = static_cast<double*>(buf);
  int i,ctr=2;
  
  ptr[0] = p_VD;
  ptr[1] = p_VQ;

  for(i=0; i < p_ngen; i++) {
    if(p_gstatus[i]) {
      ptr[ctr]   = p_delta[i];
      ptr[ctr+1] = p_dw[i];
      ctr += 2;
    }
  }
}
/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto buses 
 * @param values array containing voltage magnitude and angle
 */
void gridpack::dsimplicit::DSBus::setValues(gridpack::ComplexType *values)
{
  int i,ctr=2;
  if(p_mode == XVECTOBUS) { // Push values from X vector back onto the bus 
    p_VD = real(values[0]);
    p_VQ = real(values[1]);
    if(!p_isolated) {
      // Push the generator state variables from X onto the bus
      for(i=0; i < p_ngen; i++) {
	if(p_gstatus[i]) {
	  p_delta[i] = real(values[ctr]);
	  p_dw[i]    = real(values[ctr+1]);
	  ctr += 2;
	}
      }
    }
  } else if(p_mode == XDOTVECTOBUS) { // Push the derivatives of delta and dw onto the bus
    if(p_isolated) return;
    for(i=0; i < p_ngen; i++) {
      if(p_gstatus[i]) {
	p_deltadot[i] = real(values[ctr]);
	p_dwdot[i]    = real(values[ctr+1]);
	ctr += 2;
      }
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
bool gridpack::dsimplicit::DSBus::serialWrite(char *string,
    const int bufsize, const char *signal)
{
}

/**
 *  Simple constructor
 */
gridpack::dsimplicit::DSBranch::DSBranch(void)
{
  p_nparlines = 0;
  p_mode = INIT_X;
}

/**
 *  Simple destructor
 */
gridpack::dsimplicit::DSBranch::~DSBranch(void)
{
}

/**
 * Load values stored in DataCollection object into DSBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 *
 * Note: GridPACK stores the data for all parallel branches in the same
 *       DataCollection object. 
 */
void gridpack::dsimplicit::DSBranch::load(
    const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  double tap, shift;
  double R,X,Bc,G,B,Zm,tap2,tapr,tapi;
  double Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
  int    status,i;
  double pi=22.0/7.0;
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
void gridpack::dsimplicit::DSBranch::setMode(int mode)
{
  p_mode = mode;
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::dsimplicit::DSBranch::matrixForwardSize(int *isize, int *jsize) const
{
  // Get from and to buses
  gridpack::dsimplicit::DSBus *busf = dynamic_cast<gridpack::dsimplicit::DSBus*>(getBus1().get());
  gridpack::dsimplicit::DSBus *bust = dynamic_cast<gridpack::dsimplicit::DSBus*>(getBus2().get());
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

bool gridpack::dsimplicit::DSBranch::matrixReverseSize(int *isize, int *jsize) const
{
  // Get from and to buses
  gridpack::dsimplicit::DSBus *busf = dynamic_cast<gridpack::dsimplicit::DSBus*>(getBus1().get());
  gridpack::dsimplicit::DSBus *bust = dynamic_cast<gridpack::dsimplicit::DSBus*>(getBus2().get());
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
bool gridpack::dsimplicit::DSBranch::matrixForwardValues(ComplexType *values)
{
  int i;
  // Get from and to buses
  gridpack::dsimplicit::DSBus *busf = dynamic_cast<gridpack::dsimplicit::DSBus*>(getBus1().get());
  gridpack::dsimplicit::DSBus *bust = dynamic_cast<gridpack::dsimplicit::DSBus*>(getBus2().get());
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

bool gridpack::dsimplicit::DSBranch::matrixReverseValues(ComplexType *values)
{
  int i;
  // Get from and to buses
  gridpack::dsimplicit::DSBus *busf = dynamic_cast<gridpack::dsimplicit::DSBus*>(getBus1().get());
  gridpack::dsimplicit::DSBus *bust = dynamic_cast<gridpack::dsimplicit::DSBus*>(getBus2().get());
  // If either the from or to bus is isolated then there is no contribution to the matrix as there
  // is no flow on the branch
  if(busf->isIsolated() || bust->isIsolated()) return false;

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
bool gridpack::dsimplicit::DSBranch::serialWrite(char *string, const int
    bufsize,  const char *signal)
{
  return false;
}

bool gridpack::dsimplicit::DSBranch::getForwardSelfAdmittance(double *Gff,double *Bff)
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

bool gridpack::dsimplicit::DSBranch::getReverseSelfAdmittance(double *Gtt,double *Btt)
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

bool gridpack::dsimplicit::DSBranch::getForwardTransferAdmittance(double *Gft,double *Bft)
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

bool gridpack::dsimplicit::DSBranch::getReverseTransferAdmittance(double *Gtf,double *Btf)
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
