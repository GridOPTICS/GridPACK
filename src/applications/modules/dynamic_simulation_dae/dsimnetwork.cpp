/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsimnetwork.cpp
 * @author Shrirang Abhyankar
 * @date   02/09/19
 * @last modified by Shuangshuang Jin on 10/19/2019 to add
 *  detailed generator, exciter and governor models.
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include <iostream>
#include <dsimnetwork.hpp>
#include <gridpack/include/gridpack.hpp>
#include <gridpack/utilities/complex.hpp>
#include <constants.hpp>
#include <model_classes/classical_gen_model.hpp>
#include <model_classes/genrou.hpp>
#include <model_classes/gensal.hpp>
#include <model_classes/exdc1.hpp>
#include <model_classes/esst1a.hpp>
#include <model_classes/wsieg1.hpp>

/**
 *  Simple constructor
 */
DSimBus::DSimBus(void)
{
  p_gl = p_bl = 0.0;
  p_pl = p_ql = 0.0;
  p_ngen = p_nactivegen = 0;
  p_isolated = false;
  p_mode = NONE;
  p_VDQptr = NULL;
  p_nvar = 2;
  p_neqsgen = NULL;
  p_neqsexc= NULL;
  p_neqsgov= NULL;
  p_gen     = NULL;
}

/**
 *  Simple destructor
 */
DSimBus::~DSimBus(void)
{
  if(p_ngen) {
    for(int i=0; i < p_ngen; i++) {
      if(p_gen[i]) delete(p_gen[i]);
    }
  }
  free(p_neqsgen);
  free(p_neqsexc);
  free(p_neqsgov);
  free(p_gen);
}

/**
   Get number of equations (variables) for this bus
*/
void DSimBus::getNvar(int *nvar) const
{
  *nvar = p_nvar;
}

/**
  Set the shift value provided by TS
*/
void DSimBus::setTSshift(double shift)
{
  p_TSshift = shift;
}

/**
  *  Check if the bus is isolated. Returns true if the bus is isolated

*/
bool DSimBus::isIsolated(void) const
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
void DSimBus::getVoltagesRectangular(double *VD,double *VQ) const
{
  *VD = *p_VDQptr;
  *VQ = *(p_VDQptr+1);
}


/**
 * Load values stored in DataCollection object into DSimBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 *
 * Note: It is assumed that the network data file read has data for "solved"
 * power flow so we can do the initialization when loading the data
 */
void DSimBus::load(const
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
  printf("p_sbase = %f\n", p_sbase);

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
  //bool has_ex;
  //bool has_gv;

  // Read Generators 
  // Get number of generators incident on this bus
  data->getValue(GENERATOR_NUMBER, &p_ngen);
  if(p_ngen) {
    p_gen = (BaseGenModel**)malloc(p_ngen*sizeof(BaseGenModel*));
    p_neqsgen = (int*)malloc(p_ngen*sizeof(int));
    p_neqsexc= (int*)malloc(p_ngen*sizeof(int));
    p_neqsgov= (int*)malloc(p_ngen*sizeof(int));

    for(i=0; i < p_ngen; i++) {
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
          if (type == "EXDC1") { // SJin: newly added EXDC1 exciter model object 
            Exdc1Exc *e1Exc;
            e1Exc = new Exdc1Exc;
            boost::shared_ptr<BaseExcModel> ex;
            ex.reset(e1Exc);
            p_gen[i]->setExciter(ex);
            printf("set exciter EXDC1 successfully!\n");
          } else if (type == "ESST1A") { // SJin: newly added EXDC1 exciter model object 
            Esst1aExc *e1Exc;
            e1Exc = new Esst1aExc;
            boost::shared_ptr<BaseExcModel> ex;
            ex.reset(e1Exc);
            p_gen[i]->setExciter(ex);
            printf("set exciter ESST1a successfully!\n");
          }
        }
      } 
        
      // SJin: add governors if there's any
      has_gv = false;
      data->getValue(HAS_GOVERNOR, &has_gv, i);
      if (has_gv) {
        has_gv = true;
        if (data->getValue(GOVERNOR_MODEL, &model, i)) {
          type = util.trimQuotes(model);
          if (type == "WSIEG1") { // SJin: newly added WSIEG1 governor model object 
            Wsieg1Gov *w1Gov;
            w1Gov = new Wsieg1Gov;
            boost::shared_ptr<BaseGovModel> gv;
            gv.reset(w1Gov);
            p_gen[i]->setGovernor(gv);
            printf("set governor WSIEG1 successfully!\n");
          } 
        }
      } 

      // Read generator data stored in data collection objects
      p_gen[i]->load(data,i); // load data
      if (has_ex) p_gen[i]->getExciter()->load(data,i); // load exciter data
      if (has_gv) p_gen[i]->getGovernor()->load(data,i); // load governor model

      // Set number of equations for this generator
      p_gen[i]->vectorSize(&p_neqsgen[i]);
      if (has_ex) p_gen[i]->getExciter()->vectorSize(&p_neqsexc[i]);
      if (has_gv) p_gen[i]->getGovernor()->vectorSize(&p_neqsgov[i]);

      /* Update number of variables for this bus */
      //p_nvar += p_neqsgen[i];
      p_nvar += p_neqsgen[i]+p_neqsexc[i]+p_neqsgov[i];
      printf("p_nvar = %d (ngen:%d,nexc:%d,ngov:%d)\n",p_nvar,p_neqsgen[i],p_neqsexc[i],p_neqsgov[i]);
    }
  }

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
bool DSimBus::matrixDiagSize(int *isize, int *jsize) const
{
  *isize = p_nvar;
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
bool DSimBus::matrixDiagValues(gridpack::ComplexType *values)
{
  int i,j;
  int nvar = p_nvar;
  int VD_col_start = 0; // Starting col index for VD
  int VQ_col_start = nvar; // Starting col index for VQ
  int VD_idx=0; // Location of VD in the solution vector for this bus
  int VQ_idx=1; // Location of VQ in the solution vector for this bus
  double p_VD,p_VQ;
  int ctr=0;

  getVoltagesRectangular(&p_VD,&p_VQ);
  
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
   DSimBranch *branch = dynamic_cast<DSimBranch*>(branches[i].get());
   DSimBus *busf = dynamic_cast<DSimBus*>((branch->getBus1()).get());
   DSimBus *bust = dynamic_cast<DSimBus*>((branch->getBus2()).get());
   
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
 
 DSMode dsmode;
 for(i=0; i < p_ngen; i++) {
   if(p_gen[i]->getGenStatus()) {
     p_gen[i]->setMode(p_mode);
     p_gen[i]->setVoltage(p_VDQptr[0],p_VDQptr[1]);
     p_gen[i]->setTSshift(p_TSshift);
     //gridpack::ComplexType genval[20];
     gridpack::ComplexType genval[40]; //SJin: better make this dimension an input value to accomodate various size systems! 
     //int    nval=0,row[20],col[20];
     int    nval=0,row[40],col[40]; //SJin: better make this dimension an input value to accomodate various size systems!
     // Derivatives of gen. eqs. w.r.t. state variables
     p_gen[i]->matrixDiagEntries(&nval,row,col,genval);
     
     int var_col_start,eq_idx;
     // Insert the entries in the values array
     for(int k=0; k < nval; k++) {
       var_col_start = (2+ctr+col[k])*p_nvar;
       eq_idx        = 2 + ctr + row[k];
       values[var_col_start + eq_idx] = genval[k];
     }

     // Derivatives of generator currents w.r.t. VD, VQ
     p_gen[i]->setMode(DIG_DV);
     nval = 0;
     p_gen[i]->matrixDiagEntries(&nval,row,col,genval);
     // Insert the entries in the values array
     for(int k=0; k < nval; k++) {
       var_col_start = col[k]*p_nvar;
       eq_idx        = row[k];
       values[var_col_start + eq_idx] += genval[k]; // SJin: TBD
     }

     if(p_mode != FAULT_EVAL) {
       // Derivatives of generator equations w.r.t. VD, VQ
       p_gen[i]->setMode(DFG_DV);
       nval = 0;
       p_gen[i]->matrixDiagEntries(&nval,row,col,genval);
       // Insert the entries in the values array
       for(int k=0; k < nval; k++) {
	 var_col_start = col[k]*p_nvar;
	 eq_idx        = 2 + ctr + row[k];
	 values[var_col_start + eq_idx] = genval[k]; // SJin: TBD
       }

       // Derivatives of generator currents w.r.t. generator variables
       p_gen[i]->setMode(DIG_DX);
       nval = 0;
       p_gen[i]->matrixDiagEntries(&nval,row,col,genval);
       // Insert the entries in the values array
       for(int k=0; k < nval; k++) {
	 var_col_start = (2+ctr+col[k])*p_nvar;
	 eq_idx        = row[k];
	 values[var_col_start + eq_idx] = genval[k];
       }
     }
   }
   ctr += p_neqsgen[i];
 }  

 return true;
}

/**
 * Return size of vector block contributed by component
 * @param isize: number of vector elements
 * @return: false if network component does not contribute
 *        vector element
 */
bool DSimBus::vectorSize(int *isize) const
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
bool DSimBus::vectorValues(gridpack::ComplexType *values)
{
  gridpack::ComplexType *fbus = values;  // First two locations are for the bus current balance equations.
  gridpack::ComplexType *fgen = values + 2;
  if(p_mode == INIT_X) { /* Values go in X */
    //printf("========================================================\n");
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
	if(p_gen[i]->getGenStatus()) {
	  p_gen[i]->setVoltage(p_VDQptr[0],p_VDQptr[1]);
	  p_gen[i]->init(fgen);
	}
	// SJin: add exciter and governor vector values
	if (has_ex) {
  	  fgen += p_neqsgen[i];
	  p_gen[i]->getExciter()->setVoltage(p_VDQptr[0],p_VDQptr[1]);
    	  p_gen[i]->getExciter()->init(fgen);
  	}	
	if (has_gv) {
 	  if (has_ex) fgen += p_neqsexc[i];
 	  else fgen += p_neqsgen[i];
	  p_gen[i]->getGovernor()->init(fgen);
  	}	
      }
    }	  
  } else if(p_mode == RESIDUAL_EVAL || p_mode == FAULT_EVAL) { /* Values go in F */
    //printf("//////////////////////////////////////////////////////\n");
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
      DSimBranch *branch = dynamic_cast<DSimBranch*>(branches[i].get());
      DSimBus *busf = dynamic_cast<DSimBus*>((branch->getBus1()).get());
      DSimBus *bust = dynamic_cast<DSimBus*>((branch->getBus2()).get());
      
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
      if(p_gen[i]->getGenStatus()) {
	p_gen[i]->setMode(p_mode);
	p_gen[i]->setVoltage(p_VDQptr[0],p_VDQptr[1]);
	p_gen[i]->getCurrent(&IGD,&IGQ);
	p_gen[i]->vectorValues(fgen);
        fgen += p_neqsgen[i];
        //printf("gen:\n");
        //for (int j=2; j<p_neqsgen[i]+2;j++)
        //  printf("values[%d]=%f\n",j,values[j]);

	// SJin: add exciter and governor vector values
	if (has_ex) {
	  p_gen[i]->getExciter()->setMode(p_mode);
	  p_gen[i]->getExciter()->setVoltage(p_VDQptr[0],p_VDQptr[1]);
	  p_gen[i]->getExciter()->vectorValues(fgen);
          fgen += p_neqsexc[i];
          //printf("gen+exc:\n");
          //for (int j=0; j<p_neqsgen[i]+p_neqsexc[i]+2;j++)
          //  printf("values[%d]=%f\n",j,values[j]);
  	}	
	if (has_gv) {
	  p_gen[i]->getGovernor()->setMode(p_mode);
	  p_gen[i]->getGovernor()->vectorValues(fgen);
          fgen += p_neqsgov[i];
          //printf("gen+exc+gov:\n");
      	  //for (int j=0; j<p_neqsgen[i]+p_neqsexc[i]+p_neqsgov[i]+2;j++)
          //  printf("values[%d]=%f\n",j,values[j]);
  	}	
	
	IgenD += IGD;
	IgenQ += IGQ;
      }

      //printf("\nvectorValues:\n");
      //printf("ngen=%d,nexc=%d,ngov=%d\n",p_neqsgen[i],p_neqsexc[i],p_neqsgov[i]);
      //for (int j=0; j<p_neqsgen[i]+p_neqsexc[i]+p_neqsgov[i]+2;j++)
      //  printf("values[%d]=%f\n",j,values[j]);
      /*printf("vd,vq:\n");
      for (int j=0; j<2;j++)
        printf("values[%d]=%f\n",j,values[j]);
      printf("gen:\n");
      for (int j=2; j<p_neqsgen[i]+2;j++)
        printf("values[%d]=%f\n",j,values[j]);
      printf("exc:\n");
      for (int j=p_neqsgen[i]+2;j<p_neqsgen[i]+p_neqsexc[i]+2;j++)
        printf("values[%d]=%f\n",j,values[j]);
      printf("gov:\n");
      for (int j=p_neqsgen[i]+p_neqsexc[i]+2; j<p_neqsgen[i]+p_neqsexc[i]+p_neqsgov[i]+2;j++)
        printf("values[%d]=%f\n",j,values[j]);*/
    }

    fbus[VD_idx] = IgenQ - IshuntQ - IbrQ - IloadQ;
    fbus[VQ_idx] = IgenD - IshuntD - IbrD - IloadD;

  }
  return true;
}

/**
 * Set the model to control what matrices and vectors and built when using the
 * mapper
 * @param mode: enumerated constant for different modes
 */
void DSimBus::setMode(int mode)
{
  p_mode = (DSMode)mode;
}

int DSimBus::getXCBufSize(void)
{
  return 2*sizeof(double);
}

void DSimBus::setXCBuf(void *buf)
{
  p_VDQptr = static_cast<double*>(buf);
  
}

/**
 * Set the internal values of the voltage magnitude and phase angle. Need this
 * function to push values from vectors back onto buses 
 * @param values array containing voltage magnitude and angle
 */
void DSimBus::setValues(gridpack::ComplexType *values)
{
  //printf("........................................\n");
  int i,ctr=2;
  gridpack::ComplexType *genvals;
  DSMode mode;

  if(p_isolated) return;

  mode = p_mode;

  if(p_mode == XVECTOBUS) { // Push values from X vector back on the bus
    *p_VDQptr = real(values[0]);
    *(p_VDQptr+1) = real(values[1]);
  }
  
  genvals = values + 2;
  for(i=0; i < p_ngen; i++) {
    if(p_gen[i]->getGenStatus()) {
      p_gen[i]->setMode(mode);
      p_gen[i]->setVoltage(p_VDQptr[0],p_VDQptr[1]);
      p_gen[i]->setValues(genvals);
      genvals += p_neqsgen[i];

      // SJin: add exciter and governor vector values
      if (has_ex) {
        p_gen[i]->getExciter()->setMode(p_mode);
        p_gen[i]->getExciter()->setValues(genvals);
        genvals += p_neqsexc[i];
      }	
      if (has_gv) {
        p_gen[i]->getGovernor()->setMode(p_mode);
        p_gen[i]->getGovernor()->setValues(genvals);
        genvals += p_neqsgov[i];
      }	
    }	
 
      /*printf("setValues:\n");
      printf("ngen=%d,nexc=%d,ngov=%d\n",p_neqsgen[i],p_neqsexc[i],p_neqsgov[i]);
      for (int j=0; j<p_neqsgen[i]+p_neqsexc[i]+p_neqsgov[i]+2;j++)
        printf("values[%d]=%f\n",j,values[j]);*/
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
bool DSimBus::serialWrite(char *string,
    const int bufsize, const char *signal)
{
}

/**
 *  Simple constructor
 */
DSimBranch::DSimBranch(void)
{
  p_nparlines = 0;
  p_mode = NONE;
}

/**
 *  Simple destructor
 */
DSimBranch::~DSimBranch(void)
{
}

/**
 * Load values stored in DataCollection object into DSimBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 *
 * Note: GridPACK stores the data for all parallel branches in the same
 *       DataCollection object. 
 */
void DSimBranch::load(
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
void DSimBranch::setMode(int mode)
{
  p_mode = mode;
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool DSimBranch::matrixForwardSize(int *isize, int *jsize) const
{
  // Get from and to buses
  DSimBus *busf = dynamic_cast<DSimBus*>(getBus1().get());
  DSimBus *bust = dynamic_cast<DSimBus*>(getBus2().get());
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

bool DSimBranch::matrixReverseSize(int *isize, int *jsize) const
{
  // Get from and to buses
  DSimBus *busf = dynamic_cast<DSimBus*>(getBus1().get());
  DSimBus *bust = dynamic_cast<DSimBus*>(getBus2().get());
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
bool DSimBranch::matrixForwardValues(gridpack::ComplexType *values)
{
  int i;
  // Get from and to buses
  DSimBus *busf = dynamic_cast<DSimBus*>(getBus1().get());
  DSimBus *bust = dynamic_cast<DSimBus*>(getBus2().get());

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

bool DSimBranch::matrixReverseValues(gridpack::ComplexType *values)
{
  int i;
  // Get from and to buses
  DSimBus *busf = dynamic_cast<DSimBus*>(getBus1().get());
  DSimBus *bust = dynamic_cast<DSimBus*>(getBus2().get());
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
bool DSimBranch::serialWrite(char *string, const int
    bufsize,  const char *signal)
{
  return false;
}

bool DSimBranch::getForwardSelfAdmittance(double *Gff,double *Bff)
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

bool DSimBranch::getReverseSelfAdmittance(double *Gtt,double *Btt)
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

bool DSimBranch::getForwardTransferAdmittance(double *Gft,double *Bft)
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

bool DSimBranch::getReverseTransferAdmittance(double *Gtf,double *Btf)
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
