/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
#include <vector>
#include <iostream>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "dsf_components.hpp"
#include "lvshbl.hpp"
#include "gridpack/utilities/string_utils.hpp"


/**
 *  Simple constructor
 */
gridpack::dynamic_simulation::DSFullBus::DSFullBus(void)
{
  p_shunt_gs = 0.0;
  p_shunt_bs = 0.0;
  p_busfault = false;
  p_gfault = 0.0;
  p_bfault = 0.0;
  p_mode = YBUS;
  setReferenceBus(false);
  p_ngen = 0;
  p_negngen = 0;
  p_ngen_nodynmodel = 0;
  p_from_flag = false;
  p_Yload_change_P_flag = false;
  p_Yload_change_Q_flag = false;
  p_bConstYLoadSettoZero_P = false;
  p_bConstYLoadSettoZero_Q = false;
  p_bConstYLoadSettoValue = false;
  
  p_bconstYLoadSheddingFlag = false;
  remainConstYLoadPerc = 1.0;
  p_bscatterinjload_flag = false;
  p_bscatterinjload_flag_compensateY = false;
  p_bscatterinjloadconstcur_flag = false;
  p_to_flag = false;
  p_branchrelay_from_flag = false; 
  p_branchrelay_to_flag = false;
  p_branchtripaction_from_flag = false;
  p_branchtripaction_to_flag= false;
  p_busrelaytripflag = false; 
  p_branch = NULL;
  p_isolated = false;
  p_busvolfreq = 60.0; //renke add
  pbusvolfreq_old = 60.0; //renke add
  bcomputefreq = false;  //renke add
  p_loadimpedancer = 0.0;
  p_loadimpedancei = 0.0;
  p_scatterinjload_p = 0.0;
  p_scatterinjload_q = 0.0;
  p_scatterinjload_constcur_r = 0.0;
  p_scatterinjload_constcur_i = 0.0;
  p_ndyn_load = 0;
  p_npowerflow_load = 0;
  p_bextendedloadbus = -1;
  loadMVABase = 100.0;
  p_CmplXfmrBranch = NULL; 
  p_CmplFeederBranch = NULL;  
  p_CmplXfmrBus = NULL;
  p_CmplFeederBus = NULL; 
  p_CmplXfmr_xxf = 0.01;
  p_CmplXfmr_tap = 0.0; 
  p_pl = 0.0;
  p_ql = 0.0;
  p_relaytrippedbranch = NULL;
  p_Yload_change_r = 0.0;
  p_Yload_change_i = 0.0;

  p_line_status_change = false;
  p_gen_status_change = false;
  //p_tripactionbranch.clear();
}

/**
 *  Simple destructor
 */
gridpack::dynamic_simulation::DSFullBus::~DSFullBus(void)
{
}

  /**
     * set fault
     * @param flag => true = fault on, false otherwise
     * @param gfault => fault conductance (pu)
     * @param bfault => fault susceptance (pu)
     */
void gridpack::dynamic_simulation::DSFullBus::setFault(double gfault, double bfault)
{
  p_busfault = true;
  p_gfault = gfault;
  p_bfault = bfault;
}

/**
 * clear fault
 */
void gridpack::dynamic_simulation::DSFullBus::clearFault()
{
  p_busfault = false;
}

/**
   getGenStatus - Get the generator status 
   
   @param: ckt_id - generator id
   @return: status - generator status
   
**/
int gridpack::dynamic_simulation::DSFullBus::getGenStatus(std::string id)
{
  int gstatus=0;
  for(int i=0; i < p_ngen; i++) {
    if(p_genid[i] == id) {
      gstatus = p_gstatus[i];
      break;
    }
  }
  return gstatus;
}

/**
   getGenNum - Get the generator number given generator id
   
   @param: gen_id - generator id
   @return:  idx - the internal index for the generator that can be used to access its parameters (-1 if not found)
**/
int gridpack::dynamic_simulation::DSFullBus::getGenNum(std::string id)
{
  int gennum=-1;
  for(int i=0; i < p_ngen; i++) {
    if(p_genid[i] == id) {
      gennum = i;
      break;
    }
  }
  return gennum;
}


/**
   setGenStatus - Sets the generator status 
   
   @param: ckt_id - generator id
   @param: status - new generator status
   
   Note: This method is used to
   update the generator status and update the bus/branch
   objects. It sets up values in the bus objects
   so that incrementMatrix method called on the network Ybus
   uses these values to remove the generator admittance contributions from
   the Y-bus matrix
**/
void gridpack::dynamic_simulation::DSFullBus::setGenStatus(std::string id, int status)
{
  std::string ckt_id;
  gridpack::utility::StringUtils util;
  ckt_id = util.clean2Char(id);

  int gen_i = getGenNum(ckt_id);
  int gen_status = p_generators[gen_i]->getGenStatus();

  if(status != gen_status) {
    p_ygen = p_generators[gen_i]->NortonImpedence();
    if(!status) p_ygen = -p_ygen; // Negative sign here for removing the generator admittance in the Ybus
    p_gstatus[gen_i] = status;
    p_generators[gen_i]->SetGenServiceStatus((bool)status);
    p_gen_status_change = true;
  }
}

/**
 *  Return size of matrix block contributed by the component
 *  @param isize, jsize: number of rows and columns of matrix block
 *  @return: false if network component does not contribute matrix element
 */
bool gridpack::dynamic_simulation::DSFullBus::matrixDiagSize(int *isize, int *jsize) const
{
  if (YMBus::isIsolated()) return false;
  if (p_mode == YBUS || p_mode == YL || p_mode == PG || p_mode == YDYNLOAD || p_mode == onFY || p_mode == posFY
  || p_mode == jxd || p_mode == bus_relay || p_mode == branch_relay || p_mode == branch_trip_action
  || p_mode == bus_Yload_change_P || p_mode == bus_Yload_change_Q) {
    return YMBus::matrixDiagSize(isize,jsize);
  }  else {
    *isize = 1;
    *jsize = 1;
  }
  return true;
}

/**
 * Return the values of the matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::dynamic_simulation::DSFullBus::matrixDiagValues(ComplexType *values)
{
  if (p_mode == YBUS) {
    bool status =  YMBus::matrixDiagValues(values);
    p_ybusr = values->real();
    p_ybusi = values->imag();
    return status;
  } else if (p_mode == YL) {
    //printf("DSFullBus::matrixDiagValues, bus %d: p_pl = %f, p_ql = %f, p_voltage = %f\n", getOriginalIndex(), p_pl, p_ql, p_voltage);
    //printf("p_ybusr = %f, p_ybusi = %f\n", p_ybusr, p_ybusi);
    p_ybusr = p_ybusr+p_pl/(p_voltage*p_voltage);
    p_ybusi = p_ybusi+(-p_ql)/(p_voltage*p_voltage);
    p_loadimpedancer = p_pl/(p_voltage*p_voltage);
    p_loadimpedancei = (-p_ql)/(p_voltage*p_voltage);
    /* TBD: p_ybusr = p_ybusr+(-p_pg)/(p_voltage*p_voltage);
    p_ybusi = p_ybusi+p_qg/(p_voltage*p_voltage);*/
    //printf(".. p_ybusr = %f, p_ybusi = %f\n", p_ybusr, p_ybusi);
    gridpack::ComplexType ret(p_ybusr, p_ybusi);
    values[0] = ret;
    return true;
  } else if (p_mode == PG) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
	if(!p_gstatus[i]) continue;
	
	if (p_pg[i] < 0) {
           p_ybusr = p_ybusr+(-p_pg[i])/(p_voltage*p_voltage);
           p_ybusi = p_ybusi+p_qg[i]/(p_voltage*p_voltage);
           gridpack::ComplexType ret(p_ybusr, p_ybusi);
           values[0] = ret;
         } else {
           gridpack::ComplexType u(p_ybusr, p_ybusi);
           values[0] = u;
         }
      }
    } else {
      gridpack::ComplexType u(p_ybusr, p_ybusi);
      values[0] = u;
    }
    if (p_negngen > 0) {
      for (int i = 0; i < p_negngen; i++) {
        p_ybusr = p_ybusr+(-p_negpg[i])/(p_voltage*p_voltage);
        p_ybusi = p_ybusi+p_negqg[i]/(p_voltage*p_voltage);
        gridpack::ComplexType ret(p_ybusr, p_ybusi);
        values[0] = ret;
      }
    }
    for (int i = 0; i < p_ngen; i++) {
      if(p_gen_nodynmodel[i]) {
	p_ybusr = p_ybusr+(-p_genpg_nodynmodel[i])/(p_voltage*p_voltage);
	p_ybusi = p_ybusi+p_genqg_nodynmodel[i]/(p_voltage*p_voltage);
	gridpack::ComplexType ret(p_ybusr, p_ybusi);
	values[0] = ret;
      }
    }
    return true;
  } else if (p_mode == jxd) {
    if (p_ngen > 0) {
      for (int i = 0; i < p_ngen; i++) {
	if(!p_gstatus[i] || !p_generators.size()) continue;

        gridpack::ComplexType Y_a
          = p_generators[i]->NortonImpedence();

        p_ybusr = p_ybusr + real(Y_a);
        p_ybusi = p_ybusi + imag(Y_a);
        gridpack::ComplexType ret(p_ybusr, p_ybusi);

        values[0] = ret;
      }

    } else {
      gridpack::ComplexType u(p_ybusr, p_ybusi);
      values[0] = u;
      //return true;
    }
    //printf("idx: %d Real: %f Imag: %f\n",getOriginalIndex(),
      //real(values[0]), imag(values[0]));
    return true;
  } else if(p_mode == BUSFAULTON) {
    // Note this should be only used with
    // incrementMatrix method
    if(p_busfault) {
      /* Fault on */
      // Values to be added to Ybus
      // Note negative values since the current
      // flows out of the bus
      gridpack::ComplexType ret(-p_gfault,-p_bfault);
      values[0] = ret;
      return true;
    } else {
      return false;
    }
  } else if (p_mode == BUSFAULTOFF) {
    if(p_busfault) {
      /* Fault off */
      p_ybusr += p_gfault;
      p_ybusi += p_bfault;
      // Values to be added to Ybus
      gridpack::ComplexType ret(p_gfault,p_bfault);
      values[0] = ret;
      clearFault();
      return true;
    } else {
      return false;
    }
  } else if(p_mode == LINESTATUSCHANGE) {
    if(p_line_status_change) {
      values[0] = p_yii;
      p_line_status_change = false; // Clear line status change flag
      return true;
    } else return false;
  } else if(p_mode == GENSTATUSCHANGE) {
    if(p_gen_status_change) {
      values[0] = p_ygen;
      p_gen_status_change = false;
      return true;
    } else return false;
  } else if (p_mode == onFY) {
    if (p_from_flag) {
      //gridpack::ComplexType ret(0.0, -1.0e9);
	  double tmp1 = p_ybusr;
      double tmp2 = p_ybusi - 1.0e9;
      gridpack::ComplexType ret(tmp1, tmp2);
	  
      values[0] = ret;
      return true;
    } else {
      return false;
    }
  } else if (p_mode == bus_Yload_change_P) {
    if (p_Yload_change_P_flag) {
      //gridpack::ComplexType ret(0.0, -1.0e5);
	  double tmp1 = p_Yload_change_r/(p_voltage*p_voltage); //p_ybusr + p_Yload_change_r;
      double tmp2 = 0.0; //p_ybusi;
      gridpack::ComplexType ret(tmp1, tmp2);
	  
      values[0] = ret;
	  
	  //printf("-----Ybus change at bus idx due to load change: %d Real: %f Imag: %f\n",getOriginalIndex(), real(values[0]), imag(values[0]) );
	  
      return true;
    } else {
      return false;
    }
  } else if (p_mode == bus_Yload_change_Q) {
    if (p_Yload_change_Q_flag) {
      //gridpack::ComplexType ret(0.0, -1.0e5);
	  double tmp1 = 0.0; //p_ybusr;
      double tmp2 = - p_Yload_change_i/(p_voltage*p_voltage); //_ybusi - p_Yload_change_i;
      gridpack::ComplexType ret(tmp1, tmp2);
	  
      values[0] = ret;
      return true;
    } else {
      return false;
    }
  } else if (p_mode == posFY) { // renke modify, only for line failure with tripping at post fault stage, the bus diagnal element will change value - reducing the value of the Y of the tripped line
    if ( (p_from_flag || p_to_flag) && (p_branch!=NULL) ) {
      values[0] = dynamic_cast<gridpack::dynamic_simulation::DSFullBranch*>(p_branch)->getUpdateFactor();
      return true;
    } else {
      return false;
    }
  } else if (p_mode == bus_relay) {
      if (p_busrelaytripflag) {
      gridpack::ComplexType u(p_ybusr, p_ybusi);
      values[0] = u;
      return true;
    } else {
      return false;
    }	  
  } else if (p_mode == YDYNLOAD) {  // Dynamic load model's contribution to Y matrix
	//printf("bus %d entering YDYNLOAD mode: p_ndyn_load: %d \n", getOriginalIndex(), p_ndyn_load);
    if (p_ndyn_load>0) {
		for (int i = 0; i < p_ndyn_load; i++) {

			//printf("DSFullBus::matrixDiagValues, Bus %d here 1\n", getOriginalIndex());
			gridpack::ComplexType Y_a
				= p_loadmodels[i]->NortonImpedence();
			//printf("DSFullBus::matrixDiagValues, here 2 real(Y_a): %f imag(Y_a): %f\n",real(Y_a),imag(Y_a));

			p_ybusr = p_ybusr + real(Y_a);
			p_ybusi = p_ybusi + imag(Y_a);
			gridpack::ComplexType ret(p_ybusr, p_ybusi);
			//printf("DSFullBus::matrixDiagValues, here 3 p_ybusr: %f p_ybusi: %f\n",p_ybusr,p_ybusi);
			values[0] = ret;
		}
	}else {
      gridpack::ComplexType u(p_ybusr, p_ybusi);
      values[0] = u;
      //return true;
    }	
    return true;
  } else if (p_mode == branch_relay) {
		if (p_branchrelay_from_flag || p_branchrelay_to_flag) {
			printf("Bus %d diag element change due to branch relay trip: \n", getOriginalIndex());
			if ( p_relaytrippedbranch == NULL ){
				printf("Error: in DSFullBus::matrixDiagValues(), the branch relay does not set the branch bus correctly!\n");
				return false;
			}else{
				values[0] = dynamic_cast<gridpack::dynamic_simulation::DSFullBranch*>(p_relaytrippedbranch)->getBranchRelayTripUpdateFactor();
				printf("changed value: %f + j*%f\n", real(values[0]), imag(values[0]));
				return true;
			}
		} else {
			return false;
		}	  
  }	else if (p_mode == branch_trip_action) {
		if (p_branchtripaction_from_flag || p_branchtripaction_to_flag) {
			//printf("Bus %d diag element change due to branch relay trip: \n", getOriginalIndex());
			if ( p_vec_tripactionbranch.empty() ){
				printf("Error: in DSFullBus::matrixDiagValues(), the branch does not have the branch ckts need to trip\n");
				return false;
			}else{
				gridpack::ComplexType tmp (0.0, 0.0);
				int nbr, ibr;
				nbr = p_vec_tripactionbranch.size();
				for (ibr=0; ibr<nbr; ibr++){
					tmp += dynamic_cast<gridpack::dynamic_simulation::DSFullBranch*>(p_vec_tripactionbranch[ibr])->getBranchTripActionUpdateFactorForBus(getOriginalIndex());
				}
				values[0] = tmp;
				return true;
			} 
		}else{
			return false;
		}
  
  }	// end for  else if (p_mode == branch_trip_action)
	  
}

/**
 * Return the size of the block that this component contributes to the
 * vector
 * @param size: size of vector block
 * @return: false if component does not contribute to vector
 */
bool gridpack::dynamic_simulation::DSFullBus::vectorSize(int *size) const
{
  if (!p_isolated) {
    if (p_mode == make_INorton_full) {
      *size = 1;
    }else {
      *size = 2;
    }
    return true;
  } else {
    return false;
  }
}

/**
 * Return the values of the vector block
 * @param values: pointer to vector values
 * @return: false if network component does not contribute
 *        vector element
 */
bool gridpack::dynamic_simulation::DSFullBus::vectorValues(ComplexType *values)
{
  if (!p_isolated) {
    if (p_mode == make_INorton_full) {
      values[0] = 0;
      if (p_ngen > 0) {
        for (int i = 0; i < p_ngen; i++) {
	  if(!p_gstatus[i] || !p_generators.size()) continue;
	  values[0] += p_generators[i]->INorton(); 
        } // generator for loop
      }  // if p_ngen>0
	  
      //INorton contribution from dynamic load
      for (int i =0; i<p_ndyn_load; i++){
	values[0] += p_loadmodels[i]->INorton(); 
	
	ComplexType tmp = p_loadmodels[i]->INorton();
	//printf("DSFullBus::vectorValues, bus %d, dynamic load Inorton: %12.6f + j %12.6f \n", getOriginalIndex(),real(tmp), imag(tmp));
      }
      
	  // INorton contribution from the sheded constant Y load
	  if (p_bconstYLoadSheddingFlag){
		  
		  double inj_constY_P, inj_constY_Q;
		  gridpack::ComplexType inj_constY_S, inj_constY_cur_tmp;
		  
		  inj_constY_P = p_loadimpedancer*(p_voltage*p_voltage)*(1.0-remainConstYLoadPerc);
		  inj_constY_Q = -p_loadimpedancei*(p_voltage*p_voltage)*(1.0-remainConstYLoadPerc);
		  inj_constY_S = gridpack::ComplexType(inj_constY_P, inj_constY_Q);
		  inj_constY_cur_tmp = inj_constY_S/p_volt_full;
		  inj_constY_cur_tmp = conj(inj_constY_cur_tmp);
		  
		  values[0] += inj_constY_cur_tmp;
		  
	  }
	  	  
	  // INorton contribution from the scattered injection load
	  if (p_bscatterinjload_flag){
		  
		  double p_pl_inj_tmp, p_ql_inj_tmp;
		  gridpack::ComplexType bus_inj_S, tmp, bus_scatterload_inj_cur;
		  
		  //p_pl_inj_tmp = p_pl - p_scatterinjload_p;
		  //p_ql_inj_tmp = p_ql - p_scatterinjload_q;
		  
		  p_pl_inj_tmp = 0.0 - p_scatterinjload_p;
		  p_ql_inj_tmp = 0.0 - p_scatterinjload_q;
		  
		  bus_inj_S = gridpack::ComplexType(p_pl_inj_tmp,p_ql_inj_tmp);
		  tmp = bus_inj_S/p_volt_full;
		  bus_scatterload_inj_cur = conj(tmp);
      
		  values[0] += bus_scatterload_inj_cur;
		  /*
		  printf("rkdebugInjectionCur, DSFullBus::vectorValues, bus %d, Inorton: %12.6f + j %12.6f \n", 
					getOriginalIndex(),real(bus_scatterload_inj_cur), imag(bus_scatterload_inj_cur));
					
		  printf("rkdebugpq, DSFullBus::vectorValues, bus %d, PQ: %12.6f + j %12.6f \n", 
					getOriginalIndex(), p_pl_inj_tmp, p_ql_inj_tmp);
					
		  printf("rkdebugpq, DSFullBus::vectorValues, bus %d, voltage: %12.6f + j %12.6f \n", 
					getOriginalIndex(), real(p_volt_full), imag(p_volt_full));
		  */
  
	  }
	  
	  // INorton contribution from the scattered injection load, keep the 
	  // Y load component still at the bus, while only compenstate the difference
	  if (p_bscatterinjload_flag_compensateY){
		  
		  double p_pl_inj_tmp, p_ql_inj_tmp, inj_constY_P, inj_constY_Q;
		  gridpack::ComplexType bus_inj_S, tmp, bus_scatterload_inj_cur;
		  
		  //p_pl_inj_tmp = p_pl - p_scatterinjload_p;
		  //p_ql_inj_tmp = p_ql - p_scatterinjload_q;
		  double voltagemagtmp = abs(p_volt_full);
		  inj_constY_P = p_loadimpedancer*(voltagemagtmp*voltagemagtmp);
		  inj_constY_Q = -p_loadimpedancei*(voltagemagtmp*voltagemagtmp);
		  
		  p_pl_inj_tmp = inj_constY_P - p_scatterinjload_p;
		  p_ql_inj_tmp = inj_constY_Q - p_scatterinjload_q;
		  
		  bus_inj_S = gridpack::ComplexType(p_pl_inj_tmp,p_ql_inj_tmp);
		  tmp = bus_inj_S/p_volt_full;
		  bus_scatterload_inj_cur = conj(tmp);
      
		  values[0] += bus_scatterload_inj_cur;
		  /*
		  printf("rkdebugInjectionCur, DSFullBus::vectorValues, bus %d, Inorton: %12.6f + j %12.6f \n", 
					getOriginalIndex(),real(bus_scatterload_inj_cur), imag(bus_scatterload_inj_cur));
					
		  printf("rkdebugpq, DSFullBus::vectorValues, bus %d, PQ: %12.6f + j %12.6f \n", 
					getOriginalIndex(), p_pl_inj_tmp, p_ql_inj_tmp);
					
		  printf("rkdebugpq, DSFullBus::vectorValues, bus %d, voltage: %12.6f + j %12.6f \n", 
					getOriginalIndex(), real(p_volt_full), imag(p_volt_full));
		  */
  
	  }
	  
	  	  // INorton contribution from the scattered injection load as constant current load
	  if (p_bscatterinjloadconstcur_flag){
		  
		  gridpack::ComplexType tmp, bus_scatterload_inj_cur;
		  
		  bus_scatterload_inj_cur = gridpack::ComplexType(p_scatterinjload_constcur_r, p_scatterinjload_constcur_i);
		  
		  bus_scatterload_inj_cur = bus_scatterload_inj_cur/p_volt_full*abs(p_volt_full); //just rotate the voltage angle, keep the current magnitude unchange
		  bus_scatterload_inj_cur = conj(bus_scatterload_inj_cur);
      
		  values[0] -= bus_scatterload_inj_cur;
		  
		/*  
		  printf("rkdebugInjectionCur, DSFullBus::vectorValues, bus %d, Inorton: %12.6f + j %12.6f \n", 
					getOriginalIndex(),real(bus_scatterload_inj_cur), imag(bus_scatterload_inj_cur));
					
		  printf("rkdebugpq, DSFullBus::vectorValues, bus %d, voltage: %12.6f + j %12.6f \n", 
					getOriginalIndex(), real(p_volt_full), imag(p_volt_full));
					
		*/
		  
  
	  }
	  
      //printf("bus id = %d, values[0] = %f, %f\n", getOriginalIndex(), real(values[0]), imag(values[0])); 
      return true;
    } else {
      return false;
    }  // if mode == make Inorton
  }  // if !p_isolated
  return false;
}

void gridpack::dynamic_simulation::DSFullBus::setValues(ComplexType *values)
{
  int i;
  if (p_mode == make_INorton_full) {
	p_volt_full_old = p_volt_full; //renke add  
    p_volt_full = values[0];
  }
}

/**
 * Set values of YBus matrix. These can then be used in subsequent
 * calculations
 */
void gridpack::dynamic_simulation::DSFullBus::setYBus(void)
{
  YMBus::setYBus();
  gridpack::ComplexType ret;
  ret = YMBus::getYBus();
  p_ybusr = real(ret);
  p_ybusi = imag(ret);
}

/**
 * Add shunt
 * @param gs shunt Gshunt value
 * @param bs shunt Bshunt value
 */
void gridpack::dynamic_simulation::DSFullBus::addShunt(double gs, double bs)
{
  YMBus::addShuntValues(gs,bs);
}


void gridpack::dynamic_simulation::DSFullBus::setGeneratorObPowerBaseFlag(bool generator_observationpower_systembase)
{
  if (p_ngen > 0) {
    for (int i = 0; i < p_ngen; i++) {
      if(!p_gstatus[i] || !p_generators.size()) continue;
      p_generators[i]->setGeneratorObPowerBaseFlag(generator_observationpower_systembase);
    }
  }
}

/**
 * Set initial values of vectors for integration. 
 * These can then be used in subsequent calculations
 * @param ts time step
 */
void gridpack::dynamic_simulation::DSFullBus::initDSVect(double ts)
{
  if (p_ngen > 0) {
    for (int i = 0; i < p_ngen; i++) {
      if(!p_gstatus[i] || !p_generators.size()) continue;
      p_generators[i]->init(p_voltage,p_angle, ts);
    }
  }
  //dynamic loads
  for (int i = 0; i < p_ndyn_load; i++) {
    p_loadmodels[i]->init(p_voltage,p_angle, ts);
  }

}

/**
 * Update values for vectors in each integration time step (Predictor)
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::DSFullBus::predictor_currentInjection(bool flag)
{
  if (p_ngen == 0 && p_ndyn_load == 0) return;
  int i;
  for (i = 0; i < p_ngen; i++) {
    if(!p_gstatus[i] || !p_generators.size()) continue;
    p_generators[i]->predictor_currentInjection(flag);
  }
  
  //dynamic loads
  //printf("DSFullBus::predictor_currentInjection, bus %d, p_ndyn_load: %d", getOriginalIndex(), p_ndyn_load);
  for (i = 0; i < p_ndyn_load; i++) {

	//printf("DSFullBus::predictor_currentInjection, bus %d, entering dynamic load computation loop", getOriginalIndex());
    p_loadmodels[i]->predictor_currentInjection(flag);
  }
  
}

/**
 * Update values for vectors in each integration time step (Predictor)
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::DSFullBus::predictor(double t_inc, bool flag)
{
  if (p_ngen == 0 && p_ndyn_load == 0) return;

  int i;
  for (i = 0; i < p_ngen; i++) {
    if(!p_gstatus[i] || !p_generators.size()) continue;
    p_generators[i]->predictor(t_inc,flag);
  }
  
    //dynamic loads
  for (i = 0; i < p_ndyn_load; i++) {
    p_loadmodels[i]->predictor(t_inc,flag);
  }
}

/**
 * Update values for vectors in each integration time step (Corrector)
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::DSFullBus::corrector_currentInjection(bool flag)
{
  if (p_ngen == 0 && p_ndyn_load == 0) return;
  int i;
  for (i = 0; i < p_ngen; i++) {
    if(!p_gstatus[i] || !p_generators.size()) continue;
    p_generators[i]->corrector_currentInjection(flag);
  }
  
  //dynamic loads
  for (i = 0; i < p_ndyn_load; i++) {
	//if (!p_generators[i]->getGenStatus()) {
	//	continue;
	//}
    p_loadmodels[i]->corrector_currentInjection(flag);
  }
}

/**
 * Update values for vectors in each integration time step (Corrector)
 * @param t_inc time step increment
 * @param flag initial step if true
 */
void gridpack::dynamic_simulation::DSFullBus::corrector(double t_inc, bool flag)
{
  if (p_ngen == 0 && p_ndyn_load == 0) return;

  int i;
  for (i = 0; i < p_ngen; i++) {
    if(!p_gstatus[i] || !p_generators.size()) continue;
    p_generators[i]->corrector(t_inc,flag);
  }
  
  //dynamic loads
  for (i = 0; i < p_ndyn_load; i++) {
    p_loadmodels[i]->corrector(t_inc,flag);
  }
}

void gridpack::dynamic_simulation::DSFullBus::setWideAreaFreqforPSS(double freq){
	
  int i;
  for (i = 0; i < p_ngen; i++) {
    if(!p_gstatus[i] || !p_generators.size()) continue;
    p_generators[i]->setWideAreaFreqforPSS(freq);
  }
	
}


/**
 * Update dynamic load internal relays action
 */
void gridpack::dynamic_simulation::DSFullBus::dynamicload_post_process(double t_inc, bool flag)
{
	int i;
	
	//dynamic loads
  for (i = 0; i < p_ndyn_load; i++) {
	//if (!p_generators[i]->getGenStatus()) {
	//	continue;
	//}
    p_loadmodels[i]->dynamicload_post_process(t_inc,flag);
  }
  
}

/**
 * Get roter angle of generators
 */
double gridpack::dynamic_simulation::DSFullBus::getAngle()
{
  if (p_ngen < 0) return 0.0;
  int i;
  for (i = 0; i < p_ngen; i++) {
    if(!p_gstatus[i] || !p_generators.size()) continue;
    double angle = p_generators[i]->getAngle();
    return angle;
  }
}

/**
 * Set volt from volt_full
 */
void gridpack::dynamic_simulation::DSFullBus::setVolt(bool flag) 
{
  if (p_ngen > 0) {
    for (int i = 0; i < p_ngen; i++) {
      if(!p_gstatus[i] || !p_generators.size()) continue;
      p_generators[i]->setVoltage(p_volt_full);
    }
  }
  
  for ( int i=0; i<p_ndyn_load; i++){
	  p_loadmodels[i]->setVoltage(p_volt_full);
  }
}

/**
 * compute bus frequency and set it to dynamic load models
 */
void gridpack::dynamic_simulation::DSFullBus::updateFreq (double delta_t){
	int i;
	double dbusvoltfreq;
	if ( bcomputefreq == true ) {
		
	  computeBusVolFrequency(delta_t);
	  dbusvoltfreq = getBusVolFrequency();
	  
	  //set voltage frequency for dynamic loads
	  for (i=0; i<p_ndyn_load; i++){
	    p_loadmodels[i]->setFreq(dbusvoltfreq/60.0);
	  }
	  
	  for (int i = 0; i < p_ngen; i++) {
	    if(!p_gstatus[i] || !p_generators.size()) continue;
	    p_generators[i]->setFreq(dbusvoltfreq/60.0);
	  }
	  
	}
	
}

/**
 * Get values of YBus matrix. These can then be used in subsequent
 * calculations
 */
gridpack::ComplexType gridpack::dynamic_simulation::DSFullBus::getYBus(void)
{
  return YMBus::getYBus();
}

/**
 * Load values stored in DataCollection object into DSFullBus object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       bus that were read in when network was initialized
 */
void gridpack::dynamic_simulation::DSFullBus::load(
  const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  YMBus::load(data);
  // This function may be called more than once so clear all vectors
  p_pg.clear();
  p_qg.clear();
  p_negpg.clear();
  p_negqg.clear();
  p_gen_nodynmodel.clear();
  p_genpg_nodynmodel.clear();
  p_genqg_nodynmodel.clear();
  p_genid.clear();
  p_loadid.clear();
  p_generators.clear();
  p_loadrelays.clear();
  p_loadmodels.clear();
  p_ngen_nodynmodel = 0;

  std::string snewbustype; //renke add

  if(!data->getValue(CASE_SBASE,&p_sbase)) p_sbase = 100.0;

  // check whether the bus is the ones extended by composite load models,
  // if yes, return
  if (data->getValue(NEW_BUS_TYPE, &snewbustype)){
    if ( snewbustype=="LOW_SIDE_BUS" || snewbustype=="LOAD_BUS" ) {
      printf("This bus is a extended bus by composite load models, type: %s \n",
          snewbustype.c_str());
      if ( snewbustype=="LOW_SIDE_BUS" )
      {
        p_bextendedloadbus = 1; 
      }else{
        p_bextendedloadbus = 2; 
      }

      return; 
    }	  
  }
  if (!data->getValue("BUS_PF_VANG", &p_angle)) {
    data->getValue(BUS_VOLTAGE_ANG, &p_angle);
  }
  if (!data->getValue("BUS_PF_VMAG", &p_voltage)) {
    data->getValue(BUS_VOLTAGE_MAG, &p_voltage);
  }

  p_area = 1;
  data->getValue(BUS_AREA, &p_area);
  p_zone = 1;
  data->getValue(BUS_ZONE, &p_zone);
  double pi = 4.0*atan(1.0);
  p_angle = p_angle*pi/180.0; 

  p_shunt = true;
  bool ok = true;
  ok = data->getValue(BUS_SHUNT_GL, &p_shunt_gs);
  p_shunt = p_shunt && ok;
  ok = data->getValue(BUS_SHUNT_BL, &p_shunt_bs);
  p_shunt = p_shunt && ok;
  p_shunt_gs /= p_sbase;
  p_shunt_bs /= p_sbase; 

  // Check to see if bus is reference bus
  data->getValue(BUS_TYPE, &p_type);
  if (p_type == 3) {
    setReferenceBus(true);
  } else if (p_type == 4) {
    p_isolated = true;
  }

  bool lgen;
  int i, gstatus;
  int nrelay, irelay, relaycnt;
  std::string relay_genid;
  double pg, qg, mva, r, dstr, dtr;
  double h, d0;
  bool has_ex, has_gov, has_pss, has_plantcontroller;
  bool has_wind_aero,has_wind_dt,has_wind_tc,has_wind_pc;
  GeneratorFactory genFactory;
  RelayFactory relayFactory;
  LoadFactory loadFactory;
  p_generators.clear();
  p_negngen = 0;
  int idx;
  data->getValue(BUS_NUMBER,&idx);
  if (data->getValue(GENERATOR_NUMBER, &p_ngen)) {
    std::string genid;
    int icnt = 0;

    for (i=0; i<p_ngen; i++) {
      // set up arrays to monitor frequency
      p_previousFrequency.push_back(0.0);
      p_upIntervalStart.push_back(0.0);
      p_upStartedMonitoring.push_back(false);
      p_downIntervalStart.push_back(0.0);
      p_downStartedMonitoring.push_back(false);

      // Set up model
      int stat;
      data->getValue(GENERATOR_STAT, &stat, i);
      data->getValue(GENERATOR_PG, &pg, i);
      data->getValue(GENERATOR_QG, &qg, i);
      data->getValue(GENERATOR_ID, &genid, i);
      p_genid.push_back(genid);
      p_gstatus.push_back(stat);
      p_gen_nodynmodel.push_back(false);
      p_genpg_nodynmodel.push_back(0.0);
      p_genqg_nodynmodel.push_back(0.0);
      p_pg.push_back(pg);
      p_qg.push_back(qg);
      p_savePg.push_back(pg);
      double pmin, pmax;
      data->getValue(GENERATOR_PMIN,&pmin,i);
      data->getValue(GENERATOR_PMAX,&pmax,i);
      p_gpmin.push_back(pmin);
      p_gpmax.push_back(pmax);
      
      std::string model;
      if (data->getValue(GENERATOR_MODEL, &model, i) && pg >= 0.0) {

        BaseGeneratorModel *generator
          = genFactory.createGeneratorModel(model);
	if (model == "REGCA1" || model == "REGCB1" || model == "REGCC1" || model == "GDFORM"){
	  bcomputefreq = true;
	}
        has_ex = false;
        has_gov = false;
        has_pss = false;
	has_plantcontroller = false;
	has_wind_aero = has_wind_dt = has_wind_pc = has_wind_tc = false;
	
        data->getValue(HAS_EXCITER, &has_ex, i);
        data->getValue(HAS_GOVERNOR, &has_gov, i);
	data->getValue(HAS_PSS, &has_pss, i);
	data->getValue(HAS_PLANT_CONTROLLER, &has_plantcontroller, i);
	data->getValue(HAS_WIND_TORQUECONTROL,&has_wind_tc, i);
	data->getValue(HAS_WIND_AERODYNAMIC,&has_wind_aero, i);
	data->getValue(HAS_WIND_DRIVETRAIN,&has_wind_dt, i);
	data->getValue(HAS_WIND_PITCHCONTROL,&has_wind_pc, i);

        if (generator) {
          boost::shared_ptr<BaseGeneratorModel> basegen;
          basegen.reset(generator);
          p_generators.push_back(basegen);

          if (has_ex) {
            if (data->getValue(EXCITER_MODEL, &model, i)) {
              BaseExciterModel *exciter
                = genFactory.createExciterModel(model);
              boost::shared_ptr<BaseExciterModel> ex;
              ex.reset(exciter);
              p_generators[icnt]->setExciter(ex);
            }
          }
          if (has_gov) {
            if (data->getValue(GOVERNOR_MODEL, &model, i)) {
              BaseGovernorModel *governor
                = genFactory.createGovernorModel(model);
              boost::shared_ptr<BaseGovernorModel> gov;
              gov.reset(governor);
              p_generators[icnt]->setGovernor(gov);
            }
          }
	  if (has_pss) {
            if (data->getValue(PSS_MODEL, &model, i)) {
              BasePssModel *pssmodel
                = genFactory.createPssModel(model);
              boost::shared_ptr<BasePssModel> pss;
              pss.reset(pssmodel);
              p_generators[icnt]->setPss(pss);
            }
          }
	  if (has_plantcontroller) {
            if (data->getValue(PLANT_CONTROLLER_MODEL, &model, i)) {
              BasePlantControllerModel *plantctrlmodel
                = genFactory.createPlantControllerModel(model);
              boost::shared_ptr<BasePlantControllerModel> plantctrl;
              plantctrl.reset(plantctrlmodel);
              p_generators[icnt]->setPlantController(plantctrl);
            }
          }

	  if (has_wind_tc) {
	    if(data->getValue(WIND_TORQUECONTROL, &model, i)) {
	      BaseMechanicalModel *torquectrlmodel
		= genFactory.createMechanicalModel(model);
	      boost::shared_ptr<BaseMechanicalModel> torquectrl;
	      torquectrl.reset(torquectrlmodel);
	      p_generators[icnt]->setTorqueController(torquectrl);
	    }
	  }

	  if (has_wind_pc) {
	    if(data->getValue(WIND_PITCHCONTROL, &model, i)) {
	      BaseMechanicalModel *pitchctrlmodel
		= genFactory.createMechanicalModel(model);
	      boost::shared_ptr<BaseMechanicalModel> pitchctrl;
	      pitchctrl.reset(pitchctrlmodel);
	      p_generators[icnt]->setPitchController(pitchctrl);
	    }
	  }

	  if (has_wind_dt) {
	    if(data->getValue(WIND_DRIVETRAIN, &model, i)) {
	      BaseMechanicalModel *drivetrainmodel
		= genFactory.createMechanicalModel(model);
	      boost::shared_ptr<BaseMechanicalModel> drivetrain;
	      drivetrain.reset(drivetrainmodel);
	      p_generators[icnt]->setDriveTrainModel(drivetrain);
	    }
	  }

	  if (has_wind_aero) {
	    if(data->getValue(WIND_AERODYNAMIC, &model, i)) {
	      BaseMechanicalModel *aerodynmodel
		= genFactory.createMechanicalModel(model);
	      boost::shared_ptr<BaseMechanicalModel> aerodyn;
	      aerodyn.reset(aerodynmodel);
	      p_generators[icnt]->setAeroDynamicModel(aerodyn);
	    }
	  }



          // create relay objective associate with the generator
          // get the number of relay associate with the generator
          nrelay = 0;
          data->getValue(RELAY_NUMBER, &nrelay);
          p_generators[icnt]->ClearRelay();
          relaycnt = 0;
          if (nrelay>0) {
            for (irelay=0 ; irelay<=nrelay ; irelay++) {
              data->getValue(RELAY_GENID, &relay_genid, irelay);
              if (relay_genid == genid) {
                if (data->getValue(RELAY_MODEL, &model, irelay)) {
                  if ( model== "FRQTPAT" ) {
                    BaseRelayModel *relaymodel
                      = relayFactory.createRelayModel(model);
                    boost::shared_ptr<BaseRelayModel> relay;
                    relay.reset(relaymodel);
                    relay->load(data, irelay);
                    p_generators[icnt]->AddRelay(relay);
                    relaycnt++;
                    bcomputefreq = true;
                  }  
                }
              }
            }
          }
        }
        p_generators[icnt]->load(data,i);
        if (has_gov) p_generators[icnt]->getGovernor()->load(data,i);
        if (has_ex) p_generators[icnt]->getExciter()->load(data,i);	
	if (has_pss) p_generators[icnt]->getPss()->load(data,i);	
	if (has_plantcontroller) p_generators[icnt]->getPlantController()->load(data,i);
	if (has_wind_tc) p_generators[icnt]->getTorqueController()->load(data,i);
	if (has_wind_pc) p_generators[icnt]->getPitchController()->load(data,i);
	if (has_wind_dt) p_generators[icnt]->getDriveTrainModel()->load(data,i);
	if (has_wind_aero) p_generators[icnt]->getAeroDynamicModel()->load(data,i);

      } else if (!data->getValue(GENERATOR_MODEL, &model, i) && pg >= 0.0){ 
	// handle the generators having no dynamic model, need to convert to negative load
	BaseGeneratorModel *generator = new gridpack::dynamic_simulation::BaseGeneratorModel;
	boost::shared_ptr<BaseGeneratorModel> basegen;
	basegen.reset(generator);
	p_generators.push_back(basegen);

	p_gen_nodynmodel[i] = true;
	p_genpg_nodynmodel[i] = pg;
	p_genqg_nodynmodel[i] = qg;
	p_ngen_nodynmodel++;
      } else if (pg < 0.0) {
        p_negpg.push_back(pg);
        p_negqg.push_back(qg);
        p_negngen++;
      }
      icnt++;
    }
  }

  // add load relay (LVSHBL) assoicate with the bus, renke add
  std::string model;
  nrelay = 0;
  data->getValue(RELAY_NUMBER, &nrelay);
  p_loadrelays.clear();
  if (nrelay>0) {
    for (irelay=0 ; irelay<=nrelay ; irelay++) {
      if (data->getValue(RELAY_MODEL, &model, irelay)) {
        if (model == "LVSHBL") {
          BaseRelayModel *relaymodel
            = relayFactory.createRelayModel(model);
          boost::shared_ptr<BaseRelayModel> relay;
          relay.reset(relaymodel);
          relay->load(data, irelay);
          p_loadrelays.push_back(relay);   
        }			
      }
    }			    
  }

  // add load model
  bool bdebug_load_model = false;
  double pl, ql, ip, iq, yp, yq, totaldynReactivepower;
  p_powerflowload_p.clear();
  p_powerflowload_p_save.clear();
  p_powerflowload_q.clear();
  p_powerflowload_q_save.clear();
  p_loadid.clear();
  p_powerflowload_status.clear();
  totaldynReactivepower = 0.0;
  if (bdebug_load_model) printf("DSFullBus::load():  entering processing load model \n");
  if (data->getValue(LOAD_NUMBER, &p_npowerflow_load)) {
    std::string loadid;
    int icnt = 0;
    if (bdebug_load_model) printf("bus %d has %d power flow loads \n", idx, p_npowerflow_load);
    int istat;
    for (i=0; i<p_npowerflow_load; i++) { 
      data->getValue(LOAD_PL, &pl, i);
      data->getValue(LOAD_QL, &ql, i);
      data->getValue(LOAD_IP, &ip, i);
      data->getValue(LOAD_IQ, &iq, i);
      data->getValue(LOAD_YP, &yp, i);
      data->getValue(LOAD_YQ, &yq, i);

      /* Convert all load to constant power. Later on it will be
	 converted to whatever load model it is. This is done because
	 right now the dynamics simulation application expects the load
	 to be specified as constant power in the data file. Using unity
	 voltage magnitude for conversion
      */
      /* Combine constant P, constant I, constant Y loads */
      pl = pl + ip + yp;
      ql = ql - iq - yq;

      data->getValue(LOAD_ID, &loadid, i);
      istat = 1;
      data->getValue(LOAD_STATUS, &istat, i);
      p_powerflowload_p.push_back(pl);
      p_powerflowload_p_save.push_back(pl);
      p_powerflowload_q_save.push_back(ql);
      p_powerflowload_q.push_back(ql);
      p_loadid.push_back(loadid);  
      p_powerflowload_status.push_back(istat);

      if (bdebug_load_model) printf("%d th power flow load at bus %d: %f + j%f\n", i, idx, pl, ql);	  
      std::string model;

      if (data->getValue(LOAD_MODEL, &model, i)) {
        if (bdebug_load_model) printf("dynamic load at bus %d, model = %s \n", idx, model.c_str());
        bcomputefreq = true;
        if ( model == "CMLDBLU1" ) {
          // all the dynamic loads will be added to the extended buses
          if (bdebug_load_model) printf("pop the powerflow load with ID %s back, as the type is %s !\n",
              loadid.c_str(), model.c_str());
          p_powerflowload_p.pop_back();
          p_powerflowload_q.pop_back();
          p_loadid.pop_back();	
          p_npowerflow_load = 0;
        }else{
          BaseLoadModel *load = loadFactory.createLoadModel(model); 
          if (bdebug_load_model) printf("DSFullBus::load(): base load object created!\n");	

          if (load) {
            boost::shared_ptr<BaseLoadModel> baseload;
            baseload.reset(load);
            p_loadmodels.push_back(baseload);
            p_dynamic_loadid.push_back(loadid);

            p_loadmodels[icnt]->load(data,i, pl, ql, 0);  //last parameter 0
            //means it is not load from composite load model
            //initialize the dynamic load model to get the Qini values
            p_loadmodels[icnt]->init(p_voltage, p_angle, 0.001);
            totaldynReactivepower += p_loadmodels[icnt]->getInitReactivePower();

            icnt++;
          }  // end of if (load)	
        } // end of the judgement if the dynamic load is not CMLDBLU1
      } // end of if (data->getValue(LOAD_MODEL, &model, i))
    } // end of for (i=0; i<p_npowerflow_load; i++)
		
  } // end of if (data->getValue(LOAD_NUMBER, &p_npowerflow_load))

  p_ndyn_load = p_loadmodels.size(); 
  //p_npowerflow_load = p_powerflowload_p.size();

  //sum all the power flow load P and Q at this bus together
  p_pl = 0.0;
  p_ql = 0.0;
  for (i=0; i<p_npowerflow_load; i++){
    if(p_powerflowload_status[i]) {
      p_pl+=p_powerflowload_p[i];
      p_ql+=p_powerflowload_q[i];
    }
  }

  if (bdebug_load_model) printf(" Bus %d have %d power flow loads, total %f + j%f \n", idx, p_npowerflow_load, p_pl, p_ql);

  //get total load P and Q for all dynamic loads at this bus
  double totaldyn_p, totaldyn_q, loadtmp;
  totaldyn_p = 0.0;
  totaldyn_q = 0.0;
  for (i=0; i<p_ndyn_load; i++){
    loadtmp = p_loadmodels[i]->getDynLoadP();
    totaldyn_p+=loadtmp;
    //loadtmp = p_loadmodels[i]->getInitReactivePower();
    //totaldyn_q+=loadtmp;
  }

  totaldyn_q = totaldynReactivepower;

  //set p_pl and p_ql as the static impedance load at this bus
  p_pl-=totaldyn_p;
  p_ql-=totaldyn_q;

  if (bdebug_load_model) printf(" Bus %d have %d dynamic loads, total %f + j%f \n",
      idx, p_ndyn_load, totaldyn_p, totaldyn_q);

  p_pl /= p_sbase;
  p_ql /= p_sbase;
  
  //p_scatterinjload_p = p_pl;
  //p_scatterinjload_q = p_ql;
  

  if (bdebug_load_model) printf(" Bus %d have remaining static loads for Y-bus: p_pl: %f pu, p_ql: %f pu, \n",
      idx, p_pl, p_ql);
	  
  for (i=0; i<p_ndyn_load; i++){
    // if the dynamic load model is the first AC motor model
	p_loadmodels[i]->setSameBusStaticLoadPQ(p_pl, p_ql, p_voltage);
	
  }
   
   // renke, this is the special code to determine which bus frequency need to be updated for wide area control
  //if (getOriginalIndex() == 34 || getOriginalIndex() == 30){ //renke hardcoded
  //if (idx == 34 || idx == 30){ //renke hardcoded
	//	bcomputefreq = true;
	//  printf("--------------------!!renke debug DSFullBus::load(), Bus No.: %d set bcomputefreq as true \n", getOriginalIndex());
  //}
	  

}

/**
 * Update data collection object with current values from simulation
 * @param data: DataCollection object containing parameters for this bus
 */
void gridpack::dynamic_simulation::DSFullBus::updateData(
    boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  int i;
  std::string name;
  gridpack::ComplexType voltage = getComplexVoltage();
  double rV = real(voltage);
  double iV = imag(voltage);
  rV = sqrt(rV*rV+iV*iV);
  if (!data->setValue(BUS_VMAG_CURRENT, rV)) {
    data->addValue(BUS_VMAG_CURRENT, rV);
  }
  for (i=0; i<p_ngen; i++) {
    if (data->getValue(GENERATOR_MODEL,&name,i)) {
      p_generators[i]->updateData(data, i);
    } else {
      if (!data->setValue(GENERATOR_PG_CURRENT, p_pg[i], i)) {
        data->addValue(GENERATOR_PG_CURRENT, p_pg[i], i);
      }
      if (!data->setValue(GENERATOR_QG_CURRENT, p_qg[i], i)) {
        data->addValue(GENERATOR_QG_CURRENT, p_qg[i], i);
      }
    }
  }
  int lcnt = 0;
  for (i=0; i<p_npowerflow_load; i++) {
    if (data->getValue(LOAD_MODEL,&name,i)) {
      p_loadmodels[lcnt]->updateData(data, i);
      lcnt++;
    } else {
      if (!data->setValue(LOAD_PL_CURRENT, p_powerflowload_p[i], i)) {
        data->addValue(LOAD_PL_CURRENT, p_powerflowload_p[i], i);
      }
      if (!data->setValue(LOAD_QL_CURRENT, p_powerflowload_q[i], i)) {
        data->addValue(LOAD_QL_CURRENT, p_powerflowload_q[i], i);
      }
    }
  }
}

/**
  * set voltage for the extended buses from composite load model
  */
void gridpack::dynamic_simulation::DSFullBus::setExtendedCmplBusVoltage(
 const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  std::string snewbustype;
  int ibustype = checkExtendedLoadBus();
	
  if (ibustype == -1) return;
  
  int iorgbusno;
  if (data->getValue(BUS_NUMBER, &iorgbusno)){
	  printf("DSFullBus::setExtendedCmplBusVoltage(), Bus No.: %d, Bus Type: %d \n", iorgbusno, ibustype);
  }
  
  //get neigb bus voltage information
  std::vector<boost::shared_ptr<BaseComponent> > nghbrs;
  getNeighborBuses(nghbrs);
  
  //get neigb bus voltage information
  std::vector<boost::shared_ptr<BaseComponent> > nghbbranch;
  getNeighborBranches(nghbbranch);
  
  int i, nghbrs_size, nghbbranch_size;
  nghbrs_size= nghbrs.size();
  nghbbranch_size = nghbbranch.size();
  gridpack::dynamic_simulation::DSFullBus *bus1;
  gridpack::dynamic_simulation::DSFullBranch *branch1;
  double bus_mag, bus_ang;
  
  if ( ibustype==1) { // if the bus is the LOW_SIDE_BUS
	for ( i=0 ; i<nghbrs_size; i++ ) {
		bus1 = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(nghbrs[i].get());
		if (bus1->checkExtendedLoadBus() == -1) // if the neigb bus is the normal bus, get the value of voltage
		{
			bus_mag= bus1->getVoltage();
			bus_ang= bus1->getPhase();
			printf("   DSFullBus::setExtendedCmplBusVoltage(), find corresponding raw bus, bus no.: %d, mag: %f, ang: %f, \n", bus1->getOriginalIndex(), bus_mag, bus_ang);
			
			break;
		}
	}
	
	//set the LOW_SIDE_BUS voltage
	setVoltage(bus_mag);
	setPhase(bus_ang);
  
	for ( i=0 ; i<nghbrs_size; i++ ) {
		bus1 = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(nghbrs[i].get());
		if (bus1->checkExtendedLoadBus() == 2) // if the neigb bus is the LOAD_BUS, set the voltage
		{
			bus1->setVoltage(bus_mag);
			bus1->setPhase(bus_ang);
			p_CmplFeederBus = bus1;
			printf("   DSFullBus::setExtendedCmplBusVoltage(), find LOAD_BUS, bus no.: %d, set value, mag: %f, ang: %f, \n", bus1->getOriginalIndex(), bus_mag, bus_ang);
			break;
		}
	}
	
	for ( i=0 ; i<nghbbranch_size; i++)
	{
		branch1 = dynamic_cast<gridpack::dynamic_simulation::DSFullBranch*>(nghbbranch[i].get());
		if (branch1->checkExtendedLoadBranchType() == 1){
			setCmplXfmrPt(branch1);
			p_CmplFeederBus->setCmplXfmrPt(branch1);
			printf("   DSFullBus::setExtendedCmplBusVoltage(), find corresponding xfmr branch \n");
		}
		
		if (branch1->checkExtendedLoadBranchType() == 2){
			setCmplXfeederPt(branch1);
			p_CmplFeederBus->setCmplXfeederPt(branch1);
			printf("   DSFullBus::setExtendedCmplBusVoltage(), find corresponding feeder branch \n");
		}
	}
	
  } // end of bus type = 1, LOW_SIDE_BUS
  
  if ( ibustype==2) { // if the bus is the LOAD_BUS
	for ( i=0 ; i<nghbrs_size; i++ ) {
		bus1 = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(nghbrs[i].get());
		if (bus1->checkExtendedLoadBus() == 1) // if the neigb bus is the LOW_SIDE_BUS
		{
			p_CmplXfmrBus = bus1;
			printf("   DSFullBus::setExtendedCmplBusVoltage(), find corresponding xfmr bus \n");
			break;
		}
	}
  
  }// end of bus type = 2, LOAD_BUS
 
}

/**
 * load parameters for the extended buses from composite load model
 */
void gridpack::dynamic_simulation::DSFullBus::LoadExtendedCmplBus(
   const boost::shared_ptr<gridpack::component::DataCollection> &data)
{ 
  std::string snewbustype;
  
  // check whether the bus is the ones extended by composite load models,
  // if no, return
  if (!data->getValue(NEW_BUS_TYPE, &snewbustype)) return;
  
  if ( snewbustype!="LOW_SIDE_BUS" && snewbustype!="LOAD_BUS" ) return;

  int iorgbusno;
  if (data->getValue(BUS_NUMBER, &iorgbusno)){
	  printf("DSFullBus::LoadExtendedCmplBus(), Bus No.: %d \n", iorgbusno);
  }
  
  if(!data->getValue(CASE_SBASE,&p_sbase)) p_sbase = 100.0;
  
  //load data for LOW_SIDE_BUS
  if (snewbustype=="LOW_SIDE_BUS") { 
    p_shunt = true;
	p_shunt = p_shunt && data->getValue(LOAD_BSS, &p_shunt_bs);
	data->getValue(LOAD_MVA, &loadMVABase);
	p_shunt_bs = 0.04; //tmp code, check and remove ???
	printf("DSFullBus::LoadExtendedCmplBus(), LOW_SIDE_BUS, bss = %f, loadMVABase =%f, \n ", p_shunt_bs, loadMVABase);
	
	p_shunt_bs = p_shunt_bs*loadMVABase/p_sbase;
	setParam(BUS_SHUNT_BL, p_shunt_bs, 0);
	printf("DSFullBus::LoadExtendedCmplBus(), LOW_SIDE_BUS, p_shunt_bs = %f\n", p_shunt_bs);
	return; 
  }
	  
  //load data for LOAD_BUS
  double pi = 4.0*atan(1.0);
 
  //p_shunt = true;
  //p_shunt = p_shunt && data->getValue(LOAD_BSS, &p_shunt_bs);

  LoadFactory loadFactory;
  
  int idx, i;
  data->getValue(BUS_NUMBER,&idx);
  
  // get the total of powerflow load
  double pl, ql;
  p_powerflowload_p.clear();
  p_powerflowload_q.clear();
  p_loadid.clear();
  printf("DSFullBus::LoadExtendedCmplBus():  entering processing load model \n");
  
  /*
  if (data->getValue(LOAD_NUMBER, &p_npowerflow_load)) {
    std::string loadid;
    int icnt = 0;
    printf("bus %d has %d power flow loads \n", idx, p_npowerflow_load);
    for (i=0; i<p_npowerflow_load; i++) { 
      data->getValue(LOAD_PL, &pl, i);
      data->getValue(LOAD_QL, &ql, i);
	  data->getValue(LOAD_ID, &loadid, i);
	  p_powerflowload_p.push_back(pl);
	  p_powerflowload_q.push_back(ql);
	  p_loadid.push_back(loadid);  

      printf("%d th power flow load at bus %d: %f + j%f\n", i, idx, pl, ql);	  
      std::string model;
    }
  }
  
  //sum all the power flow load P and Q at this bus together
  p_pl = 0.0;
  p_ql = 0.0;
  for (i=0; i<p_npowerflow_load; i++){
	  p_pl+=p_powerflowload_p[i];
	  p_ql+=p_powerflowload_q[i];
  }
  */
  
  data->getValue(LOAD_PL, &p_pl);
  data->getValue(LOAD_QL, &p_ql);
  
  //get the values for initialization of the substation power flow
  int tapNum;
  double Pload_MW, Qload_MVar, vt_mag, sysMVABase, Pload_pu, Qload_pu, loadFactor, Bss;
  double Bss_pu, Xxf, Xxf_pu, Rfdr, Xfdr, Rfdr_pu, Xfdr_pu, Tfixhs, Tfixls, Vlow_mag, Vmin, Vmax;
  double tap, Tmin, Tmax, step, Imag_lowbus, Ifeeder_mag, Vload_mag, Vload_ang, PloadBus_pu, QloadBus_pu;
  double FmA, FmB, FmC, FmD;
  
  gridpack::ComplexType vt = gridpack::ComplexType(p_voltage*cos(p_angle), p_voltage*sin(p_angle)); 
  
  printf("DSFullBus::LoadExtendedCmplBus(), p_voltage: %12.6f, p_angle: %12.6f, \n", p_voltage, p_angle);
  
  sysMVABase = p_sbase;
  Pload_MW = p_pl;
  Qload_MVar =p_ql;
  vt_mag = p_voltage;
  Pload_pu  = Pload_MW/sysMVABase;        // on system mva base
  Qload_pu  = Qload_MVar/sysMVABase;  
  
  data->getValue(LOAD_MVA, &loadMVABase);
  data->getValue(LOAD_BSS, &Bss);
  data->getValue(LOAD_XXF, &Xxf);
  data->getValue(LOAD_RFDR, &Rfdr);
  data->getValue(LOAD_XFDR, &Xfdr);
  data->getValue(LOAD_TFIXHS, &Tfixhs);
  data->getValue(LOAD_TFIXLS, &Tfixls);
  data->getValue(LOAD_VMIN, &Vmin);
  data->getValue(LOAD_VMAX, &Vmax);
  data->getValue(LOAD_TMIN, &Tmin);
  data->getValue(LOAD_TMAX, &Tmax);
  data->getValue(LOAD_STEP, &step);
  data->getValue(LOAD_FMA, &FmA);
  data->getValue(LOAD_FMB, &FmB);
  data->getValue(LOAD_FMC, &FmC);
  data->getValue(LOAD_FMD, &FmD);
  
  printf("DSFullBus::LoadExtendedCmplBus(), LOAD_BUS Par: loadMVABase: %f, Pload_pu: %f, Qload_pu: %f, Bss: %f, Xxf: %f, Rfdr: %f, xfdr: %f, Tfixhs: %f, Tfixls: %f, Vmin: %f, Vmax: %f, \n", 
		loadMVABase, Pload_pu, Qload_pu, Bss, Xxf, Rfdr, Xfdr, Tfixhs, Tfixls, Vmin, Vmax);
  printf("DSFullBus::LoadExtendedCmplBus(), LOAD_BUS Par: Tmin: %f, Tmax: %f, step: %f, FmA: %f, FmB: %f, FmC: %f, FmD: %f, Vmin: %f, Vmax: %f, \n", 
		Tmin, Tmax, step, FmA, FmB, FmC, FmD, Vmin, Vmax);
  
  // calculate the mva base if CMPLDW.loadMVABase <= 0.0
          
  if (loadMVABase < 0.0) {
	  loadFactor =abs(loadMVABase);
	  loadMVABase = Pload_MW/loadFactor;
  }else if (loadMVABase == 0.0) {
	  loadFactor =0.8;
	  loadMVABase = Pload_MW/loadFactor;
  }

  Bss_pu = Bss*loadMVABase/sysMVABase;
  
  // run substation power flow
          
  Xxf_pu = Xxf*sysMVABase/loadMVABase;

  gridpack::ComplexType cplx_tmp = gridpack::ComplexType(Pload_pu, Qload_pu); 
  gridpack::ComplexType cplx_tmp_3 = cplx_tmp/vt;
  gridpack::ComplexType Ilf_pu = gridpack::ComplexType(real(cplx_tmp_3), -imag(cplx_tmp_3));
  
  printf("DSFullBus::LoadExtendedCmplBus(), Bss_pu: %12.6f, Xxf_pu: %12.6f, vt: %12.6f +j* %12.6f, Sload_pu: %12.6f +j* %12.6f, Ilf_pu: %12.6f +j* %12.6f, \n", Bss_pu, Xxf_pu, real(vt), imag(vt), real(cplx_tmp), imag(cplx_tmp), real(Ilf_pu), imag(Ilf_pu));
  
  Rfdr_pu = Rfdr*sysMVABase/loadMVABase;
  Xfdr_pu = Xfdr*sysMVABase/loadMVABase;
  
  //set the parameters of the feeder branch with the composite load model
  p_CmplFeederBranch->SetCmplFeederBranch(Rfdr_pu, Xfdr_pu);
  p_CmplFeederBranch->setParam(BRANCH_R, Rfdr_pu, 0);
  p_CmplFeederBranch->setParam(BRANCH_X, Xfdr_pu, 0);
  p_CmplFeederBranch->printDSFullBranch();
  
  // voltage behind the Xxfr equivalent impedance
  Xxf_pu = Xxf_pu*Tfixhs*Tfixhs;
  Vlow_mag = (Vmin + Vmax)/2.0;  // targeted voltage mag at the low side of the transformer

 // calculate the tap
  tap = sqrt( (vt_mag*Vlow_mag)*(vt_mag*Vlow_mag) / 
     ( (Qload_pu*Xxf_pu-vt_mag*vt_mag)*(Qload_pu*Xxf_pu-vt_mag*vt_mag) + (Xxf_pu*Pload_pu)*(Xxf_pu*Pload_pu) ));
	 
  printf("DSFullBus::LoadExtendedCmplBus(), Rfdr_pu: %12.6f, Xfdr_pu: %12.6f, Xxf_pu: %12.6f, Vlow_mag: %12.6f, tap: %12.6f, \n", Rfdr_pu, Xfdr_pu, Xxf_pu, Vlow_mag, tap); 

  // need to check if Tap is within the limit
  bool tapReachLimit = false;
  if (tap<Tmin){
	  tap = Tmin;
      tapReachLimit = true;
  }else if(tap>Tmax) {
	  tap = Tmax;
      tapReachLimit = true;
  }

  // round to the closest tap
  if(!tapReachLimit) {
	 tapNum = round((tap-1.0)/step);
     tap = 1.0+double(tapNum)*step; 
  }
  
  printf("DSFullBus::LoadExtendedCmplBus(), after round the closest tap: tap = %12.6f, \n", tap);
  
  //set the parameters of the transformer branch with the composite load model
  p_CmplXfmrBranch->SetCmplXfmrBranch(Xxf_pu, 1.0/tap);
  p_CmplXfmrBranch->setParam(BRANCH_R, 0.0, 0);
  p_CmplXfmrBranch->setParam(BRANCH_X, Xxf_pu, 0);
  p_CmplXfmrBranch->setParam(BRANCH_TAP, 1.0/tap, 0);
  p_CmplXfmrBranch->printDSFullBranch();
  
  
  // original matlab formual: volt_high =  CMPLDW.vt- 1j* CMPLDW.Xxf_pu* CMPLDW.Tfixhs^2*Ilf_pu
  cplx_tmp  = Xxf_pu*Tfixhs*Tfixhs*Ilf_pu;
  gridpack::ComplexType cplx_tmp2 = gridpack::ComplexType(-imag(cplx_tmp), real(cplx_tmp));
  gridpack::ComplexType volt_high = vt - cplx_tmp2;
  
  printf("DSFullBus::LoadExtendedCmplBus(), cplx_tmp: %12.6f + j*%12.6f , volt_high: %12.6f + j*%12.6f , \n", real(cplx_tmp), imag(cplx_tmp), real(volt_high), imag(volt_high));
  
  //voltMag_high = abs(CMPLDW.volt_low)  // cmpl check???

  gridpack::ComplexType volt_low= volt_high*Tfixls*tap/Tfixhs;
  
  printf("DSFullBus::LoadExtendedCmplBus() volt_low at XFMR bus: %f + j*%f, \n", real(volt_low), imag(volt_low));
  printf("DSFullBus::LoadExtendedCmplBus() volt_low at XFMR bus, mag: %f, angle: %f, \n", abs(volt_low), atan2(imag(volt_low), real(volt_low)));
  
  //set the voltage value of the transformer bus with the composite load model
  p_CmplXfmrBus->setVoltage(abs(volt_low));
  p_CmplXfmrBus->setPhase(atan2(imag(volt_low), real(volt_low)));
  
  // current flowing from the low voltage side into the feeder
  gridpack::ComplexType I_lowbus = Ilf_pu*Tfixhs/(tap*Tfixls);
  
  printf("DSFullBus::LoadExtendedCmplBus() I_lowbus: %f + j*%f, \n", real(I_lowbus), imag(I_lowbus));
  
  Imag_lowbus = abs(I_lowbus);

  cplx_tmp = gridpack::ComplexType(0.0, Bss_pu);
  gridpack::ComplexType Ishunt = volt_low*cplx_tmp; // Bss charging current
  
  printf("DSFullBus::LoadExtendedCmplBus() Ishunt: %f + j*%f, \n", real(Ishunt), imag(Ishunt));

  gridpack::ComplexType Ifeeder = I_lowbus-Ishunt;
  
  printf("DSFullBus::LoadExtendedCmplBus() Ifeeder: %f + j*%f, \n", real(Ifeeder), imag(Ifeeder));

  Ifeeder_mag = abs(Ifeeder);

  // Neglecting the equivalent feeder charging
  cplx_tmp = gridpack::ComplexType(Rfdr_pu, Xfdr_pu);
  gridpack::ComplexType volt_load = volt_low - cplx_tmp*Ifeeder;

  Vload_mag = abs(volt_load);
  Vload_ang = atan2(imag(volt_load), real(volt_load)); //cmpl check???
  
  printf("DSFullBus::LoadExtendedCmplBus() volt_load at load bus: %f + j*%f, \n", real(volt_load), imag(volt_load));
  printf("DSFullBus::LoadExtendedCmplBus() volt_load at load bus, mag: %f, angle: %f, \n", Vload_mag, Vload_ang);
  
  setVoltage(Vload_mag);
  setPhase(Vload_ang);
  
  gridpack::ComplexType Ifeeder_conj = gridpack::ComplexType (real(Ifeeder), -imag(Ifeeder));
  gridpack::ComplexType Sload = volt_load*Ifeeder_conj;
  PloadBus_pu = real(Sload);
  QloadBus_pu = imag(Sload);
  
  printf("DSFullBus::LoadExtendedCmplBus() PloadBus_pu: %f, QloadBus_pu: %f, \n", PloadBus_pu, QloadBus_pu);
    
  double Fel ;
  Fel = 0.0; // not supported yet, force it to zero
  
  // NOTE: If sum of load fractions FmA, FmB, FmC, FmD, Fel is <1, remainder 
  // there is static load;If sum of fractions FmA, FmB, FmC, FmD, Fel is >1, 
  // fractions are normalized to 1 and there will be no static load
  //  
  // check If sum of fractions FmA, FmB, FmC, FmD, Fel is >1
  if (FmA+FmB+FmC+FmD > 1.0) {
	  double sum = FmA+FmB+FmC+FmD;
      FmA=FmA/sum;
      FmB=FmB/sum;
      FmC=FmC/sum;
      FmD=FmD/sum;  
  }
  
  // 1) create the dynamic component models based on
  // their fractions,i.e., Fma, Fmb, Fmc, Fmd and Fel
  //
  // 2) Feed/import the data into each motor dynamic model
  // including set the initial power of each model based on
  // fractions/percentages
  
  p_loadmodels.clear();
  int icnt = 0;       
  double totalLoadRactivePower = 0.0;
  double dtempQ = 0.0;
  double Pmotor;
  
  // load data for motorw A
  if (FmA > 0.0){
	double Pmotor =  PloadBus_pu*sysMVABase*FmA;
	bcomputefreq = true;
	
	printf("dynamic load MOTORW A at bus %d\n", iorgbusno);
	BaseLoadModel *load = loadFactory.createLoadModel("MOTORW"); 
	if (load) {
       boost::shared_ptr<BaseLoadModel> baseload;
       baseload.reset(load);
       p_loadmodels.push_back(baseload);
       p_loadmodels[icnt]->load(data,0, Pmotor, 0.0, 1);  //last parameter 1 means it is load from composite load model
	   p_loadmodels[icnt]->init(Vload_mag, Vload_ang, 0.001);  //check whether we need the timestep h
	   dtempQ = p_loadmodels[icnt]->getInitReactivePower();
	   totalLoadRactivePower += dtempQ;
	   //set load factor??????
	   icnt++;
	} 
  }
  
  // load data for motorw B
  if (FmB > 0.0){
	Pmotor =  PloadBus_pu*sysMVABase*FmB;
	bcomputefreq = true;
	
	printf("dynamic load MOTORW B at bus %d\n", iorgbusno);
	BaseLoadModel *load = loadFactory.createLoadModel("MOTORW"); 
	if (load) {
       boost::shared_ptr<BaseLoadModel> baseload;
       baseload.reset(load);
       p_loadmodels.push_back(baseload);
       p_loadmodels[icnt]->load(data,1, Pmotor, 0.0, 1);  //last parameter 1 means it is load from composite load model
	   p_loadmodels[icnt]->init(Vload_mag, Vload_ang, 0.001);  //check whether we need the timestep h
	   dtempQ = p_loadmodels[icnt]->getInitReactivePower();
	   totalLoadRactivePower += dtempQ;
	   //set load factor??????
	   icnt++;
	} 
  }
  
  // load data for motorw C
  if (FmC > 0.0){
	Pmotor =  PloadBus_pu*sysMVABase*FmC;
	bcomputefreq = true;
	
	printf("dynamic load MOTORW C at bus %d\n", iorgbusno);
	BaseLoadModel *load = loadFactory.createLoadModel("MOTORW"); 
	if (load) {
       boost::shared_ptr<BaseLoadModel> baseload;
       baseload.reset(load);
       p_loadmodels.push_back(baseload);
       p_loadmodels[icnt]->load(data,2, Pmotor, 0.0, 1);  //last parameter 1 means it is load from composite load model
	   p_loadmodels[icnt]->init(Vload_mag, Vload_ang, 0.001);  //check whether we need the timestep h
	   dtempQ = p_loadmodels[icnt]->getInitReactivePower();
	   totalLoadRactivePower += dtempQ;
	   //set load factor??????
	   icnt++;
	} 
  }
  
  // load data for ac motor D
  if (FmD > 0.0){
	Pmotor =  PloadBus_pu*sysMVABase*FmD;
	bcomputefreq = true;
	
	printf("dynamic load AC Motor D at bus %d\n", iorgbusno);
	BaseLoadModel *load = loadFactory.createLoadModel("ACMTBLU1"); 
	if (load) {
       boost::shared_ptr<BaseLoadModel> baseload;
       baseload.reset(load);
       p_loadmodels.push_back(baseload);
       p_loadmodels[icnt]->load(data,3, Pmotor, 0.0, 1);  //last parameter 1 means it is load from composite load model
	   p_loadmodels[icnt]->init(Vload_mag, Vload_ang, 0.001);  //check whether we need the timestep h
	   dtempQ = p_loadmodels[icnt]->getInitReactivePower();
	   totalLoadRactivePower += dtempQ;
	   //set load factor??????
	   icnt++;
	} 
  }
  
  // All the remaining loads are modeled by static loads
  // create and initialte the static loads
  // The static load model uses the following equations;
  // P=P0(P1c*V/V0^P1e+P2c*V/V0^P2e+P3)*(1+Pf*df)
  // Q=Q0*(Q1c*V/V0^Q1e+Q2c*V/V0^Q2e+Q3)*(1+Qf*df)
  // P0=Pload*(1.-FmA-FmB-FmC-FmD-Fel)
  // Q0=P0*tan(acos(PFs))
  // P3=1-P1c-P2c
  // Q3=1-Q1c-Q2c
  double Fstatic, Pint_static, Qint_static, cmpl_pf;
  if (1.0-(FmA+FmB+FmC+FmD+Fel) > 0.0){
	Fstatic =  1.0-(FmA+FmB+FmC+FmD+Fel);
	Pint_static =  PloadBus_pu*Fstatic*sysMVABase;  //check with Qiuhua
	data->getValue(LOAD_PFS, &cmpl_pf);
	Qint_static =  Pint_static*tan(acos(cmpl_pf));
	printf("   dynamic load Static Load at bus %d\n", iorgbusno);
	printf("   DSFullBus::LoadExtendedCmplBus(), Static load model, LOAD_PFS: %f, Pint_static: %f MW, Qint_static_pu: %f MVar,\n", cmpl_pf, Pint_static, Qint_static);

	BaseLoadModel *load = loadFactory.createLoadModel("IEELBL"); 
	if (load) {
       boost::shared_ptr<BaseLoadModel> baseload;
       baseload.reset(load);
       p_loadmodels.push_back(baseload);
       p_loadmodels[icnt]->load(data,-1, Pint_static, Qint_static, 1);  //last parameter 1 means it is load from composite load model
	   p_loadmodels[icnt]->init(Vload_mag, Vload_ang, 0.001);  //check whether we need the timestep h
	   totalLoadRactivePower += Qint_static; // check with qiuhua????
	   //set load factor??????
	   icnt++;
	} 
  }
      
  p_ndyn_load = p_loadmodels.size();
  
  // caclculate the difference between the reactive part of total dynamic loads and the power flow results
  // and add the compensation var to the load bus
  double compVar = totalLoadRactivePower/sysMVABase-QloadBus_pu;  // check with qiuhua???
  //p_shunt_bs = compVar/Vload_mag/Vload_mag;
  printf("  DSFullBus::LoadExtendedCmplBus(), Bus compensation var, compVar: %f pu, \n",  compVar); 
  p_pl = 0.0;
  p_ql = -compVar;  
  printf("  DSFullBus::LoadExtendedCmplBus(), Bus Y matrix loads, p_pl: %f, p_ql: %f, \n",  p_pl, p_ql);
  	
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::dynamic_simulation::DSFullBus::setMode(int mode)
{
  if (mode == YBUS || mode == YL || mode == PG || mode == jxd || mode == YDYNLOAD) {
    YMBus::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Return the value of the voltage magnitude on this bus
 * @return: voltage magnitude
 */
double gridpack::dynamic_simulation::DSFullBus::getVoltage(void)
{
  return p_voltage;
}

/**
 * Return the value of the phase angle on this bus
 * @return: phase angle
 */
double gridpack::dynamic_simulation::DSFullBus::getPhase(void)
{
	return p_angle;
}

/**
  * Return the value of whether the bus is an extended bus due to compositeload
  * Return true if this is an extended bus due to compositeload
  */
int gridpack::dynamic_simulation::DSFullBus::checkExtendedLoadBus(void)
{
	return p_bextendedloadbus;
}

/**
 * Set the point of the related extended transformer branch of this bus, due to composite load model
 */
void gridpack::dynamic_simulation::DSFullBus::setCmplXfmrPt(gridpack::dynamic_simulation::DSFullBranch* p_CmplXfmr) 
{
	p_CmplXfmrBranch = p_CmplXfmr;
}
	
/**
 * Set the point of the related extended feeder branch of this bus, due to composite load model
 */
void gridpack::dynamic_simulation::DSFullBus::setCmplXfeederPt(gridpack::dynamic_simulation::DSFullBranch* p_CmplFeeder) 
{
	p_CmplFeederBranch = p_CmplFeeder;
}
	
/**
  * Set the value of the voltage magnitude on this bus
  */
void gridpack::dynamic_simulation::DSFullBus::setVoltage(double mag)
{
	p_voltage = mag;
}

/**
 * Set the value of the phase angle on this bus
 */
void gridpack::dynamic_simulation::DSFullBus::setPhase(double ang)
{
	p_angle = ang;
}

/**
 * Return the complex value of the voltage on this bus
 * @return: complex value of the voltage
 */
gridpack::ComplexType gridpack::dynamic_simulation::DSFullBus::getComplexVoltage(void) //renke add
{
	return p_volt_full;
}

/**
* compute the value of the voltage frequency on this bus
* @return: voltage frequency
*/
void gridpack::dynamic_simulation::DSFullBus::computeBusVolFrequency( double timestep ) //renke add
{
	const double dFREQ_SYS = 60.0;
	const double dTf = 0.1;
	const double pi = 4.0*atan(1.0);
	const double dw0 = 2.0*dFREQ_SYS*pi;
	
	double dstatex, dstatex1, ddx1, ddx2, dva_old, dva;
	
	dstatex = pbusvolfreq_old/dFREQ_SYS - 1.0;
	//dva_old = atan2(imag(p_volt_full_old), real(p_volt_full_old));
	dva_old = atan2(p_volt_full_old_imag, p_volt_full_old_real);
	dva = atan2(imag(p_volt_full), real(p_volt_full));

        //printf ("pbusvolfreq_old, %8.4f,%8.4f, %8.4f,  %8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n",
        //pbusvolfreq_old, p_busvolfreq, dva_old, dva, p_volt_full_old_real, p_volt_full_old_imag, real(p_volt_full),imag(p_volt_full)); 
      	
	//process angle changing around +180 and -180 degrees
	if ( (dva_old>170.0/180.0*pi) && (dva<0.0))  {
		dva += 2.0*pi;
	}
	else if ( dva_old<-170.0/180.0*pi && dva>0.0) {
		dva -= 2.0*pi;
	}
	
	//prediction step
	ddx1 = ((dva - dva_old) / timestep / dw0 - dstatex) / dTf;
	dstatex1  = dstatex + ddx1 * timestep;
	
	// corrective step
    ddx2 = ((dva - dva_old)/ timestep / dw0 - dstatex1) / dTf;
    dstatex   = dstatex + (ddx1 + ddx2) / 2.0 * timestep;
    
    p_busvolfreq = (dstatex+1.0)*dFREQ_SYS;
	pbusvolfreq_old  = p_busvolfreq;  
	
}

/**
 * return the value of the voltage frequency on this bus
 * @return: voltage frequency
 */
double gridpack::dynamic_simulation::DSFullBus::getBusVolFrequency(void) //renke add
{
	
	return p_busvolfreq;
}

/**
 * set the value for the bcomputefreq
 * @return: void
 */
void gridpack::dynamic_simulation::DSFullBus::setBusVolFrequencyFlag(bool flag) //renke add
{
	 //printf("---rk debug, DSFullBus::setBusVolFrequencyFlag, bus %d set flag as %d\n", getOriginalIndex(), flag);
	 bcomputefreq = flag;
}

/**
 * update the old bus voltage with this bus
 */
void gridpack::dynamic_simulation::DSFullBus::updateoldbusvoltage (void) //renke add
{
	p_volt_full_old = p_volt_full;
	pbusvolfreq_old = p_busvolfreq;
	
	p_volt_full_old_real = real(p_volt_full);
	p_volt_full_old_imag = imag(p_volt_full);
	
}

void gridpack::dynamic_simulation::DSFullBus::printbusvoltage () //renke add
{
	//printf ("busvolt_old: %8.4f + j%8.4f,  busvolt: %8.4f + j%8.4f,\n", real(p_volt_full_old), imag(p_volt_full_old), real(p_volt_full),imag(p_volt_full) );
	printf ("busvolt_old: %8.4f + j%8.4f,  busvolt: %8.4f + j%8.4f,\n", p_volt_full_old_real, p_volt_full_old_imag, real(p_volt_full),imag(p_volt_full) );
}

/**
* update the relay status associate with this bus
*/
bool gridpack::dynamic_simulation::DSFullBus::updateRelay(bool flag, double delta_t) //renke add
{
	int i, nsize, itrip, itrip_prev;
	bool bbusflag;
	double dbusvoltfreq;
	gridpack::ComplexType cbusfreq;
	boost::shared_ptr<gridpack::dynamic_simulation::BaseRelayModel> p_relay;
	std::vector<gridpack::ComplexType*> vrelayvalue;
	
	//brelayflag = false;
	bbusflag = false;
	
	//update load relays
	vrelayvalue.push_back( &p_volt_full );
	if ( !p_loadrelays.empty()) {		
	  nsize = p_loadrelays.size();
	  for ( i=0; i<nsize ; i++ ) {
	    itrip = 0;
	    itrip_prev = 0;
	    
	    p_loadrelays[i]->setMonitorVariables(vrelayvalue);
	    p_loadrelays[i]->updateRelay(delta_t);
	    p_loadrelays[i]->getTripStatus(itrip, itrip_prev);
	    printf(" DSFullBus::updateRelay LVSHBL itrip = %d, itrip_prev = %d \n", itrip, itrip_prev);
	    if ( itrip==1 && itrip_prev==0 && p_loadrelays[i]->getOperationStatus()) {
	      //set the flag
	      bbusflag = true;
	      p_loadrelays[i]->setOperationStatus(false);
	      
	      double dfrac = p_loadrelays[i]->getRelayFracPar();
	      //change the bus Y Matrix, p_ybusr, p_ybusi;
	      printf("DSFullBus::updateRelay trip load shedding bus %d: dfrac = %8.4f, p_loadimpedancer = %8.4f, p_loadimpedancei = %8.4f \n", getOriginalIndex(), dfrac, p_loadimpedancer, p_loadimpedancei);
	      
	      p_ybusr = p_ybusr-p_loadimpedancer*dfrac;		//??????check the values are passed correctly!
	      p_ybusi = p_ybusi-p_loadimpedancei*dfrac;
	    }	
	  }
	}
	
	//update generator relays
	vrelayvalue.clear();
	
	if (p_ngen != 0) {
	  
	  if ( bcomputefreq == true ) {
	    computeBusVolFrequency(delta_t);
	    dbusvoltfreq = getBusVolFrequency();
	    
	    cbusfreq = gridpack::ComplexType(dbusvoltfreq, 0.0);
	    vrelayvalue.push_back( &cbusfreq );	
	  }
		
		
	  for (i = 0; i < p_ngen; i++) {
	    if(!p_gstatus[i] || !p_generators.size()) continue;
	    int irelay, nrelay;
	    p_generators[i]->getRelayNumber(nrelay);
			
	    if (nrelay >0) {
	      for ( irelay=0 ; irelay<nrelay; irelay++ ) {
		itrip = 0;
		itrip_prev = 0;
		p_relay = p_generators[i]->getRelay(irelay);
		p_relay->setMonitorVariables(vrelayvalue);
		p_relay->updateRelay(delta_t);
		p_relay->getTripStatus(itrip, itrip_prev);
		printf(" DSFullBus::updateRelay bus frequency: %8.4f \n", dbusvoltfreq);
		printf(" DSFullBus::updateRelay FRQTPAT itrip = %d, itrip_prev = %d \n", itrip, itrip_prev);
		if ( itrip==1 && itrip_prev==0 && p_relay->getOperationStatus()) {
		  //set the flag
		  bbusflag = true;
		  p_relay->setOperationStatus(false);
		  //set the generator be out of service status, as tripped by the relay
		  p_generators[i]->SetGenServiceStatus(false);
		  
		  //change the bus Y Matrix, p_ybusr, p_ybusi;
		  gridpack::ComplexType Y_a
		    = p_generators[i]->NortonImpedence();
		  printf("DSFullBus::updateRelay tripped gen real(Y_a): %8.4f imag(Y_a): %8.4f\n",real(Y_a),imag(Y_a));
		  p_ybusr = p_ybusr - real(Y_a);
		  p_ybusi = p_ybusi - imag(Y_a);
		  
		}
	      }
	    }
	  }
	}
	
	p_busrelaytripflag = bbusflag;
	return bbusflag;
}


/**
 * Return whether or not a bus is isolated
 * @return true if bus is isolated
 */
bool gridpack::dynamic_simulation::DSFullBus::isIsolated(void) const
{
  return YMBus::isIsolated();
}

/**
 * Return the number of generators on this bus
 * @return number of generators on bus
 */
int gridpack::dynamic_simulation::DSFullBus::getNumGen(void)
{
  return p_ngen;
}

void gridpack::dynamic_simulation::DSFullBus::setIFunc(void)
{
}

void gridpack::dynamic_simulation::DSFullBus::setIJaco(void)
{
}

/**
 * Check to see if a fault event applies to this bus and set an internal
 * flag marking the bus as the "from" or "to" bus for the event
 * @param from_idx index of "from" bus for fault event
 * @param to_idx index of "to" bus for fault event
 */
void gridpack::dynamic_simulation::DSFullBus::setEvent(int from_idx, int to_idx,
  gridpack::component::BaseBranchComponent* branch_ptr)
{
  if (from_idx == getOriginalIndex()) {
    p_from_flag = true;
  } else {
    p_from_flag = false;
  }
  if (to_idx == getOriginalIndex()) {
    p_to_flag = true;
  } else {
    p_to_flag = false;
  }
  if (p_to_flag || p_from_flag) {
    p_branch = branch_ptr;
  } else {
    p_branch = NULL;
  }
}

/**
 * Clear fault event from bus
 */
void gridpack::dynamic_simulation::DSFullBus::clearEvent()
{
  p_from_flag = false;
  p_to_flag = false;
  p_branch = NULL;
}

void gridpack::dynamic_simulation::DSFullBus::setBranchTripAction(int from_idx, int to_idx,
  gridpack::component::BaseBranchComponent* branch_ptr)
{
  if (from_idx == getOriginalIndex()) {
    p_branchtripaction_from_flag = true;
  } else {
    p_branchtripaction_from_flag = false;
  }
  if (to_idx == getOriginalIndex()) {
    p_branchtripaction_to_flag = true;
  } else {
    p_branchtripaction_to_flag = false;
  }
  
  if (p_branchtripaction_to_flag || p_branchtripaction_from_flag) {
    p_vec_tripactionbranch.push_back(branch_ptr);
  } 
}

void gridpack::dynamic_simulation::DSFullBus::clearBranchTripAction()
{
	p_branchtripaction_from_flag = false;
	p_branchtripaction_to_flag = false;
	p_vec_tripactionbranch.clear();
}

void gridpack::dynamic_simulation::DSFullBus::applyConstYLoad_Change_P(double loadPChangeMW)
 {
	 p_Yload_change_P_flag = true;
	 p_Yload_change_r = loadPChangeMW/100.0;
	  
 }
 
 void gridpack::dynamic_simulation::DSFullBus::applyConstYLoad_Change_Q(double loadPChangeMVAR)
 {
	 p_Yload_change_Q_flag = true;
	 p_Yload_change_i = loadPChangeMVAR/100.0;
	  
 }
 
 bool gridpack::dynamic_simulation::DSFullBus::setConstYLoadtoZero_P( )
 {
	 if (!p_bConstYLoadSettoZero_P){
		p_Yload_change_P_flag = true;
		p_Yload_change_r = -p_pl;
		p_bConstYLoadSettoZero_P = true;
		return true;
	 }else{
		p_Yload_change_P_flag = false;
		p_Yload_change_r = 0.0; 
		return false;
	 }
	  
 }
 
   bool gridpack::dynamic_simulation::DSFullBus::setConstYLoadtoValue( double impedancer, double impedancei)
 {
	 if (!p_bConstYLoadSettoValue){
		p_Yload_change_P_flag = true;
		p_Yload_change_Q_flag = true;	
		
		double voltagemagtmp = abs(p_volt_full);
		
		gridpack::ComplexType z_a(impedancer, impedancei);
		gridpack::ComplexType Y_a;
		gridpack::ComplexType volttmp(1.0, 0.0);
		Y_a = volttmp / z_a;
		p_Yload_change_r = -p_pl + real(Y_a)*(p_voltage*p_voltage);
		p_Yload_change_i = -p_ql + imag(Y_a)*(p_voltage*p_voltage);
		
		p_bConstYLoadSettoValue = true;

		return true;		
		
	 }else{
		p_Yload_change_P_flag = false;
		p_Yload_change_r = 0.0; 
		p_Yload_change_Q_flag = false;
		p_Yload_change_i = 0.0; 
		return false;
	 }
	  
 }
 
 /*
  bool gridpack::dynamic_simulation::DSFullBus::setConstYLoadtoValue_P( double impedancer)
 {
	 if (!p_bConstYLoadSettoValue_P){
		p_Yload_change_P_flag = true;
		double voltagemagtmp = abs(p_volt_full);
		//double tmp1 = p_Yload_change_r/(p_voltage*p_voltage);
		p_Yload_change_r = -p_pl + impedancer*(p_voltage*p_voltage);
		p_bConstYLoadSettoValue_P = true;
		return true;
	 }else{
		p_Yload_change_P_flag = false;
		p_Yload_change_r = 0.0; 
		return false;
	 }
	  
 }
  */
 
 bool gridpack::dynamic_simulation::DSFullBus::setConstYLoadtoZero_Q( )
 {
	 if (!p_bConstYLoadSettoZero_Q){
		p_Yload_change_Q_flag = true;
		p_Yload_change_i = -p_ql;
		p_bConstYLoadSettoZero_Q = true;
		return true;
	 }else{
		p_Yload_change_Q_flag = false;
		p_Yload_change_i = 0.0; 
		return false;
	 }
	  
 }

 /*
  bool gridpack::dynamic_simulation::DSFullBus::setConstYLoadtoValue_Q( double impedancei)
 {
	 if (!p_bConstYLoadSettoValue_Q){
		p_Yload_change_Q_flag = true;
		double voltagemagtmp = abs(p_volt_full);
		p_Yload_change_i = -p_ql - impedancei*(p_voltage*p_voltage);
		p_bConstYLoadSettoValue_Q = true;
		return true;
	 }else{
		p_Yload_change_Q_flag = false;
		p_Yload_change_i = 0.0; 
		return false;
	 }
	  
 }
 */

void gridpack::dynamic_simulation::DSFullBus::clearConstYLoad_Change_P()
 {
	 p_Yload_change_P_flag = false;
	 p_Yload_change_r = 0.0;
	  
 }
 
 void gridpack::dynamic_simulation::DSFullBus::clearConstYLoad_Change_Q()
 {
	 p_Yload_change_Q_flag = false;
	 p_Yload_change_i = 0.0;
	  
 }

void gridpack::dynamic_simulation::DSFullBus::setBranchRelayFromBusStatus(bool sta)
{
  p_branchrelay_from_flag = sta;
}
void gridpack::dynamic_simulation::DSFullBus::setBranchRelayToBusStatus(bool sta)
{
  p_branchrelay_to_flag = sta;
}

void gridpack::dynamic_simulation::DSFullBus::setRelayTrippedbranch(gridpack::component::BaseBranchComponent* branch_ptr)
{
  p_relaytrippedbranch = branch_ptr;
}

void gridpack::dynamic_simulation::DSFullBus::clearRelayTrippedbranch()
{
  p_relaytrippedbranch = NULL;
}

bool gridpack::dynamic_simulation::DSFullBus::checkisolated()
{
  return p_isolated;
}

/**
 * Write output from buses to standard out
 * @param string (output) string with information to be printed out
 * @param bufsize size of string buffer in bytes
 * @param signal an optional character string to signal to this
 * routine what about kind of information to write
 * @return true if bus is contributing string to output, false otherwise
 */
bool gridpack::dynamic_simulation::DSFullBus::serialWrite(char *string,
    const int bufsize, const char *signal)
{
  if (p_ngen == 0 && p_ndyn_load == 0 && p_npowerflow_load == 0) return false;
  int i;
  char buf[128];
  char *ptr = string;
  int idx = getOriginalIndex();
  buf[0] = '\0';
  if (signal == NULL) {
    return false;
  } else if (!strcmp(signal,"watch_header") ||
      !strcmp(signal,"watch")) {
    if (p_ngen == 0) return false;
    int len = 0;
    bool ok;
    for (i=0; i<p_ngen; i++) {
      if(!p_generators.size()) continue;
      if (p_generators[i]->getWatch()) {
	char buf[128];
        ///printf("(DSFull::serialWrite) Got to 1\n");
        ok = p_generators[i]->serialWrite(buf,128,signal);
        ///printf("(DSFull::serialWrite) Got to 2\n");
        if (ok) {
          int slen = strlen(buf);
          if (len+slen < bufsize) sprintf(ptr,"%s",buf);
          len += slen;
          ptr += slen;
        }
      }
    }
    if (len > 0) return true;
  } else if (!strcmp(signal,"load_watch_header") ||
      !strcmp(signal,"load_watch")) {
    if (p_ndyn_load == 0) return false;
    int i;
    char buf[128];
    char *ptr = string;
    int len = 0;
    bool ok;
    for (i=0; i<p_ndyn_load; i++) {
      if (p_loadmodels[i]->getWatch()) {
        ///printf("(DSFull::serialWrite) Got to 1\n");
        ok = p_loadmodels[i]->serialWrite(buf,128,signal);
        ///printf("(DSFull::serialWrite) Got to 2\n");
        if (ok) {
          int slen = strlen(buf);
          if (len+slen < bufsize) sprintf(ptr,"%s",buf);
          len += slen;
          ptr += slen;
        }
      }
    }
    if (len > 0) return true;
  } else if (!strcmp(signal,"src_gen")) {
    if (p_source) {
      char sbuf[128];
      char *cptr = string;
      int i, len, slen = 0;
      std::string status; 
      for (i=0; i<p_ngen; i++) {
        if (p_gstatus[i]) {
          status = "  active";
        } else { 
          status = "inactive";
        }
        sprintf(sbuf,"%8d %s %s %4d %4d %14.4f %14.4f %14.4f %14.4f\n",
            getOriginalIndex(),
            p_genid[i].c_str(),status.c_str(),p_area,p_zone,p_savePg[i],
            p_pg[i],p_gpmin[i],p_gpmax[i]);
        len = strlen(sbuf);
        if (slen+len <= bufsize) {
          sprintf(cptr,"%s",sbuf);
          slen += len;
          cptr += len;
        }
      }
      if (slen>0) {
        return true;
      } else { 
        return false;
      }
    } else {
      return false;
    }
  } else if (!strcmp(signal,"sink_load")) {
    if (p_sink) {
      char sbuf[128];
      char *cptr = string;
      int i, len, slen = 0;
      std::string status;
      for (i=0; i<p_npowerflow_load; i++) {
        if (p_powerflowload_status[i]) {
          status = "  active";
        } else {
          status = "inactive";
        }
        sprintf(sbuf,"%8d %s %s %4d %4d %14.4f %14.4f %14.4f %14.4f\n",
            getOriginalIndex(),p_loadid[i].c_str(),status.c_str(),p_area,
            p_zone,p_powerflowload_p_save[i],p_powerflowload_p[i],
            p_powerflowload_q_save[i],p_powerflowload_q[i]);
        len = strlen(sbuf);
        if (slen+len <= bufsize) {
          sprintf(cptr,"%s",sbuf);
          slen += len;
          cptr += len;
        }
      }
      if (slen>0) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  } else if (strlen(signal) > 0) {
    int i;
    char buf[128];
    int len = 0;
    bool ok = true;
  //  printf("Writing for %d generators\n",p_ngen);
    for (i=0; i<p_ngen; i++) {
      if(!p_gstatus[i] || !p_generators.size()) continue;
      if (p_generators[i]->serialWrite(buf,128,signal)) {
        int slen = strlen(buf);
        if (len+slen < bufsize) sprintf(ptr,"%s",buf);
        len += slen;
        ptr += slen;
      }
    }
    if (len > 0) return true;
  }
  return false;
}

/**
 * Add constant impedance load admittance to diagonal elements of
 * Y-matrix
 */
void gridpack::dynamic_simulation::DSFullBus::addLoadAdmittance()
{
  p_ybusr = p_ybusr+p_pl/(p_voltage*p_voltage);
  p_ybusi = p_ybusi+(-p_ql)/(p_voltage*p_voltage);
  //printf("idx: %d %f %f\n", getOriginalIndex(), p_ybusr, p_ybusi);
}

/**
 * Set load on bus
 * @param pl real load
 * @param ql imaginary load
 */
void gridpack::dynamic_simulation::DSFullBus::setLoad(double pl, double ql)
{
  p_pl = pl;
  p_ql = ql;
}

/**
 * Get load on bus
 * @param pl real load
 * @param ql imaginary load
 */
void gridpack::dynamic_simulation::DSFullBus::getLoad(double *pl, double *ql)
{
  *pl = p_pl;
  *ql = p_ql;
}

/**
 * Set value of real power on individual generators
 * @param tag generator ID
 * @param value new value of real power
 * @param data data collection object associated with bus
 */
void gridpack::dynamic_simulation::DSFullBus::setGeneratorRealPower(
    std::string tag, double value, gridpack::component::DataCollection *data)
{
  int i, idx;
  idx = -1;
  for (i=0; i<p_ngen; i++) {
    if (p_genid[i] == tag) {
      idx = i;
      break;
    }
  }
  if (idx != -1) {
    p_pg[idx] = value;
    data->setValue(GENERATOR_PG,value,idx);
  } else {
    printf("DSsetGeneratorRealPower: No generator found on bus %d with id: (%s)\n",getOriginalIndex(),tag.c_str());
  }
}

/**
 * Set status variable on individual generators
 * @param tag generator ID
 * @param value new value of status
 * @param data data collection object associated with bus
 */
void gridpack::dynamic_simulation::DSFullBus::setGeneratorStatus(
    std::string tag, int status, gridpack::component::DataCollection *data)
{
  int i, idx;
  idx = -1;
  for (i=0; i<p_ngen; i++) {
    if (p_genid[i] == tag) {
      idx = i;
      break;
    }
  }
  if (idx != -1) {
    p_gstatus[idx] = status;
    data->setValue(GENERATOR_STAT,status,idx);
    if(!status) {
      p_pg[idx] = p_qg[idx] = 0.0;
      data->setValue(GENERATOR_PG, 0.0, idx);
      data->setValue(GENERATOR_QG, 0.0, idx);
    }
  } else {
    printf("DSsetGeneratorStatus: No generator found on bus %d with id: (%s)\n",getOriginalIndex(),tag.c_str());
  }
}

/**
 * Set value of real power on individual loads
 * @param tag load ID
 * @param value new value of real power
 * @param data data collection object associated with bus
 */
void gridpack::dynamic_simulation::DSFullBus::setLoadRealPower(
    std::string tag, double value, gridpack::component::DataCollection *data)
{
  /*
  int i, idx;
  idx = -1;
  for (i=0; i<p_nload; i++) {
    if (p_loadid[i] == tag) {
      idx = i;
      break;
    }
  }
  p_pl[idx] = value;
  */
  data->setValue(LOAD_PL,value);
  data->setValue(LOAD_PL,value,0);
}

#ifdef USE_FNCS
/**
 * Retrieve an opaque data item from component.
 * @param data item to retrieve from component
 * @param signal string to control behavior of routine
 * (currently ignored)
 * @return true if component is returning data item,
 * false otherwise
 */
bool gridpack::dynamic_simulation::DSFullBus::getDataItem(void *data, const char *signal)
{
  voltage_data *vdata = static_cast<voltage_data*>(data);
  vdata->busID = getOriginalIndex();
  vdata->voltage = gridpack::ComplexType(p_voltage*sin(p_angle),
      p_voltage*cos(p_angle));
}
#endif

/**
 * Set an internal parameter that specifies that the rotor speed and angle
 * for the generator corresponding to the string tag are to be printed to
 * output
 * @param tag 2-character identifier of generator
 * @param flag set to true to monitor generator
 */

void gridpack::dynamic_simulation::DSFullBus::setWatch(std::string tag, bool flag)
{
  int i;
  for (i=0; i<p_genid.size(); i++) {
    if (tag == p_genid[i]) {
      p_generators[i]->setWatch(flag);
      break;
    }
  }
}

/**
 * Return a list of watched generators
 * @return list of generator tags
 */
std::vector<std::string> gridpack::dynamic_simulation::DSFullBus::getWatchedGenerators()
{
  std::vector<std::string> ret;
  int i;
  for (i=0; i<p_ngen; i++) {
    if(!p_gstatus[i] || !p_generators.size()) continue;
    if (p_generators[i]->getWatch()) ret.push_back(p_genid[i]);
  }
  return ret;
}

/**
 * Return a vector of watched values
 * @return rotor angle and speed for all watched generators on bus
 */
std::vector<double> gridpack::dynamic_simulation::DSFullBus::getWatchedValues()
{
  std::vector<double> ret;
  int i, j;
  for (i=0; i<p_ngen; i++) {
    if(!p_generators.size()) continue;
    if (p_generators[i]->getWatch()) {
      std::vector<double> vals;
      p_generators[i]->getWatchValues(vals);
      for (j=0; j<vals.size(); j++) ret.push_back(vals[j]);
    }
  }
  return ret;
}

/**
 * Return a vector of watched values
 * @return vector of watched values for the given generator tag (id) at a bus
 */
std::vector<double> gridpack::dynamic_simulation::DSFullBus::getWatchedValues(std::string tag)
{
  std::vector<double> ret;
  int i, j;
  for (i=0; i<p_ngen; i++) {
    if(!p_generators.size()) continue;
    if (p_generators[i]->getWatch()) {
      if (p_genid[i] == tag) {
	std::vector<double> vals;
	p_generators[i]->getWatchValues(vals);
	for (j=0; j<vals.size(); j++) ret.push_back(vals[j]);
      }
    }
  }
  return ret;
}


/**
 * Return rotor speed and angle for a specific generator
 * @param idx index of generator
 * @param speed generator rotor speed
 * @param angle generator rotor angle
 * @param speed generator real power
 * @param angle generator reactive power
 */
void gridpack::dynamic_simulation::DSFullBus::getWatchedValues(int idx,
    double *speed, double *angle, double *genP, double *genQ)
{
  std::vector<double> vals;
  p_generators[idx]->getWatchValues(vals);
  *angle = vals[0];
  *speed = vals[1];
  *genP = vals[2];
  *genQ = vals[3];
  
}


/**
 * Check generators for frequency violations
 * @param start time at which monitoring begins
 * @param time current time
 * @return true if no violation has occured
 */
bool gridpack::dynamic_simulation::DSFullBus::checkFrequency(
    double start, double time)
{
  bool ret = true;
  // don't check anthing until current time exceeds start. This is to exclude
  // prefault period of simulation
  if (time >= start) {
    int i;
    for (i=0; i<p_genid.size(); i++) {
      if (p_generators[i]->getWatch()) {
        std::vector<double> vals;
        p_generators[i]->getWatchValues(vals);
        double freq = 60.0*vals[1];
        // check for upward drift
        if (!p_upStartedMonitoring[i]) {
          if (freq > 61.0) {
            p_upStartedMonitoring[i] = true;
            p_upIntervalStart[i] = time;
          }
        } else {
          if (freq - p_previousFrequency[i] > 0.0 && freq > 61.0) {
            if (time - p_upIntervalStart[i] >= 0.5) {
              ret = false;
            }
          } else {
            p_upStartedMonitoring[i] = false;
          }
        }
        // check for downward drift
        if (!p_downStartedMonitoring[i]) {
          if (freq < 59.0) {
            p_downStartedMonitoring[i] = true;
            p_downIntervalStart[i] = time;
          }
        } else {
          if (freq - p_previousFrequency[i] < 0.0 && freq < 59.0) {
            if (time - p_downIntervalStart[i] >= 0.5) {
              ret = false;
            }
          } else {
            p_downStartedMonitoring[i] = false;
          }
        }
        p_previousFrequency[i] = freq;
      }
    }
  }
  return ret;
}

/**
 * Check generators for frequency violations
 * @param limit maximum allowable frequency excursion
 * @return true if no violation has occured
 */
bool gridpack::dynamic_simulation::DSFullBus::checkFrequency(double limit)
{
  bool ret = true;
  int i;
  for (i=0; i<p_genid.size(); i++) {
    if (p_generators[i]->getWatch()) {
      std::vector<double> vals;
      p_generators[i]->getWatchValues(vals);
      double freq = 60.0*vals[1];
      if (freq > limit) ret = false;
    }
  }
  return ret;
}

/**
 * Scale value of real power on all generators
 * @param character ID for generator
 * @param value scale factor for real power
 * @param data data collection object for bus holding generators
 */
void gridpack::dynamic_simulation::DSFullBus::scaleGeneratorRealPower(
    std::string tag, double value,
    boost::shared_ptr<gridpack::component::DataCollection> data)
{
  int i;
  for (i=0; i<p_ngen; i++) {
    if (p_genid[i] == tag) {
      if (value > 0.0) {
        double excess = p_gpmax[i]-p_pg[i];
        if (excess < 0.0) {
          printf("bus: %d generator: %s excess (pt): %f (pg): %f\n",
              getOriginalIndex(),tag.c_str(),p_gpmax[i],p_pg[i]);
        }
        p_pg[i] += value*excess;
      } else {
        double slack = p_pg[i]-p_gpmin[i];
        p_pg[i] += value*slack;
      }
      data->setValue(GENERATOR_PG,p_pg[i],i);
      break;
    }
  }
}

/**
 * Scale value of power on loads
 * @param character ID for load
 * @param value scale factor for real power
 */
void gridpack::dynamic_simulation::DSFullBus::scaleLoadPower(
    std::string tag, double value)
{
  int i;
  for (i=0; i<p_npowerflow_load; i++) {
    if (p_loadid[i] == tag) {
      p_powerflowload_p[i] = value*p_powerflowload_p[i];
      p_powerflowload_q[i] = value*p_powerflowload_q[i];
      break;
    }
  }
}

/**
 * Reset real power for generators and loads back to original values
 * @param data data collection object for bus
 */
void gridpack::dynamic_simulation::DSFullBus::resetPower(
    boost::shared_ptr<gridpack::component::DataCollection> data)
{
  int i;
  for (i=0; i<p_ngen; i++) {
    p_pg[i] = p_savePg[i];
    data->setValue(GENERATOR_PG,p_pg[i],i);
  }
  for (i=0; i<p_npowerflow_load; i++) {
    p_powerflowload_p[i] = p_powerflowload_p_save[i];
    p_powerflowload_q[i] = p_powerflowload_q_save[i];
  }
}



/**
 * Get available margin for generator
 * @param tag character ID for generator
 * @param current initial generation
 * @param pmin minimum allowable generation
 * @param pmax maximum available generation
 * @param status current status of generator
 */
void gridpack::dynamic_simulation::DSFullBus::getGeneratorMargins(
    std::vector<std::string> &tag,
    std::vector<double> &current, std::vector<double> &pmin,
    std::vector<double> &pmax,std::vector<int> &status)
{
  tag.clear();
  current.clear();
  pmin.clear();
  pmax.clear();
  status.clear();
  int i;
  for (i=0; i<p_ngen; i++) {
    tag.push_back(p_genid[i]);
    current.push_back(p_savePg[i]);
    pmin.push_back(p_gpmin[i]);
    pmax.push_back(p_gpmax[i]);
    status.push_back(1);
  }
}

/**
 * Get current value of loads
 * @param tag character ID for load
 * @param current initial value of load
 * @param status current status of load
 */
void gridpack::dynamic_simulation::DSFullBus::getLoadPower(
    std::vector<std::string> &tag, std::vector<double> &pl,
    std::vector<double> &ql, std::vector<int> &status)
{
  tag.clear();
  pl.clear();
  ql.clear();
  status.clear();
  int nloads = p_powerflowload_p.size();
  int i;
  for (i=0; i<nloads; i++) {
    tag.push_back(p_loadid[i]);
    pl.push_back(p_powerflowload_p_save[i]);
    ql.push_back(p_powerflowload_q_save[i]);
    status.push_back(1);
  }
}

/**
 * Label bus as a source for real time path rating
 * @param flag identify bus as source
 */
void gridpack::dynamic_simulation::DSFullBus::setSource(bool flag)
{
  p_source = flag;
}

/**
 * Label bus as a sink for real time path rating
 * @param flag identify bus as sink
 */
void gridpack::dynamic_simulation::DSFullBus::setSink(bool flag)
{
  p_sink = flag;
}

/**
 * Store scale factor
 * @param scale factor for scaling generation or loads
 */
void gridpack::dynamic_simulation::DSFullBus::setScale(double scale)
{
  p_rtpr_scale = scale;
}

/**
 * Return pointer to generator corresponding to two-character id
 * @param id two character id labeling generator
 * @return pointer to generator device
 */
gridpack::dynamic_simulation::BaseGeneratorModel
*gridpack::dynamic_simulation::DSFullBus::getGenerator(std::string id)
{
  int i;
  for (i=0; i<p_generators.size(); i++) {
    if (id == p_genid[i]) {
      return p_generators[i].get();
    }
  }
  return NULL;
}

/**
 * Return pointer to load corresponding to two-character id
 * @param id two character id labeling load
 * @return pointer to load device
 */
gridpack::dynamic_simulation::BaseLoadModel
*gridpack::dynamic_simulation::DSFullBus::getLoad(std::string id)
{
  int i;
  for (i=0; i<p_loadmodels.size(); i++) {
    if (id == p_loadid[i]) {
      return p_loadmodels[i].get();
    }
  }
  return NULL;
}

/**
 * Return pointer to relay corresponding to two-character id
 * @param id two character id labeling relay
 * @return pointer to relay device
 */
gridpack::dynamic_simulation::BaseRelayModel
*gridpack::dynamic_simulation::DSFullBus::getRelay(std::string id)
{
}


/**
 * Get list of generator IDs
 * @return vector of generator IDs
 */
std::vector<std::string> gridpack::dynamic_simulation::DSFullBus::getGenerators()
{
  return p_genid;
}

/**
 * Get list of load IDs
 * @return vector of load IDs
 */
std::vector<std::string> gridpack::dynamic_simulation::DSFullBus::getLoads()
{
  return p_loadid;
}

/**
 * Get list of IDs for dynamic loads
 * @return vector of dynamic load IDs
 */
std::vector<std::string> gridpack::dynamic_simulation::DSFullBus::getDynamicLoads()
{
  return p_dynamic_loadid;
}

/**
 * Get area parameter for bus
 * @return bus area index
 */
int gridpack::dynamic_simulation::DSFullBus::getArea()
{
  return p_area;
}

/**
 * Get zone parameter for bus
 * @return bus zone index
 */
int gridpack::dynamic_simulation::DSFullBus::getZone()
{
  return p_zone;
}

/**
 * execute load scattering, the P and Q values of the STATIC load at certain buses will be changed to the values of 
 * the loadP and loadQ
*/
void gridpack::dynamic_simulation::DSFullBus::scatterInjectionLoad(double loadP, double loadQ){
	p_bscatterinjload_flag = true;
	p_scatterinjload_p = loadP;
	p_scatterinjload_q = loadQ;
	
}

/**
 * execute load scattering, the P and Q values of the STATIC load at certain buses will be changed to the values of 
 * the loadP and loadQ, this function keeps the Y load component of the bus still at the bus, while only compenstates the difference
*/
void gridpack::dynamic_simulation::DSFullBus::scatterInjectionLoad_compensateY(double loadP, double loadQ){
	p_bscatterinjload_flag_compensateY = true;
	p_scatterinjload_p = loadP;
	p_scatterinjload_q = loadQ;
	
}

/**
 * execute load scattering, the current values of the STATIC load at certain buses will be changed to the values of 
 * the curR and curI
*/
void gridpack::dynamic_simulation::DSFullBus::scatterInjectionLoadConstCurrent(double curR, double curI){
	p_bscatterinjloadconstcur_flag = true;
	p_scatterinjload_constcur_r = curR;
	p_scatterinjload_constcur_i = curI;
	
}

/**
 * apply load shedding for the loads in this bus
 */
void gridpack::dynamic_simulation::DSFullBus::applyLoadShedding(std::string loadid, double percentage){
	
	int iload, nload;
	nload = p_loadmodels.size();
	for (iload=0; iload<nload; iload++){
		p_loadmodels[iload]->changeLoad(percentage);
	}
}

/**
 * execute constant Y load shedding 	 
 * percentage: float load shed percentage, for example -0.2 means shed 20%
 */
void gridpack::dynamic_simulation::DSFullBus::applyConstYLoadShedding(double percentage ){
	p_bconstYLoadSheddingFlag = true;
	remainConstYLoadPerc = remainConstYLoadPerc + percentage;
	if (remainConstYLoadPerc<= 0.0){
		remainConstYLoadPerc = 0.0;
	}
	
	//printf("----------renke debug, DSFullBus::applyConstYLoadShedding, bus %d remaining load perc: %f, p_loadimpedancer: %f, p_loadimpedancei: %f\n", 
	//getOriginalIndex(), remainConstYLoadPerc, p_loadimpedancer, p_loadimpedancei);
	
}

/**
 * set the wide area control signals of the PSS of a certain generator
 * input bus_number: generator bus number
 * input bus_number: generator gen ID
 * input wideAreaControlSignal:  wide area control signal for the PSS of the generator
 */
void gridpack::dynamic_simulation::DSFullBus::setWideAreaControlSignal(std::string genid, double wideAreaControlSignal){
	
	int igen, ngen;
	gridpack::utility::StringUtils util;
	std::string clean_genid;
	ngen = p_generators.size();
	clean_genid = util.clean2Char(genid);
	
	for (igen=0; igen<ngen; igen++){
				
		//printf("----------renke debug, DSFullBus::applyGFIAdjustment, generator at bus %d with p_genid---%s---, search genid ---%s--, clean_genid --%s---\n", getOriginalIndex(), p_genid[igen].c_str(), genid.c_str(),clean_genid.c_str());
		if (clean_genid == p_genid[igen]){
			//first check whether this generator is GFI?????
			//printf("----------renke debug, DSFullBus::applyGeneratorTripping, find generator at bus %d with genid %s \n", getOriginalIndex(), genid.c_str());
			p_generators[igen]->setWideAreaFreqforPSS(wideAreaControlSignal);
		}
	}
	
}

/**
 * execute Grid Forming Inverter control parameters adjustment at this bus
 * input controlTyp: 0: GFI mp adjust; 1: GFI mq adjust; 2: GFI Pset adjust; 3: GFI Qset adjust; others: invalid
 * input bus_number: GFI gen ID
 * input newParValScaletoOrg: GFI new parameter scale factor to the very initial parameter value at the begining of dynamic simulation
 */
void gridpack::dynamic_simulation::DSFullBus::applyGFIAdjustment(int controlType, std::string genid, double newParValScaletoOrg){
	
	int igen, ngen;
	gridpack::utility::StringUtils util;
	std::string clean_genid;
	ngen = p_generators.size();
	for (igen=0; igen<ngen; igen++){
		
		clean_genid = util.clean2Char(genid);
		//printf("----------renke debug, DSFullBus::applyGFIAdjustment, generator at bus %d with p_genid---%s---, search genid ---%s--, clean_genid --%s---\n", getOriginalIndex(), p_genid[igen].c_str(), genid.c_str(),clean_genid.c_str());
		if (clean_genid == p_genid[igen]){
			//first check whether this generator is GFI?????
			//printf("----------renke debug, DSFullBus::applyGeneratorTripping, find generator at bus %d with genid %s \n", getOriginalIndex(), genid.c_str());
			p_generators[igen]->applyGeneratorParAdjustment(controlType, newParValScaletoOrg);
		}
	}
	
}

/**
 * apply generator tripping for the specific generator with genid in this bus
 */
void gridpack::dynamic_simulation::DSFullBus::applyGeneratorTripping(std::string genid){
	
	int igen, ngen;
	gridpack::utility::StringUtils util;
	std::string clean_genid;
	ngen = p_generators.size();
	for (igen=0; igen<ngen; igen++){
		clean_genid = util.clean2Char(genid);
		if (clean_genid == p_genid[igen]){
			
			//printf("----------renke debug, DSFullBus::applyGeneratorTripping, find generator at bus %d with genid %s \n", getOriginalIndex(), genid.c_str());
			p_generators[igen]->tripGenerator();
		}
	}
}

/**
 * return the fraction online from dynamic load model
 * @param idx index of dynamic load model
 * @return fraction of dynamic load that is online
 */
double gridpack::dynamic_simulation::DSFullBus::getOnlineLoadFraction(int idx)
{
  double ret = 1.0;
  if (idx >=0 && idx < p_loadmodels.size()) {
    ret = p_loadmodels[idx]->getFonline();
  } else {
    printf("No dynamic load model for index %d on bus %d\n",
        idx,getOriginalIndex());
  }
  return ret;
}
/**
 * return the total load on the bus
 * @param total_p total real power load on bus
 * @param total_q total reactive power load on bus
 */

void gridpack::dynamic_simulation::DSFullBus::getTotalLoadPower(double &total_p,
    double &total_q) const
{
  total_p = 0.0;
  total_q = 0.0;
  int i;
//#if 0
  for (i=0; i<p_npowerflow_load; i++) {
    total_p += p_powerflowload_p[i];
    total_q += p_powerflowload_q[i];
  }
  
//#else
  //total_p = p_pl;
  //total_q = p_ql;
//#endif
}

/**
 * return the power generated on the bus
 * @param total_p total active power generated on bus
 * @param total_q total reactive power generated on bus
 */
void gridpack::dynamic_simulation::DSFullBus::getTotalGeneratorPower(double &total_p,
    double &total_q) const
{
  total_p = 0.0;
  total_q = 0.0;
  int i;
  if (p_ngen != p_pg.size() || p_ngen != p_qg.size()) {
    printf("Mismatched generator sizes p_ngen: %d psize: %d qsize: %d\n",
        p_ngen,p_pg.size(),p_qg.size());
  }
  for (i=0; i<p_ngen; i++) {
    total_p += p_pg[i];
    total_q += p_qg[i];
  }
  
  total_p *= p_sbase;
  total_q *= p_sbase;
  
}

/**
 * return the real and reactive power produced by the generator indicated by
 * the tag variable
 * @param tag 2-character identifier for the generator
 * @param pg real power produced by generator
 * @param qg reactive power produced by generator
 * @return false if no generator corresponds to tag value.
 */
bool gridpack::dynamic_simulation::DSFullBus::getGeneratorPower(std::string tag,
    double &pg, double &qg) const
{
  bool ret = false;
  int i;
  pg = 0.0;
  qg = 0.0;
  for (i=0; i<p_ngen; i++) {
    if (tag == p_genid[i]) {
      ret = true;
      pg = p_pg[i]*p_sbase;
      qg = p_qg[i]*p_sbase;
      break;
    }
  }
  return ret;
}

/**
   update the diag value contributions on line status change
   @param : ybr_self - contribution for line status change
**/
void gridpack::dynamic_simulation::DSFullBus::diagValuesInsertForLineStatusChange(gridpack::ComplexType Ybr_self)
{
  p_line_status_change = true;
  p_yii = Ybr_self;
}


/**
 *  Simple constructor
 */
gridpack::dynamic_simulation::DSFullBranch::DSFullBranch(void)
{
  p_reactance.clear();
  p_resistance.clear();
  p_tap_ratio.clear();
  p_phase_shift.clear();
  p_charging.clear();
  p_shunt_admt_g1.clear();
  p_shunt_admt_b1.clear();
  p_shunt_admt_g2.clear();
  p_shunt_admt_b2.clear();
  p_xform.clear();
  p_shunt.clear();
  p_switched.clear();
  p_branch_status.clear();
  p_elems = 0;
  p_theta = 0.0;
  p_sbase = 0.0;
  p_mode = YBUS;
  p_event = false;
  p_branchrelaytripflag = false;
  p_branchactiontripflag = false;
  p_bextendedloadbranch = -1;

  p_line_status_change = false;
}

/**
 *  Simple destructor
 */
gridpack::dynamic_simulation::DSFullBranch::~DSFullBranch(void)
{
}

/**
 * Return size of off-diagonal matrix block contributed by the component
 * for the forward/reverse directions
 * @param isize, jsize: number of rows and columns of matrix block
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::dynamic_simulation::DSFullBranch::matrixForwardSize(int *isize, int *jsize) const
{
  if (p_mode == YBUS || p_mode == YL || p_mode == PG || p_mode == onFY || p_mode == posFY || p_mode == branch_trip_action
  || p_mode == jxd || p_mode == YDYNLOAD ||p_mode == bus_relay || p_mode == branch_relay || p_mode == LINESTATUSCHANGE || p_mode == GENSTATUSCHANGE) { 
    return YMBranch::matrixForwardSize(isize,jsize);
  } else {
    return false;
  }
}
bool gridpack::dynamic_simulation::DSFullBranch::matrixReverseSize(int *isize, int *jsize) const
{
  if (p_mode == YBUS || p_mode == YL || p_mode == PG || p_mode == onFY || p_mode == posFY || p_mode == branch_trip_action
  || p_mode == jxd || p_mode == YDYNLOAD || p_mode == bus_relay || p_mode == branch_relay || p_mode == LINESTATUSCHANGE) { 
    return YMBranch::matrixReverseSize(isize,jsize);
  } else {
    return false;
  }
}

/**
 * Return the values of the off-diagonal matrix block. The values are
 * returned in row-major order
 * @param values: pointer to matrix block values
 * @return: false if network component does not contribute matrix element
 */
bool gridpack::dynamic_simulation::DSFullBranch::matrixForwardValues(ComplexType *values)
{
  if (p_mode == YBUS || p_mode == YL || p_mode == PG || p_mode == jxd || p_mode == YDYNLOAD) {
	  
	//bool bstatus = YMBranch::matrixForwardValues(values);
	//gridpack::dynamic_simulation::DSFullBus *bus1 =
    //dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get());
    //gridpack::dynamic_simulation::DSFullBus *bus2 =
    //dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get());
	//if (bstatus) printf("DSFullBranch::matrixForwardValues(), Ybranch from Bus %d to Bus %d, %12.6f +j* %12.6f, \n", bus1->getOriginalIndex(), 
	//   bus2->getOriginalIndex(), real(values[0]), imag(values[0]));
	//return bstatus;  
	
    return YMBranch::matrixForwardValues(values);
  } else if(p_mode == LINESTATUSCHANGE) {
    if(p_line_status_change) {
      values[0] = p_yft;
      return true;
    } else return false;
  } else if (p_mode == posFY) {
    if (p_event) {
      values[0] = -getUpdateFactor();
      return true;
    } else {
      return false;
    }
  } else if (p_mode == branch_relay) {
	  if (p_branchrelaytripflag) {
      printf("matrix off diag forward element changes due to branch relay trip!\n");
      values[0] = -getBranchRelayTripUpdateFactor();
	  printf("changed value: %f + j*%f\n", real(values[0]), imag(values[0]));
	      
      return true;
    } else {
      return false;
    } 
  }else if (p_mode == branch_trip_action) {
		if (p_branchactiontripflag) {
			//printf("matrix off diag forward element changes due to branch relay trip!\n");
			values[0] = -getBranchTripActionUpdateFactorForward();
			//printf("changed value: %f + j*%f\n", real(values[0]), imag(values[0]));
	      
			return true;
		} else {
			return false;
		} 
  }else {
    return false;
  }
  /* if (p_mode == relay) {
 * loop through all relays and check:
 * if setrelaystatus == true
   update p_ybusr_frwd;
   update  p_ybusi_frwd;
   else 
     return false;
  */

}

bool gridpack::dynamic_simulation::DSFullBranch::matrixReverseValues(ComplexType *values)
{
  if (p_mode == YBUS || p_mode == YL || p_mode == PG || p_mode == jxd || p_mode == YDYNLOAD) {
	
	// bool bstatus = YMBranch::matrixReverseValues(values);
	// gridpack::dynamic_simulation::DSFullBus *bus1 =
    // dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get());
    // gridpack::dynamic_simulation::DSFullBus *bus2 =
    // dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get());
	// if (bstatus) printf("DSFullBranch::matrixReverseValues(), Ybranch from Bus %d to Bus %d, %12.6f +j* %12.6f, \n", bus1->getOriginalIndex(), 
	//    bus2->getOriginalIndex(), real(values[0]), imag(values[0]));
	// return bstatus;  
	  
    return YMBranch::matrixReverseValues(values);
  } else if (p_mode == LINESTATUSCHANGE) {
    if(p_line_status_change) {
      values[0] = p_ytf;
      p_line_status_change = false; // Clear line status change flag
      return true;
    } else return false;
  } else if (p_mode == posFY) {
    if (p_event) {
      values[0] = -getUpdateFactor();
      return true;
    } else {
      return false;
    }
  } else if (p_mode == branch_relay) {
	  if (p_branchrelaytripflag) {
      printf("matrix off diag reverse element changes due to branch relay trip!\n");
      values[0] = -getBranchRelayTripUpdateFactor();
	  printf("changed value: %f + j*%f\n", real(values[0]), imag(values[0]));
      return true;
    } else {
      return false;
    } 
  }else if (p_mode == branch_trip_action) {
		if (p_branchactiontripflag) {
			//printf("matrix off diag forward element changes due to branch relay trip!\n");
			values[0] = -getBranchTripActionUpdateFactorReverse();
			//printf("changed value: %f + j*%f\n", real(values[0]), imag(values[0]));
	      
			return true;
		} else {
			return false;
		} 
  }else {
    return false;
  }
}

// Calculate contributions to the admittance matrix from the branches
void gridpack::dynamic_simulation::DSFullBranch::setYBus(void)
{
  YMBranch::setYBus();
  gridpack::ComplexType ret;
  ret = YMBranch::getForwardYBus();
  p_ybusr_frwd = real(ret);
  p_ybusi_frwd = imag(ret);
  ret = YMBranch::getReverseYBus();
  p_ybusr_rvrs = real(ret);
  p_ybusi_rvrs = imag(ret);  
  // Not really a contribution to the admittance matrix but might as well
  // calculate phase angle difference between buses at each end of branch
  gridpack::dynamic_simulation::DSFullBus *bus1 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get());
  gridpack::dynamic_simulation::DSFullBus *bus2 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get());
//  if (p_xform) {
//    printf ("from %d-> to %d: p_phase_shift = %f, a = %f+%fi\n", bus1->getOriginalIndex(), bus2->getOriginalIndex(), p_phase_shift, real(a), imag(a) );
//  }
  //p_theta = bus1->getPhase() - bus2->getPhase();
  double pi = 4.0*atan(1.0);
  p_theta = (bus1->getPhase() - bus2->getPhase());
  //printf("p_phase_shift: %12.6f\n",p_phase_shift);
  //printf("p_theta: %12.6f\n",p_theta);
  //printf("p_tap_ratio: %12.6f\n",p_tap_ratio);
 
}

/**
 * Load values stored in DataCollection object into DSFullBranch object. The
 * DataCollection object will have been filled when the network was created
 * from an external configuration file
 * @param data: DataCollection object contain parameters relevant to this
 *       branch that were read in when network was initialized
 */
void gridpack::dynamic_simulation::DSFullBranch::load(
  const boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  YMBranch::load(data);

  // This function may be called more than once so clear all vectors
  p_reactance.clear();
  p_resistance.clear();
  p_tap_ratio.clear();
  p_phase_shift.clear();
  p_charging.clear();
  p_shunt_admt_g1.clear();
  p_shunt_admt_b1.clear();
  p_shunt_admt_g2.clear();
  p_shunt_admt_b2.clear();
  p_xform.clear();
  p_shunt.clear();
  p_switched.clear();
  p_branch_status.clear();
  p_linerelays.clear();
  p_relaybranchidx.clear();
  p_ckt.clear();
 
  bool bdebug_branch_load = false;
  if (bdebug_branch_load) printf("entering DSFullBranch::load() \n");

  /*
  gridpack::dynamic_simulation::DSFullBus *bus1 =
  dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get());
  printf("DSFullBranch::load() get bus 1 number\n");	
  gridpack::dynamic_simulation::DSFullBus *bus2 =
  dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get());
  printf("DSFullBranch::load() get bus 2 number\n");		

  printf("DSFullBranch::load(), Bus No.: %d to Bus No.: %d \n",
          bus1->getOriginalIndex(), bus2->getOriginalIndex());
  */

  //read line type, check this line is added by composite load model
  std::string snewbratype;

  if (data->getValue(NEW_BRANCH_TYPE, &snewbratype)){
    if ( snewbratype=="TRANSFORMER" || snewbratype=="FEEDER" ) {
      printf("This branch is a extended bus by composite load models, type: %s \n",
          snewbratype.c_str());
      if ( snewbratype=="TRANSFORMER" )
      {
        p_bextendedloadbranch = 1; 
        return;
      }else{
        p_bextendedloadbranch = 2; 
        return;
      }
    }	  
  }

  bool ok = true;
  data->getValue(BRANCH_NUM_ELEMENTS, &p_elems);
  double rvar;
  int ivar;
  bool lvar;
  double pi = 4.0*atan(1.0);
  p_active = false;
  ok = data->getValue(CASE_SBASE, &p_sbase);
  int idx;

  //renke add
  int nrelay, irelay, relaycnt; 
  std::string sckt, srelay_lineckt, smodel;
  RelayFactory relayFactory;

  nrelay = 0;
  data->getValue(RELAY_NUMBER, &nrelay); //renke add, get number of relays with this branch
  //printf ("-- component branch line relay:  -- nrelay = %d \n", nrelay);

  for (idx = 0; idx<p_elems; idx++) {
    data->getValue(BRANCH_X, &rvar, idx);
    p_reactance.push_back(rvar);
    data->getValue(BRANCH_R, &rvar, idx);
    p_resistance.push_back(rvar);
    rvar = 0.0;
    data->getValue(BRANCH_SHIFT, &rvar, idx);
    rvar = -rvar*pi/180.0;
    p_phase_shift.push_back(rvar);
    rvar = 0.0;
    data->getValue(BRANCH_TAP, &rvar, idx);
    p_tap_ratio.push_back(rvar);
    if (rvar != 0.0) {
      p_xform.push_back(true);
    } else {
      p_xform.push_back(false);
    }
    ivar = 1;
    data->getValue(BRANCH_STATUS, &ivar, idx);
    p_branch_status.push_back(ivar);
    if (ivar == 1) p_active = true;
	ok = data->getValue(BRANCH_SWITCHED, &lvar, idx);
    if (!ok) lvar = false;
    p_switched.push_back(lvar);
    bool shunt = true;
    shunt = shunt && data->getValue(BRANCH_B, &rvar, idx);
    p_charging.push_back(rvar);
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_G1, &rvar, idx);
    p_shunt_admt_g1.push_back(rvar);
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_B1, &rvar, idx);
    p_shunt_admt_b1.push_back(rvar);
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_G2, &rvar, idx);
    p_shunt_admt_g2.push_back(rvar);
    shunt = shunt && data->getValue(BRANCH_SHUNT_ADMTTNC_B2, &rvar, idx);
    p_shunt_admt_b2.push_back(rvar);
    p_shunt.push_back(shunt);

    //renke add, get line relays associated with this branch object
    if (data->getValue(BRANCH_CKT, &sckt, idx)) {
      p_ckt.push_back(sckt);
      //printf("branch %d element, ckt %s \n", idx, sckt.c_str());
    }

    if (nrelay>0) {
      for (irelay=0 ; irelay<nrelay ; irelay++) {
        data->getValue(RELAY_ID, &srelay_lineckt, irelay);
        data->getValue(RELAY_MODEL, &smodel, irelay);
        //printf("branch relay irelay: %d, model: %s \n",
        //         irelay, smodel.c_str());
        //printf("branch relay irelay: %d, ckt: %s \n",
        //         irelay, srelay_lineckt.c_str());
        if (srelay_lineckt == sckt && smodel == "DISTR1") {
          printf("find a distr1 relay with ckt %s \n", sckt.c_str());
          p_relaybranchidx.push_back(idx);

          BaseRelayModel *relaymodel
            = relayFactory.createRelayModel(smodel);
          boost::shared_ptr<BaseRelayModel> relay;
          relay.reset(relaymodel);
          relay->load(data, irelay);
          p_linerelays.push_back(relay);	
        }  
      }
    }
  }
}

/**
 * Evaluate branch flows for the to and from bus on the branch
 */
void gridpack::dynamic_simulation::DSFullBranch::evaluateBranchFlow()
{
  int i;
  double pi = 4.0*atan(1.0);

  gridpack::dynamic_simulation::DSFullBus *bus1 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get());
  gridpack::dynamic_simulation::DSFullBus *bus2 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get());

  //get bus voltages
  p_branchfrombusvolt = bus1->getComplexVoltage();
  p_branchtobusvolt = bus2->getComplexVoltage();

  // printf ("Branch volts bus1, %8.4f+%8.4fj,  bus 2, %8.4f+%8.4fj,\n", real(p_branchfrombusvolt), 
  // 		imag(p_branchfrombusvolt), real(p_branchtobusvolt),imag(p_branchtobusvolt) );

  // define branch from and to bus p + jq
  gridpack::ComplexType c_branchfrombuspq, c_branchtobuspq;

  if (!p_branchfrombuspq.empty()){
    p_branchfrombuspq.clear();
  }

  if (!p_branchtobuspq.empty()){
    p_branchtobuspq.clear();
  }

  for ( i=0 ; i<p_elems ; i++ ) {

    // printf ("Branch %d impedance, %8.4f+%8.4fj \n", i, p_resistance[i], p_reactance[i]); 
    // printf ("Branch %d shunt susceptance, %8.4fj \n", i, p_charging[i]); 
    // printf ("Branch %d phase shift, %8.4f \n", i, p_phase_shift[i]); 
    // printf ("Branch %d tap ratio, %8.4f \n", i, p_tap_ratio[i]); 

    if (p_xform[i]) {
      // For transformer, we don't consider charging susceptance
      // TODO: In PSS/E DataFormats.pdf, phase shift value is in degree. 
      // However, in getTransformer inside this script file, it seems cos and sin take degree as input. Need to confirm.

      // calculate complex turn ratio, this part of code is referred to getTransformer
      gridpack::ComplexType a(cos(p_phase_shift[i]),sin(p_phase_shift[i]));
      gridpack::ComplexType c_xformsecvolt;
      a = p_tap_ratio[i]*a;

      // calculate secondary side voltage of idea transformer
      c_xformsecvolt = p_branchfrombusvolt / a;

      // calculate secondary side current and primary side current
      gridpack::ComplexType c_xformpricur, c_xformseccur, c_xformz;
      c_xformz = gridpack::ComplexType(p_resistance[i], p_reactance[i]);
      c_xformseccur = (c_xformsecvolt - p_branchtobusvolt)/c_xformz;
      c_xformpricur = c_xformseccur / a; // assume pri and sec currents of transformer are in the same direction

      // calculate primary and secondary side complex power
      c_branchfrombuspq = p_branchfrombusvolt * conj(c_xformpricur);
      c_branchtobuspq = p_branchtobusvolt * conj(c_xformseccur);

      // store from and to buses complex power into vectors
      p_branchfrombuspq.push_back(c_branchfrombuspq); // active power P and reactive power Q at "from" bus (P_from + j Q_from)
      p_branchtobuspq.push_back(c_branchtobuspq); // active power P and reactive power Q at "to" bus (P_to + j Q_to)

    } else {
      // For transmission line, we consider charging susceptance
      gridpack::ComplexType c_Z, c_branchcurr;
      gridpack::ComplexType c_Ysh, c_frombusshuntcurr, c_tobusshuntcurr;

      // calculate line impedance and charging susceptance (on either end of the line)
      c_Z = gridpack::ComplexType(p_resistance[i], p_reactance[i]);
      c_Ysh = gridpack::ComplexType(0, p_charging[i]/2.0); // shunt susceptance at either end of a branch

      // calculate line current
      c_branchcurr = (p_branchfrombusvolt - p_branchtobusvolt)/c_Z;
  
      // printf ("Branch %d element current, %8.4f+%8.4fj \n", i, real(c_branchcurr), 
      // 	imag(c_branchcurr) );

      // calculate currents flowing through charging susceptance on from and to ends
      c_frombusshuntcurr = c_Ysh * p_branchfrombusvolt;
      c_tobusshuntcurr = c_Ysh * p_branchtobusvolt;

      // calculate the p + jq at from and to buses
      c_branchfrombuspq = p_branchfrombusvolt * conj(c_branchcurr + c_frombusshuntcurr);
      c_branchtobuspq = p_branchtobusvolt * conj(c_branchcurr - c_tobusshuntcurr);

      // store from and to buses complex power into vectors
      p_branchfrombuspq.push_back(c_branchfrombuspq); // active power P and reactive power Q at "from" bus (P_from + j Q_from)
      p_branchtobuspq.push_back(c_branchtobuspq); // active power P and reactive power Q at "to" bus (P_to + j Q_to)

      // printf ("Branch %d element from bus apparent power, %8.4f+%8.4fj \n", i, real(p_branchfrombuspq[i]), 
      // 	imag(p_branchfrombuspq[i]) ); 

      // printf ("Branch %d element to bus apparent power, %8.4f+%8.4fj \n", i, real(p_branchtobuspq[i]), 
      // 	imag(p_branchtobuspq[i]) );
    }

  }
}

/**
 * Update data collection object with current values from simulation
 * @param data: DataCollection object containing parameters for this branch
 */
void gridpack::dynamic_simulation::DSFullBranch::updateData(
    boost::shared_ptr<gridpack::component::DataCollection> &data)
{
  int i;
  evaluateBranchFlow();
  for (i=0; i<p_elems; i++) {
    // Here, it does not need to differentiate transformers or lines for storing the variables.
    // Treating both transformers and lines as branches
    double pf = real(p_branchfrombuspq[i]);
    double qf = imag(p_branchfrombuspq[i]);
    double pt = real(p_branchtobuspq[i]);
    double qt = imag(p_branchtobuspq[i]);
    if (!data->setValue(BRANCH_FROM_P_CURRENT, pf, i)) {
      data->addValue(BRANCH_FROM_P_CURRENT, pf, i);
    }
    if (!data->setValue(BRANCH_FROM_Q_CURRENT, qf, i)) {
      data->addValue(BRANCH_FROM_Q_CURRENT, qf, i);
    }
    if (!data->setValue(BRANCH_TO_P_CURRENT, pt, i)) {
      data->addValue(BRANCH_TO_P_CURRENT, pt, i);
    }
    if (!data->setValue(BRANCH_TO_Q_CURRENT, qt, i)) {
      data->addValue(BRANCH_TO_Q_CURRENT, qt, i);
    }
    
  }
}

/**
 * update branch current
 */
void gridpack::dynamic_simulation::DSFullBranch::updateBranchCurrent() //renke add
{
	int i;
	double dbranchR, dbranchX;
	
	gridpack::dynamic_simulation::DSFullBus *bus1 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get());
    gridpack::dynamic_simulation::DSFullBus *bus2 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get());
	
	//get bus voltages
	gridpack::ComplexType c_Z, c_branchcurr;
	p_branchfrombusvolt = bus1->getComplexVoltage();
	p_branchtobusvolt = bus2->getComplexVoltage();
	
	//printf ("Branch volts bus1, %8.4f+%8.4fj,  bus 2, %8.4f+%8.4fj,\n", real(p_branchfrombusvolt), 
	//		imag(p_branchfrombusvolt), real(p_branchtobusvolt),imag(p_branchtobusvolt) );
	
	if (!p_branchcurrent.empty()){
		p_branchcurrent.clear();
	}
	
	for ( i=0 ; i<p_elems ; i++ ) {
		dbranchR = p_resistance[i];
		dbranchX = p_reactance[i];
		c_Z = gridpack::ComplexType(dbranchR, dbranchX);
		c_branchcurr = (p_branchfrombusvolt - p_branchtobusvolt)/c_Z;
		p_branchcurrent.push_back(c_branchcurr);
		//printf ("Branch %d element current, %8.4f+%8.4fj \n", real(p_branchcurrent[i]), 
		//	imag(p_branchcurrent[i]), i );
	}
	
}

/**
* update the relay status associate with this branch
*/
bool gridpack::dynamic_simulation::DSFullBranch::updateRelay(bool flag, double delta_t) //renke add
{
	int irelay, nrelay, itrip, itrip_prev, ibranch;
	double dbusvoltfreq;
	bool bbranchflag;
	gridpack::ComplexType cbusfreq;
	boost::shared_ptr<gridpack::dynamic_simulation::BaseRelayModel> p_relay;
	std::vector<gridpack::ComplexType*> vrelayvalue;
	bbranchflag = false;
	p_newtripbranchcktidx.clear();
	gridpack::dynamic_simulation::DSFullBus *bus1 =
      dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get());
    gridpack::dynamic_simulation::DSFullBus *bus2 =
      dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get());
	bus1->setBranchRelayFromBusStatus(false);
	bus1->setBranchRelayToBusStatus(false);
	bus2->setBranchRelayFromBusStatus(false);
	bus2->setBranchRelayToBusStatus(false);
	bus1->clearRelayTrippedbranch();
	bus2->clearRelayTrippedbranch();
	int idxbus1 = getBus1OriginalIndex();
    int idxbus2 = getBus2OriginalIndex();
	
	//update line relays
	
	if (!p_linerelays.empty()) {
		nrelay = p_linerelays.size();
		updateBranchCurrent();
		
		for ( irelay=0 ; irelay<nrelay ; irelay++ ){
			
			itrip = 0;
			itrip_prev = 0;
			vrelayvalue.clear();
			vrelayvalue.push_back( &p_branchfrombusvolt ); //make sure the volt is at the from bus
			ibranch = p_relaybranchidx[irelay];
			vrelayvalue.push_back( &(p_branchcurrent[ibranch]) );
			p_linerelays[irelay]->setMonitorVariables(vrelayvalue);
			p_linerelays[irelay]->updateRelay(delta_t);
			p_linerelays[irelay]->getTripStatus( itrip, itrip_prev );
			/*
			printf(" from bus volt = %8.4f + j*%8.4f ; branch current = %3.6f + j*%3.6f\n", 
				real(p_branchfrombusvolt), imag(p_branchfrombusvolt), 
				real(p_branchcurrent[ibranch]), imag(p_branchcurrent[ibranch]));
			*/
			printf(" DSFullBranch::updateRelay DISTR1 itrip = %d, itrip_prev = %d \n", itrip, itrip_prev);
			if ( itrip==1 && itrip_prev==0 && p_linerelays[irelay]->getOperationStatus()) {
				bbranchflag = true;
				p_linerelays[irelay]->setOperationStatus(false);
								
				bus1->setBranchRelayFromBusStatus(true);
				bus2->setBranchRelayToBusStatus(true);
				
				// add more code handling the branch related Y matrix?? Shuangshuang tbd
				printf(" DSFullBranch::updateRelay DISTR1 trip!!!  \n");
                                p_branch_status[ibranch] = 0;
								p_newtripbranchcktidx.push_back(ibranch);
								bus1->setRelayTrippedbranch(this);
								bus2->setRelayTrippedbranch(this);
								
                                //setLineStatus(p_ckt[ibranch], false);
			}
		}
	}
	
	p_branchrelaytripflag = bbranchflag;
	return bbranchflag;
}

/**
 * Set the mode to control what matrices and vectors are built when using
 * the mapper
 * @param mode: enumerated constant for different modes
 */
void gridpack::dynamic_simulation::DSFullBranch::setMode(int mode)
{
  if (mode == YBUS || mode == YL || mode == PG || mode == jxd || mode == YDYNLOAD) {
    YMBranch::setMode(gridpack::ymatrix::YBus);
  }
  p_mode = mode;
}

/**
 * Return the complex admittance of the branch
 * @return: complex addmittance of branch
 */
gridpack::ComplexType gridpack::dynamic_simulation::DSFullBranch::getAdmittance(void)
{
  int i;
  gridpack::ComplexType ret(0.0,0.0);
  for (i=0; i<p_elems; i++) {
    gridpack::ComplexType tmp(p_resistance[i], p_reactance[i]);
    if (!p_xform[i] && p_branch_status[i] == 1) {
      tmp = -1.0/tmp;
    } else {
      tmp = gridpack::ComplexType(0.0,0.0);
    }
    ret += tmp;
  }
  return ret;
}

/**
 * Return transformer contribution from the branch to the calling
 * bus
 * @param bus: pointer to the bus making the call
 * @return: contribution to Y matrix from branch
 */
gridpack::ComplexType
gridpack::dynamic_simulation::DSFullBranch::getTransformer(gridpack::dynamic_simulation::DSFullBus *bus)
{
  int i;
  gridpack::ComplexType ret(0.0,0.0);
  for (i=0; i<p_elems; i++) {
    gridpack::ComplexType tmp(p_resistance[i],p_reactance[i]);
    gridpack::ComplexType tmpB(0.0,0.5*p_charging[i]);
    if (p_xform[i] && p_branch_status[i] == 1) {
      tmp = -1.0/tmp;
      tmp = tmp - tmpB;
      gridpack::ComplexType a(cos(p_phase_shift[i]),sin(p_phase_shift[i]));
      a = p_tap_ratio[i]*a;
      if (bus == getBus1().get()) {
        tmp = tmp/(conj(a)*a);
      } else if (bus == getBus2().get()) {
        // tmp is unchanged
      }
    } else {
      tmp = gridpack::ComplexType(0.0,0.0);
    }
    ret += tmp;
  }
  return ret;
}

/**
 * Return the contribution to a bus from shunts
 * @param bus: pointer to the bus making the call
 * @return: contribution to Y matrix from shunts associated with branches
 */
gridpack::ComplexType
gridpack::dynamic_simulation::DSFullBranch::getShunt(gridpack::dynamic_simulation::DSFullBus *bus)
{
  double retr, reti;
  retr = 0.0;
  reti = 0.0;
  int i;
  for (i=0; i<p_elems; i++) {
    double tmpr, tmpi;
    if (p_shunt[i] && p_branch_status[i] == 1) {
      tmpr = 0.0;
      tmpi = 0.0;
      if (!p_xform[i]) {
        tmpi = 0.5*p_charging[i];
        tmpr = 0.0;
      }
      // HACK: pointer comparison, maybe could handle this better
      if (bus == getBus1().get()) {
        tmpr += p_shunt_admt_g1[i];
        tmpi += p_shunt_admt_b1[i];
      } else if (bus == getBus2().get()) {        tmpr += p_shunt_admt_g2[i];        tmpi += p_shunt_admt_b2[i];
      } else {
        // TODO: Some kind of error
      }
    } else { 
      tmpr = 0.0;
      tmpi = 0.0;
    }
    retr += tmpr;
    reti += tmpi;
  }
  return gridpack::ComplexType(retr,reti);
}

gridpack::ComplexType
gridpack::dynamic_simulation::DSFullBranch::getPosfy11YbusUpdateFactor(int sw2_2, int sw3_2)
{ 
  double retr, reti;
  int i;
  gridpack::dynamic_simulation::DSFullBus *bus1 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get());
  gridpack::dynamic_simulation::DSFullBus *bus2 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get());
  if (bus1->getOriginalIndex() == sw2_2+1 && bus2->getOriginalIndex() == sw3_2+1) {
    for (i=0; i<p_elems; i++) {
      gridpack::ComplexType myValue(p_resistance[i], p_reactance[i]);
	  // tbd, have not consider the transformer ratio and line shunt capacitance
      myValue = 1.0 / myValue;
      //printf("myValue = %f+%fi\n", real(myValue), imag(myValue));
      //printf("%f %f\n", p_resistance, p_reactance);
      retr = real(myValue);
      reti = imag(myValue);
      return gridpack::ComplexType(retr, reti);
    }
  } else {
    return gridpack::ComplexType(-999.0, -999.0); // return a dummy value
  }
}

gridpack::ComplexType 
gridpack::dynamic_simulation::DSFullBranch::getUpdateFactor()
{ 
  int i;
  gridpack::ComplexType ret(0.0,0.0);
  for (i=0; i<p_elems; i++) {
    gridpack::ComplexType tmp(p_resistance[i], p_reactance[i]);
    tmp = -1.0 / tmp;
    ret += tmp;
  }
  return ret;
}

//renke add, update the contributions from the branch relay trip
gridpack::ComplexType 
gridpack::dynamic_simulation::DSFullBranch::getBranchRelayTripUpdateFactor()
{ 
  int i, idx, ntripbranch;
  gridpack::ComplexType ret(0.0,0.0);
  
  if (!p_newtripbranchcktidx.empty()) {
	  ntripbranch = p_newtripbranchcktidx.size();
	  printf("DSFullBranch::getBranchRelayTripUpdateFactor ntripbranch = %d\n", ntripbranch);
	  for (i=0; i<ntripbranch; i++ ){
		  idx = p_newtripbranchcktidx[i];
		  printf("DSFullBranch::getBranchRelayTripUpdateFactor idx = %d\n", idx);
		  gridpack::ComplexType tmp(p_resistance[idx], p_reactance[idx]);
		  // tbd, have not consider the transformer ratio and line shunt capacitance
          tmp = -1.0 / tmp;
          ret += tmp;
	  }
  }

  return ret;
}

/**
 * Return the updating factor that will be applied to the ybus matrix at
 * the branch trip action
 * @return: value of update factor
*/
gridpack::ComplexType 
gridpack::dynamic_simulation::DSFullBranch::getBranchTripActionUpdateFactorForBus(int busNo)
{ 
	int ibr, idx, ntripbranch;
	gridpack::ComplexType ret(0.0,0.0);
  
	if (!p_vec_tripaction_branchcktidx.empty()) {
		ntripbranch = p_vec_tripaction_branchcktidx.size();
		//printf("DSFullBranch::getBranchRelayTripUpdateFactor ntripbranch = %d\n", ntripbranch);
		for (ibr=0; ibr<ntripbranch; ibr++ ){
			idx = p_vec_tripaction_branchcktidx[ibr];
			//printf("DSFullBranch::getBranchRelayTripUpdateFactor idx = %d\n", idx);
			if (p_branch_status[idx]) {
				
				gridpack::ComplexType tmp(p_resistance[idx],p_reactance[idx]);
				gridpack::ComplexType tmpB(0.0,0.5*p_charging[idx]);
				tmp = -1.0/tmp;
				tmp = tmp - tmpB;
				
				if (p_xform[idx]){
					gridpack::ComplexType a(cos(p_phase_shift[idx]),sin(p_phase_shift[idx]));
					a = p_tap_ratio[idx]*a;
					if ( (!p_switched[idx] && busNo == dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get())->getOriginalIndex()) ||
						(p_switched[idx] && busNo == dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get())->getOriginalIndex()) ) {
						tmp = tmp/(conj(a)*a);
					} else {
						// tmp is unchanged
					}
	
				}

				ret += tmp; 
			}

		}
	}
		
	return ret;
}

gridpack::ComplexType 
gridpack::dynamic_simulation::DSFullBranch::getBranchTripActionUpdateFactorForward()
{ 
  int ibr, idx, ntripbranch;
  gridpack::ComplexType ret(0.0,0.0);
  
  if (!p_vec_tripaction_branchcktidx.empty()) {
	  ntripbranch = p_vec_tripaction_branchcktidx.size();
	  //printf("DSFullBranch::getBranchRelayTripUpdateFactor ntripbranch = %d\n", ntripbranch);
	  for (ibr=0; ibr<ntripbranch; ibr++ ){
		  idx = p_vec_tripaction_branchcktidx[ibr];
		  //printf("DSFullBranch::getBranchRelayTripUpdateFactor idx = %d\n", idx);
		  if (p_branch_status[idx]) {
			gridpack::ComplexType tmp(p_resistance[idx], p_reactance[idx]);		  
			tmp = -1.0 / tmp;
			// tbd, have not consider the transformer ratio and line shunt capacitance
			gridpack::ComplexType a(cos(p_phase_shift[idx]),sin(p_phase_shift[idx]));
			if (p_xform[idx]) a = p_tap_ratio[idx]*a;
			if (p_switched[idx]) a = conj(a);
			if (p_xform[idx]) {
				tmp = tmp/conj(a);
				
			}
			
			ret += tmp;
		  }
 
	  }
  }

  return ret;
}

gridpack::ComplexType 
gridpack::dynamic_simulation::DSFullBranch::getBranchTripActionUpdateFactorReverse()
{ 
  int ibr, idx, ntripbranch;
  gridpack::ComplexType ret(0.0,0.0);
  
  if (!p_vec_tripaction_branchcktidx.empty()) {
	  ntripbranch = p_vec_tripaction_branchcktidx.size();
	  //printf("DSFullBranch::getBranchRelayTripUpdateFactor ntripbranch = %d\n", ntripbranch);
	  for (ibr=0; ibr<ntripbranch; ibr++ ){
		  idx = p_vec_tripaction_branchcktidx[ibr];
		  //printf("DSFullBranch::getBranchRelayTripUpdateFactor idx = %d\n", idx);
		  if (p_branch_status[idx]) {
			gridpack::ComplexType tmp(p_resistance[idx], p_reactance[idx]);		  
			tmp = -1.0 / tmp;
			// tbd, have not consider the transformer ratio and line shunt capacitance
			gridpack::ComplexType a(cos(p_phase_shift[idx]),sin(p_phase_shift[idx]));
			if (p_xform[idx]) a = p_tap_ratio[idx]*a;
			if (p_switched[idx]) a = conj(a);
			if (p_xform[idx]) {
				tmp = tmp/a;
				
			}
			
			ret += tmp;
		  }
 
	  }
  }

  return ret;
}

/**
 * Clear fault event from branch
 */
void gridpack::dynamic_simulation::DSFullBranch::clearEvent()
{
  p_event = false;
}

/**
 * Check to see if an event applies to this branch and set appropriate internal
 * parameters
 * @param event a struct containing parameters that describe a fault event in
 * a dyanamic simulation
 */
void gridpack::dynamic_simulation::DSFullBranch::setEvent(
    const gridpack::dynamic_simulation::Event &event)
{
  int idx1 = getBus1OriginalIndex();
  int idx2 = getBus2OriginalIndex();
  
  /*
  if (event.isGenerator) {
    if (event.from_idx == event.to_idx) {
      if (event.from_idx == idx1) {
        dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
          (getBus1().get())->setEvent(idx1,idx1,this);
      } else if (event.from_idx == idx2) {
        dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
          (getBus2().get())->setEvent(idx2,idx2,this);
      }
    }
  }
  */
  
  if (event.isLine) {
    // Check to see if event refers to this bus
    if (idx1 == event.from_idx && idx2 == event.to_idx) {
      p_event = true;
    } else {
      p_event = false;
    }
    if (p_event) {
      dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (getBus1().get())->setEvent(idx1,idx2,this);
      dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (getBus2().get())->setEvent(idx1,idx2,this);
    }
  }
}

bool gridpack::dynamic_simulation::DSFullBranch::setBranchTripAction(
    const std::string ckt_tag)
{
	int ickt;
	//printf("----renke debug load shed, DSFullBranch::setBranchTripAction, ckt_tag: ,%s, \n", ckt_tag.c_str());
	for (ickt=0; ickt<p_elems; ickt++) {
		
		//printf("----renke debug load shed, DSFullBranch::setBranchTripAction, ickt: %d, p_ckt[ickt]: ,%s, status: %d \n", ickt, p_ckt[ickt].c_str(), p_branch_status[ickt]);
		
		if (ckt_tag == p_ckt[ickt]  && p_branch_status[ickt] == 1){
			//printf("----renke debug load shed, DSFullBranch::setBranchTripAction, ckt_tag: ,%s, \n", ckt_tag.c_str());
			//printf("----renke debug load shed, DSFullBranch::setBranchTripAction, ickt: %d, p_ckt[ickt]: ,%s, status: %d \n", ickt, p_ckt[ickt].c_str(), p_branch_status[ickt]);
			//printf("----renke debug load shed, DSFullBranch::setBranchTripAction, find the circuit of the branch!!!! \n");
			
			p_branchactiontripflag = true;
	
			p_vec_tripaction_branchcktidx.push_back(ickt);
	
			int idx1 = getBus1OriginalIndex();
			int idx2 = getBus2OriginalIndex();
			dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
				(getBus1().get())->setBranchTripAction(idx1,idx2,this);
			dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
				(getBus2().get())->setBranchTripAction(idx1,idx2,this);
				
			return true;
		}
	}
	
	return false;

}

bool gridpack::dynamic_simulation::DSFullBranch::setBranchTripAction( )
{
	int ickt;
	for (ickt=0; ickt<p_elems; ickt++) {
		if (!p_xform[ickt] && p_branch_status[ickt] == 1){
			p_branchactiontripflag = true;
	
			p_vec_tripaction_branchcktidx.push_back(ickt);
	
			int idx1 = getBus1OriginalIndex();
			int idx2 = getBus2OriginalIndex();
			
			//printf("---------testing, DSFullBranch::setBranchTripAction(), Tripped Branch from Bus %d to Bus %d, ckt idx: %d \n", idx1, idx2, ickt);
			dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
				(getBus1().get())->setBranchTripAction(idx1,idx2,this);
			dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
				(getBus2().get())->setBranchTripAction(idx1,idx2,this);
				
			return true;
		}
	}
	
	return false;

}

void gridpack::dynamic_simulation::DSFullBranch::clearBranchTripAction()
{
	p_branchactiontripflag = false;
	
	int ibr, nbr, idx;
	nbr = p_vec_tripaction_branchcktidx.size();
	for (ibr=0 ; ibr<nbr ; ibr++){
		idx = p_vec_tripaction_branchcktidx[ibr];
		p_branch_status[idx] = 0;
	}
	
	p_vec_tripaction_branchcktidx.clear();
	
	int idx1 = getBus1OriginalIndex();
    int idx2 = getBus2OriginalIndex();
	dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (getBus1().get())->clearBranchTripAction();
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>
        (getBus2().get())->clearBranchTripAction();
}

/**
 * Set parameters of the transformer branch due to composite load model
 */
void gridpack::dynamic_simulation::DSFullBranch::SetCmplXfmrBranch(double dx, double dtap)
{
	p_elems = 1;
	p_active = true;
	p_sbase = 100.0;
	p_reactance.push_back(dx);
	p_resistance.push_back(0.0);
	p_phase_shift.push_back(0.0);
	p_tap_ratio.push_back(dtap);
	p_xform.push_back(true);
	p_branch_status.push_back(1);
	p_shunt.push_back(false);
	p_ckt.push_back("1"); 
}

/**
 * Set parameters of the feeder branch due to composite load model
 */
void gridpack::dynamic_simulation::DSFullBranch::SetCmplFeederBranch(double dr, double dx)
{
	p_elems = 1;
	p_active = true;
	p_sbase = 100.0;
	p_reactance.push_back(dx);
	p_resistance.push_back(dr);
	p_phase_shift.push_back(0.0);
	p_tap_ratio.push_back(1.0);
	p_xform.push_back(false);
	p_branch_status.push_back(1);
	p_shunt.push_back(false);
	p_ckt.push_back("1"); 
}
/*
 * print the content of the DSFullBranch
 */
void gridpack::dynamic_simulation::DSFullBranch::printDSFullBranch()
{
	int i, branchelem;
	gridpack::dynamic_simulation::DSFullBus *bus1 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get());
  gridpack::dynamic_simulation::DSFullBus *bus2 =
    dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get());
	
	branchelem = p_reactance.size();
	printf("DSFullBranch::printDSFullBranch(), Branch from Bus %d (isolated: %d) to Bus %d (isolated: %d), branchelem: %d \n", bus1->getOriginalIndex(), bus1->checkisolated(), 
		bus2->getOriginalIndex(), bus2->checkisolated(), branchelem);
	for ( i=0; i<p_elems; i++ ){
		printf("DSFullBranch::printDSFullBranch(), %d-th elem: cktid: %s, p_active: %d, p_reactance: %f, p_resistance: %f, p_phase_shift: %f, p_tap_ratio: %f, p_xform: %d, p_branch_status: %d, p_shunt: %d \n", 
                       i, p_ckt[i].c_str(), p_active, p_reactance[i], p_resistance[i], p_phase_shift[i], p_tap_ratio[i], static_cast<int>(p_xform[i]), p_branch_status[i], static_cast<int>(p_shunt[i]));
	}	
}

/**
 * check the type of the extended load branch type variable: p_bextendedloadbranch
 */
int gridpack::dynamic_simulation::DSFullBranch::checkExtendedLoadBranchType(void)
{
	return p_bextendedloadbranch;
}

/**
 * Return contributions to Y-matrix from a specific transmission element
 * @param tag character string for transmission element
 * @param Yjj contribution at "to" bus
 * @param Yji contribution for ji_th Y-matrix element
 */
void gridpack::dynamic_simulation::DSFullBranch::getRvrsLineElements(const std::string tag,
   gridpack::ComplexType *Yjj, gridpack::ComplexType *Yji)
{
  YMBranch::getRvrsLineElements(tag, Yjj, Yji);
}

/**
 * Return contributions to Y-matrix from a specific transmission element
 * @param tag character string for transmission element
 * @param Yii contribution at "from" bus
 * @param Yij contribution for ij_th Y-matrix element
 */
void gridpack::dynamic_simulation::DSFullBranch::getLineElements(const std::string tag,
   gridpack::ComplexType *Yii, gridpack::ComplexType *Yij)
{
  YMBranch::getLineElements(tag, Yii, Yij);
}

/**
   setLineStatus - Sets the line status and updates the associated
   branch and bus objects. 
   
   @param: ckt_id - circuit id
   @param: status - new line status
   
   Note: This method is used to
   update the branch status and update the bus/branch
   objects. It sets up values in the bus and branch objects
   so that incrementMatrix method called on the network Ybus
   uses these values to remove the branch contributions from
   the Y-bus matrix
**/
void gridpack::dynamic_simulation::DSFullBranch::setLineStatus(std::string id, int status)
{
  std::string ckt_id;
  gridpack::utility::StringUtils util;
  ckt_id = util.clean2Char(id);

  bool line_status = getLineStatus(ckt_id);
  if(line_status != (bool)status) {
    
    gridpack::ComplexType yii,yij,yji,yjj;
    getLineElements(ckt_id,&yii,&yij);
    getRvrsLineElements(ckt_id,&yjj,&yji);
    gridpack::dynamic_simulation::DSFullBus *bus1
      = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus1().get());
    gridpack::dynamic_simulation::DSFullBus *bus2
      = dynamic_cast<gridpack::dynamic_simulation::DSFullBus*>(getBus2().get());

    if(!status) {
      /* Negative sign for removing these values */
      p_yft = -yij;
      p_ytf = -yji;
      bus1->diagValuesInsertForLineStatusChange(-yii);
      bus2->diagValuesInsertForLineStatusChange(-yjj);
    } else {
      p_yft = yij;
      p_ytf = yji;
      bus1->diagValuesInsertForLineStatusChange(yii);
      bus2->diagValuesInsertForLineStatusChange(yjj);
    }

    p_line_status_change = true;

    YMBranch::setLineStatus(ckt_id,status);
  }
}
  


