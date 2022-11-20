/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: May 10, 2021
 *      Author: Renke Huang
 */
#ifndef REECA1_HPP	
#define REECA1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Reeca1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Reeca1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Reeca1Parser()
    {
    }

    /**
     * Extract data from _data_struct and store it in data collection object
     * @param data_struct data struct object
     * @param data data collection object
     * @param gen_id index of generator
     */
    void extract(_data_struct &data_struct,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
      // HAS_EXCITER
      if (!data->getValue(HAS_EXCITER,&bval,g_id)) {
        data->addValue(HAS_EXCITER, true, g_id);
      } else {
        data->setValue(HAS_EXCITER, true, g_id);
      }
	  
      // EXCITER_MODEL              "MODEL"        string
      std::string stmp;
      if (!data->getValue(EXCITER_MODEL,&stmp,g_id)) {
        data->addValue(EXCITER_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(EXCITER_MODEL, data_struct.model, g_id);
      }
	  
	  int ival;
      if (!data->getValue(GENERATOR_REECA_IREG,&ival,g_id)) {
        data->addValue(GENERATOR_REECA_IREG, data_struct.reeca1_ireg, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IREG, data_struct.reeca1_ireg, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_PFFLAG,&ival,g_id)) {
        data->addValue(GENERATOR_REECA_PFFLAG, data_struct.reeca1_pfflag, g_id);
      } else {
        data->setValue(GENERATOR_REECA_PFFLAG, data_struct.reeca1_pfflag, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VFLAG,&ival,g_id)) {
        data->addValue(GENERATOR_REECA_VFLAG, data_struct.reeca1_vflag, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VFLAG, data_struct.reeca1_vflag, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_QFLAG,&ival,g_id)) {
        data->addValue(GENERATOR_REECA_QFLAG, data_struct.reeca1_qflag, g_id);
      } else {
        data->setValue(GENERATOR_REECA_QFLAG, data_struct.reeca1_qflag, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_PFLAG,&ival,g_id)) {
        data->addValue(GENERATOR_REECA_PFLAG, data_struct.reeca1_pflag, g_id);
      } else {
        data->setValue(GENERATOR_REECA_PFLAG, data_struct.reeca1_pflag, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_PQFLAG,&ival,g_id)) {
        data->addValue(GENERATOR_REECA_PQFLAG, data_struct.reeca1_pqflag, g_id);
      } else {
        data->setValue(GENERATOR_REECA_PQFLAG, data_struct.reeca1_pqflag, g_id);
      }
	  
	  //J parameters start here

      if (!data->getValue(GENERATOR_REECA_VDIP,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VDIP,data_struct.reeca1_vdip, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VDIP, data_struct.reeca1_vdip, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VUP,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VUP,data_struct.reeca1_vup, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VUP, data_struct.reeca1_vup, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_TRV,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_TRV,data_struct.reeca1_trv, g_id);
      } else {
        data->setValue(GENERATOR_REECA_TRV, data_struct.reeca1_trv, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_DBD1,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_DBD1,data_struct.reeca1_dbd1, g_id);
      } else {
        data->setValue(GENERATOR_REECA_DBD1, data_struct.reeca1_dbd1, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_DBD2,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_DBD2,data_struct.reeca1_dbd2, g_id);
      } else {
        data->setValue(GENERATOR_REECA_DBD2, data_struct.reeca1_dbd2, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_KQV,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_KQV,data_struct.reeca1_kqv, g_id);
      } else {
        data->setValue(GENERATOR_REECA_KQV, data_struct.reeca1_kqv, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_IQH1,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_IQH1,data_struct.reeca1_lqh1, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IQH1, data_struct.reeca1_lqh1, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_IQL1,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_IQL1,data_struct.reeca1_lql1, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IQL1, data_struct.reeca1_lql1, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VREF0,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VREF0,data_struct.reeca1_vref0, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VREF0, data_struct.reeca1_vref0, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_IQFRZ,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_IQFRZ,data_struct.reeca1_lqfrz, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IQFRZ, data_struct.reeca1_lqfrz, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_THLD,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_THLD,data_struct.reeca1_thld, g_id);
      } else {
        data->setValue(GENERATOR_REECA_THLD, data_struct.reeca1_thld, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_THLD2,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_THLD2,data_struct.reeca1_thld2, g_id);
      } else {
        data->setValue(GENERATOR_REECA_THLD2, data_struct.reeca1_thld2, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_TP,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_TP,data_struct.reeca1_tp, g_id);
      } else {
        data->setValue(GENERATOR_REECA_TP, data_struct.reeca1_tp, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_QMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_QMAX,data_struct.reeca1_qmax, g_id);
      } else {
        data->setValue(GENERATOR_REECA_QMAX, data_struct.reeca1_qmax, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_QMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_QMIN,data_struct.repca1_qmin, g_id);
      } else {
        data->setValue(GENERATOR_REECA_QMIN, data_struct.repca1_qmin, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VMAX,data_struct.reeca1_vmax, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VMAX, data_struct.reeca1_vmax, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VMIN,data_struct.reeca1_vmin, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VMIN, data_struct.reeca1_vmin, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_KQP,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_KQP,data_struct.reeca1_kqp, g_id);
      } else {
        data->setValue(GENERATOR_REECA_KQP, data_struct.reeca1_kqp, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_KQI,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_KQI,data_struct.reeca1_kqi, g_id);
      } else {
        data->setValue(GENERATOR_REECA_KQI, data_struct.reeca1_kqi, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_KVP,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_KVP,data_struct.reeca1_kvp, g_id);
      } else {
        data->setValue(GENERATOR_REECA_KVP, data_struct.reeca1_kvp, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_KVI,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_KVI,data_struct.reeca1_kvi, g_id);
      } else {
        data->setValue(GENERATOR_REECA_KVI, data_struct.reeca1_kvi, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VBIAS,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VBIAS,data_struct.reeca1_vbias, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VBIAS, data_struct.reeca1_vbias, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_TIQ,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_TIQ,data_struct.reeca1_tiq, g_id);
      } else {
        data->setValue(GENERATOR_REECA_TIQ, data_struct.reeca1_tiq, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_DPMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_DPMAX,data_struct.reeca1_dpmax, g_id);
      } else {
        data->setValue(GENERATOR_REECA_DPMAX, data_struct.reeca1_dpmax, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_DPMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_DPMIN,data_struct.reeca1_dpmin, g_id);
      } else {
        data->setValue(GENERATOR_REECA_DPMIN, data_struct.reeca1_dpmin, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_PMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_PMAX,data_struct.reeca1_pmax, g_id);
      } else {
        data->setValue(GENERATOR_REECA_PMAX, data_struct.reeca1_pmax, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_PMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_PMIN,data_struct.reeca1_pmin, g_id);
      } else {
        data->setValue(GENERATOR_REECA_PMIN, data_struct.reeca1_pmin, g_id);
      }
	  
	  
	  if (!data->getValue(GENERATOR_REECA_IMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_IMAX,data_struct.reeca1_imax, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IMAX, data_struct.reeca1_imax, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_TPORD,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_TPORD,data_struct.reeca1_tpord, g_id);
      } else {
        data->setValue(GENERATOR_REECA_TPORD, data_struct.reeca1_tpord, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VQ1,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VQ1,data_struct.reeca1_vq1, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VQ1, data_struct.reeca1_vq1, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_IQ1,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_IQ1,data_struct.reeca1_iq1, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IQ1, data_struct.reeca1_iq1, g_id);
      }
	  
	  
	  if (!data->getValue(GENERATOR_REECA_VQ2,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VQ2,data_struct.reeca1_vq2, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VQ2, data_struct.reeca1_vq2, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_IQ2,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_IQ2,data_struct.reeca1_iq2, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IQ2, data_struct.reeca1_iq2, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VQ3,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VQ3,data_struct.reeca1_vq3, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VQ3, data_struct.reeca1_vq3, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_IQ3,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_IQ3,data_struct.reeca1_iq3, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IQ3, data_struct.reeca1_iq3, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VQ4,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VQ4,data_struct.reeca1_vq4, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VQ4, data_struct.reeca1_vq4, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_IQ4,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_IQ4,data_struct.reeca1_iq4, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IQ4, data_struct.reeca1_iq4, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VP1,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VP1,data_struct.reeca1_vp1, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VP1, data_struct.reeca1_vp1, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_IP1,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_IP1,data_struct.reeca1_ip1, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IP1, data_struct.reeca1_ip1, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VP2,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VP2,data_struct.reeca1_vp2, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VP2, data_struct.reeca1_vp2, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_IP2,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_IP2,data_struct.reeca1_ip2, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IP2, data_struct.reeca1_ip2, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VP3,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VP3,data_struct.reeca1_vp3, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VP3, data_struct.reeca1_vp3, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_IP3,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_IP3,data_struct.reeca1_ip3, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IP3, data_struct.reeca1_ip3, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_VP4,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_VP4,data_struct.reeca1_vp4, g_id);
      } else {
        data->setValue(GENERATOR_REECA_VP4, data_struct.reeca1_vp4, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REECA_IP4,&rval,g_id)) {
        data->addValue(GENERATOR_REECA_IP4,data_struct.reeca1_ip4, g_id);
      } else {
        data->setValue(GENERATOR_REECA_IP4, data_struct.reeca1_ip4, g_id);
      }


    }

    /**
     * Parser list of strings and store results in data collection object
     * @param split_line list of tokens from .dyr file
     * @param data data collection object
     * @param gen_id index of generator
     */
    void parse(std::vector<std::string> &split_line,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
      int nstr = split_line.size();
      // HAS_EXCITER
      if (!data->getValue(HAS_EXCITER,&bval,g_id)) {
        data->addValue(HAS_EXCITER, true, g_id);
      } else {
        data->setValue(HAS_EXCITER, true, g_id);
      }
	  
      // EXCITER_MODEL              "MODEL"                  string
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(EXCITER_MODEL, &stmp, g_id)) {
        data->addValue(EXCITER_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(EXCITER_MODEL, model.c_str(), g_id);
      }

      int ival;
      if (nstr > 3) {
        if (!data->getValue(GENERATOR_REECA_IREG,&ival,g_id)) {
          data->addValue(GENERATOR_REECA_IREG, atoi(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IREG, atoi(split_line[3].c_str()), g_id);
        }
      } 

      if (nstr > 4) {
        if (!data->getValue(GENERATOR_REECA_PFFLAG,&ival,g_id)) {
          data->addValue(GENERATOR_REECA_PFFLAG, atoi(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_PFFLAG, atoi(split_line[4].c_str()), g_id);
        }
      } 


      if (nstr > 5) {
        if (!data->getValue(GENERATOR_REECA_VFLAG,&ival,g_id)) {
          data->addValue(GENERATOR_REECA_VFLAG, atoi(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VFLAG, atoi(split_line[5].c_str()), g_id);
        }
      } 

      if (nstr > 6) {

        if (!data->getValue(GENERATOR_REECA_QFLAG,&ival,g_id)) {
          data->addValue(GENERATOR_REECA_QFLAG, atoi(split_line[6].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_QFLAG, atoi(split_line[6].c_str()), g_id);
        }
		
      } 

      if (nstr > 7) {
        if (!data->getValue(GENERATOR_REECA_PFLAG,&ival,g_id)) {
          data->addValue(GENERATOR_REECA_PFLAG, atoi(split_line[7].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_PFLAG, atoi(split_line[7].c_str()), g_id);
        }
      }


      if (nstr > 8) {
        if (!data->getValue(GENERATOR_REECA_PQFLAG,&ival,g_id)) {
          data->addValue(GENERATOR_REECA_PQFLAG, atoi(split_line[8].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_PQFLAG, atoi(split_line[8].c_str()), g_id);
        }
      } 


// start J parameters here

      if (nstr > 9) {
        if (!data->getValue(GENERATOR_REECA_VDIP,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VDIP, atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VDIP, atof(split_line[9].c_str()), g_id);
        }
      } 

      if (nstr > 10) {
        if (!data->getValue(GENERATOR_REECA_VUP,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VUP, atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VUP, atof(split_line[10].c_str()), g_id);
        }
      } 


      if (nstr > 11) {
        if (!data->getValue(GENERATOR_REECA_TRV,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_TRV, atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_TRV, atof(split_line[11].c_str()), g_id);
        }
      } 


      if (nstr > 12) {
        if (!data->getValue(GENERATOR_REECA_DBD1,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_DBD1, atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_DBD1, atof(split_line[12].c_str()), g_id);
        }
      } 

      if (nstr > 13) {
        if (!data->getValue(GENERATOR_REECA_DBD2,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_DBD2,
              atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_DBD2,
              atof(split_line[13].c_str()), g_id);
        }
      } 


      if (nstr > 14) {
        if (!data->getValue(GENERATOR_REECA_KQV,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_KQV, atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_KQV, atof(split_line[14].c_str()), g_id);
        }
      } 


      if (nstr > 15) {
        if (!data->getValue(GENERATOR_REECA_IQH1,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_IQH1, atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IQH1, atof(split_line[15].c_str()), g_id);
        }
      } 
	  

      if (nstr > 16) {
        if (!data->getValue(GENERATOR_REECA_IQL1,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_IQL1, atof(split_line[16].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IQL1, atof(split_line[16].c_str()), g_id);
        }
      } 
	  
	  if (nstr > 17) {
        if (!data->getValue(GENERATOR_REECA_VREF0,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VREF0, atof(split_line[17].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VREF0, atof(split_line[17].c_str()), g_id);
        }
      } 
	  
	  if (nstr > 18  ) {
        if (!data->getValue(GENERATOR_REECA_IQFRZ,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_IQFRZ, atof(split_line[ 18].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IQFRZ, atof(split_line[ 18].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  19 ) {
        if (!data->getValue(GENERATOR_REECA_THLD,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_THLD, atof(split_line[ 19].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_THLD, atof(split_line[ 19].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  20 ) {
        if (!data->getValue(GENERATOR_REECA_THLD2,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_THLD2, atof(split_line[ 20].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_THLD2, atof(split_line[ 20].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  21 ) {
        if (!data->getValue(GENERATOR_REECA_TP,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_TP, atof(split_line[ 21].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_TP, atof(split_line[ 21].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  22 ) {
        if (!data->getValue(GENERATOR_REECA_QMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_QMAX, atof(split_line[ 22].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_QMAX, atof(split_line[ 22].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  23 ) {
        if (!data->getValue(GENERATOR_REECA_QMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_QMIN, atof(split_line[ 23].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_QMIN, atof(split_line[ 23].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  24 ) {
        if (!data->getValue(GENERATOR_REECA_VMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VMAX, atof(split_line[ 24].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VMAX, atof(split_line[ 24].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  25 ) {
        if (!data->getValue(GENERATOR_REECA_VMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VMIN, atof(split_line[ 25].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VMIN, atof(split_line[ 25].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  26 ) {
        if (!data->getValue(GENERATOR_REECA_KQP,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_KQP, atof(split_line[ 26].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_KQP, atof(split_line[ 26].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  27 ) {
        if (!data->getValue(GENERATOR_REECA_KQI,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_KQI, atof(split_line[ 27].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_KQI, atof(split_line[ 27].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  28 ) {
        if (!data->getValue(GENERATOR_REECA_KVP,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_KVP, atof(split_line[ 28].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_KVP, atof(split_line[ 28].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  29 ) {
        if (!data->getValue(GENERATOR_REECA_KVI,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_KVI, atof(split_line[ 29].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_KVI, atof(split_line[ 29].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  30 ) {
        if (!data->getValue(GENERATOR_REECA_VBIAS,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VBIAS, atof(split_line[ 30].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VBIAS, atof(split_line[ 30].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  31 ) {
        if (!data->getValue(GENERATOR_REECA_TIQ,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_TIQ, atof(split_line[ 31].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_TIQ, atof(split_line[ 32].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  32 ) {
        if (!data->getValue(GENERATOR_REECA_DPMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_DPMAX, atof(split_line[ 32].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_DPMAX, atof(split_line[ 32].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  33 ) {
        if (!data->getValue(GENERATOR_REECA_DPMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_DPMIN, atof(split_line[ 33].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_DPMIN, atof(split_line[ 33].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  34 ) {
        if (!data->getValue(GENERATOR_REECA_PMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_PMAX, atof(split_line[ 34].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_PMAX, atof(split_line[ 34].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   35) {
        if (!data->getValue(GENERATOR_REECA_PMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_PMIN, atof(split_line[ 35].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_PMIN, atof(split_line[ 35].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   36) {
        if (!data->getValue(GENERATOR_REECA_IMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_IMAX, atof(split_line[ 36].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IMAX, atof(split_line[ 36].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   37) {
        if (!data->getValue(GENERATOR_REECA_TPORD,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_TPORD, atof(split_line[ 37].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_TPORD, atof(split_line[ 37].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   38) {
        if (!data->getValue(GENERATOR_REECA_VQ1,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VQ1, atof(split_line[ 38 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VQ1, atof(split_line[ 38 ].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   39) {
        if (!data->getValue(GENERATOR_REECA_IQ1,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_IQ1, atof(split_line[ 39 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IQ1, atof(split_line[ 39 ].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   40) {
        if (!data->getValue(GENERATOR_REECA_VQ2,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VQ2, atof(split_line[ 40 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VQ2, atof(split_line[ 40 ].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   41) {
        if (!data->getValue(GENERATOR_REECA_IQ2,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_IQ2, atof(split_line[ 41 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IQ2, atof(split_line[ 41 ].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   42) {
        if (!data->getValue(GENERATOR_REECA_VQ3,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VQ3, atof(split_line[ 42 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VQ3, atof(split_line[ 42 ].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   43) {
        if (!data->getValue(GENERATOR_REECA_IQ3,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_IQ3, atof(split_line[ 43 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IQ3, atof(split_line[ 43 ].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   44) {
        if (!data->getValue(GENERATOR_REECA_VQ4,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VQ4, atof(split_line[ 44 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VQ4, atof(split_line[ 44 ].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   45) {
        if (!data->getValue(GENERATOR_REECA_IQ4,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_IQ4, atof(split_line[ 45 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IQ4, atof(split_line[  45].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   46) {
        if (!data->getValue(GENERATOR_REECA_VP1,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VP1, atof(split_line[ 46 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VP1, atof(split_line[ 46 ].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   47) {
        if (!data->getValue(GENERATOR_REECA_IP1,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_IP1, atof(split_line[ 47 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IP1, atof(split_line[  47].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   48) {
        if (!data->getValue(GENERATOR_REECA_VP2,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VP2, atof(split_line[ 48 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VP2, atof(split_line[ 48 ].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   49) {
        if (!data->getValue(GENERATOR_REECA_IP2,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_IP2, atof(split_line[ 49 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IP2, atof(split_line[  49].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   50) {
        if (!data->getValue(GENERATOR_REECA_VP3,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VP3, atof(split_line[ 50 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VP3, atof(split_line[  50].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   51) {
        if (!data->getValue(GENERATOR_REECA_IP3,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_IP3, atof(split_line[ 51 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IP3, atof(split_line[  51].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   52) {
        if (!data->getValue(GENERATOR_REECA_VP4,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_VP4, atof(split_line[  52].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_VP4, atof(split_line[ 52 ].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   53) {
        if (!data->getValue(GENERATOR_REECA_IP4,&rval,g_id)) {
          data->addValue(GENERATOR_REECA_IP4, atof(split_line[ 53 ].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REECA_IP4, atof(split_line[ 53 ].c_str()), g_id);
        }
      } 
	  
    }

    /**
     * Parse list of strings store results in data_struct object
     * @param split_line list of tokens from .dyr file
     * @param data data struct that stores information from file
     */
    void store(std::vector<std::string> &split_line,_data_struct &data)
    {
      // GENERATOR_BUSNUMBER               "I"                   integer
      int o_idx;
      o_idx = atoi(split_line[0].c_str());
      data.bus_id = o_idx;

      // Clean up 2 character tag for generator ID
      gridpack::utility::StringUtils util;
      std::string tag = util.clean2Char(split_line[2]);
      strcpy(data.gen_id, tag.c_str());

      std::string sval;

      sval = util.trimQuotes(split_line[1]);
      util.toUpper(sval);

      // GENERATOR_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size();
     
      if (nstr > 3) {
        data.reeca1_ireg = atoi(split_line[3].c_str());
      } 

     
      if (nstr > 4) {
        data.reeca1_pfflag = atoi(split_line[4].c_str());
      } 

    
      if (nstr > 5) {
        data.reeca1_vflag = atoi(split_line[5].c_str());
      } 

                             
      if (nstr > 6) {
        data.reeca1_qflag = atoi(split_line[6].c_str());
      } 


      if (nstr > 7) {
        data.reeca1_pflag = atoi(split_line[7].c_str());
      }


      if (nstr > 8) {
        data.reeca1_pqflag = atoi(split_line[8].c_str());
      } 

// J parameters start here
      if (nstr > 9) {
        data.reeca1_vdip = atof(split_line[9].c_str());
      } 


      if (nstr > 10) {
        data.reeca1_vup = atof(split_line[10].c_str());
      } 

      if (nstr > 11) {
        data.reeca1_trv = atof(split_line[11].c_str());
      } 


      if (nstr > 12) {
        data.reeca1_dbd1 = atof(split_line[12].c_str());
      } 


      if (nstr > 13) {
        data.reeca1_dbd2 = atof(split_line[13].c_str());
      } 


      if (nstr > 14) {
        data.reeca1_kqv = atof(split_line[14].c_str());
      } 


      if (nstr > 15) {
        data.reeca1_lqh1 = atof(split_line[15].c_str());
      } 
	  
	  if (nstr > 16) {
        data.reeca1_lql1 = atof(split_line[16].c_str());
      } 
	  
	  if (nstr > 17) {
        data.reeca1_vref0 = atof(split_line[17].c_str());
      } 
	  
	  	  if (nstr >  18) {
        data.reeca1_lqfrz = atof(split_line[ 18 ].c_str());
      } 
	  
	  if (nstr >  19) {
        data.reeca1_thld = atof(split_line[ 19 ].c_str());
      } 
	  
	  if (nstr >  20) {
        data.reeca1_thld2 = atof(split_line[ 20 ].c_str());
      } 
	  
	  if (nstr >  21) {
        data.reeca1_tp = atof(split_line[ 21 ].c_str());
      } 
	  
	  if (nstr >  22) {
        data.reeca1_qmax = atof(split_line[ 22 ].c_str());
      } 
	  
	  if (nstr >  23) {
        data.reeca1_qmin = atof(split_line[ 23 ].c_str());
      } 
	  
	  if (nstr >  24) {
        data.reeca1_vmax = atof(split_line[ 24 ].c_str());
      } 
	  
	  if (nstr >  25) {
        data.reeca1_vmin = atof(split_line[  25].c_str());
      } 
	  
	  if (nstr >  26) {
        data.reeca1_kqp = atof(split_line[ 26 ].c_str());
      } 
	  
	  if (nstr >  27) {
        data.reeca1_kqi = atof(split_line[  27].c_str());
      } 
	  
	  if (nstr >  28) {
        data.reeca1_kvp = atof(split_line[ 28 ].c_str());
      } 
	  
	  if (nstr >  29) {
        data.reeca1_kvi = atof(split_line[ 29 ].c_str());
      } 
	  
	  if (nstr >  30) {
        data.reeca1_vbias = atof(split_line[  30].c_str());
      } 
	  
	  if (nstr >  31) {
        data.reeca1_tiq = atof(split_line[  31].c_str());
      } 
	  
	  if (nstr >  32) {
        data.reeca1_dpmax = atof(split_line[  32].c_str());
      } 
	  
	  if (nstr >  33) {
        data.reeca1_dpmin = atof(split_line[  33].c_str());
      } 
	  
	  if (nstr >  34) {
        data.reeca1_pmax = atof(split_line[  34].c_str());
      } 
	  
	  if (nstr >  35) {
        data.reeca1_pmin = atof(split_line[  35].c_str());
      } 
	  
	  if (nstr >  36) {
        data.reeca1_imax = atof(split_line[ 36 ].c_str());
      } 
	  
	  if (nstr >  37 ) {
        data.reeca1_tpord = atof(split_line[ 37].c_str());
      } 
	  
	  if (nstr >   38) {
        data.reeca1_vq1 = atof(split_line[ 38 ].c_str());
      } 
	  
	  if (nstr >   39) {
        data.reeca1_iq1 = atof(split_line[ 39 ].c_str());
      } 
	  
	  if (nstr >   40) {
        data.reeca1_vq2 = atof(split_line[ 40].c_str());
      } 
	  
	  if (nstr >   41) {
        data.reeca1_iq2 = atof(split_line[ 41].c_str());
      } 
	  
	  if (nstr >   42) {
        data.reeca1_vq3 = atof(split_line[ 42].c_str());
      } 
	  
	  if (nstr >   43) {
        data.reeca1_iq3 = atof(split_line[ 43].c_str());
      } 
	  
	  if (nstr >   44) {
        data.reeca1_vq4 = atof(split_line[ 44].c_str());
      } 
	  
	  if (nstr >   45) {
        data.reeca1_iq4 = atof(split_line[ 45].c_str());
      } 
	  
	  if (nstr >   46) {
        data.reeca1_vp1 = atof(split_line[ 46].c_str());
      } 
	  
	  if (nstr >   47) {
        data.reeca1_ip1 = atof(split_line[ 47].c_str());
      } 
	  
	  if (nstr >   48) {
        data.reeca1_vp2 = atof(split_line[ 48].c_str());
      } 
	  
	  if (nstr >   49) {
        data.reeca1_ip2 = atof(split_line[ 49].c_str());
      } 
	  
	  if (nstr >  50 ) {
        data.reeca1_vp3 = atof(split_line[ 50].c_str());
      } 
	  
	  if (nstr >  51 ) {
        data.reeca1_ip3 = atof(split_line[ 51].c_str());
      } 
	  
	  if (nstr >  52 ) {
        data.reeca1_vp4 = atof(split_line[ 52].c_str());
      } 
	  
	  if (nstr >  53 ) {
        data.reeca1_ip4 = atof(split_line[ 53].c_str());
      } 
	  
	  
    }
	
};
}  // parser
}  // gridpack
#endif
