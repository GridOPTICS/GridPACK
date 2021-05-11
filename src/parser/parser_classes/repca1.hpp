/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: May 10, 2021
 *      Author: Renke Huang
 */
#ifndef REPCA1_HPP	
#define REPCA1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Repca1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Repca1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Repca1Parser()
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
      // GENERATOR_MODEL              "MODEL"        string
      std::string stmp;
      if (!data->getValue(GENERATOR_MODEL,&stmp,g_id)) {
        data->addValue(GENERATOR_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(GENERATOR_MODEL, data_struct.model, g_id);
      }
	  
	  int ival;
      if (!data->getValue(GENERATOR_REPCA_IREG,&ival,g_id)) {
        data->addValue(GENERATOR_REPCA_IREG, data_struct.repca1_ireg, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_IREG, data_struct.repca1_ireg, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_BRCH_BUS_FROM,&ival,g_id)) {
        data->addValue(GENERATOR_REPCA_BRCH_BUS_FROM, data_struct.repca1_brh_bus_from, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_BRCH_BUS_FROM, data_struct.repca1_brh_bus_from, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_BRCH_BUS_TO,&ival,g_id)) {
        data->addValue(GENERATOR_REPCA_BRCH_BUS_TO, data_struct.repca1_brh_bus_to, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_BRCH_BUS_TO, data_struct.repca1_brh_bus_to, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_BRCH_CKT,&stmp,g_id)) {
        data->addValue(GENERATOR_REPCA_BRCH_CKT, data_struct.repca1_brh_ckt, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_BRCH_CKT, data_struct.repca1_brh_ckt, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_VC_FLAG,&ival,g_id)) {
        data->addValue(GENERATOR_REPCA_VC_FLAG, data_struct.repca1_vcflag, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_VC_FLAG, data_struct.repca1_vcflag, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_REF_FLAG,&ival,g_id)) {
        data->addValue(GENERATOR_REPCA_REF_FLAG, data_struct.repca1_refflag, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_REF_FLAG, data_struct.repca1_refflag, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_F_FLAG,&ival,g_id)) {
        data->addValue(GENERATOR_REPCA_F_FLAG, data_struct.repca1_fflag, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_F_FLAG, data_struct.repca1_fflag, g_id);
      }
	  
	  //J parameters start here

      if (!data->getValue(GENERATOR_REPCA_TFLTR,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_TFLTR,data_struct.repca1_tfltr, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_TFLTR, data_struct.repca1_tfltr, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_KP,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_KP,data_struct.repca1_kp, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_KP, data_struct.repca1_kp, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_KI,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_KI,data_struct.repca1_ki, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_KI, data_struct.repca1_ki, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_TFT,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_TFT,data_struct.repca1_tft, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_TFT, data_struct.repca1_tft, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_TFV,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_TFV,data_struct.repca1_tfv, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_TFV, data_struct.repca1_tfv, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_VFRZ,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_VFRZ,data_struct.repca1_vfrz, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_VFRZ, data_struct.repca1_vfrz, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_RC,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_RC,data_struct.repca1_rc, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_RC, data_struct.repca1_rc, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_XC,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_XC,data_struct.repca1_xc, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_XC, data_struct.repca1_xc, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_KC,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_KC,data_struct.repca1_kc, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_KC, data_struct.repca1_kc, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_EMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_EMAX,data_struct.repca1_emax, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_EMAX, data_struct.repca1_emax, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_EMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_EMIN,data_struct.repca1_emin, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_EMIN, data_struct.repca1_emin, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_DBD1,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_DBD1,data_struct.repca1_dbd1, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_DBD1, data_struct.repca1_dbd1, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_DBD2,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_DBD2,data_struct.repca1_dbd2, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_DBD2, data_struct.repca1_dbd2, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_QMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_QMAX,data_struct.repca1_qmax, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_QMAX, data_struct.repca1_qmax, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_QMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_QMIN,data_struct.repca1_qmin, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_QMIN, data_struct.repca1_qmin, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_KPG,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_KPG,data_struct.repca1_kpg, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_KPG, data_struct.repca1_kpg, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_KIG,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_KIG,data_struct.repca1_kig, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_KIG, data_struct.repca1_kig, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_TP,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_TP,data_struct.repca1_tp, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_TP, data_struct.repca1_tp, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_FDBD1,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_FDBD1,data_struct.repca1_fdbd1, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_FDBD1, data_struct.repca1_fdbd1, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_FDBD2,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_FDBD2,data_struct.repca1_fdbd2, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_FDBD2, data_struct.repca1_fdbd2, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_FEMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_FEMAX,data_struct.repca1_femax, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_FEMAX, data_struct.repca1_femax, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_FEMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_FEMIN,data_struct.repca1_femin, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_FEMIN, data_struct.repca1_femin, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_PMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_PMAX,data_struct.repca1_pmax, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_PMAX, data_struct.repca1_pmax, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_PMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_PMIN,data_struct.repca1_pmin, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_PMIN, data_struct.repca1_pmin, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_TG,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_TG,data_struct.repca1_tg, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_TG, data_struct.repca1_tg, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_DDN,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_DDN,data_struct.repca1_ddn, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_DDN, data_struct.repca1_ddn, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REPCA_DUP,&rval,g_id)) {
        data->addValue(GENERATOR_REPCA_DUP,data_struct.repca1_dup, g_id);
      } else {
        data->setValue(GENERATOR_REPCA_DUP, data_struct.repca1_dup, g_id);
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
      int nstr = split_line.size();
      // GENERATOR_MODEL              "MODEL"                  string
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(GENERATOR_MODEL, &stmp, g_id)) {
        data->addValue(GENERATOR_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(GENERATOR_MODEL, model.c_str(), g_id);
      }

      int ival;
      if (nstr > 3) {
        if (!data->getValue(GENERATOR_REPCA_IREG,&ival,g_id)) {
          data->addValue(GENERATOR_REPCA_IREG, atoi(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_IREG, atoi(split_line[3].c_str()), g_id);
        }
      } 

      if (nstr > 4) {
        if (!data->getValue(GENERATOR_REPCA_BRCH_BUS_FROM,&ival,g_id)) {
          data->addValue(GENERATOR_REPCA_BRCH_BUS_FROM, atoi(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_BRCH_BUS_FROM, atoi(split_line[4].c_str()), g_id);
        }
      } 


      if (nstr > 5) {
        if (!data->getValue(GENERATOR_REPCA_BRCH_BUS_TO,&ival,g_id)) {
          data->addValue(GENERATOR_REPCA_BRCH_BUS_TO, atoi(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_BRCH_BUS_TO, atoi(split_line[5].c_str()), g_id);
        }
      } 


      if (nstr > 6) {
		/*
        if (!data->getValue(GENERATOR_REPCA_BRCH_CKT,&ival,g_id)) {
          data->addValue(GENERATOR_REPCA_BRCH_CKT, atoi(split_line[6].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_BRCH_CKT, atoi(split_line[6].c_str()), g_id);
        }
		*/

	    std::string idtmp = util.trimQuotes(split_line[6]);
        util.toUpper(idtmp);
        if (!data->getValue(GENERATOR_REPCA_BRCH_CKT, &stmp, g_id)) {
          data->addValue(GENERATOR_REPCA_BRCH_CKT, idtmp.c_str(), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_BRCH_CKT, idtmp.c_str(), g_id);
        }
		
      } 

      if (nstr > 7) {
        if (!data->getValue(GENERATOR_REPCA_VC_FLAG,&ival,g_id)) {
          data->addValue(GENERATOR_REPCA_VC_FLAG, atoi(split_line[7].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_VC_FLAG, atoi(split_line[7].c_str()), g_id);
        }
      }


      if (nstr > 8) {
        if (!data->getValue(GENERATOR_REPCA_REF_FLAG,&ival,g_id)) {
          data->addValue(GENERATOR_REPCA_REF_FLAG, atoi(split_line[8].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_REF_FLAG, atoi(split_line[8].c_str()), g_id);
        }
      } 


      if (nstr > 9) {
        if (!data->getValue(GENERATOR_REPCA_F_FLAG,&ival,g_id)) {
          data->addValue(GENERATOR_REPCA_F_FLAG, atoi(split_line[9].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_F_FLAG, atoi(split_line[9].c_str()), g_id);
        }
      } 

// start J parameters here

      if (nstr > 10) {
        if (!data->getValue(GENERATOR_REPCA_TFLTR,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_TFLTR, atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_TFLTR, atof(split_line[10].c_str()), g_id);
        }
      } 


      if (nstr > 11) {
        if (!data->getValue(GENERATOR_REPCA_KP,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_KP, atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_KP, atof(split_line[11].c_str()), g_id);
        }
      } 


      if (nstr > 12) {
        if (!data->getValue(GENERATOR_REPCA_KI,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_KI, atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_KI, atof(split_line[12].c_str()), g_id);
        }
      } 

      if (nstr > 13) {
        if (!data->getValue(GENERATOR_REPCA_TFT,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_TFT,
              atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_TFT,
              atof(split_line[13].c_str()), g_id);
        }
      } 


      if (nstr > 14) {
        if (!data->getValue(GENERATOR_REPCA_TFV,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_TFV, atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_TFV, atof(split_line[14].c_str()), g_id);
        }
      } 


      if (nstr > 15) {
        if (!data->getValue(GENERATOR_REPCA_VFRZ,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_VFRZ, atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_VFRZ, atof(split_line[15].c_str()), g_id);
        }
      } 
	  

      if (nstr > 16) {
        if (!data->getValue(GENERATOR_REPCA_RC,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_RC, atof(split_line[16].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_RC, atof(split_line[16].c_str()), g_id);
        }
      } 
	  
	  if (nstr > 17) {
        if (!data->getValue(GENERATOR_REPCA_XC,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_XC, atof(split_line[17].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_XC, atof(split_line[17].c_str()), g_id);
        }
      } 
	  
	  if (nstr > 18  ) {
        if (!data->getValue(GENERATOR_REPCA_KC,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_KC, atof(split_line[ 18].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_KC, atof(split_line[ 18].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  19 ) {
        if (!data->getValue(GENERATOR_REPCA_EMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_EMAX, atof(split_line[ 19].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_EMAX, atof(split_line[ 19].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  20 ) {
        if (!data->getValue(GENERATOR_REPCA_EMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_EMIN, atof(split_line[ 20].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_EMIN, atof(split_line[ 20].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  21 ) {
        if (!data->getValue(GENERATOR_REPCA_DBD1,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_DBD1, atof(split_line[ 21].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_DBD1, atof(split_line[ 21].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  22 ) {
        if (!data->getValue(GENERATOR_REPCA_DBD2,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_DBD2, atof(split_line[ 22].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_DBD2, atof(split_line[ 22].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  23 ) {
        if (!data->getValue(GENERATOR_REPCA_QMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_QMAX, atof(split_line[ 23].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_QMAX, atof(split_line[ 23].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  24 ) {
        if (!data->getValue(GENERATOR_REPCA_QMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_QMIN, atof(split_line[ 24].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_QMIN, atof(split_line[ 24].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  25 ) {
        if (!data->getValue(GENERATOR_REPCA_KPG,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_KPG, atof(split_line[ 25].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_KPG, atof(split_line[ 25].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  26 ) {
        if (!data->getValue(GENERATOR_REPCA_KIG,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_KIG, atof(split_line[ 26].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_KIG, atof(split_line[ 26].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  27 ) {
        if (!data->getValue(GENERATOR_REPCA_TP,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_TP, atof(split_line[ 27].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_TP, atof(split_line[ 27].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  28 ) {
        if (!data->getValue(GENERATOR_REPCA_FDBD1,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_FDBD1, atof(split_line[ 28].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_FDBD1, atof(split_line[ 28].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  29 ) {
        if (!data->getValue(GENERATOR_REPCA_FDBD2,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_FDBD2, atof(split_line[ 29].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_FDBD2, atof(split_line[ 29].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  30 ) {
        if (!data->getValue(GENERATOR_REPCA_FEMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_FEMAX, atof(split_line[ 30].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_FEMAX, atof(split_line[ 30].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  31 ) {
        if (!data->getValue(GENERATOR_REPCA_FEMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_FEMIN, atof(split_line[ 31].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_FEMIN, atof(split_line[ 32].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  32 ) {
        if (!data->getValue(GENERATOR_REPCA_PMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_PMAX, atof(split_line[ 32].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_PMAX, atof(split_line[ 32].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  33 ) {
        if (!data->getValue(GENERATOR_REPCA_PMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_PMIN, atof(split_line[ 33].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_PMIN, atof(split_line[ 33].c_str()), g_id);
        }
      } 
	  
	  if (nstr >  34 ) {
        if (!data->getValue(GENERATOR_REPCA_TG,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_TG, atof(split_line[ 34].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_TG, atof(split_line[ 34].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   35) {
        if (!data->getValue(GENERATOR_REPCA_DDN,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_DDN, atof(split_line[ 35].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_DDN, atof(split_line[ 35].c_str()), g_id);
        }
      } 
	  
	  if (nstr >   36) {
        if (!data->getValue(GENERATOR_REPCA_DUP,&rval,g_id)) {
          data->addValue(GENERATOR_REPCA_DUP, atof(split_line[ 36].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REPCA_DUP, atof(split_line[ 36].c_str()), g_id);
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
        data.repca1_ireg = atoi(split_line[3].c_str());
      } 

     
      if (nstr > 4) {
        data.repca1_brh_bus_from = atoi(split_line[4].c_str());
      } 

    
      if (nstr > 5) {
        data.repca1_brh_bus_to = atoi(split_line[5].c_str());
      } 

                             
      if (nstr > 6) {
		tag = util.clean2Char(split_line[6]);
        strcpy(data.repca1_brh_ckt, tag.c_str());
      } 


      if (nstr > 7) {
        data.repca1_vcflag = atoi(split_line[7].c_str());
      }


      if (nstr > 8) {
        data.repca1_refflag = atoi(split_line[8].c_str());
      } 


      if (nstr > 9) {
        data.repca1_fflag = atoi(split_line[9].c_str());
      } 

/// start J parameters here

      if (nstr > 10) {
        data.repca1_tfltr = atof(split_line[10].c_str());
      } 

      if (nstr > 11) {
        data.repca1_kp = atof(split_line[11].c_str());
      } 


      if (nstr > 12) {
        data.repca1_ki = atof(split_line[12].c_str());
      } 


      if (nstr > 13) {
        data.repca1_tft = atof(split_line[13].c_str());
      } 


      if (nstr > 14) {
        data.repca1_tfv = atof(split_line[14].c_str());
      } 


      if (nstr > 15) {
        data.repca1_vfrz = atof(split_line[15].c_str());
      } 
	  
	  if (nstr > 16) {
        data.repca1_rc = atof(split_line[16].c_str());
      } 
	  
	  if (nstr > 17) {
        data.repca1_xc = atof(split_line[17].c_str());
      } 
	  
	  if (nstr >  18) {
        data.repca1_kc = atof(split_line[ 18 ].c_str());
      } 
	  
	  if (nstr >  19) {
        data.repca1_emax = atof(split_line[ 19 ].c_str());
      } 
	  
	  if (nstr >  20) {
        data.repca1_emin = atof(split_line[ 20 ].c_str());
      } 
	  
	  if (nstr >  21) {
        data.repca1_dbd1 = atof(split_line[ 21 ].c_str());
      } 
	  
	  if (nstr >  22) {
        data.repca1_dbd2 = atof(split_line[ 22 ].c_str());
      } 
	  
	  if (nstr >  23) {
        data.repca1_qmax = atof(split_line[ 23 ].c_str());
      } 
	  
	  if (nstr >  24) {
        data.repca1_qmin = atof(split_line[ 24 ].c_str());
      } 
	  
	  if (nstr >  25) {
        data.repca1_kpg = atof(split_line[  25].c_str());
      } 
	  
	  if (nstr >  26) {
        data.repca1_kig = atof(split_line[ 26 ].c_str());
      } 
	  
	  if (nstr >  27) {
        data.repca1_tp = atof(split_line[  27].c_str());
      } 
	  
	  if (nstr >  28) {
        data.repca1_fdbd1 = atof(split_line[ 28 ].c_str());
      } 
	  
	  if (nstr >  29) {
        data.repca1_fdbd2 = atof(split_line[ 29 ].c_str());
      } 
	  
	  if (nstr >  30) {
        data.repca1_femax = atof(split_line[  30].c_str());
      } 
	  
	  if (nstr >  31) {
        data.repca1_femin = atof(split_line[  31].c_str());
      } 
	  
	  if (nstr >  32) {
        data.repca1_pmax = atof(split_line[  32].c_str());
      } 
	  
	  if (nstr >  33) {
        data.repca1_pmin = atof(split_line[  33].c_str());
      } 
	  
	  if (nstr >  34) {
        data.repca1_tg = atof(split_line[  34].c_str());
      } 
	  
	  if (nstr >  35) {
        data.repca1_ddn = atof(split_line[  35].c_str());
      } 
	  
	  if (nstr >  36) {
        data.repca1_dup = atof(split_line[ 36 ].c_str());
      } 

 
    }
	
};
}  // parser
}  // gridpack
#endif
