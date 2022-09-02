/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: May 10, 2021
 *      Author: Renke Huang
 */
#ifndef REGCA1_HPP
#define REGCA1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Regca1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Regca1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Regca1Parser()
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
      if (!data->getValue(GENERATOR_REGCA_LVPLSW,&ival,g_id)) {
        data->addValue(GENERATOR_REGCA_LVPLSW, data_struct.regca1_lvplsw, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_LVPLSW, data_struct.regca1_lvplsw, g_id);
      }

      if (!data->getValue(GENERATOR_REGCA_TG,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_TG,data_struct.regca1_tg, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_TG, data_struct.regca1_tg, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REGCA_RRPWR,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_RRPWR,data_struct.regca1_rrpwr, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_RRPWR, data_struct.regca1_rrpwr, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REGCA_BRKPT,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_BRKPT,data_struct.regca1_brkpt, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_BRKPT, data_struct.regca1_brkpt, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REGCA_ZEROX,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_ZEROX,data_struct.regca1_zerox, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_ZEROX, data_struct.regca1_zerox, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REGCA_LVPL1,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_LVPL1,data_struct.regca1_lvpl1, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_LVPL1, data_struct.regca1_lvpl1, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REGCA_VOLIM,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_VOLIM,data_struct.regca1_volim, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_VOLIM, data_struct.regca1_volim, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REGCA_LVPNT1,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_LVPNT1,data_struct.regca1_lvpnt1, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_LVPNT1, data_struct.regca1_lvpnt1, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REGCA_LVPNT0,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_LVPNT0,data_struct.regca1_lvpnt0, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_LVPNT0, data_struct.regca1_lvpnt0, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REGCA_LOLIM,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_LOLIM,data_struct.regca1_lolim, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_LOLIM, data_struct.regca1_lolim, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REGCA_TFLTR,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_TFLTR,data_struct.regca1_tfltr, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_TFLTR, data_struct.regca1_tfltr, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REGCA_KHV,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_KHV,data_struct.regca1_khv, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_KHV, data_struct.regca1_khv, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REGCA_LQRMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_LQRMAX,data_struct.regca1_lqrmax, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_LQRMAX, data_struct.regca1_lqrmax, g_id);
      }
	  
	  if (!data->getValue(GENERATOR_REGCA_LQRMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_LQRMIN,data_struct.regca1_lqrmin, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_LQRMIN, data_struct.regca1_lqrmin, g_id);
      }  
	  
	  if (!data->getValue(GENERATOR_REGCA_ACCEL,&rval,g_id)) {
        data->addValue(GENERATOR_REGCA_ACCEL,data_struct.regca1_accel, g_id);
      } else {
        data->setValue(GENERATOR_REGCA_ACCEL, data_struct.regca1_accel, g_id);
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
        if (!data->getValue(GENERATOR_REGCA_LVPLSW,&ival,g_id)) {
          data->addValue(GENERATOR_REGCA_LVPLSW, atoi(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_LVPLSW, atoi(split_line[3].c_str()), g_id);
        }
      } 

      if (nstr > 4) {
        if (!data->getValue(GENERATOR_REGCA_TG,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_TG, atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_TG, atof(split_line[4].c_str()), g_id);
        }
      } 


      if (nstr > 5) {
        if (!data->getValue(GENERATOR_REGCA_RRPWR,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_RRPWR, atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_RRPWR, atof(split_line[5].c_str()), g_id);
        }
      } 


      if (nstr > 6) {
        if (!data->getValue(GENERATOR_REGCA_BRKPT,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_BRKPT, atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_BRKPT, atof(split_line[6].c_str()), g_id);
        }
      } 


      if (nstr > 7) {
        if (!data->getValue(GENERATOR_REGCA_ZEROX,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_ZEROX, atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_ZEROX, atof(split_line[7].c_str()), g_id);
        }
      }


      if (nstr > 8) {
        if (!data->getValue(GENERATOR_REGCA_LVPL1,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_LVPL1, atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_LVPL1, atof(split_line[8].c_str()), g_id);
        }
      } 


      if (nstr > 9) {
        if (!data->getValue(GENERATOR_REGCA_VOLIM,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_VOLIM, atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_VOLIM, atof(split_line[9].c_str()), g_id);
        }
      } 


      if (nstr > 10) {
        if (!data->getValue(GENERATOR_REGCA_LVPNT1,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_LVPNT1, atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_LVPNT1, atof(split_line[10].c_str()), g_id);
        }
      } 


      if (nstr > 11) {
        if (!data->getValue(GENERATOR_REGCA_LVPNT0,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_LVPNT0, atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_LVPNT0, atof(split_line[11].c_str()), g_id);
        }
      } 


      if (nstr > 12) {
        if (!data->getValue(GENERATOR_REGCA_LOLIM,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_LOLIM, atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_LOLIM, atof(split_line[12].c_str()), g_id);
        }
      } 

      if (nstr > 13) {
        if (!data->getValue(GENERATOR_REGCA_TFLTR,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_TFLTR,
              atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_TFLTR,
              atof(split_line[13].c_str()), g_id);
        }
      } 


      if (nstr > 14) {
        if (!data->getValue(GENERATOR_REGCA_KHV,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_KHV, atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_KHV, atof(split_line[14].c_str()), g_id);
        }
      } 


      if (nstr > 15) {
        if (!data->getValue(GENERATOR_REGCA_LQRMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_LQRMAX, atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_LQRMAX, atof(split_line[15].c_str()), g_id);
        }
      } 
	  

      if (nstr > 16) {
        if (!data->getValue(GENERATOR_REGCA_LQRMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_LQRMIN, atof(split_line[16].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_LQRMIN, atof(split_line[16].c_str()), g_id);
        }
      } 
	  
	  if (nstr > 17) {
        if (!data->getValue(GENERATOR_REGCA_ACCEL,&rval,g_id)) {
          data->addValue(GENERATOR_REGCA_ACCEL, atof(split_line[17].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGCA_ACCEL, atof(split_line[17].c_str()), g_id);
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
        data.regca1_lvplsw = atoi(split_line[3].c_str());
      } 

     
      if (nstr > 4) {
        data.regca1_tg = atof(split_line[4].c_str());
      } 

    
      if (nstr > 5) {
        data.regca1_rrpwr = atof(split_line[5].c_str());
      } 

                             
      if (nstr > 6) {
        data.regca1_brkpt = atof(split_line[6].c_str());
      } 


      if (nstr > 7) {
        data.regca1_zerox = atof(split_line[7].c_str());
      }


      if (nstr > 8) {
        data.regca1_lvpl1 = atof(split_line[8].c_str());
      } 


      if (nstr > 9) {
        data.regca1_volim = atof(split_line[9].c_str());
      } 


      if (nstr > 10) {
        data.regca1_lvpnt1 = atof(split_line[10].c_str());
      } 

      if (nstr > 11) {
        data.regca1_lvpnt0 = atof(split_line[11].c_str());
      } 


      if (nstr > 12) {
        data.regca1_lolim = atof(split_line[12].c_str());
      } 


      if (nstr > 13) {
        data.regca1_tfltr = atof(split_line[13].c_str());
      } 


      if (nstr > 14) {
        data.regca1_khv = atof(split_line[14].c_str());
      } 


      if (nstr > 15) {
        data.regca1_lqrmax = atof(split_line[15].c_str());
      } 
	  
	  if (nstr > 16) {
        data.regca1_lqrmin = atof(split_line[16].c_str());
      } 
	  
	  if (nstr > 17) {
        data.regca1_accel = atof(split_line[17].c_str());
      } 
    }
	
};
}  // parser
}  // gridpack
#endif
