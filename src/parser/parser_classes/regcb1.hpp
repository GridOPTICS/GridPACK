/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: Dec 31, 2022
 *      Author: Shrirang Abhyankar
 */
#ifndef REGCB1_HPP
#define REGCB1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Regcb1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Regcb1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Regcb1Parser()
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
      if (!data->getValue(GENERATOR_REGC_RATEFLAG,&ival,g_id)) {
        data->addValue(GENERATOR_REGC_RATEFLAG, data_struct.regc_rateflag, g_id);
      } else {
        data->setValue(GENERATOR_REGC_RATEFLAG, data_struct.regc_rateflag, g_id);
      }

      if (!data->getValue(GENERATOR_REGC_PQFLAG,&ival,g_id)) {
        data->addValue(GENERATOR_REGC_PQFLAG, data_struct.regc_pqflag, g_id);
      } else {
        data->setValue(GENERATOR_REGC_PQFLAG, data_struct.regc_pqflag, g_id);
      }

      if (!data->getValue(GENERATOR_REGC_TG,&rval,g_id)) {
        data->addValue(GENERATOR_REGC_TG,data_struct.regc_tg, g_id);
      } else {
        data->setValue(GENERATOR_REGC_TG, data_struct.regc_tg, g_id);
      }

      if (!data->getValue(GENERATOR_REGC_TFLTR,&rval,g_id)) {
        data->addValue(GENERATOR_REGC_TFLTR,data_struct.regc_tfltr, g_id);
      } else {
        data->setValue(GENERATOR_REGC_TFLTR, data_struct.regc_tfltr, g_id);
      }

      if (!data->getValue(GENERATOR_REGC_IQRMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REGC_IQRMAX,data_struct.regc_iqrmax, g_id);
      } else {
        data->setValue(GENERATOR_REGC_IQRMAX, data_struct.regc_iqrmax, g_id);
      }
	  
      if (!data->getValue(GENERATOR_REGC_IQRMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REGC_IQRMIN,data_struct.regc_iqrmin, g_id);
      } else {
        data->setValue(GENERATOR_REGC_IQRMIN, data_struct.regc_iqrmin, g_id);
      }  

      if (!data->getValue(GENERATOR_REGC_RRPWR,&rval,g_id)) {
        data->addValue(GENERATOR_REGC_RRPWR,data_struct.regc_rrpwr, g_id);
      } else {
        data->setValue(GENERATOR_REGC_RRPWR, data_struct.regc_rrpwr, g_id);
      }

      if (!data->getValue(GENERATOR_REGC_TE,&rval,g_id)) {
        data->addValue(GENERATOR_REGC_TE,data_struct.regc_te, g_id);
      } else {
        data->setValue(GENERATOR_REGC_TE, data_struct.regc_te, g_id);
      }

      if (!data->getValue(GENERATOR_REGC_IMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REGC_IMAX,data_struct.regc_imax, g_id);
      } else {
        data->setValue(GENERATOR_REGC_IMAX, data_struct.regc_imax, g_id);
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
        if (!data->getValue(GENERATOR_REGC_RATEFLAG,&ival,g_id)) {
          data->addValue(GENERATOR_REGC_RATEFLAG, atoi(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGC_RATEFLAG, atoi(split_line[3].c_str()), g_id);
        }
      } 

      if (nstr > 4) {
        if (!data->getValue(GENERATOR_REGC_PQFLAG,&ival,g_id)) {
          data->addValue(GENERATOR_REGC_PQFLAG, atoi(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGC_PQFLAG, atoi(split_line[4].c_str()), g_id);
        }
      } 

      if (nstr > 5) {
        if (!data->getValue(GENERATOR_REGC_TG,&rval,g_id)) {
          data->addValue(GENERATOR_REGC_TG, atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGC_TG, atof(split_line[5].c_str()), g_id);
        }
      } 


      if (nstr > 6) {
        if (!data->getValue(GENERATOR_REGC_TFLTR,&rval,g_id)) {
          data->addValue(GENERATOR_REGC_TFLTR, atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGC_TFLTR, atof(split_line[6].c_str()), g_id);
        }
      } 

      if (nstr > 7) {
        if (!data->getValue(GENERATOR_REGC_IQRMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REGC_IQRMAX, atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGC_IQRMAX, atof(split_line[7].c_str()), g_id);
        }
      } 
	  

      if (nstr > 8) {
        if (!data->getValue(GENERATOR_REGC_IQRMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REGC_IQRMIN, atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGC_IQRMIN, atof(split_line[8].c_str()), g_id);
        }
      } 

      if (nstr > 9) {
        if (!data->getValue(GENERATOR_REGC_RRPWR,&rval,g_id)) {
          data->addValue(GENERATOR_REGC_RRPWR, atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGC_RRPWR, atof(split_line[9].c_str()), g_id);
        }
      }

      if (nstr > 10) {
        if (!data->getValue(GENERATOR_REGC_TE,&rval,g_id)) {
          data->addValue(GENERATOR_REGC_TE, atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGC_TE, atof(split_line[10].c_str()), g_id);
        }
      } 

      if (nstr > 11) {
        if (!data->getValue(GENERATOR_REGC_IMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REGC_IMAX, atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGC_IMAX, atof(split_line[11].c_str()), g_id);
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
        data.regc_rateflag = atoi(split_line[3].c_str());
      } 

      if (nstr > 4) {
        data.regc_pqflag = atoi(split_line[4].c_str());
      } 

      if (nstr > 5) {
        data.regc_tg = atof(split_line[5].c_str());
      } 

      if (nstr > 6) {
        data.regc_tfltr = atof(split_line[6].c_str());
      } 

      if (nstr > 7) {
        data.regc_iqrmax = atof(split_line[7].c_str());
      } 
	  
      if (nstr > 8) {
        data.regc_iqrmin = atof(split_line[8].c_str());
      } 
    
      if (nstr > 9) {
        data.regc_rrpwr = atof(split_line[9].c_str());
      } 

      if (nstr > 10) {
        data.regc_te = atof(split_line[10].c_str());
      } 
	  
      if (nstr > 11) {
        data.regc_imax = atof(split_line[11].c_str());
      } 
    }
	
};
}  // parser
}  // gridpack
#endif
