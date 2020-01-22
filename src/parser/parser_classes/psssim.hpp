/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: July 31, 2019 
 *      Author: Shuangshuang Jin
 */
#ifndef PSSSIM_HPP
#define PSSSIM_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class PsssimParser
{
  public:
    /**
     * Constructor
     */
    explicit PsssimParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~PsssimParser()
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

      // HAS_PSS
      if (!data->getValue(HAS_PSS,&bval,g_id)) {
        data->addValue(HAS_PSS, true, g_id);
      } else {
        data->setValue(HAS_PSS, true, g_id);
      }

      // PSSSIM_MODEL
      std::string stmp;
      if (!data->getValue(PSS_MODEL, &stmp, g_id)) {
        data->addValue(PSS_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(PSS_MODEL, data_struct.model, g_id);
      }

      // PSSSIM_INPUTTYPE 
      int ival;
      if (!data->getValue(PSSSIM_INPUTTYPE,&ival,g_id)) {
        data->addValue(PSSSIM_INPUTTYPE, data_struct.psssim_inputtype, g_id);
      } else {
        data->setValue(PSSSIM_INPUTTYPE, data_struct.psssim_inputtype, g_id);
      }

      // PSSSIM_BUS1 
      if (!data->getValue(PSSSIM_BUS1,&ival,g_id)) {
        data->addValue(PSSSIM_BUS1, data_struct.psssim_bus1, g_id);
      } else {
        data->setValue(PSSSIM_BUS1, data_struct.psssim_bus1, g_id);
      }

      // PSSSIM_BUS2 
      if (!data->getValue(PSSSIM_BUS2,&ival,g_id)) {
        data->addValue(PSSSIM_BUS2, data_struct.psssim_bus2, g_id);
      } else {
        data->setValue(PSSSIM_BUS2, data_struct.psssim_bus2, g_id);
      }

      // PSSSIM_BUS3 
      if (!data->getValue(PSSSIM_BUS3,&ival,g_id)) {
        data->addValue(PSSSIM_BUS3, data_struct.psssim_bus3, g_id);
      } else {
        data->setValue(PSSSIM_BUS3, data_struct.psssim_bus3, g_id);
      }

      // PSSSIM_BUS4 
      if (!data->getValue(PSSSIM_BUS4,&ival,g_id)) {
        data->addValue(PSSSIM_BUS4, data_struct.psssim_bus4, g_id);
      } else {
        data->setValue(PSSSIM_BUS4, data_struct.psssim_bus4, g_id);
      }

      // PSSSIM_BUS5 
      if (!data->getValue(PSSSIM_BUS5,&ival,g_id)) {
        data->addValue(PSSSIM_BUS5, data_struct.psssim_bus5, g_id);
      } else {
        data->setValue(PSSSIM_BUS5, data_struct.psssim_bus5, g_id);
      }

      // PSSSIM_BUS6 
      if (!data->getValue(PSSSIM_BUS6,&ival,g_id)) {
        data->addValue(PSSSIM_BUS6, data_struct.psssim_bus6, g_id);
      } else {
        data->setValue(PSSSIM_BUS6, data_struct.psssim_bus6, g_id);
      }

      // PSSSIM_GAINK
      if (!data->getValue(PSSSIM_GAINK,&rval,g_id)) {
        data->addValue(PSSSIM_GAINK, data_struct.psssim_gaink, g_id);
      } else {
        data->setValue(PSSSIM_GAINK, data_struct.psssim_gaink, g_id);
      }

      // PSSSIM_TW
      if (!data->getValue(PSSSIM_TW,&rval,g_id)) {
        data->addValue(PSSSIM_TW, data_struct.psssim_tw, g_id);
      } else {
        data->setValue(PSSSIM_TW, data_struct.psssim_tw, g_id);
      }

      // PSSSIM_T1
      if (!data->getValue(PSSSIM_T1,&rval,g_id)) {
        data->addValue(PSSSIM_T1, data_struct.psssim_t1, g_id);
      } else {
        data->setValue(PSSSIM_T1, data_struct.psssim_t1, g_id);
      }

      // PSSSIM_T2
      if (!data->getValue(PSSSIM_T2,&rval,g_id)) {
        data->addValue(PSSSIM_T2, data_struct.psssim_t2, g_id);
      } else {
        data->setValue(PSSSIM_T2, data_struct.psssim_t2, g_id);
      }

      // PSSSIM_T3
      if (!data->getValue(PSSSIM_T3,&rval,g_id)) {
        data->addValue(PSSSIM_T3, data_struct.psssim_t3, g_id);
      } else {
        data->setValue(PSSSIM_T3, data_struct.psssim_t3, g_id);
      }

      // PSSSIM_T4
      if (!data->getValue(PSSSIM_T4,&rval,g_id)) {
        data->addValue(PSSSIM_T4, data_struct.psssim_t4, g_id);
      } else {
        data->setValue(PSSSIM_T4, data_struct.psssim_t4, g_id);
      }

      // PSSSIM_MAXOUT
      if (!data->getValue(PSSSIM_MAXOUT,&rval,g_id)) {
        data->addValue(PSSSIM_MAXOUT, data_struct.psssim_maxout, g_id);
      } else {
        data->setValue(PSSSIM_MAXOUT, data_struct.psssim_maxout, g_id);
      }

      // PSSSIM_MINOUT
      if (!data->getValue(PSSSIM_MINOUT,&rval,g_id)) {
        data->addValue(PSSSIM_MINOUT, data_struct.psssim_minout, g_id);
      } else {
        data->setValue(PSSSIM_MINOUT, data_struct.psssim_minout, g_id);
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
      // HAS_PSS
      if (!data->getValue(HAS_PSS,&bval,g_id)) {
        data->addValue(HAS_PSS, true, g_id);
      } else {
        data->setValue(HAS_PSS, true, g_id);
      }

      // PSSSIM_MODEL
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(PSS_MODEL,&stmp,g_id)) {
        data->addValue(PSS_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(PSS_MODEL, model.c_str(), g_id);
      }

      // PSSSIM_INPUTTYPE 
      int ival;
      if (nstr > 3) {
        if (!data->getValue(PSSSIM_INPUTTYPE,&ival,g_id)) {
          data->addValue(PSSSIM_INPUTTYPE,
              atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_INPUTTYPE,
              atof(split_line[3].c_str()), g_id);
        }
      } 

      // PSSSIM_BUS1
      if (nstr > 4) {
        if (!data->getValue(PSSSIM_BUS1,&ival,g_id)) {
          data->addValue(PSSSIM_BUS1,
              atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_BUS1,
              atof(split_line[4].c_str()), g_id);
        }
      } 

      // PSSSIM_BUS2 
      if (nstr > 5) {
        if (!data->getValue(PSSSIM_BUS2,&ival,g_id)) {
          data->addValue(PSSSIM_BUS2,
              atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_BUS2,
              atof(split_line[5].c_str()), g_id);
        }
      } 

      // PSSSIM_BUS3
      if (nstr > 6) {
        if (!data->getValue(PSSSIM_BUS3,&ival,g_id)) {
          data->addValue(PSSSIM_BUS3,
              atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_BUS3,
              atof(split_line[6].c_str()), g_id);
        }
      } 

      // PSSSIM_BUS4
      if (nstr > 7) {
        if (!data->getValue(PSSSIM_BUS4,&ival,g_id)) {
          data->addValue(PSSSIM_BUS4,
              atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_BUS4,
              atof(split_line[7].c_str()), g_id);
        }
      } 

      // PSSSIM_BUS5
      if (nstr > 8) {
        if (!data->getValue(PSSSIM_BUS5,&ival,g_id)) {
          data->addValue(PSSSIM_BUS5,
              atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_BUS5,
              atof(split_line[8].c_str()), g_id);
        }
      } 

      // PSSSIM_BUS6
      if (nstr > 9) {
        if (!data->getValue(PSSSIM_BUS6,&ival,g_id)) {
          data->addValue(PSSSIM_BUS6,
              atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_BUS6,
              atof(split_line[9].c_str()), g_id);
        }
      }

      // PSSSIM_GAINK
      if (nstr > 10) {
        if (!data->getValue(PSSSIM_GAINK,&rval,g_id)) {
          data->addValue(PSSSIM_GAINK,
              atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_GAINK,
              atof(split_line[10].c_str()), g_id);
        }
      } 

      // PSSSIM_TW
      if (nstr > 11) {
        if (!data->getValue(PSSSIM_TW,&rval,g_id)) {
          data->addValue(PSSSIM_TW,
              atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_TW,
              atof(split_line[11].c_str()), g_id);
        }
      } 

      // PSSSIM_T1
      if (nstr > 12) {
        if (!data->getValue(PSSSIM_T1,&rval,g_id)) {
          data->addValue(PSSSIM_T1,
              atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_T1,
              atof(split_line[12].c_str()), g_id);
        }
      } 

      // PSSSIM_T2
      if (nstr > 13) {
        if (!data->getValue(PSSSIM_T2,&rval,g_id)) {
          data->addValue(PSSSIM_T2,
              atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_T2,
              atof(split_line[13].c_str()), g_id);
        }
      } 

      // PSSSIM_T3
      if (nstr > 14) {
        if (!data->getValue(PSSSIM_T3,&rval,g_id)) {
          data->addValue(PSSSIM_T3,
              atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_T3,
              atof(split_line[14].c_str()), g_id);
        }
      } 

      // PSSSIM_T4
      if (nstr > 15) {
        if (!data->getValue(PSSSIM_T4,&rval,g_id)) {
          data->addValue(PSSSIM_T4,
              atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_T4,
              atof(split_line[15].c_str()), g_id);
        }
      } 

      // PSSSIM_MAXOUT
      if (nstr > 16) {
        if (!data->getValue(PSSSIM_MAXOUT,&rval,g_id)) {
          data->addValue(PSSSIM_MAXOUT,
              atof(split_line[16].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_MAXOUT,
              atof(split_line[16].c_str()), g_id);
        }
      } 

      // PSSSIM_MINOUT
      if (nstr > 17) {
        if (!data->getValue(PSSSIM_MINOUT,&rval,g_id)) {
          data->addValue(PSSSIM_MINOUT,
              atof(split_line[17].c_str()), g_id);
        } else {
          data->setValue(PSSSIM_MINOUT,
              atof(split_line[17].c_str()), g_id);
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
      // PSSSIM_BUSNUMBER               "I"                   integer
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

      // PSSSIM_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size();
      // PSSSIM_INPUTTYPE
      if (nstr > 3) {
        data.psssim_inputtype= atof(split_line[3].c_str());
      }

      // PSSSIM_BUS1
      if (nstr > 4) {
        data.psssim_bus1= atof(split_line[4].c_str());
      }

      // PSSSIM_BUS2
      if (nstr > 5) {
        data.psssim_bus2 = atof(split_line[5].c_str());
      }

      // PSSSIM_BUS3
      if (nstr > 6) {
        data.psssim_bus3 = atof(split_line[6].c_str());
      }

      // PSSSIM_BUS4
      if (nstr > 7) {
        data.psssim_bus4 = atof(split_line[7].c_str());
      }

      // PSSSIM_BUS5
      if (nstr > 8) {
        data.psssim_bus5 = atof(split_line[8].c_str());
      }

      // PSSSIM_BUS6
      if (nstr > 9) {
        data.psssim_bus6 = atof(split_line[9].c_str());
      }

      // PSSSIM_GAINK
      if (nstr > 10) {
        data.psssim_gaink = atof(split_line[10].c_str());
      }

      // PSSSIM_TW
      if (nstr > 11) {
        data.psssim_tw = atof(split_line[11].c_str());
      }

      // PSSSIM_T1
      if (nstr > 12) {
        data.psssim_t1 = atof(split_line[12].c_str());
      }

      // PSSSIM_T2
      if (nstr > 13) {
        data.psssim_t2 = atof(split_line[13].c_str());
      }

      // PSSSIM_T3
      if (nstr > 14) {
        data.psssim_t3 = atof(split_line[14].c_str());
      }

      // PSSSIM_T4
      if (nstr > 15) {
        data.psssim_t4 = atof(split_line[15].c_str());
      }

      // PSSSIM_MAXOUT
      if (nstr > 16) {
        data.psssim_maxout = atof(split_line[16].c_str());
      }

      // PSSSIM_MINOUT
      if (nstr > 17) {
        data.psssim_minout = atof(split_line[17].c_str());
      }

    }
};
}  // parser
}  // gridpack
#endif
