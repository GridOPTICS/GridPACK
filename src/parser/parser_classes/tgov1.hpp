/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: June 18, 2020
 *      Author: Bruce Palmer
 */
#ifndef TGOV1_HPP
#define TGOV1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Tgov1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Tgov1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Tgov1Parser()
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
      // HAS_GOVERNOR
      if (!data->getValue(HAS_GOVERNOR,&bval,g_id)) {
        data->addValue(HAS_GOVERNOR, true, g_id);
      } else {
        data->setValue(HAS_GOVERNOR, true, g_id);
      }

      // GOVERNOR_NAME
      std::string stmp;
      if (!data->getValue(GOVERNOR_MODEL, &stmp, g_id)) {
        data->addValue(GOVERNOR_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(GOVERNOR_MODEL, data_struct.model, g_id);
      }

      // GOVERNOR_R
      if (!data->getValue(GOVERNOR_R,&rval,g_id)) {
        data->addValue(GOVERNOR_R, data_struct.gv_r, g_id);
      } else {
        data->setValue(GOVERNOR_R, data_struct.gv_r, g_id);
      }

      // GOVERNOR_T1
      if (!data->getValue(GOVERNOR_T1,&rval,g_id)) {
        data->addValue(GOVERNOR_T1, data_struct.gv_t1, g_id);
      } else {
        data->setValue(GOVERNOR_T1, data_struct.gv_t1, g_id);
      }

      // GOVERNOR_VMAX
      if (!data->getValue(GOVERNOR_VMAX,&rval,g_id)) {
        data->addValue(GOVERNOR_VMAX, data_struct.vmax, g_id);
      } else {
        data->setValue(GOVERNOR_VMAX, data_struct.vmax, g_id);
      }

      // GOVERNOR_VMIN
      if (!data->getValue(GOVERNOR_VMIN,&rval,g_id)) {
        data->addValue(GOVERNOR_VMIN, data_struct.vmin, g_id);
      } else {
        data->setValue(GOVERNOR_VMIN, data_struct.vmin, g_id);
      }

      // GOVERNOR_T2
      if (!data->getValue(GOVERNOR_T2,&rval,g_id)) {
        data->addValue(GOVERNOR_T2, data_struct.gv_t2, g_id);
      } else {
        data->setValue(GOVERNOR_T2, data_struct.gv_t2, g_id);
      }

      // GOVERNOR_T3
      if (!data->getValue(GOVERNOR_T3,&rval,g_id)) {
        data->addValue(GOVERNOR_T3, data_struct.gv_t3, g_id);
      } else {
        data->setValue(GOVERNOR_T3, data_struct.gv_t3, g_id);
      }

      // GOVERNOR_DT
      if (!data->getValue(GOVERNOR_DT,&rval,g_id)) {
        data->addValue(GOVERNOR_DT, data_struct.gv_dt, g_id);
      } else {
        data->setValue(GOVERNOR_DT, data_struct.gv_dt, g_id);
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
      bool bval;
      // HAS_GOVERNOR
      if (!data->getValue(HAS_GOVERNOR,&bval,g_id)) {
        data->addValue(HAS_GOVERNOR, true, g_id);
      } else {
        data->setValue(HAS_GOVERNOR, true, g_id);
      }

      // GOVERNOR_MODEL
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(GOVERNOR_MODEL,&stmp,g_id)) {
        data->addValue(GOVERNOR_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(GOVERNOR_MODEL, model.c_str(), g_id);
      }

      // GOVERNOR_R
      if (nstr > 3) {
        if (!data->getValue(GOVERNOR_R,&rval,g_id)) {
          data->addValue(GOVERNOR_R, atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_R, atof(split_line[3].c_str()), g_id);
        }
      } 

      // GOVERNOR_T1
      if (nstr > 4) {
        if (!data->getValue(GOVERNOR_T1,&rval,g_id)) {
          data->addValue(GOVERNOR_T1, atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_T1, atof(split_line[4].c_str()), g_id);
        }
      } 

      // GOVERNOR_VMAX
      if (nstr > 5) {
        if (!data->getValue(GOVERNOR_VMAX,&rval,g_id)) {
          data->addValue(GOVERNOR_VMAX, atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_VMAX, atof(split_line[5].c_str()), g_id);
        }
      } 

      // GOVERNOR_VMIN
      if (nstr > 6) {
        if (!data->getValue(GOVERNOR_VMIN,&rval,g_id)) {
          data->addValue(GOVERNOR_VMIN, atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_VMIN, atof(split_line[6].c_str()), g_id);
        }
      } 

      // GOVERNOR_T2
      if (nstr > 7) {
        if (!data->getValue(GOVERNOR_T2,&rval,g_id)) {
          data->addValue(GOVERNOR_T2, atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_T2, atof(split_line[7].c_str()), g_id);
        }
      } 

      // GOVERNOR_T3
      if (nstr > 8) {
        if (!data->getValue(GOVERNOR_T3,&rval,g_id)) {
          data->addValue(GOVERNOR_T3, atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_T3, atof(split_line[8].c_str()), g_id);
        }
      } 

      // GOVERNOR_DT
      if (nstr > 9) {
        if (!data->getValue(GOVERNOR_DT,&rval,g_id)) {
          data->addValue(GOVERNOR_DT, atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_DT, atof(split_line[9].c_str()), g_id);
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
      // GOVERNOR_BUSNUMBER               "I"                   integer
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

      // GOVERNOR_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size();
      // GOVERNOR_R
      if (nstr > 3) {
        data.gv_r = atof(split_line[3].c_str());
      }

      // GOVERNOR_T1
      if (nstr > 4) {
        data.gv_t1 = atof(split_line[4].c_str());
      }

      // GOVERNOR_VMAX
      if (nstr > 5) {
        data.vmax = atof(split_line[5].c_str());
      }

      // GOVERNOR_VMIN
      if (nstr > 6) {
        data.vmin = atof(split_line[6].c_str());
      }

      // GOVERNOR_T2
      if (nstr > 7) {
        data.gv_t2 = atof(split_line[7].c_str());
      }

      // GOVERNOR_T3
      if (nstr > 8) {
        data.gv_t3 = atof(split_line[8].c_str());
      }

      // GOVERNOR_DT
      if (nstr > 9) {
        data.gv_dt = atof(split_line[9].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
