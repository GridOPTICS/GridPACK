/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: June 17, 2016
 *      Author: Bruce Palmer
 */
#ifndef SEXS_HPP
#define SEXS_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class SexsParser
{
  public:
    /**
     * Constructor
     */
    explicit SexsParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~SexsParser()
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

      // EXCITER_MODEL
      std::string stmp;
      if (!data->getValue(EXCITER_MODEL, &stmp, g_id)) {
        data->addValue(EXCITER_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(EXCITER_MODEL, data_struct.model, g_id);
      }

      // EXCITER_TA_OVER_TB
      if (!data->getValue(EXCITER_TA_OVER_TB,&rval,g_id)) {
        data->addValue(EXCITER_TA_OVER_TB, data_struct.ex_ta_over_tb, g_id);
      } else {
        data->setValue(EXCITER_TA_OVER_TB, data_struct.ex_ta_over_tb, g_id);
      }

      // EXCITER_TB
      if (!data->getValue(EXCITER_TB,&rval,g_id)) {
        data->addValue(EXCITER_TB, data_struct.ex_tb, g_id);
      } else {
        data->setValue(EXCITER_TB, data_struct.ex_tb, g_id);
      }

      // EXCITER_K
      if (!data->getValue(EXCITER_K,&rval,g_id)) {
        data->addValue(EXCITER_K, data_struct.ex_k, g_id);
      } else {
        data->setValue(EXCITER_K, data_struct.ex_k, g_id);
      }

      // EXCITER_TE
      if (!data->getValue(EXCITER_TE,&rval,g_id)) {
        data->addValue(EXCITER_TE, data_struct.ex_te, g_id);
      } else {
        data->setValue(EXCITER_TE, data_struct.ex_te, g_id);
      }


      // EXCITER_EMIN
      if (!data->getValue(EXCITER_EMIN,&rval,g_id)) {
        data->addValue(EXCITER_EMIN, data_struct.ex_emin, g_id);
      } else {
        data->setValue(EXCITER_EMIN, data_struct.ex_emin, g_id);
      }

      // EXCITER_EMAX
      if (!data->getValue(EXCITER_EMAX,&rval,g_id)) {
        data->addValue(EXCITER_EMAX, data_struct.ex_emax, g_id);
      } else {
        data->setValue(EXCITER_EMAX, data_struct.ex_emax, g_id);
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

      // EXCITER_MODEL
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(EXCITER_MODEL,&stmp,g_id)) {
        data->addValue(EXCITER_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(EXCITER_MODEL, model.c_str(), g_id);
      }

      // EXCITER_TA_OVER_TB
      if (nstr > 3) {
        if (!data->getValue(EXCITER_TA_OVER_TB,&rval,g_id)) {
          data->addValue(EXCITER_TA_OVER_TB,
              atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TA_OVER_TB,
              atof(split_line[3].c_str()), g_id);
        }
      } 

      // EXCITER_TB
      if (nstr > 4) {
        if (!data->getValue(EXCITER_TB,&rval,g_id)) {
          data->addValue(EXCITER_TB,
              atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TB,
              atof(split_line[4].c_str()), g_id);
        }
      } 

      // EXCITER_K
      if (nstr > 5) {
        if (!data->getValue(EXCITER_K,&rval,g_id)) {
          data->addValue(EXCITER_K,
              atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(EXCITER_K,
              atof(split_line[5].c_str()), g_id);
        }
      }

      // EXCITER_TE
      if (nstr > 6) {
        if (!data->getValue(EXCITER_TE,&rval,g_id)) {
          data->addValue(EXCITER_TE,
              atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TE,
              atof(split_line[6].c_str()), g_id);
        }
      }

      // EXCITER_EMIN
      if (nstr > 7) {
        if (!data->getValue(EXCITER_EMIN,&rval,g_id)) {
          data->addValue(EXCITER_EMIN,
              atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(EXCITER_EMIN,
              atof(split_line[7].c_str()), g_id);
        }
      }

      // EXCITER_EMAX
      if (nstr > 8) {
        if (!data->getValue(EXCITER_EMAX,&rval,g_id)) {
          data->addValue(EXCITER_EMAX,
              atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(EXCITER_EMAX,
              atof(split_line[8].c_str()), g_id);
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
      // EXCITER_BUSNUMBER               "I"                   integer
      int o_idx;
      o_idx = atoi(split_line[0].c_str());
      data.bus_id = o_idx;

      // Clean up 2 character tag for generator ID
      gridpack::utility::StringUtils util;
      std::string tag = util.clean2Char(split_line[2]);
      strcpy(data.gen_id, tag.c_str());

      std::string sval;
      double rval;
      int ival;

      sval = util.trimQuotes(split_line[1]);
      util.toUpper(sval);

      // EXCITER_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size();
      // EXCITER_TA_OVER_TB
      if (nstr > 3) {
        data.ex_ta_over_tb = atof(split_line[3].c_str());
      }

      // EXCITER_TB
      if (nstr > 4) {
        data.ex_tb = atof(split_line[4].c_str());
      }

      // EXCITER_K
      if (nstr > 5) {
        data.ex_k = atof(split_line[5].c_str());
      }

      // EXCITER_TE
      if (nstr > 6) {
        data.ex_te = atof(split_line[6].c_str());
      }

      // EXCITER_EMIN
      if (nstr > 7) {
        data.ex_emin = atof(split_line[7].c_str());
      }

      // EXCITER_EMAX
      if (nstr > 8) {
        data.ex_emax = atof(split_line[8].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
