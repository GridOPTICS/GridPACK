/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: November 2, 2022
 *      Author: Bruce Palmer
 */
#ifndef WTARA1_HPP
#define WTARA1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Wtara1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Wtara1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Wtara1Parser()
    {
    }

    /**
     * Extract data from _data_struct and store it in data collection object
     * @param data_struct data struct object
     * @param data data collection object
     * @param gen_id index of user model
     */
    void extract(_data_struct &data_struct,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
      int ival;
      // HAS_WIND
      if (!data->getValue(HAS_WIND,&bval,g_id)) {
        data->addValue(HAS_WIND, true, g_id);
      } else {
        data->setValue(HAS_WIND, true, g_id);
      }

      // WIND_MODEL
      std::string stmp;
      if (!data->getValue(WIND_MODEL, &stmp, g_id)) {
        data->addValue(WIND_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(WIND_MODEL, data_struct.model, g_id);
      }

      // WIND_ID
      if (!data->getValue(WIND_ID, &stmp, g_id)) {
        data->addValue(WIND_ID, data_struct.gen_id, g_id);
      } else {
        data->setValue(WIND_ID, data_struct.gen_id, g_id);
      }

      // WIND_KA
      if (!data->getValue(WIND_KA,&rval,g_id)) {
        data->addValue(WIND_KA, data_struct.wind_ka, g_id);
      } else {
        data->setValue(WIND_KA, data_struct.wind_ka, g_id);
      }

      // WIND_THETA
      if (!data->getValue(WIND_THETA,&rval,g_id)) {
        data->addValue(WIND_THETA, data_struct.wind_theta, g_id);
      } else {
        data->setValue(WIND_THETA, data_struct.wind_theta, g_id);
      }
    }

    /**
     * Parser list of strings and store results in data collection object
     * @param split_line list of tokens from .dyr file
     * @param data data collection object
     * @param gen_id index of user model
     */
    void parse(std::vector<std::string> &split_line,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      int ival;
      int nstr = split_line.size();
      bool bval;
      // HAS_WIND
      if (!data->getValue(HAS_WIND,&bval,g_id)) {
        data->addValue(HAS_WIND, true, g_id);
      } else {
        data->setValue(HAS_WIND, true, g_id);
      }

      // WIND_NAME
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(WIND_MODEL,&stmp,g_id)) {
        data->addValue(WIND_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(WIND_MODEL, model.c_str(), g_id);
      }

      // WIND_ID
      std::string tag;
      tag = util.clean2Char(split_line[2]);
      if (!data->getValue(WIND_ID,&stmp,g_id)) {
        data->addValue(WIND_ID, tag.c_str(), g_id);
      } else {
        data->setValue(WIND_ID, tag.c_str(), g_id);
      }

      // Use counter to keep track of additional parameters
      if (nstr > 3) {
        if (!data->getValue(WIND_KA,&rval,g_id)) {
          data->addValue(WIND_KA, atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(WIND_KA, atof(split_line[3].c_str()), g_id);
        }
      }

      if (nstr > 4) {
        if (!data->getValue(WIND_THETA,&rval,g_id)) {
          data->addValue(WIND_THETA, atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(WIND_THETA, atof(split_line[4].c_str()), g_id);
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
      // WIND_BUSNUMBER               "I"                   integer
      int o_idx;
      int ncnt, ni, nc, ns, nv;
      o_idx = atoi(split_line[0].c_str());
      data.bus_id = o_idx;

      // Clean up 2 character tag for user model ID
      gridpack::utility::StringUtils util;
      std::string tag = util.clean2Char(split_line[2]);
      strcpy(data.gen_id, tag.c_str());

      std::string sval;
      sval = util.trimQuotes(split_line[1]);
      util.toUpper(sval);

      // WIND_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size();

      // Wind model parameters

      if (nstr > 3) {
        data.wind_ka = atof(split_line[3].c_str());
      }

      if (nstr > 4) {
        data.wind_theta = atof(split_line[4].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
