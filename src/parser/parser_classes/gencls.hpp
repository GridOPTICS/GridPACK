/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: June 16, 2016
 *      Author: Bruce Palmer
 */
#ifndef GENCLS_HPP
#define GENCLS_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class GenclsParser
{
  public:
    /**
     * Constructor
     */
    explicit GenclsParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~GenclsParser()
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
      std::string stmp;
      // GENERATOR_MODEL              "MODEL"        string
      if (!data->getValue(GENERATOR_MODEL,&stmp,g_id)) {
        data->addValue(GENERATOR_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(GENERATOR_MODEL, data_struct.model, g_id);
      }

      // GENERATOR_INERTIA_CONSTANT_H                float
      if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H,&rval,g_id)) {
        data->addValue(GENERATOR_INERTIA_CONSTANT_H,
            data_struct.inertia, g_id);
      } else {
        data->setValue(GENERATOR_INERTIA_CONSTANT_H,
            data_struct.inertia, g_id);
      }

      // GENERATOR_DAMPING_COEFFICIENT_0             float
      if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&rval,g_id)) {
        data->addValue(GENERATOR_DAMPING_COEFFICIENT_0,
            data_struct.damping, g_id);
      } else {
        data->setValue(GENERATOR_DAMPING_COEFFICIENT_0,
            data_struct.damping, g_id);
      }

      // GENERATOR_TRANSIENT_REACTANCE               float
      if (!data->getValue(GENERATOR_TRANSIENT_REACTANCE,&rval,g_id)) {
        data->addValue(GENERATOR_TRANSIENT_REACTANCE,
            data_struct.reactance, g_id);
      } else {
        data->setValue(GENERATOR_TRANSIENT_REACTANCE,
            data_struct.reactance, g_id);
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
      if (!data->getValue(GENERATOR_MODEL,&stmp,g_id)) {
        data->addValue(GENERATOR_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(GENERATOR_MODEL, model.c_str(), g_id);
      }

      // GENERATOR_INERTIA_CONSTANT_H                           float
      if (nstr > 3) {
        if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H,&rval,g_id)) {
          data->addValue(GENERATOR_INERTIA_CONSTANT_H,
              atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_INERTIA_CONSTANT_H,
              atof(split_line[3].c_str()), g_id);
        }
      } 

      // GENERATOR_DAMPING_COEFFICIENT_0                           float
      if (nstr > 4) {
        if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&rval,g_id)) {
          data->addValue(GENERATOR_DAMPING_COEFFICIENT_0,
              atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_DAMPING_COEFFICIENT_0,
              atof(split_line[4].c_str()), g_id);
        }
      }
      // GENERATOR_TRANSIENT_REACTANCE               float
      if (nstr > 5) {
        if (!data->getValue(GENERATOR_TRANSIENT_REACTANCE,&rval,g_id)) {
          data->addValue(GENERATOR_TRANSIENT_REACTANCE,
              atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_TRANSIENT_REACTANCE,
              atof(split_line[5].c_str()), g_id);
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
      // GENERATOR_INERTIA_CONSTANT_H                           float
      if (nstr > 3) {
        data.inertia = atof(split_line[3].c_str());
      } 

      // GENERATOR_DAMPING_COEFFICIENT_0                           float
      if (nstr > 4) {
        data.damping = atof(split_line[4].c_str());
      }

      // GENERATOR_TRANSIENT_REACTANCE                             float
      if (nstr > 5) {
        data.reactance = atof(split_line[5].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
