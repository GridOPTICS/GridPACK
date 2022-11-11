/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: November 2, 2022
 *      Author: Bruce Palmer
 */
#ifndef WTDTA1_HPP
#define WTDTA1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Wtdta1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Wtdta1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Wtdta1Parser()
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
      // HAS_WIND_DRIVETRAIN
      if (!data->getValue(HAS_WIND_DRIVETRAIN,&bval,g_id)) {
        data->addValue(HAS_WIND_DRIVETRAIN, true, g_id);
      } else {
        data->setValue(HAS_WIND_DRIVETRAIN, true, g_id);
      }

      // WIND_DRIVETRAIN
      std::string stmp;
      if (!data->getValue(WIND_DRIVETRAIN, &stmp, g_id)) {
        data->addValue(WIND_DRIVETRAIN, data_struct.model, g_id);
      } else {
        data->setValue(WIND_DRIVETRAIN, data_struct.model, g_id);
      }

      // WIND_ID
      if (!data->getValue(WIND_ID, &stmp, g_id)) {
        data->addValue(WIND_ID, data_struct.gen_id, g_id);
      } else {
        data->setValue(WIND_ID, data_struct.gen_id, g_id);
      }

      // WIND_H
      if (!data->getValue(WIND_DT_H,&rval,g_id)) {
        data->addValue(WIND_DT_H, data_struct.wind_h, g_id);
      } else {
        data->setValue(WIND_DT_H, data_struct.wind_h, g_id);
      }

      // WIND_DAMP
      if (!data->getValue(WIND_DT_DAMP,&rval,g_id)) {
        data->addValue(WIND_DT_DAMP, data_struct.wind_damp, g_id);
      } else {
        data->setValue(WIND_DT_DAMP, data_struct.wind_damp, g_id);
      }

      // WIND_HFRAC
      if (!data->getValue(WIND_DT_HFRAC,&rval,g_id)) {
        data->addValue(WIND_DT_HFRAC, data_struct.wind_hfrac, g_id);
      } else {
        data->setValue(WIND_DT_HFRAC, data_struct.wind_hfrac, g_id);
      }

      // WIND_FREQ1
      if (!data->getValue(WIND_DT_FREQ1,&rval,g_id)) {
        data->addValue(WIND_DT_FREQ1, data_struct.wind_freq1, g_id);
      } else {
        data->setValue(WIND_DT_FREQ1, data_struct.wind_freq1, g_id);
      }

      // WIND_DSHAFT
      if (!data->getValue(WIND_DT_DSHAFT,&rval,g_id)) {
        data->addValue(WIND_DT_DSHAFT, data_struct.wind_dshaft, g_id);
      } else {
        data->setValue(WIND_DT_DSHAFT, data_struct.wind_dshaft, g_id);
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
      // HAS_WIND_DRIVETRAIN
      if (!data->getValue(HAS_WIND_DRIVETRAIN,&bval,g_id)) {
        data->addValue(HAS_WIND_DRIVETRAIN, true, g_id);
      } else {
        data->setValue(HAS_WIND_DRIVETRAIN, true, g_id);
      }

      // WIND_NAME
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(WIND_DRIVETRAIN,&stmp,g_id)) {
        data->addValue(WIND_DRIVETRAIN, model.c_str(), g_id);
      } else {
        data->setValue(WIND_DRIVETRAIN, model.c_str(), g_id);
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
        if (!data->getValue(WIND_DT_H,&rval,g_id)) {
          data->addValue(WIND_DT_H, atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(WIND_DT_H, atof(split_line[3].c_str()), g_id);
        }
      }

      if (nstr > 4) {
        if (!data->getValue(WIND_DT_DAMP,&rval,g_id)) {
          data->addValue(WIND_DT_DAMP, atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(WIND_DT_DAMP, atof(split_line[4].c_str()), g_id);
        }
      }

      if (nstr > 5) {
        if (!data->getValue(WIND_DT_HFRAC,&rval,g_id)) {
          data->addValue(WIND_DT_HFRAC, atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(WIND_DT_HFRAC, atof(split_line[5].c_str()), g_id);
        }
      }

      if (nstr > 6) {
        if (!data->getValue(WIND_DT_FREQ1,&rval,g_id)) {
          data->addValue(WIND_DT_FREQ1, atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(WIND_DT_FREQ1, atof(split_line[6].c_str()), g_id);
        }
      }

      if (nstr > 7) {
        if (!data->getValue(WIND_DT_DSHAFT,&rval,g_id)) {
          data->addValue(WIND_DT_DSHAFT, atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(WIND_DT_DSHAFT, atof(split_line[7].c_str()), g_id);
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

      // WIND_DRIVETRAIN              "MODEL"                  integer
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size();

      // Wind model parameters

      if (nstr > 3) {
        data.wind_h = atof(split_line[3].c_str());
      }

      if (nstr > 4) {
        data.wind_damp = atof(split_line[4].c_str());
      }

      if (nstr > 5) {
        data.wind_hfrac = atof(split_line[5].c_str());
      }

      if (nstr > 6) {
        data.wind_freq1 = atof(split_line[6].c_str());
      }

      if (nstr > 7) {
        data.wind_dshaft = atof(split_line[7].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
