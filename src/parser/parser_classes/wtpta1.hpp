/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: November 2, 2022
 *      Author: Bruce Palmer
 */
#ifndef WTPTA1_HPP
#define WTPTA1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Wtpta1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Wtpta1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Wtpta1Parser()
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
      // HAS_WIND_PITCHCONTROL
      if (!data->getValue(HAS_WIND_PITCHCONTROL,&bval,g_id)) {
        data->addValue(HAS_WIND_PITCHCONTROL, true, g_id);
      } else {
        data->setValue(HAS_WIND_PITCHCONTROL, true, g_id);
      }

      // WIND_PITCHCONTROL
      std::string stmp;
      if (!data->getValue(WIND_PITCHCONTROL, &stmp, g_id)) {
        data->addValue(WIND_PITCHCONTROL, data_struct.model, g_id);
      } else {
        data->setValue(WIND_PITCHCONTROL, data_struct.model, g_id);
      }

      // WIND_ID
      if (!data->getValue(WIND_ID, &stmp, g_id)) {
        data->addValue(WIND_ID, data_struct.gen_id, g_id);
      } else {
        data->setValue(WIND_ID, data_struct.gen_id, g_id);
      }

      // WIND_KIW
      if (!data->getValue(WIND_KIW,&rval,g_id)) {
        data->addValue(WIND_KIW, data_struct.wind_kiw, g_id);
      } else {
        data->setValue(WIND_KIW, data_struct.wind_kiw, g_id);
      }

      // WIND_KPW
      if (!data->getValue(WIND_KPW,&rval,g_id)) {
        data->addValue(WIND_KPW, data_struct.wind_kpw, g_id);
      } else {
        data->setValue(WIND_KPW, data_struct.wind_kpw, g_id);
      }

      // WIND_KIC
      if (!data->getValue(WIND_KIC,&rval,g_id)) {
        data->addValue(WIND_KIC, data_struct.wind_kic, g_id);
      } else {
        data->setValue(WIND_KIC, data_struct.wind_kic, g_id);
      }

      // WIND_KPC
      if (!data->getValue(WIND_KPC,&rval,g_id)) {
        data->addValue(WIND_KPC, data_struct.wind_kpc, g_id);
      } else {
        data->setValue(WIND_KPC, data_struct.wind_kpc, g_id);
      }

      // WIND_KCC
      if (!data->getValue(WIND_KCC,&rval,g_id)) {
        data->addValue(WIND_KCC, data_struct.wind_kcc, g_id);
      } else {
        data->setValue(WIND_KCC, data_struct.wind_kcc, g_id);
      }

      // WIND_BR_TP
      if (!data->getValue(WIND_BR_TP,&rval,g_id)) {
        data->addValue(WIND_BR_TP, data_struct.wind_tp, g_id);
      } else {
        data->setValue(WIND_BR_TP, data_struct.wind_tp, g_id);
      }

      // WIND_TETAMAX
      if (!data->getValue(WIND_TETAMAX,&rval,g_id)) {
        data->addValue(WIND_TETAMAX, data_struct.wind_tetamax, g_id);
      } else {
        data->setValue(WIND_TETAMAX, data_struct.wind_tetamax, g_id);
      }

      // WIND_TETAMIN
      if (!data->getValue(WIND_TETAMIN,&rval,g_id)) {
        data->addValue(WIND_TETAMIN, data_struct.wind_tetamin, g_id);
      } else {
        data->setValue(WIND_TETAMIN, data_struct.wind_tetamin, g_id);
      }

      // WIND_RTETAMAX
      if (!data->getValue(WIND_RTETAMAX,&rval,g_id)) {
        data->addValue(WIND_RTETAMAX, data_struct.wind_rtetamax, g_id);
      } else {
        data->setValue(WIND_RTETAMAX, data_struct.wind_rtetamax, g_id);
      }

      // WIND_RTETAMIN
      if (!data->getValue(WIND_RTETAMIN,&rval,g_id)) {
        data->addValue(WIND_RTETAMIN, data_struct.wind_rtetamin, g_id);
      } else {
        data->setValue(WIND_RTETAMIN, data_struct.wind_rtetamin, g_id);
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
      // HAS_WIND_PITCHCONTROL
      if (!data->getValue(HAS_WIND_PITCHCONTROL,&bval,g_id)) {
        data->addValue(HAS_WIND_PITCHCONTROL, true, g_id);
      } else {
        data->setValue(HAS_WIND_PITCHCONTROL, true, g_id);
      }

      // WIND_NAME
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(WIND_PITCHCONTROL,&stmp,g_id)) {
        data->addValue(WIND_PITCHCONTROL, model.c_str(), g_id);
      } else {
        data->setValue(WIND_PITCHCONTROL, model.c_str(), g_id);
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
        if (!data->getValue(WIND_KIW,&rval,g_id)) {
          data->addValue(WIND_KIW, atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(WIND_KIW, atof(split_line[3].c_str()), g_id);
        }
      }

      if (nstr > 4) {
        if (!data->getValue(WIND_KPW,&rval,g_id)) {
          data->addValue(WIND_KPW, atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(WIND_KPW, atof(split_line[4].c_str()), g_id);
        }
      }

      if (nstr > 5) {
        if (!data->getValue(WIND_KIC,&rval,g_id)) {
          data->addValue(WIND_KIC, atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(WIND_KIC, atof(split_line[5].c_str()), g_id);
        }
      }

      if (nstr > 6) {
        if (!data->getValue(WIND_KPC,&rval,g_id)) {
          data->addValue(WIND_KPC, atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(WIND_KPC, atof(split_line[6].c_str()), g_id);
        }
      }

      if (nstr > 7) {
        if (!data->getValue(WIND_KCC,&rval,g_id)) {
          data->addValue(WIND_KCC, atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(WIND_KCC, atof(split_line[7].c_str()), g_id);
        }
      }

      if (nstr > 8) {
        if (!data->getValue(WIND_BR_TP,&rval,g_id)) {
          data->addValue(WIND_BR_TP, atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(WIND_BR_TP, atof(split_line[8].c_str()), g_id);
        }
      }

      if (nstr > 9) {
        if (!data->getValue(WIND_TETAMAX,&rval,g_id)) {
          data->addValue(WIND_TETAMAX, atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(WIND_TETAMAX, atof(split_line[9].c_str()), g_id);
        }
      }

      if (nstr > 10) {
        if (!data->getValue(WIND_TETAMIN,&rval,g_id)) {
          data->addValue(WIND_TETAMIN, atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(WIND_TETAMIN, atof(split_line[10].c_str()), g_id);
        }
      }

      if (nstr > 11) {
        if (!data->getValue(WIND_RTETAMAX,&rval,g_id)) {
          data->addValue(WIND_RTETAMAX, atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(WIND_RTETAMAX, atof(split_line[11].c_str()), g_id);
        }
      }

      if (nstr > 12) {
        if (!data->getValue(WIND_RTETAMIN,&rval,g_id)) {
          data->addValue(WIND_RTETAMIN, atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(WIND_RTETAMIN, atof(split_line[12].c_str()), g_id);
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

      // WIND_PITCHCONTROL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size();

      // Wind model parameters

      if (nstr > 3) {
        data.wind_kiw = atof(split_line[3].c_str());
      }

      if (nstr > 4) {
        data.wind_kpw = atof(split_line[4].c_str());
      }

      if (nstr > 5) {
        data.wind_kic = atof(split_line[5].c_str());
      }

      if (nstr > 6) {
        data.wind_kpc = atof(split_line[6].c_str());
      }

      if (nstr > 7) {
        data.wind_kcc = atof(split_line[7].c_str());
      }

      if (nstr > 8) {
        data.wind_tp = atof(split_line[8].c_str());
      }

      if (nstr > 9) {
        data.wind_tetamax = atof(split_line[9].c_str());
      }

      if (nstr > 10) {
        data.wind_tetamin = atof(split_line[10].c_str());
      }

      if (nstr > 11) {
        data.wind_rtetamax = atof(split_line[11].c_str());
      }

      if (nstr > 12) {
        data.wind_rtetamin = atof(split_line[12].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
