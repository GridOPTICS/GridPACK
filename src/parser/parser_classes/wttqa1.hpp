/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: November 2, 2022
 *      Author: Bruce Palmer
 */
#ifndef WTTQA1_HPP
#define WTTQA1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Wttqa1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Wttqa1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Wttqa1Parser()
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
      // HAS_WIND_TORQUECONTROL
      if (!data->getValue(HAS_WIND_TORQUECONTROL,&bval,g_id)) {
        data->addValue(HAS_WIND_TORQUECONTROL, true, g_id);
      } else {
        data->setValue(HAS_WIND_TORQUECONTROL, true, g_id);
      }

      // WIND_TORQUECONTROL
      std::string stmp;
      if (!data->getValue(WIND_TORQUECONTROL, &stmp, g_id)) {
        data->addValue(WIND_TORQUECONTROL, data_struct.model, g_id);
      } else {
        data->setValue(WIND_TORQUECONTROL, data_struct.model, g_id);
      }

      // WIND_ID
      if (!data->getValue(WIND_ID, &stmp, g_id)) {
        data->addValue(WIND_ID, data_struct.gen_id, g_id);
      } else {
        data->setValue(WIND_ID, data_struct.gen_id, g_id);
      }

      // WIND_TFLAG
      if (!data->getValue(WIND_TC_TFLAG,&ival,g_id)) {
        data->addValue(WIND_TC_TFLAG, data_struct.wind_tflag, g_id);
      } else {
        data->setValue(WIND_TC_TFLAG, data_struct.wind_tflag, g_id);
      }

      // WIND_KPP
      if (!data->getValue(WIND_TC_KPP,&rval,g_id)) {
        data->addValue(WIND_TC_KPP, data_struct.wind_kpp, g_id);
      } else {
        data->setValue(WIND_TC_KPP, data_struct.wind_kpp, g_id);
      }

      // WIND_KIP
      if (!data->getValue(WIND_TC_KIP,&rval,g_id)) {
        data->addValue(WIND_TC_KIP, data_struct.wind_kip, g_id);
      } else {
        data->setValue(WIND_TC_KIP, data_struct.wind_kip, g_id);
      }

      // WIND_PF_TP
      if (!data->getValue(WIND_TC_TP,&rval,g_id)) {
        data->addValue(WIND_TC_TP, data_struct.wind_tp, g_id);
      } else {
        data->setValue(WIND_TC_TP, data_struct.wind_tp, g_id);
      }

      // WIND_TWREF
      if (!data->getValue(WIND_TC_TWREF,&rval,g_id)) {
        data->addValue(WIND_TC_TWREF, data_struct.wind_twref, g_id);
      } else {
        data->setValue(WIND_TC_TWREF, data_struct.wind_twref, g_id);
      }

      // WIND_TEMAX
      if (!data->getValue(WIND_TC_TEMAX,&rval,g_id)) {
        data->addValue(WIND_TC_TEMAX, data_struct.wind_temax, g_id);
      } else {
        data->setValue(WIND_TC_TEMAX, data_struct.wind_temax, g_id);
      }

      // WIND_TEMIN
      if (!data->getValue(WIND_TC_TEMIN,&rval,g_id)) {
        data->addValue(WIND_TC_TEMIN, data_struct.wind_temin, g_id);
      } else {
        data->setValue(WIND_TC_TEMIN, data_struct.wind_temin, g_id);
      }

      // WIND_P1
      if (!data->getValue(WIND_TC_P1,&rval,g_id)) {
        data->addValue(WIND_TC_P1, data_struct.wind_p1, g_id);
      } else {
        data->setValue(WIND_TC_P1, data_struct.wind_p1, g_id);
      }

      // WIND_SPD1
      if (!data->getValue(WIND_TC_SPD1,&rval,g_id)) {
        data->addValue(WIND_TC_SPD1, data_struct.wind_spd1, g_id);
      } else {
        data->setValue(WIND_TC_SPD1, data_struct.wind_spd1, g_id);
      }

      // WIND_P2
      if (!data->getValue(WIND_TC_P2,&rval,g_id)) {
        data->addValue(WIND_TC_P2, data_struct.wind_p2, g_id);
      } else {
        data->setValue(WIND_TC_P2, data_struct.wind_p2, g_id);
      }

      // WIND_SPD2
      if (!data->getValue(WIND_TC_SPD2,&rval,g_id)) {
        data->addValue(WIND_TC_SPD2, data_struct.wind_spd2, g_id);
      } else {
        data->setValue(WIND_TC_SPD2, data_struct.wind_spd2, g_id);
      }

      // WIND_P3
      if (!data->getValue(WIND_TC_P3,&rval,g_id)) {
        data->addValue(WIND_TC_P3, data_struct.wind_p3, g_id);
      } else {
        data->setValue(WIND_TC_P3, data_struct.wind_p3, g_id);
      }

      // WIND_SPD3
      if (!data->getValue(WIND_TC_SPD3,&rval,g_id)) {
        data->addValue(WIND_TC_SPD3, data_struct.wind_spd3, g_id);
      } else {
        data->setValue(WIND_TC_SPD3, data_struct.wind_spd3, g_id);
      }

      // WIND_P4
      if (!data->getValue(WIND_TC_P4,&rval,g_id)) {
        data->addValue(WIND_TC_P4, data_struct.wind_p4, g_id);
      } else {
        data->setValue(WIND_TC_P4, data_struct.wind_p4, g_id);
      }

      // WIND_SPD4
      if (!data->getValue(WIND_TC_SPD4,&rval,g_id)) {
        data->addValue(WIND_TC_SPD4, data_struct.wind_spd4, g_id);
      } else {
        data->setValue(WIND_TC_SPD4, data_struct.wind_spd4, g_id);
      }

      // WIND_TRATE
      if (!data->getValue(WIND_TC_TRATE,&rval,g_id)) {
        data->addValue(WIND_TC_TRATE, data_struct.wind_trate, g_id);
      } else {
        data->setValue(WIND_TC_TRATE, data_struct.wind_trate, g_id);
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
      // HAS_WIND_TORQUECONTROL
      if (!data->getValue(HAS_WIND_TORQUECONTROL,&bval,g_id)) {
        data->addValue(HAS_WIND_TORQUECONTROL, true, g_id);
      } else {
        data->setValue(HAS_WIND_TORQUECONTROL, true, g_id);
      }

      // WIND_NAME
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(WIND_TORQUECONTROL,&stmp,g_id)) {
        data->addValue(WIND_TORQUECONTROL, model.c_str(), g_id);
      } else {
        data->setValue(WIND_TORQUECONTROL, model.c_str(), g_id);
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
        if (!data->getValue(WIND_TC_TFLAG,&ival,g_id)) {
          data->addValue(WIND_TC_TFLAG, atoi(split_line[3].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_TFLAG, atoi(split_line[3].c_str()), g_id);
        }
      }

      if (nstr > 4) {
        if (!data->getValue(WIND_TC_KPP,&rval,g_id)) {
          data->addValue(WIND_TC_KPP, atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_KPP, atof(split_line[4].c_str()), g_id);
        }
      }

      if (nstr > 5) {
        if (!data->getValue(WIND_TC_KIP,&rval,g_id)) {
          data->addValue(WIND_TC_KIP, atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_KIP, atof(split_line[5].c_str()), g_id);
        }
      }

      if (nstr > 6) {
        if (!data->getValue(WIND_TC_TP,&rval,g_id)) {
          data->addValue(WIND_TC_TP, atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_TP, atof(split_line[6].c_str()), g_id);
        }
      }

      if (nstr > 7) {
        if (!data->getValue(WIND_TC_TWREF,&rval,g_id)) {
          data->addValue(WIND_TC_TWREF, atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_TWREF, atof(split_line[7].c_str()), g_id);
        }
      }

      if (nstr > 8) {
        if (!data->getValue(WIND_TC_TEMAX,&rval,g_id)) {
          data->addValue(WIND_TC_TEMAX, atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_TEMAX, atof(split_line[8].c_str()), g_id);
        }
      }

      if (nstr > 9) {
        if (!data->getValue(WIND_TC_TEMIN,&rval,g_id)) {
          data->addValue(WIND_TC_TEMIN, atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_TEMIN, atof(split_line[9].c_str()), g_id);
        }
      }

      if (nstr > 10) {
        if (!data->getValue(WIND_TC_P1,&rval,g_id)) {
          data->addValue(WIND_TC_P1, atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_P1, atof(split_line[10].c_str()), g_id);
        }
      }

      if (nstr > 11) {
        if (!data->getValue(WIND_TC_SPD1,&rval,g_id)) {
          data->addValue(WIND_TC_SPD1, atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_SPD1, atof(split_line[11].c_str()), g_id);
        }
      }

      if (nstr > 12) {
        if (!data->getValue(WIND_TC_P2,&rval,g_id)) {
          data->addValue(WIND_TC_P2, atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_P2, atof(split_line[12].c_str()), g_id);
        }
      }

      if (nstr > 13) {
        if (!data->getValue(WIND_TC_SPD2,&rval,g_id)) {
          data->addValue(WIND_TC_SPD2, atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_SPD2, atof(split_line[13].c_str()), g_id);
        }
      }

      if (nstr > 14) {
        if (!data->getValue(WIND_TC_P3,&rval,g_id)) {
          data->addValue(WIND_TC_P3, atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_P3, atof(split_line[14].c_str()), g_id);
        }
      }

      if (nstr > 15) {
        if (!data->getValue(WIND_TC_SPD3,&rval,g_id)) {
          data->addValue(WIND_TC_SPD3, atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_SPD3, atof(split_line[15].c_str()), g_id);
        }
      }

      if (nstr > 16) {
        if (!data->getValue(WIND_TC_P4,&rval,g_id)) {
          data->addValue(WIND_TC_P4, atof(split_line[16].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_P4, atof(split_line[16].c_str()), g_id);
        }
      }

      if (nstr > 17) {
        if (!data->getValue(WIND_TC_SPD4,&rval,g_id)) {
          data->addValue(WIND_TC_SPD4, atof(split_line[17].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_SPD4, atof(split_line[17].c_str()), g_id);
        }
      }

      if (nstr > 18) {
        if (!data->getValue(WIND_TC_TRATE,&rval,g_id)) {
          data->addValue(WIND_TC_TRATE, atof(split_line[18].c_str()), g_id);
        } else {
          data->setValue(WIND_TC_TRATE, atof(split_line[18].c_str()), g_id);
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

      // WIND_TORQUECONTROL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size();

      // Wind model parameters
      if (nstr > 3) {
        data.wind_tflag = atoi(split_line[3].c_str());
      }

      if (nstr > 4) {
        data.wind_kpp = atof(split_line[4].c_str());
      }

      if (nstr > 5) {
        data.wind_kip = atof(split_line[5].c_str());
      }

      if (nstr > 6) {
        data.wind_tp = atof(split_line[6].c_str());
      }

      if (nstr > 7) {
        data.wind_twref = atof(split_line[7].c_str());
      }

      if (nstr > 8) {
        data.wind_temax = atof(split_line[8].c_str());
      }

      if (nstr > 9) {
        data.wind_temin = atof(split_line[9].c_str());
      }

      if (nstr > 10) {
        data.wind_p1 = atof(split_line[10].c_str());
      }

      if (nstr > 11) {
        data.wind_spd1 = atof(split_line[11].c_str());
      }

      if (nstr > 12) {
        data.wind_p2 = atof(split_line[12].c_str());
      }

      if (nstr > 13) {
        data.wind_spd2 = atof(split_line[13].c_str());
      }

      if (nstr > 14) {
        data.wind_p3 = atof(split_line[14].c_str());
      }

      if (nstr > 15) {
        data.wind_spd3 = atof(split_line[15].c_str());
      }

      if (nstr > 16) {
        data.wind_p4 = atof(split_line[16].c_str());
      }

      if (nstr > 17) {
        data.wind_spd4 = atof(split_line[17].c_str());
      }

      if (nstr > 18) {
        data.wind_trate = atof(split_line[18].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
