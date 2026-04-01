/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: 2026-03-27
 *      Author: Yousu Chen
 *
 *  SCRX — Simple Controllable Rectifier Exciter
 *
 *  DYR format:
 *    I, 'SCRX', ID, TATB, TB, K, TE, EMIN, EMAX, CSWITCH, RC_RFD
 *
 *  Indices:
 *    [3] = TATB (= TA/TB)
 *    [4] = TB
 *    [5] = K
 *    [6] = TE
 *    [7] = EMIN
 *    [8] = EMAX
 *    [9] = CSWITCH  (0=bus-fed, 1=solid-fed)
 *    [10] = RC_RFD   (Rc/Rfd ratio, >=0; 0=no field current limit)
 */
#ifndef SCRX_HPP
#define SCRX_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class ScrxParser
{
  public:
    /**
     * Constructor
     */
    explicit ScrxParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~ScrxParser()
    {
    }

    /**
     * Extract data from _data_struct and store it in data collection object
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

      // EXCITER_SWITCH (CSWITCH: 0=bus-fed, 1=solid-fed)
      if (!data->getValue(EXCITER_SWITCH,&rval,g_id)) {
        data->addValue(EXCITER_SWITCH, data_struct.rswitch, g_id);
      } else {
        data->setValue(EXCITER_SWITCH, data_struct.rswitch, g_id);
      }
    }

    /**
     * Parser list of strings and store results in data collection object
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

      // EXCITER_TA_OVER_TB [3]
      if (nstr > 3) {
        if (!data->getValue(EXCITER_TA_OVER_TB,&rval,g_id)) {
          data->addValue(EXCITER_TA_OVER_TB, atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TA_OVER_TB, atof(split_line[3].c_str()), g_id);
        }
      }

      // EXCITER_TB [4]
      if (nstr > 4) {
        if (!data->getValue(EXCITER_TB,&rval,g_id)) {
          data->addValue(EXCITER_TB, atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TB, atof(split_line[4].c_str()), g_id);
        }
      }

      // EXCITER_K [5]
      if (nstr > 5) {
        if (!data->getValue(EXCITER_K,&rval,g_id)) {
          data->addValue(EXCITER_K, atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(EXCITER_K, atof(split_line[5].c_str()), g_id);
        }
      }

      // EXCITER_TE [6]
      if (nstr > 6) {
        if (!data->getValue(EXCITER_TE,&rval,g_id)) {
          data->addValue(EXCITER_TE, atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TE, atof(split_line[6].c_str()), g_id);
        }
      }

      // EXCITER_EMIN [7]
      if (nstr > 7) {
        if (!data->getValue(EXCITER_EMIN,&rval,g_id)) {
          data->addValue(EXCITER_EMIN, atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(EXCITER_EMIN, atof(split_line[7].c_str()), g_id);
        }
      } else {
        if (!data->getValue(EXCITER_EMIN,&rval,g_id)) {
          data->addValue(EXCITER_EMIN, 0.0, g_id);
        }
      }

      // EXCITER_EMAX [8]
      if (nstr > 8) {
        if (!data->getValue(EXCITER_EMAX,&rval,g_id)) {
          data->addValue(EXCITER_EMAX, atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(EXCITER_EMAX, atof(split_line[8].c_str()), g_id);
        }
      } else {
        if (!data->getValue(EXCITER_EMAX,&rval,g_id)) {
          data->addValue(EXCITER_EMAX, 999.0, g_id);
        }
      }

      // EXCITER_SWITCH (CSWITCH) [9]
      if (nstr > 9) {
        if (!data->getValue(EXCITER_SWITCH,&rval,g_id)) {
          data->addValue(EXCITER_SWITCH, atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(EXCITER_SWITCH, atof(split_line[9].c_str()), g_id);
        }
      } else {
        // Default: solid-fed (CSWITCH=1), behaves like SEXS
        if (!data->getValue(EXCITER_SWITCH,&rval,g_id)) {
          data->addValue(EXCITER_SWITCH, 1.0, g_id);
        }
      }

      // EXCITER_RC_RFD (Rc/Rfd) [10]
      if (nstr > 10) {
        if (!data->getValue(EXCITER_RC_RFD,&rval,g_id)) {
          data->addValue(EXCITER_RC_RFD, atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(EXCITER_RC_RFD, atof(split_line[10].c_str()), g_id);
        }
      } else {
        if (!data->getValue(EXCITER_RC_RFD,&rval,g_id)) {
          data->addValue(EXCITER_RC_RFD, 0.0, g_id);
        }
      }
    }

    /**
     * Parse list of strings store results in data_struct object
     */
    void store(std::vector<std::string> &split_line, _data_struct &data)
    {
      int o_idx = atoi(split_line[0].c_str());
      data.bus_id = o_idx;

      gridpack::utility::StringUtils util;
      std::string tag = util.clean2Char(split_line[2]);
      strcpy(data.gen_id, tag.c_str());

      std::string sval = util.trimQuotes(split_line[1]);
      util.toUpper(sval);
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size();

      // TATB [3]
      if (nstr > 3) data.ex_ta_over_tb = atof(split_line[3].c_str());

      // TB [4]
      if (nstr > 4) data.ex_tb = atof(split_line[4].c_str());

      // K [5]
      if (nstr > 5) data.ex_k = atof(split_line[5].c_str());

      // TE [6]
      if (nstr > 6) data.ex_te = atof(split_line[6].c_str());

      // EMIN [7]
      if (nstr > 7) {
        data.ex_emin = atof(split_line[7].c_str());
      } else {
        data.ex_emin = 0.0;
      }

      // EMAX [8]
      if (nstr > 8) {
        data.ex_emax = atof(split_line[8].c_str());
      } else {
        data.ex_emax = 999.0;
      }

      // CSWITCH [9] → rswitch (1=solid-fed default)
      if (nstr > 9) {
        data.rswitch = atof(split_line[9].c_str());
      } else {
        data.rswitch = 1.0;
      }

      // RC_RFD [10] → ex_kc (Rc/Rfd ratio, 0=no limit)
      if (nstr > 10) {
        data.ex_kc = atof(split_line[10].c_str());
      } else {
        data.ex_kc = 0.0;
      }
    }
};
}  // parser
}  // gridpack
#endif
