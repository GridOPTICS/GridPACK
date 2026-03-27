/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: March 2026
 *      Author: Yousu Chen
 *
 *  EXST1 exciter parser.
 *  DYR format: I, 'EXST1', ID, TR, VIMAX, VIMIN, TC, TB, KA, TA, VRMAX, VRMIN, KC, KF, TF
 *  Indices:    0    1      2   3    4      5      6   7   8   9   10     11    12  13  14
 */
#ifndef EXST1_HPP
#define EXST1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Exst1Parser
{
  public:
    explicit Exst1Parser() {}
    virtual ~Exst1Parser() {}

    /**
     * Extract data from _data_struct and store it in data collection object
     */
    void extract(_data_struct &data_struct,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
      if (!data->getValue(HAS_EXCITER,&bval,g_id)) {
        data->addValue(HAS_EXCITER, true, g_id);
      } else {
        data->setValue(HAS_EXCITER, true, g_id);
      }
      std::string stmp;
      if (!data->getValue(EXCITER_MODEL, &stmp, g_id)) {
        data->addValue(EXCITER_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(EXCITER_MODEL, data_struct.model, g_id);
      }
      // TR
      if (!data->getValue(EXCITER_TR,&rval,g_id)) {
        data->addValue(EXCITER_TR, data_struct.ex_tr, g_id);
      } else {
        data->setValue(EXCITER_TR, data_struct.ex_tr, g_id);
      }
      // VIMAX
      if (!data->getValue(EXCITER_VIMAX,&rval,g_id)) {
        data->addValue(EXCITER_VIMAX, data_struct.vimax, g_id);
      } else {
        data->setValue(EXCITER_VIMAX, data_struct.vimax, g_id);
      }
      // VIMIN
      if (!data->getValue(EXCITER_VIMIN,&rval,g_id)) {
        data->addValue(EXCITER_VIMIN, data_struct.vimin, g_id);
      } else {
        data->setValue(EXCITER_VIMIN, data_struct.vimin, g_id);
      }
      // TC
      if (!data->getValue(EXCITER_TC,&rval,g_id)) {
        data->addValue(EXCITER_TC, data_struct.ex_tc, g_id);
      } else {
        data->setValue(EXCITER_TC, data_struct.ex_tc, g_id);
      }
      // TB
      if (!data->getValue(EXCITER_TB,&rval,g_id)) {
        data->addValue(EXCITER_TB, data_struct.ex_tb, g_id);
      } else {
        data->setValue(EXCITER_TB, data_struct.ex_tb, g_id);
      }
      // KA
      if (!data->getValue(EXCITER_KA,&rval,g_id)) {
        data->addValue(EXCITER_KA, data_struct.ex_ka, g_id);
      } else {
        data->setValue(EXCITER_KA, data_struct.ex_ka, g_id);
      }
      // TA
      if (!data->getValue(EXCITER_TA,&rval,g_id)) {
        data->addValue(EXCITER_TA, data_struct.ex_ta, g_id);
      } else {
        data->setValue(EXCITER_TA, data_struct.ex_ta, g_id);
      }
      // VRMAX
      if (!data->getValue(EXCITER_VRMAX,&rval,g_id)) {
        data->addValue(EXCITER_VRMAX, data_struct.vrmax, g_id);
      } else {
        data->setValue(EXCITER_VRMAX, data_struct.vrmax, g_id);
      }
      // VRMIN
      if (!data->getValue(EXCITER_VRMIN,&rval,g_id)) {
        data->addValue(EXCITER_VRMIN, data_struct.vrmin, g_id);
      } else {
        data->setValue(EXCITER_VRMIN, data_struct.vrmin, g_id);
      }
      // KC
      if (!data->getValue(EXCITER_KC,&rval,g_id)) {
        data->addValue(EXCITER_KC, data_struct.ex_kc, g_id);
      } else {
        data->setValue(EXCITER_KC, data_struct.ex_kc, g_id);
      }
      // KF
      if (!data->getValue(EXCITER_KF,&rval,g_id)) {
        data->addValue(EXCITER_KF, data_struct.ex_kf, g_id);
      } else {
        data->setValue(EXCITER_KF, data_struct.ex_kf, g_id);
      }
      // TF (stored in TF1 slot)
      if (!data->getValue(EXCITER_TF1,&rval,g_id)) {
        data->addValue(EXCITER_TF1, data_struct.ex_tf, g_id);
      } else {
        data->setValue(EXCITER_TF1, data_struct.ex_tf, g_id);
      }
    }

    /**
     * Parse list of strings and store results in data collection object
     */
    void parse(std::vector<std::string> &split_line,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
      int nstr = split_line.size();
      if (!data->getValue(HAS_EXCITER,&bval,g_id)) {
        data->addValue(HAS_EXCITER, true, g_id);
      } else {
        data->setValue(HAS_EXCITER, true, g_id);
      }
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(EXCITER_MODEL,&stmp,g_id)) {
        data->addValue(EXCITER_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(EXCITER_MODEL, model.c_str(), g_id);
      }
      // [3] TR
      if (nstr > 3) {
        if (!data->getValue(EXCITER_TR,&rval,g_id)) {
          data->addValue(EXCITER_TR, atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TR, atof(split_line[3].c_str()), g_id);
        }
      }
      // [4] VIMAX
      if (nstr > 4) {
        if (!data->getValue(EXCITER_VIMAX,&rval,g_id)) {
          data->addValue(EXCITER_VIMAX, atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VIMAX, atof(split_line[4].c_str()), g_id);
        }
      }
      // [5] VIMIN
      if (nstr > 5) {
        if (!data->getValue(EXCITER_VIMIN,&rval,g_id)) {
          data->addValue(EXCITER_VIMIN, atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VIMIN, atof(split_line[5].c_str()), g_id);
        }
      }
      // [6] TC
      if (nstr > 6) {
        if (!data->getValue(EXCITER_TC,&rval,g_id)) {
          data->addValue(EXCITER_TC, atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TC, atof(split_line[6].c_str()), g_id);
        }
      }
      // [7] TB
      if (nstr > 7) {
        if (!data->getValue(EXCITER_TB,&rval,g_id)) {
          data->addValue(EXCITER_TB, atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TB, atof(split_line[7].c_str()), g_id);
        }
      }
      // [8] KA
      if (nstr > 8) {
        if (!data->getValue(EXCITER_KA,&rval,g_id)) {
          data->addValue(EXCITER_KA, atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KA, atof(split_line[8].c_str()), g_id);
        }
      }
      // [9] TA
      if (nstr > 9) {
        if (!data->getValue(EXCITER_TA,&rval,g_id)) {
          data->addValue(EXCITER_TA, atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TA, atof(split_line[9].c_str()), g_id);
        }
      }
      // [10] VRMAX
      if (nstr > 10) {
        if (!data->getValue(EXCITER_VRMAX,&rval,g_id)) {
          data->addValue(EXCITER_VRMAX, atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VRMAX, atof(split_line[10].c_str()), g_id);
        }
      }
      // [11] VRMIN
      if (nstr > 11) {
        if (!data->getValue(EXCITER_VRMIN,&rval,g_id)) {
          data->addValue(EXCITER_VRMIN, atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VRMIN, atof(split_line[11].c_str()), g_id);
        }
      }
      // [12] KC
      if (nstr > 12) {
        if (!data->getValue(EXCITER_KC,&rval,g_id)) {
          data->addValue(EXCITER_KC, atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KC, atof(split_line[12].c_str()), g_id);
        }
      }
      // [13] KF
      if (nstr > 13) {
        if (!data->getValue(EXCITER_KF,&rval,g_id)) {
          data->addValue(EXCITER_KF, atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KF, atof(split_line[13].c_str()), g_id);
        }
      }
      // [14] TF (stored in TF1 slot)
      if (nstr > 14) {
        if (!data->getValue(EXCITER_TF1,&rval,g_id)) {
          data->addValue(EXCITER_TF1, atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TF1, atof(split_line[14].c_str()), g_id);
        }
      }
    }

    /**
     * Parse list of strings and store results in data_struct object
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
      if (nstr > 3)  data.ex_tr  = atof(split_line[3].c_str());
      if (nstr > 4)  data.vimax  = atof(split_line[4].c_str());
      if (nstr > 5)  data.vimin  = atof(split_line[5].c_str());
      if (nstr > 6)  data.ex_tc  = atof(split_line[6].c_str());
      if (nstr > 7)  data.ex_tb  = atof(split_line[7].c_str());
      if (nstr > 8)  data.ex_ka  = atof(split_line[8].c_str());
      if (nstr > 9)  data.ex_ta  = atof(split_line[9].c_str());
      if (nstr > 10) data.vrmax  = atof(split_line[10].c_str());
      if (nstr > 11) data.vrmin  = atof(split_line[11].c_str());
      if (nstr > 12) data.ex_kc  = atof(split_line[12].c_str());
      if (nstr > 13) data.ex_kf  = atof(split_line[13].c_str());
      if (nstr > 14) data.ex_tf  = atof(split_line[14].c_str());
    }
};
}  // parser
}  // gridpack
#endif
