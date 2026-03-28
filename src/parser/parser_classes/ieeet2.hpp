/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: 2026-03-27
 *      Author: Yousu Chen
 *
 *  IEEET2 exciter parser.
 *  DYR format: I, 'IEEET2', ID, TR, KA, TA, VRMAX, VRMIN, KE, TE, KF, TF1, TF2,
 *                           E1, SE1, E2, SE2
 *  Indices:    0    1       2   3   4   5   6      7      8   9   10  11   12
 *                          13  14  15  16
 *
 *  Same as IEEET1 except:
 *  - [12] = TF2  (second rate-feedback time constant, not SWITCH)
 *  - VR limits scale with terminal voltage Vt (handled in model, not parser)
 */
#ifndef IEEET2_HPP
#define IEEET2_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Ieeet2Parser
{
  public:
    explicit Ieeet2Parser() {}
    virtual ~Ieeet2Parser() {}

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
#define E2EXT(KEY, FIELD) \
      if (!data->getValue(KEY,&rval,g_id)) { \
        data->addValue(KEY, data_struct.FIELD, g_id); \
      } else { \
        data->setValue(KEY, data_struct.FIELD, g_id); \
      }
      E2EXT(EXCITER_TR,   ex_tr)
      E2EXT(EXCITER_KA,   ex_ka)
      E2EXT(EXCITER_TA,   ex_ta)
      E2EXT(EXCITER_VRMAX, vrmax)
      E2EXT(EXCITER_VRMIN, vrmin)
      E2EXT(EXCITER_KE,   ex_ke)
      E2EXT(EXCITER_TE,   ex_te)
      E2EXT(EXCITER_KF,   ex_kf)
      E2EXT(EXCITER_TF1,  tf1)
      E2EXT(EXCITER_TF2,  tf2)
      E2EXT(EXCITER_E1,   ex_e1)
      E2EXT(EXCITER_SE1,  se1)
      E2EXT(EXCITER_E2,   ex_e2)
      E2EXT(EXCITER_SE2,  se2)
#undef E2EXT
    }

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
#define E2PRS(IDX, KEY) \
      if (nstr > IDX) { \
        if (!data->getValue(KEY,&rval,g_id)) { \
          data->addValue(KEY, atof(split_line[IDX].c_str()), g_id); \
        } else { \
          data->setValue(KEY, atof(split_line[IDX].c_str()), g_id); \
        } \
      }
      E2PRS( 3, EXCITER_TR)
      E2PRS( 4, EXCITER_KA)
      E2PRS( 5, EXCITER_TA)
      E2PRS( 6, EXCITER_VRMAX)
      E2PRS( 7, EXCITER_VRMIN)
      E2PRS( 8, EXCITER_KE)
      E2PRS( 9, EXCITER_TE)
      E2PRS(10, EXCITER_KF)
      E2PRS(11, EXCITER_TF1)
      E2PRS(12, EXCITER_TF2)
      E2PRS(13, EXCITER_E1)
      E2PRS(14, EXCITER_SE1)
      E2PRS(15, EXCITER_E2)
      E2PRS(16, EXCITER_SE2)
#undef E2PRS
    }

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
      if (nstr >  3) data.ex_tr  = atof(split_line[3].c_str());
      if (nstr >  4) data.ex_ka  = atof(split_line[4].c_str());
      if (nstr >  5) data.ex_ta  = atof(split_line[5].c_str());
      if (nstr >  6) data.vrmax  = atof(split_line[6].c_str());
      if (nstr >  7) data.vrmin  = atof(split_line[7].c_str());
      if (nstr >  8) data.ex_ke  = atof(split_line[8].c_str());
      if (nstr >  9) data.ex_te  = atof(split_line[9].c_str());
      if (nstr > 10) data.ex_kf  = atof(split_line[10].c_str());
      if (nstr > 11) data.tf1    = atof(split_line[11].c_str());
      if (nstr > 12) data.tf2    = atof(split_line[12].c_str());
      if (nstr > 13) data.ex_e1  = atof(split_line[13].c_str());
      if (nstr > 14) data.se1    = atof(split_line[14].c_str());
      if (nstr > 15) data.ex_e2  = atof(split_line[15].c_str());
      if (nstr > 16) data.se2    = atof(split_line[16].c_str());
    }
};
}  // parser
}  // gridpack
#endif
