/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: 2026-03-27
 *      Author: Yousu Chen
 *
 *  STAB2A PSS parser.
 *  DYR format: I, 'STAB2A', ID, KT, T, T1, T2, T3, T4, H1, H2
 *  Indices:    0    1       2   3   4   5   6   7   8   9   10
 *
 *  KT  — gain
 *  T   — washout time constant
 *  T1  — lead-lag 1 lead TC
 *  T2  — lead-lag 1 lag TC
 *  T3  — lead-lag 2 lead TC
 *  T4  — lead-lag 2 lag TC
 *  H1  — maximum output
 *  H2  — minimum output
 */
#ifndef STAB2A_HPP
#define STAB2A_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Stab2aParser
{
  public:
    explicit Stab2aParser() {}
    virtual ~Stab2aParser() {}

    void extract(_data_struct &data_struct,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
      if (!data->getValue(HAS_PSS,&bval,g_id)) {
        data->addValue(HAS_PSS, true, g_id);
      } else {
        data->setValue(HAS_PSS, true, g_id);
      }
      std::string stmp;
      if (!data->getValue(PSS_MODEL, &stmp, g_id)) {
        data->addValue(PSS_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(PSS_MODEL, data_struct.model, g_id);
      }
#define S2EXT(KEY, FIELD) \
      if (!data->getValue(KEY,&rval,g_id)) { \
        data->addValue(KEY, data_struct.FIELD, g_id); \
      } else { \
        data->setValue(KEY, data_struct.FIELD, g_id); \
      }
      S2EXT(STAB2A_KT, stab2a_kt)
      S2EXT(STAB2A_T,  stab2a_t)
      S2EXT(STAB2A_T1, stab2a_t1)
      S2EXT(STAB2A_T2, stab2a_t2)
      S2EXT(STAB2A_T3, stab2a_t3)
      S2EXT(STAB2A_T4, stab2a_t4)
      S2EXT(STAB2A_H1, stab2a_h1)
      S2EXT(STAB2A_H2, stab2a_h2)
#undef S2EXT
    }

    void parse(std::vector<std::string> &split_line,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
      int nstr = split_line.size();
      if (!data->getValue(HAS_PSS,&bval,g_id)) {
        data->addValue(HAS_PSS, true, g_id);
      } else {
        data->setValue(HAS_PSS, true, g_id);
      }
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(PSS_MODEL,&stmp,g_id)) {
        data->addValue(PSS_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(PSS_MODEL, model.c_str(), g_id);
      }
#define S2PRS(IDX, KEY) \
      if (nstr > IDX) { \
        if (!data->getValue(KEY,&rval,g_id)) { \
          data->addValue(KEY, atof(split_line[IDX].c_str()), g_id); \
        } else { \
          data->setValue(KEY, atof(split_line[IDX].c_str()), g_id); \
        } \
      }
      S2PRS( 3, STAB2A_KT)
      S2PRS( 4, STAB2A_T)
      S2PRS( 5, STAB2A_T1)
      S2PRS( 6, STAB2A_T2)
      S2PRS( 7, STAB2A_T3)
      S2PRS( 8, STAB2A_T4)
      S2PRS( 9, STAB2A_H1)
      S2PRS(10, STAB2A_H2)
#undef S2PRS
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
      if (nstr >  3) data.stab2a_kt = atof(split_line[3].c_str());
      if (nstr >  4) data.stab2a_t  = atof(split_line[4].c_str());
      if (nstr >  5) data.stab2a_t1 = atof(split_line[5].c_str());
      if (nstr >  6) data.stab2a_t2 = atof(split_line[6].c_str());
      if (nstr >  7) data.stab2a_t3 = atof(split_line[7].c_str());
      if (nstr >  8) data.stab2a_t4 = atof(split_line[8].c_str());
      if (nstr >  9) data.stab2a_h1 = atof(split_line[9].c_str());
      if (nstr > 10) data.stab2a_h2 = atof(split_line[10].c_str());
    }
};
}  // parser
}  // gridpack
#endif
