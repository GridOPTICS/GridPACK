/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: March 2026
 *      Author: Yousu Chen
 *
 *  IEEEG1 steam turbine governor parser.
 *  DYR format: I, 'IEEEG1', ID, K, T1, T2, T3, UO, UC, PMAX, PMIN,
 *                           T4, K1, K2, T5, K3, K4, T6, K5, K6, T7, K7, K8
 *  Indices:    0    1       2   3   4   5   6   7   8   9    10
 *                          11  12  13  14  15  16  17  18  19  20  21  22
 *
 *  IEEEG1 is the base steam governor. It has no deadband (Db1=0) and no
 *  nonlinear valve (NGV = identity). The Wsieg1Model handles both cases
 *  correctly when these extra parameters are absent.
 */
#ifndef IEEEG1_HPP
#define IEEEG1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Ieeeg1Parser
{
  public:
    explicit Ieeeg1Parser() {}
    virtual ~Ieeeg1Parser() {}

    /**
     * Extract data from _data_struct and store it in data collection object
     */
    void extract(_data_struct &data_struct,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
      int ival;
      if (!data->getValue(HAS_GOVERNOR,&bval,g_id)) {
        data->addValue(HAS_GOVERNOR, true, g_id);
      } else {
        data->setValue(HAS_GOVERNOR, true, g_id);
      }
      std::string stmp;
      if (!data->getValue(GOVERNOR_MODEL, &stmp, g_id)) {
        data->addValue(GOVERNOR_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(GOVERNOR_MODEL, data_struct.model, g_id);
      }
#define GEXTRACT(KEY, FIELD) \
      if (!data->getValue(KEY,&rval,g_id)) { \
        data->addValue(KEY, data_struct.FIELD, g_id); \
      } else { \
        data->setValue(KEY, data_struct.FIELD, g_id); \
      }
      GEXTRACT(GOVERNOR_K,    gv_k)
      GEXTRACT(GOVERNOR_T1,   gv_t1)
      GEXTRACT(GOVERNOR_T2,   gv_t2)
      GEXTRACT(GOVERNOR_T3,   gv_t3)
      GEXTRACT(GOVERNOR_UO,   gv_uo)
      GEXTRACT(GOVERNOR_UC,   gv_uc)
      GEXTRACT(GOVERNOR_PMAX, pmax)
      GEXTRACT(GOVERNOR_PMIN, pmin)
      GEXTRACT(GOVERNOR_T4,   gv_t4)
      GEXTRACT(GOVERNOR_K1,   gv_k1)
      GEXTRACT(GOVERNOR_K2,   gv_k2)
      GEXTRACT(GOVERNOR_T5,   gv_t5)
      GEXTRACT(GOVERNOR_K3,   gv_k3)
      GEXTRACT(GOVERNOR_K4,   gv_k4)
      GEXTRACT(GOVERNOR_T6,   gv_t6)
      GEXTRACT(GOVERNOR_K5,   gv_k5)
      GEXTRACT(GOVERNOR_K6,   gv_k6)
      GEXTRACT(GOVERNOR_T7,   gv_t7)
      GEXTRACT(GOVERNOR_K7,   gv_k7)
      GEXTRACT(GOVERNOR_K8,   gv_k8)
#undef GEXTRACT
    }

    /**
     * Parse list of strings and store results in data collection object
     * DYR format: I, 'IEEEG1', ID, JBUS, M, K, T1, T2, T3, UO, UC, PMAX, PMIN,
     *                          T4, K1, K2, T5, K3, K4, T6, K5, K6, T7, K7, K8
     * DYR indices: [3]=JBUS, [4]=M, [5]=K, [6]=T1, [7]=T2, [8]=T3, [9]=UO, [10]=UC,
     *              [11]=PMAX, [12]=PMIN, [13]=T4, [14]=K1, [15]=K2,
     *              [16]=T5, [17]=K3, [18]=K4, [19]=T6, [20]=K5, [21]=K6,
     *              [22]=T7, [23]=K7, [24]=K8
     */
    void parse(std::vector<std::string> &split_line,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
      int ival;
      int nstr = split_line.size();
      if (!data->getValue(HAS_GOVERNOR,&bval,g_id)) {
        data->addValue(HAS_GOVERNOR, true, g_id);
      } else {
        data->setValue(HAS_GOVERNOR, true, g_id);
      }
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(GOVERNOR_MODEL,&stmp,g_id)) {
        data->addValue(GOVERNOR_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(GOVERNOR_MODEL, model.c_str(), g_id);
      }
      if (nstr > 3) {
        if (!data->getValue(GOVERNOR_JBUS,&ival,g_id))
          data->addValue(GOVERNOR_JBUS, atoi(split_line[3].c_str()), g_id);
        else
          data->setValue(GOVERNOR_JBUS, atoi(split_line[3].c_str()), g_id);
      }
      if (nstr > 4) {
        if (!data->getValue(GOVERNOR_M,&ival,g_id))
          data->addValue(GOVERNOR_M, atoi(split_line[4].c_str()), g_id);
        else
          data->setValue(GOVERNOR_M, atoi(split_line[4].c_str()), g_id);
      }
#define GPARSE(IDX, KEY) \
      if (nstr > IDX) { \
        if (!data->getValue(KEY,&rval,g_id)) { \
          data->addValue(KEY, atof(split_line[IDX].c_str()), g_id); \
        } else { \
          data->setValue(KEY, atof(split_line[IDX].c_str()), g_id); \
        } \
      }
      GPARSE( 5, GOVERNOR_K)
      GPARSE( 6, GOVERNOR_T1)
      GPARSE( 7, GOVERNOR_T2)
      GPARSE( 8, GOVERNOR_T3)
      GPARSE( 9, GOVERNOR_UO)
      GPARSE(10, GOVERNOR_UC)
      GPARSE(11, GOVERNOR_PMAX)
      GPARSE(12, GOVERNOR_PMIN)
      GPARSE(13, GOVERNOR_T4)
      GPARSE(14, GOVERNOR_K1)
      GPARSE(15, GOVERNOR_K2)
      GPARSE(16, GOVERNOR_T5)
      GPARSE(17, GOVERNOR_K3)
      GPARSE(18, GOVERNOR_K4)
      GPARSE(19, GOVERNOR_T6)
      GPARSE(20, GOVERNOR_K5)
      GPARSE(21, GOVERNOR_K6)
      GPARSE(22, GOVERNOR_T7)
      GPARSE(23, GOVERNOR_K7)
      GPARSE(24, GOVERNOR_K8)
#undef GPARSE
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
      if (nstr >  3) data.jbus  = atoi(split_line[3].c_str());
      if (nstr >  4) data.gv_m  = atoi(split_line[4].c_str());
      if (nstr >  5) data.gv_k  = atof(split_line[5].c_str());
      if (nstr >  6) data.gv_t1 = atof(split_line[6].c_str());
      if (nstr >  7) data.gv_t2 = atof(split_line[7].c_str());
      if (nstr >  8) data.gv_t3 = atof(split_line[8].c_str());
      if (nstr >  9) data.gv_uo = atof(split_line[9].c_str());
      if (nstr > 10) data.gv_uc = atof(split_line[10].c_str());
      if (nstr > 11) data.pmax  = atof(split_line[11].c_str());
      if (nstr > 12) data.pmin  = atof(split_line[12].c_str());
      if (nstr > 13) data.gv_t4 = atof(split_line[13].c_str());
      if (nstr > 14) data.gv_k1 = atof(split_line[14].c_str());
      if (nstr > 15) data.gv_k2 = atof(split_line[15].c_str());
      if (nstr > 16) data.gv_t5 = atof(split_line[16].c_str());
      if (nstr > 17) data.gv_k3 = atof(split_line[17].c_str());
      if (nstr > 18) data.gv_k4 = atof(split_line[18].c_str());
      if (nstr > 19) data.gv_t6 = atof(split_line[19].c_str());
      if (nstr > 20) data.gv_k5 = atof(split_line[20].c_str());
      if (nstr > 21) data.gv_k6 = atof(split_line[21].c_str());
      if (nstr > 22) data.gv_t7 = atof(split_line[22].c_str());
      if (nstr > 23) data.gv_k7 = atof(split_line[23].c_str());
      if (nstr > 24) data.gv_k8 = atof(split_line[24].c_str());
    }
};
}  // parser
}  // gridpack
#endif
