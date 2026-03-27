/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: March 2026
 *      Author: Yousu Chen
 *
 *  IEESGO steam turbine governor parser.
 *  DYR format: I, 'IEESGO', ID, T1, T2, T3, T4, T5, T6, K1, K2, K3, PMAX, PMIN /
 *  Indices:    [3]=T1, [4]=T2, [5]=T3, [6]=T4, [7]=T5, [8]=T6,
 *              [9]=K1, [10]=K2, [11]=K3, [12]=PMAX, [13]=PMIN
 */
#ifndef IEESGO_HPP
#define IEESGO_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class IeesgoParser
{
  public:
    explicit IeesgoParser() {}
    virtual ~IeesgoParser() {}

    void extract(_data_struct &data_struct,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
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
      GEXTRACT(GOVERNOR_T1,   gv_t1)
      GEXTRACT(GOVERNOR_T2,   gv_t2)
      GEXTRACT(GOVERNOR_T3,   gv_t3)
      GEXTRACT(GOVERNOR_T4,   gv_t4)
      GEXTRACT(GOVERNOR_T5,   gv_t5)
      GEXTRACT(GOVERNOR_T6,   gv_t6)
      GEXTRACT(GOVERNOR_K1,   gv_k1)
      GEXTRACT(GOVERNOR_K2,   gv_k2)
      GEXTRACT(GOVERNOR_K3,   gv_k3)
      GEXTRACT(GOVERNOR_PMAX, pmax)
      GEXTRACT(GOVERNOR_PMIN, pmin)
#undef GEXTRACT
    }

    /**
     * Parse list of strings and store results in data collection object
     * DYR format: I, 'IEESGO', ID, T1, T2, T3, T4, T5, T6, K1, K2, K3, PMAX, PMIN
     * Indices: [3]=T1, [4]=T2, [5]=T3, [6]=T4, [7]=T5, [8]=T6,
     *          [9]=K1, [10]=K2, [11]=K3, [12]=PMAX, [13]=PMIN
     */
    void parse(std::vector<std::string> &split_line,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
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
#define GPARSE(IDX, KEY) \
      if (nstr > IDX) { \
        if (!data->getValue(KEY,&rval,g_id)) { \
          data->addValue(KEY, atof(split_line[IDX].c_str()), g_id); \
        } else { \
          data->setValue(KEY, atof(split_line[IDX].c_str()), g_id); \
        } \
      }
      GPARSE( 3, GOVERNOR_T1)
      GPARSE( 4, GOVERNOR_T2)
      GPARSE( 5, GOVERNOR_T3)
      GPARSE( 6, GOVERNOR_T4)
      GPARSE( 7, GOVERNOR_T5)
      GPARSE( 8, GOVERNOR_T6)
      GPARSE( 9, GOVERNOR_K1)
      GPARSE(10, GOVERNOR_K2)
      GPARSE(11, GOVERNOR_K3)
      GPARSE(12, GOVERNOR_PMAX)
      GPARSE(13, GOVERNOR_PMIN)
#undef GPARSE
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
      if (nstr >  3) data.gv_t1 = atof(split_line[3].c_str());
      if (nstr >  4) data.gv_t2 = atof(split_line[4].c_str());
      if (nstr >  5) data.gv_t3 = atof(split_line[5].c_str());
      if (nstr >  6) data.gv_t4 = atof(split_line[6].c_str());
      if (nstr >  7) data.gv_t5 = atof(split_line[7].c_str());
      if (nstr >  8) data.gv_t6 = atof(split_line[8].c_str());
      if (nstr >  9) data.gv_k1 = atof(split_line[9].c_str());
      if (nstr > 10) data.gv_k2 = atof(split_line[10].c_str());
      if (nstr > 11) data.gv_k3 = atof(split_line[11].c_str());
      if (nstr > 12) data.pmax  = atof(split_line[12].c_str());
      if (nstr > 13) data.pmin  = atof(split_line[13].c_str());
    }
};
}  // parser
}  // gridpack
#endif
