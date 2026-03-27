/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: March 2026
 *      Author: Yousu Chen
 *
 *  ST2CUT PSS parser (two-input stabilizer).
 *  DYR format: I, 'ST2CUT', ID, MODE, BUSR, MODE2, BUSR2, K1, K2, T1, T2,
 *                           T3, T4, T5, T6, T7, T8, T9, T10,
 *                           LSMAX, LSMIN, VCU, VCL
 *  Indices:    0    1       2   3     4    5     6    7   8   9   10
 *                          11  12  13  14  15  16  17  18
 *                          19    20   21  22
 */
#ifndef ST2CUT_HPP
#define ST2CUT_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class St2cutParser
{
  public:
    explicit St2cutParser() {}
    virtual ~St2cutParser() {}

    void extract(_data_struct &data_struct,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
      int ival;
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
#define SEXT_I(KEY, FIELD) \
      if (!data->getValue(KEY,&ival,g_id)) { \
        data->addValue(KEY, data_struct.FIELD, g_id); \
      } else { \
        data->setValue(KEY, data_struct.FIELD, g_id); \
      }
#define SEXT_D(KEY, FIELD) \
      if (!data->getValue(KEY,&rval,g_id)) { \
        data->addValue(KEY, data_struct.FIELD, g_id); \
      } else { \
        data->setValue(KEY, data_struct.FIELD, g_id); \
      }
      SEXT_I(ST2CUT_MODE,  st2cut_mode)
      SEXT_I(ST2CUT_MODE2, st2cut_mode2)
      SEXT_D(ST2CUT_K1,    st2cut_k1)
      SEXT_D(ST2CUT_K2,    st2cut_k2)
      SEXT_D(ST2CUT_T1,    st2cut_t1)
      SEXT_D(ST2CUT_T2,    st2cut_t2)
      SEXT_D(ST2CUT_T3,    st2cut_t3)
      SEXT_D(ST2CUT_T4,    st2cut_t4)
      SEXT_D(ST2CUT_T5,    st2cut_t5)
      SEXT_D(ST2CUT_T6,    st2cut_t6)
      SEXT_D(ST2CUT_T7,    st2cut_t7)
      SEXT_D(ST2CUT_T8,    st2cut_t8)
      SEXT_D(ST2CUT_T9,    st2cut_t9)
      SEXT_D(ST2CUT_T10,   st2cut_t10)
      SEXT_D(ST2CUT_LSMAX, st2cut_lsmax)
      SEXT_D(ST2CUT_LSMIN, st2cut_lsmin)
      SEXT_D(ST2CUT_VCU,   st2cut_vcu)
      SEXT_D(ST2CUT_VCL,   st2cut_vcl)
#undef SEXT_I
#undef SEXT_D
    }

    /**
     * [3]=MODE, [4]=BUSR(skip), [5]=MODE2, [6]=BUSR2(skip),
     * [7]=K1, [8]=K2, [9]=T1, [10]=T2,
     * [11]=T3, [12]=T4, [13]=T5, [14]=T6, [15]=T7, [16]=T8, [17]=T9, [18]=T10,
     * [19]=LSMAX, [20]=LSMIN, [21]=VCU, [22]=VCL
     */
    void parse(std::vector<std::string> &split_line,
        gridpack::component::DataCollection *data, int g_id)
    {
      double rval;
      bool bval;
      int ival;
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
#define SPRS_I(IDX, KEY) \
      if (nstr > IDX) { \
        if (!data->getValue(KEY,&ival,g_id)) { \
          data->addValue(KEY, atoi(split_line[IDX].c_str()), g_id); \
        } else { \
          data->setValue(KEY, atoi(split_line[IDX].c_str()), g_id); \
        } \
      }
#define SPRS_D(IDX, KEY) \
      if (nstr > IDX) { \
        if (!data->getValue(KEY,&rval,g_id)) { \
          data->addValue(KEY, atof(split_line[IDX].c_str()), g_id); \
        } else { \
          data->setValue(KEY, atof(split_line[IDX].c_str()), g_id); \
        } \
      }
      SPRS_I( 3, ST2CUT_MODE)
      // [4] = BUSR1, skip
      SPRS_I( 5, ST2CUT_MODE2)
      // [6] = BUSR2, skip
      SPRS_D( 7, ST2CUT_K1)
      SPRS_D( 8, ST2CUT_K2)
      SPRS_D( 9, ST2CUT_T1)
      SPRS_D(10, ST2CUT_T2)
      SPRS_D(11, ST2CUT_T3)
      SPRS_D(12, ST2CUT_T4)
      SPRS_D(13, ST2CUT_T5)
      SPRS_D(14, ST2CUT_T6)
      SPRS_D(15, ST2CUT_T7)
      SPRS_D(16, ST2CUT_T8)
      SPRS_D(17, ST2CUT_T9)
      SPRS_D(18, ST2CUT_T10)
      SPRS_D(19, ST2CUT_LSMAX)
      SPRS_D(20, ST2CUT_LSMIN)
      SPRS_D(21, ST2CUT_VCU)
      SPRS_D(22, ST2CUT_VCL)
#undef SPRS_I
#undef SPRS_D
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
      if (nstr >  3) data.st2cut_mode  = atoi(split_line[3].c_str());
      // [4] = BUSR1, skip
      if (nstr >  5) data.st2cut_mode2 = atoi(split_line[5].c_str());
      // [6] = BUSR2, skip
      if (nstr >  7) data.st2cut_k1    = atof(split_line[7].c_str());
      if (nstr >  8) data.st2cut_k2    = atof(split_line[8].c_str());
      if (nstr >  9) data.st2cut_t1    = atof(split_line[9].c_str());
      if (nstr > 10) data.st2cut_t2    = atof(split_line[10].c_str());
      if (nstr > 11) data.st2cut_t3    = atof(split_line[11].c_str());
      if (nstr > 12) data.st2cut_t4    = atof(split_line[12].c_str());
      if (nstr > 13) data.st2cut_t5    = atof(split_line[13].c_str());
      if (nstr > 14) data.st2cut_t6    = atof(split_line[14].c_str());
      if (nstr > 15) data.st2cut_t7    = atof(split_line[15].c_str());
      if (nstr > 16) data.st2cut_t8    = atof(split_line[16].c_str());
      if (nstr > 17) data.st2cut_t9    = atof(split_line[17].c_str());
      if (nstr > 18) data.st2cut_t10   = atof(split_line[18].c_str());
      if (nstr > 19) data.st2cut_lsmax = atof(split_line[19].c_str());
      if (nstr > 20) data.st2cut_lsmin = atof(split_line[20].c_str());
      if (nstr > 21) data.st2cut_vcu   = atof(split_line[21].c_str());
      if (nstr > 22) data.st2cut_vcl   = atof(split_line[22].c_str());
    }
};
}  // parser
}  // gridpack
#endif
