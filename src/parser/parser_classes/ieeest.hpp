/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: March 2026
 *      Author: Yousu Chen
 *
 *  IEEEST PSS parser.
 *  DYR format: I, 'IEEEST', ID, MODE, BUSR, A1, A2, A3, A4, A5, A6,
 *                           T1, T2, T3, T4, T5, T6, KS, LSMAX, LSMIN, VCU, VCL
 *  Indices:    0    1       2   3     4    5   6   7   8   9   10
 *                          11  12  13  14  15  16  17  18     19   20   21
 */
#ifndef IEEEST_HPP
#define IEEEST_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class IeeeStParser
{
  public:
    explicit IeeeStParser() {}
    virtual ~IeeeStParser() {}

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
#define IEXT_I(KEY, FIELD) \
      if (!data->getValue(KEY,&ival,g_id)) { \
        data->addValue(KEY, data_struct.FIELD, g_id); \
      } else { \
        data->setValue(KEY, data_struct.FIELD, g_id); \
      }
#define IEXT_D(KEY, FIELD) \
      if (!data->getValue(KEY,&rval,g_id)) { \
        data->addValue(KEY, data_struct.FIELD, g_id); \
      } else { \
        data->setValue(KEY, data_struct.FIELD, g_id); \
      }
      IEXT_I(IEEEST_MODE,  ieeest_mode)
      IEXT_D(IEEEST_A1,    ieeest_a1)
      IEXT_D(IEEEST_A2,    ieeest_a2)
      IEXT_D(IEEEST_A3,    ieeest_a3)
      IEXT_D(IEEEST_A4,    ieeest_a4)
      IEXT_D(IEEEST_A5,    ieeest_a5)
      IEXT_D(IEEEST_A6,    ieeest_a6)
      IEXT_D(IEEEST_T1,    ieeest_t1)
      IEXT_D(IEEEST_T2,    ieeest_t2)
      IEXT_D(IEEEST_T3,    ieeest_t3)
      IEXT_D(IEEEST_T4,    ieeest_t4)
      IEXT_D(IEEEST_T5,    ieeest_t5)
      IEXT_D(IEEEST_T6,    ieeest_t6)
      IEXT_D(IEEEST_KS,    ieeest_ks)
      IEXT_D(IEEEST_LSMAX, ieeest_lsmax)
      IEXT_D(IEEEST_LSMIN, ieeest_lsmin)
      IEXT_D(IEEEST_VCU,   ieeest_vcu)
      IEXT_D(IEEEST_VCL,   ieeest_vcl)
#undef IEXT_I
#undef IEXT_D
    }

    /**
     * [3]=MODE, [4]=BUSR(skip), [5]=A1, [6]=A2, [7]=A3, [8]=A4, [9]=A5, [10]=A6,
     * [11]=T1, [12]=T2, [13]=T3, [14]=T4, [15]=T5, [16]=T6,
     * [17]=KS, [18]=LSMAX, [19]=LSMIN, [20]=VCU, [21]=VCL
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
#define IPRS_I(IDX, KEY) \
      if (nstr > IDX) { \
        if (!data->getValue(KEY,&ival,g_id)) { \
          data->addValue(KEY, atoi(split_line[IDX].c_str()), g_id); \
        } else { \
          data->setValue(KEY, atoi(split_line[IDX].c_str()), g_id); \
        } \
      }
#define IPRS_D(IDX, KEY) \
      if (nstr > IDX) { \
        if (!data->getValue(KEY,&rval,g_id)) { \
          data->addValue(KEY, atof(split_line[IDX].c_str()), g_id); \
        } else { \
          data->setValue(KEY, atof(split_line[IDX].c_str()), g_id); \
        } \
      }
      IPRS_I( 3, IEEEST_MODE)
      // [4] = BUSR, skip
      IPRS_D( 5, IEEEST_A1)
      IPRS_D( 6, IEEEST_A2)
      IPRS_D( 7, IEEEST_A3)
      IPRS_D( 8, IEEEST_A4)
      IPRS_D( 9, IEEEST_A5)
      IPRS_D(10, IEEEST_A6)
      IPRS_D(11, IEEEST_T1)
      IPRS_D(12, IEEEST_T2)
      IPRS_D(13, IEEEST_T3)
      IPRS_D(14, IEEEST_T4)
      IPRS_D(15, IEEEST_T5)
      IPRS_D(16, IEEEST_T6)
      IPRS_D(17, IEEEST_KS)
      IPRS_D(18, IEEEST_LSMAX)
      IPRS_D(19, IEEEST_LSMIN)
      IPRS_D(20, IEEEST_VCU)
      IPRS_D(21, IEEEST_VCL)
#undef IPRS_I
#undef IPRS_D
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
      if (nstr >  3) data.ieeest_mode  = atoi(split_line[3].c_str());
      // [4] = BUSR, skip
      if (nstr >  5) data.ieeest_a1    = atof(split_line[5].c_str());
      if (nstr >  6) data.ieeest_a2    = atof(split_line[6].c_str());
      if (nstr >  7) data.ieeest_a3    = atof(split_line[7].c_str());
      if (nstr >  8) data.ieeest_a4    = atof(split_line[8].c_str());
      if (nstr >  9) data.ieeest_a5    = atof(split_line[9].c_str());
      if (nstr > 10) data.ieeest_a6    = atof(split_line[10].c_str());
      if (nstr > 11) data.ieeest_t1    = atof(split_line[11].c_str());
      if (nstr > 12) data.ieeest_t2    = atof(split_line[12].c_str());
      if (nstr > 13) data.ieeest_t3    = atof(split_line[13].c_str());
      if (nstr > 14) data.ieeest_t4    = atof(split_line[14].c_str());
      if (nstr > 15) data.ieeest_t5    = atof(split_line[15].c_str());
      if (nstr > 16) data.ieeest_t6    = atof(split_line[16].c_str());
      if (nstr > 17) data.ieeest_ks    = atof(split_line[17].c_str());
      if (nstr > 18) data.ieeest_lsmax = atof(split_line[18].c_str());
      if (nstr > 19) data.ieeest_lsmin = atof(split_line[19].c_str());
      if (nstr > 20) data.ieeest_vcu   = atof(split_line[20].c_str());
      if (nstr > 21) data.ieeest_vcl   = atof(split_line[21].c_str());
    }
};
}  // parser
}  // gridpack
#endif
