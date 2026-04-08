/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: 2026-03-27
 *      Author: Yousu Chen
 *
 *  ESST3A exciter parser.
 *  DYR format: I, 'ESST3A', ID, TR, VIMAX, VIMIN, KM, TC, TB, KA, TA, VRMAX, VRMIN,
 *                           KG, KP, KI, VBMAX, KC, XL, VGMAX, THETAP, TM, VMMAX, VMMIN
 *  Indices:    0     1      2   3    4      5     6   7   8   9   10    11     12
 *                           13  14  15    16    17  18    19     20    21    22    23
 */
#ifndef ESST3A_HPP
#define ESST3A_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Esst3aParser
{
  public:
    explicit Esst3aParser() {}
    virtual ~Esst3aParser() {}

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
      if (!data->getValue(EXCITER_MODEL,&stmp,g_id)) {
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
      // KM
      if (!data->getValue(EXCITER_KM,&rval,g_id)) {
        data->addValue(EXCITER_KM, data_struct.ex_km, g_id);
      } else {
        data->setValue(EXCITER_KM, data_struct.ex_km, g_id);
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
      // KG
      if (!data->getValue(EXCITER_KG,&rval,g_id)) {
        data->addValue(EXCITER_KG, data_struct.ex_kg, g_id);
      } else {
        data->setValue(EXCITER_KG, data_struct.ex_kg, g_id);
      }
      // KP
      if (!data->getValue(EXCITER_KP,&rval,g_id)) {
        data->addValue(EXCITER_KP, data_struct.ex_kp, g_id);
      } else {
        data->setValue(EXCITER_KP, data_struct.ex_kp, g_id);
      }
      // KI
      if (!data->getValue(EXCITER_KI,&rval,g_id)) {
        data->addValue(EXCITER_KI, data_struct.ex_ki, g_id);
      } else {
        data->setValue(EXCITER_KI, data_struct.ex_ki, g_id);
      }
      // VBMAX
      if (!data->getValue(EXCITER_VBMAX,&rval,g_id)) {
        data->addValue(EXCITER_VBMAX, data_struct.vbmax, g_id);
      } else {
        data->setValue(EXCITER_VBMAX, data_struct.vbmax, g_id);
      }
      // KC
      if (!data->getValue(EXCITER_KC,&rval,g_id)) {
        data->addValue(EXCITER_KC, data_struct.ex_kc, g_id);
      } else {
        data->setValue(EXCITER_KC, data_struct.ex_kc, g_id);
      }
      // XL
      if (!data->getValue(EXCITER_XL,&rval,g_id)) {
        data->addValue(EXCITER_XL, data_struct.ex_xl, g_id);
      } else {
        data->setValue(EXCITER_XL, data_struct.ex_xl, g_id);
      }
      // VGMAX
      if (!data->getValue(EXCITER_VGMAX,&rval,g_id)) {
        data->addValue(EXCITER_VGMAX, data_struct.vgmax, g_id);
      } else {
        data->setValue(EXCITER_VGMAX, data_struct.vgmax, g_id);
      }
      // THETAP
      if (!data->getValue(EXCITER_THETAP,&rval,g_id)) {
        data->addValue(EXCITER_THETAP, data_struct.thetap, g_id);
      } else {
        data->setValue(EXCITER_THETAP, data_struct.thetap, g_id);
      }
      // TM
      if (!data->getValue(EXCITER_TM,&rval,g_id)) {
        data->addValue(EXCITER_TM, data_struct.ex_tm, g_id);
      } else {
        data->setValue(EXCITER_TM, data_struct.ex_tm, g_id);
      }
      // VMMAX
      if (!data->getValue(EXCITER_VMMAX,&rval,g_id)) {
        data->addValue(EXCITER_VMMAX, data_struct.vmmax, g_id);
      } else {
        data->setValue(EXCITER_VMMAX, data_struct.vmmax, g_id);
      }
      // VMMIN
      if (!data->getValue(EXCITER_VMMIN,&rval,g_id)) {
        data->addValue(EXCITER_VMMIN, data_struct.vmmin, g_id);
      } else {
        data->setValue(EXCITER_VMMIN, data_struct.vmmin, g_id);
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
        if (!data->getValue(EXCITER_TR,&rval,g_id))
          data->addValue(EXCITER_TR, atof(split_line[3].c_str()), g_id);
        else
          data->setValue(EXCITER_TR, atof(split_line[3].c_str()), g_id);
      }
      // [4] VIMAX
      if (nstr > 4) {
        if (!data->getValue(EXCITER_VIMAX,&rval,g_id))
          data->addValue(EXCITER_VIMAX, atof(split_line[4].c_str()), g_id);
        else
          data->setValue(EXCITER_VIMAX, atof(split_line[4].c_str()), g_id);
      }
      // [5] VIMIN
      if (nstr > 5) {
        if (!data->getValue(EXCITER_VIMIN,&rval,g_id))
          data->addValue(EXCITER_VIMIN, atof(split_line[5].c_str()), g_id);
        else
          data->setValue(EXCITER_VIMIN, atof(split_line[5].c_str()), g_id);
      }
      // [6] KM
      if (nstr > 6) {
        if (!data->getValue(EXCITER_KM,&rval,g_id))
          data->addValue(EXCITER_KM, atof(split_line[6].c_str()), g_id);
        else
          data->setValue(EXCITER_KM, atof(split_line[6].c_str()), g_id);
      }
      // [7] TC
      if (nstr > 7) {
        if (!data->getValue(EXCITER_TC,&rval,g_id))
          data->addValue(EXCITER_TC, atof(split_line[7].c_str()), g_id);
        else
          data->setValue(EXCITER_TC, atof(split_line[7].c_str()), g_id);
      }
      // [8] TB
      if (nstr > 8) {
        if (!data->getValue(EXCITER_TB,&rval,g_id))
          data->addValue(EXCITER_TB, atof(split_line[8].c_str()), g_id);
        else
          data->setValue(EXCITER_TB, atof(split_line[8].c_str()), g_id);
      }
      // [9] KA
      if (nstr > 9) {
        if (!data->getValue(EXCITER_KA,&rval,g_id))
          data->addValue(EXCITER_KA, atof(split_line[9].c_str()), g_id);
        else
          data->setValue(EXCITER_KA, atof(split_line[9].c_str()), g_id);
      }
      // [10] TA
      if (nstr > 10) {
        if (!data->getValue(EXCITER_TA,&rval,g_id))
          data->addValue(EXCITER_TA, atof(split_line[10].c_str()), g_id);
        else
          data->setValue(EXCITER_TA, atof(split_line[10].c_str()), g_id);
      }
      // [11] VRMAX
      if (nstr > 11) {
        if (!data->getValue(EXCITER_VRMAX,&rval,g_id))
          data->addValue(EXCITER_VRMAX, atof(split_line[11].c_str()), g_id);
        else
          data->setValue(EXCITER_VRMAX, atof(split_line[11].c_str()), g_id);
      }
      // [12] VRMIN
      if (nstr > 12) {
        if (!data->getValue(EXCITER_VRMIN,&rval,g_id))
          data->addValue(EXCITER_VRMIN, atof(split_line[12].c_str()), g_id);
        else
          data->setValue(EXCITER_VRMIN, atof(split_line[12].c_str()), g_id);
      }
      // [13] KG
      if (nstr > 13) {
        if (!data->getValue(EXCITER_KG,&rval,g_id))
          data->addValue(EXCITER_KG, atof(split_line[13].c_str()), g_id);
        else
          data->setValue(EXCITER_KG, atof(split_line[13].c_str()), g_id);
      }
      // [14] KP
      if (nstr > 14) {
        if (!data->getValue(EXCITER_KP,&rval,g_id))
          data->addValue(EXCITER_KP, atof(split_line[14].c_str()), g_id);
        else
          data->setValue(EXCITER_KP, atof(split_line[14].c_str()), g_id);
      }
      // [15] KI
      if (nstr > 15) {
        if (!data->getValue(EXCITER_KI,&rval,g_id))
          data->addValue(EXCITER_KI, atof(split_line[15].c_str()), g_id);
        else
          data->setValue(EXCITER_KI, atof(split_line[15].c_str()), g_id);
      }
      // [16] VBMAX
      if (nstr > 16) {
        if (!data->getValue(EXCITER_VBMAX,&rval,g_id))
          data->addValue(EXCITER_VBMAX, atof(split_line[16].c_str()), g_id);
        else
          data->setValue(EXCITER_VBMAX, atof(split_line[16].c_str()), g_id);
      }
      // [17] KC
      if (nstr > 17) {
        if (!data->getValue(EXCITER_KC,&rval,g_id))
          data->addValue(EXCITER_KC, atof(split_line[17].c_str()), g_id);
        else
          data->setValue(EXCITER_KC, atof(split_line[17].c_str()), g_id);
      }
      // [18] XL
      if (nstr > 18) {
        if (!data->getValue(EXCITER_XL,&rval,g_id))
          data->addValue(EXCITER_XL, atof(split_line[18].c_str()), g_id);
        else
          data->setValue(EXCITER_XL, atof(split_line[18].c_str()), g_id);
      }
      // [19] VGMAX
      if (nstr > 19) {
        if (!data->getValue(EXCITER_VGMAX,&rval,g_id))
          data->addValue(EXCITER_VGMAX, atof(split_line[19].c_str()), g_id);
        else
          data->setValue(EXCITER_VGMAX, atof(split_line[19].c_str()), g_id);
      }
      // [20] THETAP
      if (nstr > 20) {
        if (!data->getValue(EXCITER_THETAP,&rval,g_id))
          data->addValue(EXCITER_THETAP, atof(split_line[20].c_str()), g_id);
        else
          data->setValue(EXCITER_THETAP, atof(split_line[20].c_str()), g_id);
      }
      // [21] TM
      if (nstr > 21) {
        if (!data->getValue(EXCITER_TM,&rval,g_id))
          data->addValue(EXCITER_TM, atof(split_line[21].c_str()), g_id);
        else
          data->setValue(EXCITER_TM, atof(split_line[21].c_str()), g_id);
      }
      // [22] VMMAX
      if (nstr > 22) {
        if (!data->getValue(EXCITER_VMMAX,&rval,g_id))
          data->addValue(EXCITER_VMMAX, atof(split_line[22].c_str()), g_id);
        else
          data->setValue(EXCITER_VMMAX, atof(split_line[22].c_str()), g_id);
      }
      // [23] VMMIN
      if (nstr > 23) {
        if (!data->getValue(EXCITER_VMMIN,&rval,g_id))
          data->addValue(EXCITER_VMMIN, atof(split_line[23].c_str()), g_id);
        else
          data->setValue(EXCITER_VMMIN, atof(split_line[23].c_str()), g_id);
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
      if (nstr > 3)  data.ex_tr   = atof(split_line[3].c_str());
      if (nstr > 4)  data.vimax   = atof(split_line[4].c_str());
      if (nstr > 5)  data.vimin   = atof(split_line[5].c_str());
      if (nstr > 6)  data.ex_km   = atof(split_line[6].c_str());
      if (nstr > 7)  data.ex_tc   = atof(split_line[7].c_str());
      if (nstr > 8)  data.ex_tb   = atof(split_line[8].c_str());
      if (nstr > 9)  data.ex_ka   = atof(split_line[9].c_str());
      if (nstr > 10) data.ex_ta   = atof(split_line[10].c_str());
      if (nstr > 11) data.vrmax   = atof(split_line[11].c_str());
      if (nstr > 12) data.vrmin   = atof(split_line[12].c_str());
      if (nstr > 13) data.ex_kg   = atof(split_line[13].c_str());
      if (nstr > 14) data.ex_kp   = atof(split_line[14].c_str());
      if (nstr > 15) data.ex_ki   = atof(split_line[15].c_str());
      if (nstr > 16) data.vbmax   = atof(split_line[16].c_str());
      if (nstr > 17) data.ex_kc   = atof(split_line[17].c_str());
      if (nstr > 18) data.ex_xl   = atof(split_line[18].c_str());
      if (nstr > 19) data.vgmax   = atof(split_line[19].c_str());
      if (nstr > 20) data.thetap  = atof(split_line[20].c_str());
      if (nstr > 21) data.ex_tm   = atof(split_line[21].c_str());
      if (nstr > 22) data.vmmax   = atof(split_line[22].c_str());
      if (nstr > 23) data.vmmin   = atof(split_line[23].c_str());
    }
};
}  // parser
}  // gridpack
#endif
