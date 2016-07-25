/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: June 20, 2016
 *      Author: Bruce Palmer
 */
#ifndef ESST1A_HPP
#define ESST1A_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Esst1aParser
{
  public:
    /**
     * Constructor
     */
    explicit Esst1aParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Esst1aParser()
    {
    }

    /**
     * Extract data from _data_struct and store it in data collection object
     * @param data_struct data struct object
     * @param data data collection object
     * @param gen_id index of generator
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

      // EXCITER_UEL
      if (!data->getValue(EXCITER_UEL,&rval,g_id)) {
        data->addValue(EXCITER_UEL, data_struct.uel, g_id);
      } else {
        data->setValue(EXCITER_UEL, data_struct.uel, g_id);
      }

      // EXCITER_VOS
      if (!data->getValue(EXCITER_VOS,&rval,g_id)) {
        data->addValue(EXCITER_VOS, data_struct.vos, g_id);
      } else {
        data->setValue(EXCITER_VOS, data_struct.vos, g_id);
      }

      // EXCITER_TR
      if (!data->getValue(EXCITER_TR,&rval,g_id)) {
        data->addValue(EXCITER_TR, data_struct.ex_tr, g_id);
      } else {
        data->setValue(EXCITER_TR, data_struct.ex_tr, g_id);
      }

      // EXCITER_VIMAX
      if (!data->getValue(EXCITER_VIMAX,&rval,g_id)) {
        data->addValue(EXCITER_VIMAX, data_struct.vimax, g_id);
      } else {
        data->setValue(EXCITER_VIMAX, data_struct.vimax, g_id);
      }

      // EXCITER_VIMIN
      if (!data->getValue(EXCITER_VIMIN,&rval,g_id)) {
        data->addValue(EXCITER_VIMIN, data_struct.vimin, g_id);
      } else {
        data->setValue(EXCITER_VIMIN, data_struct.vimin, g_id);
      }

      // EXCITER_TC
      if (!data->getValue(EXCITER_TC,&rval,g_id)) {
        data->addValue(EXCITER_TC, data_struct.ex_tc, g_id);
      } else {
        data->setValue(EXCITER_TC, data_struct.ex_tc, g_id);
      }

      // EXCITER_TB
      if (!data->getValue(EXCITER_TB,&rval,g_id)) {
        data->addValue(EXCITER_TB, data_struct.ex_tb, g_id);
      } else {
        data->setValue(EXCITER_TB, data_struct.ex_tb, g_id);
      }

      // EXCITER_TC1
      if (!data->getValue(EXCITER_TC1,&rval,g_id)) {
        data->addValue(EXCITER_TC1, data_struct.tc1, g_id);
      } else {
        data->setValue(EXCITER_TC1, data_struct.tc1, g_id);
      }

      // EXCITER_TB1
      if (!data->getValue(EXCITER_TB1,&rval,g_id)) {
        data->addValue(EXCITER_TB1, data_struct.tb1, g_id);
      } else {
        data->setValue(EXCITER_TB1, data_struct.tb1, g_id);
      }

      // EXCITER_KA
      if (!data->getValue(EXCITER_KA,&rval,g_id)) {
        data->addValue(EXCITER_KA, data_struct.ex_ka, g_id);
      } else {
        data->setValue(EXCITER_KA, data_struct.ex_ka, g_id);
      }

      // EXCITER_TA
      if (!data->getValue(EXCITER_TA,&rval,g_id)) {
        data->addValue(EXCITER_TA, data_struct.ex_ta, g_id);
      } else {
        data->setValue(EXCITER_TA, data_struct.ex_ta, g_id);
      }

      // EXCITER_VAMAX
      if (!data->getValue(EXCITER_VAMAX,&rval,g_id)) {
        data->addValue(EXCITER_VAMAX, data_struct.vamax, g_id);
      } else {
        data->setValue(EXCITER_VAMAX, data_struct.vamax, g_id);
      }

      // EXCITER_VAMIN
      if (!data->getValue(EXCITER_VAMIN,&rval,g_id)) {
        data->addValue(EXCITER_VAMIN, data_struct.vamin, g_id);
      } else {
        data->setValue(EXCITER_VAMIN, data_struct.vamin, g_id);
      }

      // EXCITER_VRMAX
      if (!data->getValue(EXCITER_VRMAX,&rval,g_id)) {
        data->addValue(EXCITER_VRMAX, data_struct.vrmax, g_id);
      } else {
        data->setValue(EXCITER_VRMAX, data_struct.vrmax, g_id);
      }

      // EXCITER_VRMIN
      if (!data->getValue(EXCITER_VRMIN,&rval,g_id)) {
        data->addValue(EXCITER_VRMIN, data_struct.vrmin, g_id);
      } else {
        data->setValue(EXCITER_VRMIN, data_struct.vrmin, g_id);
      }

      // EXCITER_KC
      if (!data->getValue(EXCITER_KC,&rval,g_id)) {
        data->addValue(EXCITER_KC, data_struct.ex_kc, g_id);
      } else {
        data->setValue(EXCITER_KC, data_struct.ex_kc, g_id);
      }

      // EXCITER_KF
      if (!data->getValue(EXCITER_KF,&rval,g_id)) {
        data->addValue(EXCITER_KF, data_struct.ex_kf, g_id);
      } else {
        data->setValue(EXCITER_KF, data_struct.ex_kf, g_id);
      }

      // EXCITER_KLR
      if (!data->getValue(EXCITER_KLR,&rval,g_id)) {
        data->addValue(EXCITER_KLR, data_struct.klr, g_id);
      } else {
        data->setValue(EXCITER_KLR, data_struct.klr, g_id);
      }

      // EXCITER_ILR
      if (!data->getValue(EXCITER_ILR,&rval,g_id)) {
        data->addValue(EXCITER_ILR, data_struct.ilr, g_id);
      } else {
        data->setValue(EXCITER_ILR, data_struct.ilr, g_id);
      }

      // EXCITER_TF
      if (!data->getValue(EXCITER_TF,&rval,g_id)) {
        data->addValue(EXCITER_TF, data_struct.ex_tf, g_id);
      } else {
        data->setValue(EXCITER_TF, data_struct.ex_tf, g_id);
      }
    }

    /**
     * Parser list of strings and store results in data collection object
     * @param split_line list of tokens from .dyr file
     * @param data data collection object
     * @param gen_id index of generator
     */
    void parse(std::vector<std::string> &split_line,
        gridpack::component::DataCollection *data,
        int g_id)
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

      // EXCITER_UEL
      if (nstr > 3) {
        if (!data->getValue(EXCITER_UEL,&rval,g_id)) {
          data->addValue(EXCITER_UEL,
              atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(EXCITER_UEL,
              atof(split_line[3].c_str()), g_id);
        }
      } 

      // EXCITER_VOS
      if (nstr > 4) {
        if (!data->getValue(EXCITER_VOS,&rval,g_id)) {
          data->addValue(EXCITER_VOS,
              atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VOS,
              atof(split_line[4].c_str()), g_id);
        }
      } 

      // EXCITER_TR
      if (nstr > 5) {
        if (!data->getValue(EXCITER_TR,&rval,g_id)) {
          data->addValue(EXCITER_TR,
              atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TR,
              atof(split_line[5].c_str()), g_id);
        }
      }

      // EXCITER_VIMAX
      if (nstr > 6) {
        if (!data->getValue(EXCITER_VIMAX,&rval,g_id)) {
          data->addValue(EXCITER_VIMAX,
              atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VIMAX,
              atof(split_line[6].c_str()), g_id);
        }
      }

      // EXCITER_VIMIN
      if (nstr > 7) {
        if (!data->getValue(EXCITER_VIMIN,&rval,g_id)) {
          data->addValue(EXCITER_VIMIN,
              atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VIMIN,
              atof(split_line[7].c_str()), g_id);
        }
      }

      // EXCITER_TC
      if (nstr > 8) {
        if (!data->getValue(EXCITER_TC,&rval,g_id)) {
          data->addValue(EXCITER_TC,
              atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TC,
              atof(split_line[8].c_str()), g_id);
        }
      } 

      // EXCITER_TB
      if (nstr > 9) {
        if (!data->getValue(EXCITER_TB,&rval,g_id)) {
          data->addValue(EXCITER_TB,
              atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TB,
              atof(split_line[9].c_str()), g_id);
        }
      } 

      // EXCITER_TC1
      if (nstr > 10) {
        if (!data->getValue(EXCITER_TC1,&rval,g_id)) {
          data->addValue(EXCITER_TC1,
              atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TC1,
              atof(split_line[10].c_str()), g_id);
        }
      } 

      // EXCITER_TB1
      if (nstr > 11) {
        if (!data->getValue(EXCITER_TB1,&rval,g_id)) {
          data->addValue(EXCITER_TB1,
              atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TB1,
              atof(split_line[11].c_str()), g_id);
        }
      } 

      // EXCITER_KA
      if (nstr > 12) {
        if (!data->getValue(EXCITER_KA,&rval,g_id)) {
          data->addValue(EXCITER_KA,
              atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KA,
              atof(split_line[12].c_str()), g_id);
        }
      } 

      // EXCITER_TA
      if (nstr > 13) {
        if (!data->getValue(EXCITER_TA,&rval,g_id)) {
          data->addValue(EXCITER_TA,
              atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TA,
              atof(split_line[13].c_str()), g_id);
        }
      } 

      // EXCITER_VAMAX
      if (nstr > 14) {
        if (!data->getValue(EXCITER_VAMAX,&rval,g_id)) {
          data->addValue(EXCITER_VAMAX,
              atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VAMAX,
              atof(split_line[14].c_str()), g_id);
        }
      } 

      // EXCITER_VAMIN
      if (nstr > 15) {
        if (!data->getValue(EXCITER_VAMIN,&rval,g_id)) {
          data->addValue(EXCITER_VAMIN,
              atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VAMIN,
              atof(split_line[15].c_str()), g_id);
        }
      } 

      // EXCITER_VRMAX
      if (nstr > 16) {
        if (!data->getValue(EXCITER_VRMAX,&rval,g_id)) {
          data->addValue(EXCITER_VRMAX,
              atof(split_line[16].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VRMAX,
              atof(split_line[16].c_str()), g_id);
        }
      } 

      // EXCITER_VRMIN
      if (nstr > 17) {
        if (!data->getValue(EXCITER_VRMAX,&rval,g_id)) {
          data->addValue(EXCITER_VRMAX,
              atof(split_line[17].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VRMAX,
              atof(split_line[17].c_str()), g_id);
        }
      } 

      // EXCITER_KC
      if (nstr > 18) {
        if (!data->getValue(EXCITER_KC,&rval,g_id)) {
          data->addValue(EXCITER_KC,
              atof(split_line[18].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KC,
              atof(split_line[18].c_str()), g_id);
        }
      } 

      // EXCITER_KF
      if (nstr > 19) {
        if (!data->getValue(EXCITER_KF,&rval,g_id)) {
          data->addValue(EXCITER_KF,
              atof(split_line[19].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KF,
              atof(split_line[19].c_str()), g_id);
        }
      } 

      // EXCITER_TF
      if (nstr > 20) {
        if (!data->getValue(EXCITER_TF,&rval,g_id)) {
          data->addValue(EXCITER_TF,
              atof(split_line[20].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TF,
              atof(split_line[20].c_str()), g_id);
        }
      } 

      // EXCITER_KLR
      if (nstr > 21) {
        if (!data->getValue(EXCITER_KLR,&rval,g_id)) {
          data->addValue(EXCITER_KLR,
              atof(split_line[21].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KLR,
              atof(split_line[21].c_str()), g_id);
        }
      } 

      // EXCITER_ILR
      if (nstr > 22) {
        if (!data->getValue(EXCITER_ILR,&rval,g_id)) {
          data->addValue(EXCITER_ILR,
              atof(split_line[22].c_str()), g_id);
        } else {
          data->setValue(EXCITER_ILR,
              atof(split_line[22].c_str()), g_id);
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
      // EXCITOR_BUSNUMBER               "I"                   integer
      int o_idx;
      o_idx = atoi(split_line[0].c_str());
      data.bus_id = o_idx;

      // Clean up 2 character tag for generator ID
      gridpack::utility::StringUtils util;
      std::string tag = util.clean2Char(split_line[2]);
      strcpy(data.gen_id, tag.c_str());

      std::string sval;

      sval = util.trimQuotes(split_line[1]);
      util.toUpper(sval);

      // EXCITER_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());
      int nstr = split_line.size();

      // EXCITER_UEL
      if (nstr > 3) {
        data.uel = atof(split_line[3].c_str());
      }

      // EXCITER_VOS
      if (nstr > 4) {
        data.vos = atof(split_line[4].c_str());
      }

      // EXCITER_TR
      if (nstr > 5) {
        data.ex_tr = atof(split_line[5].c_str());
      }

      // EXCITER_VIMAX
      if (nstr > 6) {
        data.vimax = atof(split_line[6].c_str());
      }

      // EXCITER_VIMIN
      if (nstr > 7) {
        data.vimin = atof(split_line[7].c_str());
      }

      // EXCITER_TC
      if (nstr > 8) {
        data.ex_tc = atof(split_line[8].c_str());
      }

      // EXCITER_TB
      if (nstr > 9) {
        data.ex_tb = atof(split_line[9].c_str());
      }

      // EXCITER_TC1
      if (nstr > 10) {
        data.tc1 = atof(split_line[10].c_str());
      }

      // EXCITER_TB1
      if (nstr > 11) {
        data.tb1 = atof(split_line[11].c_str());
      }

      // EXCITER_KA
      if (nstr > 12) {
        data.ex_ka = atof(split_line[12].c_str());
      }

      // EXCITER_TA
      if (nstr > 13) {
        data.ex_ta = atof(split_line[13].c_str());
      }

      // EXCITER_VAMAX
      if (nstr > 14) {
        data.vamax = atof(split_line[14].c_str());
      }

      // EXCITER_VAMIN
      if (nstr > 15) {
        data.vamin = atof(split_line[15].c_str());
      }

      // EXCITER_VRMAX
      if (nstr > 16) {
        data.vrmax = atof(split_line[16].c_str());
      }

      // EXCITER_VRMIN
      if (nstr > 17) {
        data.vrmin = atof(split_line[17].c_str());
      }

      // EXCITER_KC
      if (nstr > 18) {
        data.ex_kc = atof(split_line[18].c_str());
      }

      // EXCITER_KF
      if (nstr > 19) {
        data.ex_kf = atof(split_line[19].c_str());
      }

      // EXCITER_TF
      if (nstr > 20) {
        data.ex_tf = atof(split_line[20].c_str());
      }

      // EXCITER_KLR
      if (nstr > 21) {
        data.klr = atof(split_line[21].c_str());
      }

      // EXCITER_ILR
      if (nstr > 22) {
        data.ilr = atof(split_line[22].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
