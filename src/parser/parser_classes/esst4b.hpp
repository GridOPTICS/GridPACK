/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: June 20, 2016
 *      Author: Bruce Palmer
 */
#ifndef ESST4B_HPP
#define ESST4B_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Esst4bParser
{
  public:
    /**
     * Constructor
     */
    explicit Esst4bParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Esst4bParser()
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

      // EXCITER_TR
      if (!data->getValue(EXCITER_TR,&rval,g_id)) {
        data->addValue(EXCITER_TR, data_struct.ex_tr, g_id);
      } else {
        data->setValue(EXCITER_TR, data_struct.ex_tr, g_id);
      }

      // EXCITER_KPR
      if (!data->getValue(EXCITER_KPR,&rval,g_id)) {
        data->addValue(EXCITER_KPR, data_struct.kpr, g_id);
      } else {
        data->setValue(EXCITER_KPR, data_struct.kpr, g_id);
      }

      // EXCITER_KIR
      if (!data->getValue(EXCITER_KIR,&rval,g_id)) {
        data->addValue(EXCITER_KIR, data_struct.kir, g_id);
      } else {
        data->setValue(EXCITER_KIR, data_struct.kir, g_id);
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

      // EXCITER_TA
      if (!data->getValue(EXCITER_TA,&rval,g_id)) {
        data->addValue(EXCITER_TA, data_struct.ex_ta, g_id);
      } else {
        data->setValue(EXCITER_TA, data_struct.ex_ta, g_id);
      }

      // EXCITER_KPM
      if (!data->getValue(EXCITER_KPM,&rval,g_id)) {
        data->addValue(EXCITER_KPM, data_struct.kpm, g_id);
      } else {
        data->setValue(EXCITER_KPM, data_struct.kpm, g_id);
      }

      // EXCITER_KIM
      if (!data->getValue(EXCITER_KIM,&rval,g_id)) {
        data->addValue(EXCITER_KIM, data_struct.kim, g_id);
      } else {
        data->setValue(EXCITER_KIM, data_struct.kim, g_id);
      }

      // EXCITER_VMMAX
      if (!data->getValue(EXCITER_VMMAX,&rval,g_id)) {
        data->addValue(EXCITER_VMMAX, data_struct.vmmax, g_id);
      } else {
        data->setValue(EXCITER_VMMAX, data_struct.vmmax, g_id);
      }

      // EXCITER_VMMIN
      if (!data->getValue(EXCITER_VMMIN,&rval,g_id)) {
        data->addValue(EXCITER_VMMIN, data_struct.vmmin, g_id);
      } else {
        data->setValue(EXCITER_VMMIN, data_struct.vmmin, g_id);
      }

      // EXCITER_KG
      if (!data->getValue(EXCITER_KG,&rval,g_id)) {
        data->addValue(EXCITER_KG, data_struct.ex_kg, g_id);
      } else {
        data->setValue(EXCITER_KG, data_struct.ex_kg, g_id);
      }

      // EXCITER_KP
      if (!data->getValue(EXCITER_KP,&rval,g_id)) {
        data->addValue(EXCITER_KP, data_struct.ex_kp, g_id);
      } else {
        data->setValue(EXCITER_KP, data_struct.ex_kp, g_id);
      }

      // EXCITER_KI
      if (!data->getValue(EXCITER_KI,&rval,g_id)) {
        data->addValue(EXCITER_KI, data_struct.ex_ki, g_id);
      } else {
        data->setValue(EXCITER_KI, data_struct.ex_ki, g_id);
      }

      // EXCITER_VBMAX
      if (!data->getValue(EXCITER_VBMAX,&rval,g_id)) {
        data->addValue(EXCITER_VBMAX, data_struct.vbmax, g_id);
      } else {
        data->setValue(EXCITER_VBMAX, data_struct.vbmax, g_id);
      }

      // EXCITER_KC
      if (!data->getValue(EXCITER_KC,&rval,g_id)) {
        data->addValue(EXCITER_KC, data_struct.ex_kc, g_id);
      } else {
        data->setValue(EXCITER_KC, data_struct.ex_kc, g_id);
      }

      // EXCITER_XL
      if (!data->getValue(EXCITER_XL,&rval,g_id)) {
        data->addValue(EXCITER_XL, data_struct.ex_xl, g_id);
      } else {
        data->setValue(EXCITER_XL, data_struct.ex_xl, g_id);
      }

      // EXCITER_THETAP
      if (!data->getValue(EXCITER_THETAP,&rval,g_id)) {
        data->addValue(EXCITER_THETAP, data_struct.thetap, g_id);
      } else {
        data->setValue(EXCITER_THETAP, data_struct.thetap, g_id);
      }
    }

    /**
     * Parser list of strings and store results in data collection object
     * @param split_line list of tokens from .dyr file
     * @param data data collection object
     * @param gen_id index of generator
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

      // EXCITER_TR
      if (nstr > 3) {
        if (!data->getValue(EXCITER_TR,&rval,g_id)) {
          data->addValue(EXCITER_TR,
              atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TR,
              atof(split_line[3].c_str()), g_id);
        }
      } 

      // EXCITER_KPR
      if (nstr > 4) {
        if (!data->getValue(EXCITER_KPR,&rval,g_id)) {
          data->addValue(EXCITER_KPR,
              atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KPR,
              atof(split_line[4].c_str()), g_id);
        }
      } 

      // EXCITER_KIR
      if (nstr > 5) {
        if (!data->getValue(EXCITER_KIR,&rval,g_id)) {
          data->addValue(EXCITER_KIR,
              atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KIR,
              atof(split_line[5].c_str()), g_id);
        }
      }

      // EXCITER_VRMAX
      if (nstr > 6) {
        if (!data->getValue(EXCITER_VRMAX,&rval,g_id)) {
          data->addValue(EXCITER_VRMAX,
              atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VRMAX,
              atof(split_line[6].c_str()), g_id);
        }
      }

      // EXCITER_VRMIN
      if (nstr > 7) {
        if (!data->getValue(EXCITER_VRMIN,&rval,g_id)) {
          data->addValue(EXCITER_VRMIN,
              atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VRMIN,
              atof(split_line[7].c_str()), g_id);
        }
      }

      // EXCITER_TA
      if (nstr > 8) {
        if (!data->getValue(EXCITER_TA,&rval,g_id)) {
          data->addValue(EXCITER_TA,
              atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(EXCITER_TA,
              atof(split_line[8].c_str()), g_id);
        }
      } 

      // EXCITER_KPM
      if (nstr > 9) {
        if (!data->getValue(EXCITER_KPM,&rval,g_id)) {
          data->addValue(EXCITER_KPM,
              atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KPM,
              atof(split_line[9].c_str()), g_id);
        }
      } 

      // EXCITER_KIM
      if (nstr > 10) {
        if (!data->getValue(EXCITER_KIM,&rval,g_id)) {
          data->addValue(EXCITER_KIM,
              atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KIM,
              atof(split_line[10].c_str()), g_id);
        }
      } 

      // EXCITER_VMMAX
      if (nstr > 11) {
        if (!data->getValue(EXCITER_VMMAX,&rval,g_id)) {
          data->addValue(EXCITER_VMMAX,
              atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VMMAX,
              atof(split_line[11].c_str()), g_id);
        }
      } 

      // EXCITER_VMMIN
      if (nstr > 12) {
        if (!data->getValue(EXCITER_VMMIN,&rval,g_id)) {
          data->addValue(EXCITER_VMMIN,
              atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VMMIN,
              atof(split_line[12].c_str()), g_id);
        }
      } 

      // EXCITER_KG
      if (nstr > 13) {
        if (!data->getValue(EXCITER_KG,&rval,g_id)) {
          data->addValue(EXCITER_KG,
              atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KG,
              atof(split_line[13].c_str()), g_id);
        }
      } 

      // EXCITER_KP
      if (nstr > 14) {
        if (!data->getValue(EXCITER_KP,&rval,g_id)) {
          data->addValue(EXCITER_KP,
              atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KP,
              atof(split_line[14].c_str()), g_id);
        }
      } 

      // EXCITER_KI
      if (nstr > 15) {
        if (!data->getValue(EXCITER_KI,&rval,g_id)) {
          data->addValue(EXCITER_KI,
              atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KI,
              atof(split_line[15].c_str()), g_id);
        }
      } 

      // EXCITER_VBMAX
      if (nstr > 16) {
        if (!data->getValue(EXCITER_VBMAX,&rval,g_id)) {
          data->addValue(EXCITER_VBMAX,
              atof(split_line[16].c_str()), g_id);
        } else {
          data->setValue(EXCITER_VBMAX,
              atof(split_line[16].c_str()), g_id);
        }
      } 

      // EXCITER_KC
      if (nstr > 17) {
        if (!data->getValue(EXCITER_KC,&rval,g_id)) {
          data->addValue(EXCITER_KC,
              atof(split_line[17].c_str()), g_id);
        } else {
          data->setValue(EXCITER_KC,
              atof(split_line[17].c_str()), g_id);
        }
      } 

      // EXCITER_XL
      if (nstr > 18) {
        if (!data->getValue(EXCITER_XL,&rval,g_id)) {
          data->addValue(EXCITER_XL,
              atof(split_line[18].c_str()), g_id);
        } else {
          data->setValue(EXCITER_XL,
              atof(split_line[18].c_str()), g_id);
        }
      } 

      // EXCITER_THETAP
      if (nstr > 19) {
        if (!data->getValue(EXCITER_THETAP,&rval,g_id)) {
          data->addValue(EXCITER_THETAP,
              atof(split_line[19].c_str()), g_id);
        } else {
          data->setValue(EXCITER_THETAP,
              atof(split_line[19].c_str()), g_id);
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
      // EXCITER_BUSNUMBER               "I"                   integer
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
      // EXCITER_TR
      if (nstr > 3) {
        data.ex_tr = atof(split_line[3].c_str());
      }

      // EXCITER_KPR
      if (nstr > 4) {
        data.kpr = atof(split_line[4].c_str());
      }

      // EXCITER_KIR
      if (nstr > 5) {
        data.kir = atof(split_line[5].c_str());
      }

      // EXCITER_VRMAX
      if (nstr > 6) {
        data.vrmax = atof(split_line[6].c_str());
      }

      // EXCITER_VRMIN
      if (nstr > 7) {
        data.vrmin = atof(split_line[7].c_str());
      }

      // EXCITER_TA
      if (nstr > 8) {
        data.ex_ta = atof(split_line[8].c_str());
      }

      // EXCITER_KPM
      if (nstr > 9) {
        data.kpm = atof(split_line[9].c_str());
      }

      // EXCITER_KIM
      if (nstr > 10) {
        data.kim = atof(split_line[10].c_str());
      }

      // EXCITER_VMMAX
      if (nstr > 11) {
        data.vmmax = atof(split_line[11].c_str());
      }

      // EXCITER_VMMIN
      if (nstr > 12) {
        data.vmmin = atof(split_line[12].c_str());
      }

      // EXCITER_KG
      if (nstr > 13) {
        data.ex_kg = atof(split_line[13].c_str());
      }

      // EXCITER_KP
      if (nstr > 14) {
        data.ex_kp = atof(split_line[14].c_str());
      }

      // EXCITER_KI
      if (nstr > 15) {
        data.ex_ki = atof(split_line[15].c_str());
      }

      // EXCITER_VBMAX
      if (nstr > 16) {
        data.vbmax = atof(split_line[16].c_str());
      }

      // EXCITER_KC
      if (nstr > 17) {
        data.ex_kc = atof(split_line[17].c_str());
      }

      // EXCITER_XL
      if (nstr > 18) {
        data.ex_xl = atof(split_line[18].c_str());
      }

      // EXCITER_THETAP
      if (nstr > 19) {
        data.thetap = atof(split_line[19].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
