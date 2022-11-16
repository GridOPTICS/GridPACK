/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: November 7, 2022
 *      Author: Shrirang Abhyankar
 */
#ifndef HYGOV_HPP
#define HYGOV_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class HygovParser
{
  public:
    /**
     * Constructor
     */
    explicit HygovParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~HygovParser()
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
      // HAS_GOVERNOR
      if (!data->getValue(HAS_GOVERNOR,&bval,g_id)) {
        data->addValue(HAS_GOVERNOR, true, g_id);
      } else {
        data->setValue(HAS_GOVERNOR, true, g_id);
      }

      // GOVERNOR_NAME
      std::string stmp;
      if (!data->getValue(GOVERNOR_MODEL, &stmp, g_id)) {
        data->addValue(GOVERNOR_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(GOVERNOR_MODEL, data_struct.model, g_id);
      }

      // GOVERNOR_R
      if (!data->getValue(GOVERNOR_R,&rval,g_id)) {
        data->addValue(GOVERNOR_R, data_struct.gv_r, g_id);
      } else {
        data->setValue(GOVERNOR_R, data_struct.gv_r, g_id);
      }

      // GOVERNOR_r
      if (!data->getValue(GOVERNOR_r,&rval,g_id)) {
        data->addValue(GOVERNOR_r, data_struct.gv_r2, g_id);
      } else {
        data->setValue(GOVERNOR_r, data_struct.gv_r2, g_id);
      }


      // GOVERNOR_TR
      if (!data->getValue(GOVERNOR_TR,&rval,g_id)) {
        data->addValue(GOVERNOR_TR, data_struct.gv_tr, g_id);
      } else {
        data->setValue(GOVERNOR_TR, data_struct.gv_tr, g_id);
      }

      // GOVERNOR_TF
      if (!data->getValue(GOVERNOR_TF,&rval,g_id)) {
        data->addValue(GOVERNOR_TF, data_struct.gv_tf, g_id);
      } else {
        data->setValue(GOVERNOR_TF, data_struct.gv_tf, g_id);
      }

      // GOVERNOR_TG
      if (!data->getValue(GOVERNOR_TG,&rval,g_id)) {
        data->addValue(GOVERNOR_TG, data_struct.gv_tg, g_id);
      } else {
        data->setValue(GOVERNOR_TG, data_struct.gv_tg, g_id);
      }

      // GOVERNOR_VELM
      if (!data->getValue(GOVERNOR_VELM,&rval,g_id)) {
        data->addValue(GOVERNOR_VELM, data_struct.gv_velm, g_id);
      } else {
        data->setValue(GOVERNOR_VELM, data_struct.gv_velm, g_id);
      }

      // GOVERNOR_GMAX
      if (!data->getValue(GOVERNOR_GMAX,&rval,g_id)) {
        data->addValue(GOVERNOR_GMAX, data_struct.gv_gmax, g_id);
      } else {
        data->setValue(GOVERNOR_GMAX, data_struct.gv_gmax, g_id);
      }

      // GOVERNOR_GMIN
      if (!data->getValue(GOVERNOR_GMIN,&rval,g_id)) {
        data->addValue(GOVERNOR_GMIN, data_struct.gv_gmin, g_id);
      } else {
        data->setValue(GOVERNOR_GMIN, data_struct.gv_gmin, g_id);
      }

      // GOVERNOR_TW
      if (!data->getValue(GOVERNOR_TW,&rval,g_id)) {
        data->addValue(GOVERNOR_TW, data_struct.gv_tw, g_id);
      } else {
        data->setValue(GOVERNOR_TW, data_struct.gv_tw, g_id);
      }

      // GOVERNOR_AT
      if (!data->getValue(GOVERNOR_AT,&rval,g_id)) {
        data->addValue(GOVERNOR_AT, data_struct.gv_at, g_id);
      } else {
        data->setValue(GOVERNOR_AT, data_struct.gv_at, g_id);
      }

      // GOVERNOR_DT
      if (!data->getValue(GOVERNOR_DT,&rval,g_id)) {
        data->addValue(GOVERNOR_DT, data_struct.gv_dt, g_id);
      } else {
        data->setValue(GOVERNOR_DT, data_struct.gv_dt, g_id);
      }

      // GOVERNOR_QNL
      if (!data->getValue(GOVERNOR_QNL,&rval,g_id)) {
        data->addValue(GOVERNOR_QNL, data_struct.gv_qnl, g_id);
      } else {
        data->setValue(GOVERNOR_QNL, data_struct.gv_qnl, g_id);
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
      int nstr = split_line.size();
      bool bval;
      // HAS_GOVERNOR
      if (!data->getValue(HAS_GOVERNOR,&bval,g_id)) {
        data->addValue(HAS_GOVERNOR, true, g_id);
      } else {
        data->setValue(HAS_GOVERNOR, true, g_id);
      }

      // GOVERNOR_MODEL
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(GOVERNOR_MODEL,&stmp,g_id)) {
        data->addValue(GOVERNOR_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(GOVERNOR_MODEL, model.c_str(), g_id);
      }

      // GOVERNOR_R
      if (nstr > 3) {
        if (!data->getValue(GOVERNOR_R,&rval,g_id)) {
          data->addValue(GOVERNOR_R, atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_R, atof(split_line[3].c_str()), g_id);
        }
      } 

      // GOVERNOR_r
      if (nstr > 4) {
        if (!data->getValue(GOVERNOR_r,&rval,g_id)) {
          data->addValue(GOVERNOR_r, atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_r, atof(split_line[4].c_str()), g_id);
        }
      } 

      // GOVERNOR_TR
      if (nstr > 5) {
        if (!data->getValue(GOVERNOR_TR,&rval,g_id)) {
          data->addValue(GOVERNOR_TR, atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TR, atof(split_line[5].c_str()), g_id);
        }
      } 

      // GOVERNOR_TF
      if (nstr > 6) {
        if (!data->getValue(GOVERNOR_TF,&rval,g_id)) {
          data->addValue(GOVERNOR_TF, atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TF, atof(split_line[6].c_str()), g_id);
        }
      } 

      // GOVERNOR_TG
      if (nstr > 7) {
        if (!data->getValue(GOVERNOR_TG,&rval,g_id)) {
          data->addValue(GOVERNOR_TG, atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TG, atof(split_line[7].c_str()), g_id);
        }
      } 

      // GOVERNOR_VELM
      if (nstr > 8) {
        if (!data->getValue(GOVERNOR_VELM,&rval,g_id)) {
          data->addValue(GOVERNOR_VELM, atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_VELM, atof(split_line[8].c_str()), g_id);
        }
      } 

      // GOVERNOR_GMAX
      if (nstr > 9) {
        if (!data->getValue(GOVERNOR_GMAX,&rval,g_id)) {
          data->addValue(GOVERNOR_GMAX, atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_GMAX, atof(split_line[9].c_str()), g_id);
        }
      } 

      // GOVERNOR_GMIN
      if (nstr > 10) {
        if (!data->getValue(GOVERNOR_GMIN,&rval,g_id)) {
          data->addValue(GOVERNOR_GMIN, atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_GMIN, atof(split_line[10].c_str()), g_id);
        }
      } 

      // GOVERNOR_TW
      if (nstr > 11) {
        if (!data->getValue(GOVERNOR_TW,&rval,g_id)) {
          data->addValue(GOVERNOR_TW, atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TW, atof(split_line[11].c_str()), g_id);
        }
      }

      // GOVERNOR_AT
      if (nstr > 12) {
        if (!data->getValue(GOVERNOR_AT,&rval,g_id)) {
          data->addValue(GOVERNOR_AT, atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_AT, atof(split_line[12].c_str()), g_id);
        }
      }

      // GOVERNOR_DT
      if (nstr > 13) {
        if (!data->getValue(GOVERNOR_DT,&rval,g_id)) {
          data->addValue(GOVERNOR_DT, atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_DT, atof(split_line[13].c_str()), g_id);
        }
      }

      // GOVERNOR_QNL
      if (nstr > 14) {
        if (!data->getValue(GOVERNOR_QNL,&rval,g_id)) {
          data->addValue(GOVERNOR_QNL, atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_QNL, atof(split_line[14].c_str()), g_id);
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
      // GOVERNOR_BUSNUMBER               "I"                   integer
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

      // GOVERNOR_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size();
      // GOVERNOR_R
      if (nstr > 3) {
        data.gv_r = atof(split_line[3].c_str());
      }

      // GOVERNOR_r
      if (nstr > 4) {
        data.gv_r2 = atof(split_line[4].c_str());
      }

      // GOVERNOR_TR
      if (nstr > 5) {
        data.gv_tr = atof(split_line[5].c_str());
      }

      // GOVERNOR_TF
      if (nstr > 6) {
        data.gv_tf = atof(split_line[6].c_str());
      }

      // GOVERNOR_TG
      if (nstr > 7) {
        data.gv_tg = atof(split_line[7].c_str());
      }

      // GOVERNOR_VELM
      if (nstr > 8) {
        data.gv_velm = atof(split_line[8].c_str());
      }

      // GOVERNOR_GMAX
      if (nstr > 9) {
        data.gv_gmax = atof(split_line[9].c_str());
      }

      // GOVERNOR_GMIN
      if (nstr > 10) {
        data.gv_gmin = atof(split_line[10].c_str());
      }

      // GOVERNOR_TW
      if (nstr > 11) {
        data.gv_tw = atof(split_line[11].c_str());
      }

      // GOVERNOR_AT
      if (nstr > 12) {
        data.gv_at = atof(split_line[12].c_str());
      }

      // GOVERNOR_DT
      if (nstr > 13) {
        data.gv_dt = atof(split_line[13].c_str());
      }

      // GOVERNOR_QNL
      if (nstr > 14) {
        data.gv_qnl = atof(split_line[14].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
