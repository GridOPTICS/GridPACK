/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: June 16, 2016
 *      Author: Bruce Palmer
 */
#ifndef WSHYGP_HPP
#define WSHYGP_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class WshygpParser
{
  public:
    /**
     * Constructor
     */
    explicit WshygpParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~WshygpParser()
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

      // GOVERNOR_DB1
      if (!data->getValue(GOVERNOR_DB1,&rval,g_id)) {
        data->addValue(GOVERNOR_DB1, data_struct.db1, g_id);
      } else {
        data->setValue(GOVERNOR_DB1, data_struct.db1, g_id);
      }

      // GOVERNOR_ERR
      if (!data->getValue(GOVERNOR_ERR,&rval,g_id)) {
        data->addValue(GOVERNOR_ERR, data_struct.err, g_id);
      } else {
        data->setValue(GOVERNOR_ERR, data_struct.err, g_id);
      }

      // GOVERNOR_TD
      if (!data->getValue(GOVERNOR_TD,&rval,g_id)) {
        data->addValue(GOVERNOR_TD, data_struct.gv_td, g_id);
      } else {
        data->setValue(GOVERNOR_TD, data_struct.gv_td, g_id);
      }

      // GOVERNOR_KI
      if (!data->getValue(GOVERNOR_KI,&rval,g_id)) {
        data->addValue(GOVERNOR_KI, data_struct.gv_ki, g_id);
      } else {
        data->setValue(GOVERNOR_KI, data_struct.gv_ki, g_id);
      }

      // GOVERNOR_TF
      if (!data->getValue(GOVERNOR_TF,&rval,g_id)) {
        data->addValue(GOVERNOR_TF, data_struct.gv_tf, g_id);
      } else {
        data->setValue(GOVERNOR_TF, data_struct.gv_tf, g_id);
      }

      // GOVERNOR_KD
      if (!data->getValue(GOVERNOR_KD,&rval,g_id)) {
        data->addValue(GOVERNOR_KD, data_struct.gv_kd, g_id);
      } else {
        data->setValue(GOVERNOR_KD, data_struct.gv_kd, g_id);
      }

      // GOVERNOR_KP
      if (!data->getValue(GOVERNOR_KP,&rval,g_id)) {
        data->addValue(GOVERNOR_KP, data_struct.gv_kp, g_id);
      } else {
        data->setValue(GOVERNOR_KP, data_struct.gv_kp, g_id);
      }

      // GOVERNOR_R
      if (!data->getValue(GOVERNOR_R,&rval,g_id)) {
        data->addValue(GOVERNOR_R, data_struct.gv_r, g_id);
      } else {
        data->setValue(GOVERNOR_R, data_struct.gv_r, g_id);
      }

      // GOVERNOR_TT
      if (!data->getValue(GOVERNOR_TT,&rval,g_id)) {
        data->addValue(GOVERNOR_TT, data_struct.gv_tt, g_id);
      } else {
        data->setValue(GOVERNOR_TT, data_struct.gv_tt, g_id);
      }

      // GOVERNOR_KG
      if (!data->getValue(GOVERNOR_KG,&rval,g_id)) {
        data->addValue(GOVERNOR_KG, data_struct.gv_kg, g_id);
      } else {
        data->setValue(GOVERNOR_KG, data_struct.gv_kg, g_id);
      }

      // GOVERNOR_TP
      if (!data->getValue(GOVERNOR_TP,&rval,g_id)) {
        data->addValue(GOVERNOR_TP, data_struct.gv_tp, g_id);
      } else {
        data->setValue(GOVERNOR_TP, data_struct.gv_tp, g_id);
      }

      // GOVERNOR_VELOPEN
      if (!data->getValue(GOVERNOR_VELOPEN,&rval,g_id)) {
        data->addValue(GOVERNOR_VELOPEN, data_struct.velopen, g_id);
      } else {
        data->setValue(GOVERNOR_VELOPEN, data_struct.velopen, g_id);
      }

      // GOVERNOR_VELCLOSE
      if (!data->getValue(GOVERNOR_VELCLOSE,&rval,g_id)) {
        data->addValue(GOVERNOR_VELCLOSE, data_struct.velclose, g_id);
      } else {
        data->setValue(GOVERNOR_VELCLOSE, data_struct.velclose, g_id);
      }

      // GOVERNOR_PMAX
      if (!data->getValue(GOVERNOR_PMAX,&rval,g_id)) {
        data->addValue(GOVERNOR_PMAX, data_struct.pmax, g_id);
      } else {
        data->setValue(GOVERNOR_PMAX, data_struct.pmax, g_id);
      }

      // GOVERNOR_PMIN
      if (!data->getValue(GOVERNOR_PMIN,&rval,g_id)) {
        data->addValue(GOVERNOR_PMIN, data_struct.pmin, g_id);
      } else {
        data->setValue(GOVERNOR_PMIN, data_struct.pmin, g_id);
      }

      // GOVERNOR_DB2
      if (!data->getValue(GOVERNOR_DB2,&rval,g_id)) {
        data->addValue(GOVERNOR_DB2, data_struct.db2, g_id);
      } else {
        data->setValue(GOVERNOR_DB2, data_struct.db2, g_id);
      }

      // GOVERNOR_GV1
      if (!data->getValue(GOVERNOR_GV1,&rval,g_id)) {
        data->addValue(GOVERNOR_GV1, data_struct.gv1, g_id);
      } else {
        data->setValue(GOVERNOR_GV1, data_struct.gv1, g_id);
      }

      // GOVERNOR_PGV1
      if (!data->getValue(GOVERNOR_PGV1,&rval,g_id)) {
        data->addValue(GOVERNOR_PGV1, data_struct.pgv1, g_id);
      } else {
        data->setValue(GOVERNOR_PGV1, data_struct.pgv1, g_id);
      }

      // GOVERNOR_GV2
      if (!data->getValue(GOVERNOR_GV2,&rval,g_id)) {
        data->addValue(GOVERNOR_GV2, data_struct.gv2, g_id);
      } else {
        data->setValue(GOVERNOR_GV2, data_struct.gv2, g_id);
      }

      // GOVERNOR_PGV2
      if (!data->getValue(GOVERNOR_PGV2,&rval,g_id)) {
        data->addValue(GOVERNOR_PGV2, data_struct.pgv2, g_id);
      } else {
        data->setValue(GOVERNOR_PGV2, data_struct.pgv2, g_id);
      }

      // GOVERNOR_GV3
      if (!data->getValue(GOVERNOR_GV3,&rval,g_id)) {
        data->addValue(GOVERNOR_GV3, data_struct.gv3, g_id);
      } else {
        data->setValue(GOVERNOR_GV3, data_struct.gv3, g_id);
      }

      // GOVERNOR_PGV3
      if (!data->getValue(GOVERNOR_PGV3,&rval,g_id)) {
        data->addValue(GOVERNOR_PGV3, data_struct.pgv3, g_id);
      } else {
        data->setValue(GOVERNOR_PGV3, data_struct.pgv3, g_id);
      }

      // GOVERNOR_GV4
      if (!data->getValue(GOVERNOR_GV4,&rval,g_id)) {
        data->addValue(GOVERNOR_GV4, data_struct.gv4, g_id);
      } else {
        data->setValue(GOVERNOR_GV4, data_struct.gv4, g_id);
      }

      // GOVERNOR_PGV4
      if (!data->getValue(GOVERNOR_PGV4,&rval,g_id)) {
        data->addValue(GOVERNOR_PGV4, data_struct.pgv4, g_id);
      } else {
        data->setValue(GOVERNOR_PGV4, data_struct.pgv4, g_id);
      }

      // GOVERNOR_GV5
      if (!data->getValue(GOVERNOR_GV5,&rval,g_id)) {
        data->addValue(GOVERNOR_GV5, data_struct.gv5, g_id);
      } else {
        data->setValue(GOVERNOR_GV5, data_struct.gv5, g_id);
      }

      // GOVERNOR_PGV5
      if (!data->getValue(GOVERNOR_PGV5,&rval,g_id)) {
        data->addValue(GOVERNOR_PGV5, data_struct.pgv5, g_id);
      } else {
        data->setValue(GOVERNOR_PGV5, data_struct.pgv5, g_id);
      }

      // GOVERNOR_ATURB
      if (!data->getValue(GOVERNOR_ATURB,&rval,g_id)) {
        data->addValue(GOVERNOR_ATURB, data_struct.aturb, g_id);
      } else {
        data->setValue(GOVERNOR_ATURB, data_struct.aturb, g_id);
      }

      // GOVERNOR_BTURB
      if (!data->getValue(GOVERNOR_BTURB,&rval,g_id)) {
        data->addValue(GOVERNOR_BTURB, data_struct.bturb, g_id);
      } else {
        data->setValue(GOVERNOR_BTURB, data_struct.bturb, g_id);
      }

      // GOVERNOR_TTURB
      if (!data->getValue(GOVERNOR_TTURB,&rval,g_id)) {
        data->addValue(GOVERNOR_TTURB, data_struct.tturb, g_id);
      } else {
        data->setValue(GOVERNOR_TTURB, data_struct.tturb, g_id);
      }

      // GOVERNOR_TRATE
      if (!data->getValue(GOVERNOR_TRATE,&rval,g_id)) {
        data->addValue(GOVERNOR_TRATE, data_struct.trate, g_id);
      } else {
        data->setValue(GOVERNOR_TRATE, data_struct.trate, g_id);
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

      // GOVERNOR_DB1
      if (nstr > 3) {
        if (!data->getValue(GOVERNOR_DB1,&rval,g_id)) {
          data->addValue(GOVERNOR_DB1,
              atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_DB1,
              atof(split_line[3].c_str()), g_id);
        }
      } 

      // GOVERNOR_ERR
      if (nstr > 4) {
        if (!data->getValue(GOVERNOR_ERR,&rval,g_id)) {
          data->addValue(GOVERNOR_ERR,
              atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_ERR,
              atof(split_line[4].c_str()), g_id);
        }
      } 

      // GOVERNOR_TD
      if (nstr > 5) {
        if (!data->getValue(GOVERNOR_TD,&rval,g_id)) {
          data->addValue(GOVERNOR_TD,
              atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TD,
              atof(split_line[5].c_str()), g_id);
        }
      } 

      // GOVERNOR_KI
      if (nstr > 6) {
        if (!data->getValue(GOVERNOR_KI,&rval,g_id)) {
          data->addValue(GOVERNOR_KI,
              atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_KI,
              atof(split_line[6].c_str()), g_id);
        }
      } 

      // GOVERNOR_TF
      if (nstr > 7) {
        if (!data->getValue(GOVERNOR_TF,&rval,g_id)) {
          data->addValue(GOVERNOR_TF,
              atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TF,
              atof(split_line[7].c_str()), g_id);
        }
      } 

      // GOVERNOR_KD
      if (nstr > 8) {
        if (!data->getValue(GOVERNOR_KD,&rval,g_id)) {
          data->addValue(GOVERNOR_KD,
              atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_KD,
              atof(split_line[8].c_str()), g_id);
        }
      } 

      // GOVERNOR_KP
      if (nstr > 9) {
        if (!data->getValue(GOVERNOR_KP,&rval,g_id)) {
          data->addValue(GOVERNOR_KP,
              atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_KP,
              atof(split_line[9].c_str()), g_id);
        }
      }

      // GOVERNOR_R
      if (nstr > 10) {
        if (!data->getValue(GOVERNOR_R,&rval,g_id)) {
          data->addValue(GOVERNOR_R,
              atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_R,
              atof(split_line[10].c_str()), g_id);
        }
      } 

      // GOVERNOR_TT
      if (nstr > 11) {
        if (!data->getValue(GOVERNOR_TT,&rval,g_id)) {
          data->addValue(GOVERNOR_TT,
              atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TT,
              atof(split_line[11].c_str()), g_id);
        }
      } 

      // GOVERNOR_KG
      if (nstr > 12) {
        if (!data->getValue(GOVERNOR_KG,&rval,g_id)) {
          data->addValue(GOVERNOR_KG,
              atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_KG,
              atof(split_line[12].c_str()), g_id);
        }
      } 

      // GOVERNOR_TP
      if (nstr > 13) {
        if (!data->getValue(GOVERNOR_TP,&rval,g_id)) {
          data->addValue(GOVERNOR_TP,
              atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TP,
              atof(split_line[13].c_str()), g_id);
        }
      } 

      // GOVERNOR_VELOPEN
      if (nstr > 14) {
        if (!data->getValue(GOVERNOR_VELOPEN,&rval,g_id)) {
          data->addValue(GOVERNOR_VELOPEN,
              atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_VELOPEN,
              atof(split_line[14].c_str()), g_id);
        }
      } 

      // GOVERNOR_VELCLOSE
      if (nstr > 15) {
        if (!data->getValue(GOVERNOR_VELCLOSE,&rval,g_id)) {
          data->addValue(GOVERNOR_VELCLOSE,
              atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_VELCLOSE,
              atof(split_line[15].c_str()), g_id);
        }
      } 

      // GOVERNOR_PMAX
      if (nstr > 16) {
        if (!data->getValue(GOVERNOR_PMAX,&rval,g_id)) {
          data->addValue(GOVERNOR_PMAX,
              atof(split_line[16].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PMAX,
              atof(split_line[16].c_str()), g_id);
        }
      } 

      // GOVERNOR_PMIN
      if (nstr > 17) {
        if (!data->getValue(GOVERNOR_PMIN,&rval,g_id)) {
          data->addValue(GOVERNOR_PMIN,
              atof(split_line[17].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PMIN,
              atof(split_line[17].c_str()), g_id);
        }
      } 

      // GOVERNOR_DB2
      if (nstr > 18) {
        if (!data->getValue(GOVERNOR_DB2,&rval,g_id)) {
          data->addValue(GOVERNOR_DB2,
              atof(split_line[18].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_DB2,
              atof(split_line[18].c_str()), g_id);
        }
      } 

      // GOVERNOR_GV1
      if (nstr > 19) {
        if (!data->getValue(GOVERNOR_GV1,&rval,g_id)) {
          data->addValue(GOVERNOR_GV1,
              atof(split_line[19].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_GV1,
              atof(split_line[19].c_str()), g_id);
        }
      } 

      // GOVERNOR_PGV1
      if (nstr > 20) {
        if (!data->getValue(GOVERNOR_PGV1,&rval,g_id)) {
          data->addValue(GOVERNOR_PGV1,
              atof(split_line[20].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PGV1,
              atof(split_line[20].c_str()), g_id);
        }
      } 

      // GOVERNOR_GV2
      if (nstr > 21) {
        if (!data->getValue(GOVERNOR_GV2,&rval,g_id)) {
          data->addValue(GOVERNOR_GV2,
              atof(split_line[21].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_GV2,
              atof(split_line[21].c_str()), g_id);
        }
      } 

      // GOVERNOR_PGV2
      if (nstr > 22) {
        if (!data->getValue(GOVERNOR_PGV2,&rval,g_id)) {
          data->addValue(GOVERNOR_PGV2,
              atof(split_line[22].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PGV2,
              atof(split_line[22].c_str()), g_id);
        }
      } 

      // GOVERNOR_GV3
      if (nstr > 23) {
        if (!data->getValue(GOVERNOR_GV3,&rval,g_id)) {
          data->addValue(GOVERNOR_GV3,
              atof(split_line[23].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_GV3,
              atof(split_line[23].c_str()), g_id);
        }
      } 

      // GOVERNOR_PGV3
      if (nstr > 24) {
        if (!data->getValue(GOVERNOR_PGV3,&rval,g_id)) {
          data->addValue(GOVERNOR_PGV3,
              atof(split_line[24].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PGV3,
              atof(split_line[24].c_str()), g_id);
        }
      } 

      // GOVERNOR_GV4
      if (nstr > 25) {
        if (!data->getValue(GOVERNOR_GV4,&rval,g_id)) {
          data->addValue(GOVERNOR_GV4,
              atof(split_line[25].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_GV4,
              atof(split_line[25].c_str()), g_id);
        }
      } 

      // GOVERNOR_PGV4
      if (nstr > 26) {
        if (!data->getValue(GOVERNOR_PGV4,&rval,g_id)) {
          data->addValue(GOVERNOR_PGV4,
              atof(split_line[26].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PGV4,
              atof(split_line[26].c_str()), g_id);
        }
      } 

      // GOVERNOR_GV5
      if (nstr > 27) {
        if (!data->getValue(GOVERNOR_GV5,&rval,g_id)) {
          data->addValue(GOVERNOR_GV5,
              atof(split_line[27].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_GV5,
              atof(split_line[27].c_str()), g_id);
        }
      } 

      // GOVERNOR_PGV5
      if (nstr > 28) {
        if (!data->getValue(GOVERNOR_PGV5,&rval,g_id)) {
          data->addValue(GOVERNOR_PGV5,
              atof(split_line[28].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PGV5,
              atof(split_line[28].c_str()), g_id);
        }
      } 

      // GOVERNOR_ATURB
      if (nstr > 29) {
        if (!data->getValue(GOVERNOR_ATURB,&rval,g_id)) {
          data->addValue(GOVERNOR_ATURB,
              atof(split_line[29].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_ATURB,
              atof(split_line[29].c_str()), g_id);
        }
      } 

      // GOVERNOR_BTURB
      if (nstr > 30) {
        if (!data->getValue(GOVERNOR_BTURB,&rval,g_id)) {
          data->addValue(GOVERNOR_BTURB,
              atof(split_line[30].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_BTURB,
              atof(split_line[30].c_str()), g_id);
        }
      } 

      // GOVERNOR_TTURB
      if (nstr > 31) {
        if (!data->getValue(GOVERNOR_TTURB,&rval,g_id)) {
          data->addValue(GOVERNOR_TTURB,
              atof(split_line[31].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TTURB,
              atof(split_line[31].c_str()), g_id);
        }
      } 

      // GOVERNOR_TRATE
      if (nstr > 32) {
        if (!data->getValue(GOVERNOR_TRATE,&rval,g_id)) {
          data->addValue(GOVERNOR_TRATE,
              atof(split_line[32].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TRATE,
              atof(split_line[32].c_str()), g_id);
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
      // GOVERNOR_DB1
      if (nstr > 3) {
        data.db1 = atof(split_line[3].c_str());
      }

      // GOVERNOR_ERR
      if (nstr > 4) {
        data.err = atof(split_line[4].c_str());
      }

      // GOVERNOR_TD
      if (nstr > 5) {
        data.gv_td = atof(split_line[5].c_str());
      }

      // GOVERNOR_KL
      if (nstr > 6) {
        data.gv_ki = atof(split_line[6].c_str());
      }

      // GOVERNOR_TF
      if (nstr > 7) {
        data.gv_tf = atof(split_line[7].c_str());
      }

      // GOVERNOR_KD
      if (nstr > 8) {
        data.gv_kd = atof(split_line[8].c_str());
      }

      // GOVERNOR_KP
      if (nstr > 9) {
        data.gv_kp = atof(split_line[9].c_str());
      }

      // GOVERNOR_R
      if (nstr > 10) {
        data.gv_r = atof(split_line[10].c_str());
      }

      // GOVERNOR_TT
      if (nstr > 11) {
        data.gv_tt = atof(split_line[11].c_str());
      }

      // GOVERNOR_KG
      if (nstr > 12) {
        data.gv_kg = atof(split_line[12].c_str());
      }

      // GOVERNOR_TP
      if (nstr > 13) {
        data.gv_tp = atof(split_line[13].c_str());
      }

      // GOVERNOR_VELOPEN
      if (nstr > 14) {
        data.velopen = atof(split_line[14].c_str());
      }

      // GOVERNOR_VELCLOSE
      if (nstr > 15) {
        data.velclose = atof(split_line[15].c_str());
      }

      // GOVERNOR_PMAX
      if (nstr > 16) {
        data.pmax = atof(split_line[16].c_str());
      }

      // GOVERNOR_PMIN
      if (nstr > 17) {
        data.pmin = atof(split_line[17].c_str());
      }

      // GOVERNOR_DB2
      if (nstr > 18) {
        data.db2 = atof(split_line[18].c_str());
      }

      // GOVERNOR_GV1
      if (nstr > 19) {
        data.gv1 = atof(split_line[19].c_str());
      }

      // GOVERNOR_PGV1
      if (nstr > 20) {
        data.pgv1 = atof(split_line[20].c_str());
      }

      // GOVERNOR_GV2
      if (nstr > 21) {
        data.gv2 = atof(split_line[21].c_str());
      }

      // GOVERNOR_PGV2
      if (nstr > 22) {
        data.pgv2 = atof(split_line[22].c_str());
      }

      // GOVERNOR_GV3
      if (nstr > 23) {
        data.gv3 = atof(split_line[23].c_str());
      }

      // GOVERNOR_PGV3
      if (nstr > 24) {
        data.pgv3 = atof(split_line[24].c_str());
      }

      // GOVERNOR_GV4
      if (nstr > 25) {
        data.gv4 = atof(split_line[25].c_str());
      }

      // GOVERNOR_PGV4
      if (nstr > 26) {
        data.pgv4 = atof(split_line[26].c_str());
      }

      // GOVERNOR_GV5
      if (nstr > 27) {
        data.gv5 = atof(split_line[27].c_str());
      }

      // GOVERNOR_PGV5
      if (nstr > 28) {
        data.pgv5 = atof(split_line[28].c_str());
      }

      // GOVERNOR_ATURB
      if (nstr > 29) {
        data.aturb = atof(split_line[29].c_str());
      }

      // GOVERNOR_BTURB
      if (nstr > 30) {
        data.bturb = atof(split_line[30].c_str());
      }

      // GOVERNOR_TTURB
      if (nstr > 31) {
        data.tturb = atof(split_line[31].c_str());
      }

      // GOVERNOR_TRATE
      if (nstr > 32) {
        data.trate = atof(split_line[32].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
