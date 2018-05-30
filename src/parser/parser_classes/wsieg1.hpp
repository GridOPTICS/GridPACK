/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: June 16, 2016
 *      Author: Bruce Palmer
 */
#ifndef WSIEG1_HPP
#define WSIEG1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Wsieg1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Wsieg1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Wsieg1Parser()
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
      int ival;
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

      // GOVERNOR_JBUS
      if (!data->getValue(GOVERNOR_JBUS,&ival,g_id)) {
        data->addValue(GOVERNOR_JBUS, data_struct.jbus, g_id);
      } else {
        data->setValue(GOVERNOR_JBUS, data_struct.jbus, g_id);
      }

      // GOVERNOR_M
      if (!data->getValue(GOVERNOR_M,&ival,g_id)) {
        data->addValue(GOVERNOR_M, data_struct.gv_m, g_id);
      } else {
        data->setValue(GOVERNOR_M, data_struct.gv_m, g_id);
      }

      // GOVERNOR_K
      if (!data->getValue(GOVERNOR_K,&rval,g_id)) {
        data->addValue(GOVERNOR_K, data_struct.gv_k, g_id);
      } else {
        data->setValue(GOVERNOR_K, data_struct.gv_k, g_id);
      }

      // GOVERNOR_T1
      if (!data->getValue(GOVERNOR_T1,&rval,g_id)) {
        data->addValue(GOVERNOR_T1, data_struct.gv_t1, g_id);
      } else {
        data->setValue(GOVERNOR_T1, data_struct.gv_t1, g_id);
      }

      // GOVERNOR_T2
      if (!data->getValue(GOVERNOR_T2,&rval,g_id)) {
        data->addValue(GOVERNOR_T2, data_struct.gv_t2, g_id);
      } else {
        data->setValue(GOVERNOR_T2, data_struct.gv_t2, g_id);
      }

      // GOVERNOR_T3
      if (!data->getValue(GOVERNOR_T3,&rval,g_id)) {
        data->addValue(GOVERNOR_T3, data_struct.gv_t2, g_id);
      } else {
        data->setValue(GOVERNOR_T3, data_struct.gv_t3, g_id);
      }

      // GOVERNOR_UO
      if (!data->getValue(GOVERNOR_UO,&rval,g_id)) {
        data->addValue(GOVERNOR_UO, data_struct.gv_uo, g_id);
      } else {
        data->setValue(GOVERNOR_UO, data_struct.gv_uo, g_id);
      }

      // GOVERNOR_UC
      if (!data->getValue(GOVERNOR_UC,&rval,g_id)) {
        data->addValue(GOVERNOR_UC, data_struct.gv_uc, g_id);
      } else {
        data->setValue(GOVERNOR_UC, data_struct.gv_uc, g_id);
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

      // GOVERNOR_T4
      if (!data->getValue(GOVERNOR_T4,&rval,g_id)) {
        data->addValue(GOVERNOR_T4, data_struct.gv_t4, g_id);
      } else {
        data->setValue(GOVERNOR_T4, data_struct.gv_t4, g_id);
      }

      // GOVERNOR_K1
      if (!data->getValue(GOVERNOR_K1,&rval,g_id)) {
        data->addValue(GOVERNOR_K1, data_struct.gv_k1, g_id);
      } else {
        data->setValue(GOVERNOR_K1, data_struct.gv_k1, g_id);
      }

      // GOVERNOR_K2
      if (!data->getValue(GOVERNOR_K2,&rval,g_id)) {
        data->addValue(GOVERNOR_K2, data_struct.gv_k2, g_id);
      } else {
        data->setValue(GOVERNOR_K2, data_struct.gv_k2, g_id);
      }

      // GOVERNOR_T5
      if (!data->getValue(GOVERNOR_T5,&rval,g_id)) {
        data->addValue(GOVERNOR_T5, data_struct.gv_t5, g_id);
      } else {
        data->setValue(GOVERNOR_T5, data_struct.gv_t5, g_id);
      }

      // GOVERNOR_K3
      if (!data->getValue(GOVERNOR_K3,&rval,g_id)) {
        data->addValue(GOVERNOR_K3, data_struct.gv_k3, g_id);
      } else {
        data->setValue(GOVERNOR_K3, data_struct.gv_k3, g_id);
      }

      // GOVERNOR_K4
      if (!data->getValue(GOVERNOR_K4,&rval,g_id)) {
        data->addValue(GOVERNOR_K4, data_struct.gv_k4, g_id);
      } else {
        data->setValue(GOVERNOR_K4, data_struct.gv_k4, g_id);
      }

      // GOVERNOR_T6
      if (!data->getValue(GOVERNOR_T6,&rval,g_id)) {
        data->addValue(GOVERNOR_T6, data_struct.gv_t6, g_id);
      } else {
        data->setValue(GOVERNOR_T6, data_struct.gv_t6, g_id);
      }

      // GOVERNOR_K5
      if (!data->getValue(GOVERNOR_K5,&rval,g_id)) {
        data->addValue(GOVERNOR_K5, data_struct.gv_k5, g_id);
      } else {
        data->setValue(GOVERNOR_K5, data_struct.gv_k5, g_id);
      }

      // GOVERNOR_K6
      if (!data->getValue(GOVERNOR_K6,&rval,g_id)) {
        data->addValue(GOVERNOR_K6, data_struct.gv_k6, g_id);
      } else {
        data->setValue(GOVERNOR_K6, data_struct.gv_k6, g_id);
      }

      // GOVERNOR_T7
      if (!data->getValue(GOVERNOR_T7,&rval,g_id)) {
        data->addValue(GOVERNOR_T7, data_struct.gv_t7, g_id);
      } else {
        data->setValue(GOVERNOR_T7, data_struct.gv_t7, g_id);
      }

      // GOVERNOR_K7
      if (!data->getValue(GOVERNOR_K7,&rval,g_id)) {
        data->addValue(GOVERNOR_K7, data_struct.gv_k7, g_id);
      } else {
        data->setValue(GOVERNOR_K7, data_struct.gv_k7, g_id);
      }

      // GOVERNOR_K8
      if (!data->getValue(GOVERNOR_K8,&rval,g_id)) {
        data->addValue(GOVERNOR_K8, data_struct.gv_k8, g_id);
      } else {
        data->setValue(GOVERNOR_K8, data_struct.gv_k8, g_id);
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

      // GOVERNOR_IBLOCK
      if (!data->getValue(GOVERNOR_IBLOCK,&ival,g_id)) {
        data->addValue(GOVERNOR_IBLOCK, data_struct.iblock, g_id);
      } else {
        data->setValue(GOVERNOR_IBLOCK, data_struct.iblock, g_id);
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
      int ival;
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

      // GOVERNOR_JBUS
      if (nstr > 3) {
        if (!data->getValue(GOVERNOR_JBUS,&ival,g_id)) {
          data->addValue(GOVERNOR_JBUS,
              atoi(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_JBUS,
              atoi(split_line[3].c_str()), g_id);
        }
      } 

      // GOVERNOR_M
      if (nstr > 4) {
        if (!data->getValue(GOVERNOR_M,&ival,g_id)) {
          data->addValue(GOVERNOR_M,
              atoi(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_M,
              atoi(split_line[4].c_str()), g_id);
        }
      } 

      // GOVERNOR_K
      if (nstr > 5) {
        if (!data->getValue(GOVERNOR_K,&rval,g_id)) {
          data->addValue(GOVERNOR_K,
              atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_K,
              atof(split_line[5].c_str()), g_id);
        }
      } 

      // GOVERNOR_T1
      if (nstr > 6) {
        if (!data->getValue(GOVERNOR_T1,&rval,g_id)) {
          data->addValue(GOVERNOR_T1,
              atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_T1,
              atof(split_line[6].c_str()), g_id);
        }
      } 

      // GOVERNOR_T2
      if (nstr > 7) {
        if (!data->getValue(GOVERNOR_T2,&rval,g_id)) {
          data->addValue(GOVERNOR_T2,
              atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_T2,
              atof(split_line[7].c_str()), g_id);
        }
      } 

      // GOVERNOR_T3
      if (nstr > 8) {
        if (!data->getValue(GOVERNOR_T3,&rval,g_id)) {
          data->addValue(GOVERNOR_T3,
              atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_T3,
              atof(split_line[8].c_str()), g_id);
        }
      } 

      // GOVERNOR_UO
      if (nstr > 9) {
        if (!data->getValue(GOVERNOR_UO,&rval,g_id)) {
          data->addValue(GOVERNOR_UO,
              atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_UO,
              atof(split_line[9].c_str()), g_id);
        }
      }

      // GOVERNOR_UC
      if (nstr > 10) {
        if (!data->getValue(GOVERNOR_UC,&rval,g_id)) {
          data->addValue(GOVERNOR_UC,
              atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_UC,
              atof(split_line[10].c_str()), g_id);
        }
      } 

      // GOVERNOR_PMAX
      if (nstr > 11) {
        if (!data->getValue(GOVERNOR_PMAX,&rval,g_id)) {
          data->addValue(GOVERNOR_PMAX,
              atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PMAX,
              atof(split_line[11].c_str()), g_id);
        }
      } 

      // GOVERNOR_PMIN
      if (nstr > 12) {
        if (!data->getValue(GOVERNOR_PMIN,&rval,g_id)) {
          data->addValue(GOVERNOR_PMIN,
              atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PMIN,
              atof(split_line[12].c_str()), g_id);
        }
      } 

      // GOVERNOR_T4
      if (nstr > 13) {
        if (!data->getValue(GOVERNOR_T4,&rval,g_id)) {
          data->addValue(GOVERNOR_T4,
              atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_T4,
              atof(split_line[13].c_str()), g_id);
        }
      } 

      // GOVERNOR_K1
      if (nstr > 14) {
        if (!data->getValue(GOVERNOR_K1,&rval,g_id)) {
          data->addValue(GOVERNOR_K1,
              atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_K1,
              atof(split_line[14].c_str()), g_id);
        }
      } 

      // GOVERNOR_K2
      if (nstr > 15) {
        if (!data->getValue(GOVERNOR_K2,&rval,g_id)) {
          data->addValue(GOVERNOR_K2,
              atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_K2,
              atof(split_line[15].c_str()), g_id);
        }
      } 

      // GOVERNOR_T5
      if (nstr > 16) {
        if (!data->getValue(GOVERNOR_T5,&rval,g_id)) {
          data->addValue(GOVERNOR_T5,
              atof(split_line[16].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_T5,
              atof(split_line[16].c_str()), g_id);
        }
      } 

      // GOVERNOR_K3
      if (nstr > 17) {
        if (!data->getValue(GOVERNOR_K3,&rval,g_id)) {
          data->addValue(GOVERNOR_K3,
              atof(split_line[17].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_K3,
              atof(split_line[17].c_str()), g_id);
        }
      } 

      // GOVERNOR_K4
      if (nstr > 18) {
        if (!data->getValue(GOVERNOR_K4,&rval,g_id)) {
          data->addValue(GOVERNOR_K4,
              atof(split_line[18].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_K4,
              atof(split_line[18].c_str()), g_id);
        }
      } 

      // GOVERNOR_T6
      if (nstr > 19) {
        if (!data->getValue(GOVERNOR_T6,&rval,g_id)) {
          data->addValue(GOVERNOR_T6,
              atof(split_line[19].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_T6,
              atof(split_line[19].c_str()), g_id);
        }
      } 

      // GOVERNOR_K5
      if (nstr > 20) {
        if (!data->getValue(GOVERNOR_K5,&rval,g_id)) {
          data->addValue(GOVERNOR_K5,
              atof(split_line[20].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_K5,
              atof(split_line[20].c_str()), g_id);
        }
      } 

      // GOVERNOR_K6
      if (nstr > 21) {
        if (!data->getValue(GOVERNOR_K6,&rval,g_id)) {
          data->addValue(GOVERNOR_K6,
              atof(split_line[21].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_K6,
              atof(split_line[21].c_str()), g_id);
        }
      } 

      // GOVERNOR_T7
      if (nstr > 22) {
        if (!data->getValue(GOVERNOR_T7,&rval,g_id)) {
          data->addValue(GOVERNOR_T7,
              atof(split_line[22].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_T7,
              atof(split_line[22].c_str()), g_id);
        }
      } 

      // GOVERNOR_K7
      if (nstr > 23) {
        if (!data->getValue(GOVERNOR_K7,&rval,g_id)) {
          data->addValue(GOVERNOR_K7,
              atof(split_line[23].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_K7,
              atof(split_line[23].c_str()), g_id);
        }
      } 

      // GOVERNOR_K8
      if (nstr > 24) {
        if (!data->getValue(GOVERNOR_K8,&rval,g_id)) {
          data->addValue(GOVERNOR_K8,
              atof(split_line[24].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_K8,
              atof(split_line[24].c_str()), g_id);
        }
      } 

      // GOVERNOR_DB1
      if (nstr > 25) {
        if (!data->getValue(GOVERNOR_DB1,&rval,g_id)) {
          data->addValue(GOVERNOR_DB1,
              atof(split_line[25].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_DB1,
              atof(split_line[25].c_str()), g_id);
        }
      } 

      // GOVERNOR_ERR
      if (nstr > 26) {
        if (!data->getValue(GOVERNOR_ERR,&rval,g_id)) {
          data->addValue(GOVERNOR_ERR,
              atof(split_line[26].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_ERR,
              atof(split_line[26].c_str()), g_id);
        }
      } 

      // GOVERNOR_DB2
      if (nstr > 27) {
        if (!data->getValue(GOVERNOR_DB2,&rval,g_id)) {
          data->addValue(GOVERNOR_DB2,
              atof(split_line[27].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_DB2,
              atof(split_line[27].c_str()), g_id);
        }
      } 

      // GOVERNOR_GV1
      if (nstr > 28) {
        if (!data->getValue(GOVERNOR_GV1,&rval,g_id)) {
          data->addValue(GOVERNOR_GV1,
              atof(split_line[28].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_GV1,
              atof(split_line[28].c_str()), g_id);
        }
      } 

      // GOVERNOR_PGV1
      if (nstr > 29) {
        if (!data->getValue(GOVERNOR_PGV1,&rval,g_id)) {
          data->addValue(GOVERNOR_PGV1,
              atof(split_line[29].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PGV1,
              atof(split_line[29].c_str()), g_id);
        }
      } 

      // GOVERNOR_GV2
      if (nstr > 30) {
        if (!data->getValue(GOVERNOR_GV2,&rval,g_id)) {
          data->addValue(GOVERNOR_GV2,
              atof(split_line[30].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_GV2,
              atof(split_line[30].c_str()), g_id);
        }
      } 

      // GOVERNOR_PGV2
      if (nstr > 31) {
        if (!data->getValue(GOVERNOR_PGV2,&rval,g_id)) {
          data->addValue(GOVERNOR_PGV2,
              atof(split_line[31].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PGV2,
              atof(split_line[31].c_str()), g_id);
        }
      } 

      // GOVERNOR_GV3
      if (nstr > 32) {
        if (!data->getValue(GOVERNOR_GV3,&rval,g_id)) {
          data->addValue(GOVERNOR_GV3,
              atof(split_line[32].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_GV3,
              atof(split_line[32].c_str()), g_id);
        }
      } 

      // GOVERNOR_PGV3
      if (nstr > 33) {
        if (!data->getValue(GOVERNOR_PGV3,&rval,g_id)) {
          data->addValue(GOVERNOR_PGV3,
              atof(split_line[33].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PGV3,
              atof(split_line[33].c_str()), g_id);
        }
      } 

      // GOVERNOR_GV4
      if (nstr > 34) {
        if (!data->getValue(GOVERNOR_GV4,&rval,g_id)) {
          data->addValue(GOVERNOR_GV4,
              atof(split_line[34].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_GV4,
              atof(split_line[34].c_str()), g_id);
        }
      } 

      // GOVERNOR_PGV4
      if (nstr > 35) {
        if (!data->getValue(GOVERNOR_PGV4,&rval,g_id)) {
          data->addValue(GOVERNOR_PGV4,
              atof(split_line[35].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PGV4,
              atof(split_line[35].c_str()), g_id);
        }
      } 

      // GOVERNOR_GV5
      if (nstr > 36) {
        if (!data->getValue(GOVERNOR_GV5,&rval,g_id)) {
          data->addValue(GOVERNOR_GV5,
              atof(split_line[36].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_GV4,
              atof(split_line[36].c_str()), g_id);
        }
      } 

      // GOVERNOR_PGV5
      if (nstr > 37) {
        if (!data->getValue(GOVERNOR_PGV5,&rval,g_id)) {
          data->addValue(GOVERNOR_PGV5,
              atof(split_line[37].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_PGV5,
              atof(split_line[37].c_str()), g_id);
        }
      } 

      // GOVERNOR_IBLOCK
      if (nstr > 38) {
        if (!data->getValue(GOVERNOR_IBLOCK,&ival,g_id)) {
          data->addValue(GOVERNOR_PGV5,
              atoi(split_line[38].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_IBLOCK,
              atoi(split_line[38].c_str()), g_id);
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
      // GOVERNOR_JBUS
      if (nstr > 3) {
        data.jbus = atoi(split_line[3].c_str());
      }

      // GOVERNOR_M
      if (nstr > 4) {
        data.gv_m = atoi(split_line[4].c_str());
      }

      // GOVERNOR_K
      if (nstr > 5) {
        data.gv_k = atof(split_line[5].c_str());
      }

      // GOVERNOR_T1
      if (nstr > 6) {
        data.gv_t1 = atof(split_line[6].c_str());
      }

      // GOVERNOR_T2
      if (nstr > 7) {
        data.gv_t2 = atof(split_line[7].c_str());
      }

      // GOVERNOR_T3
      if (nstr > 8) {
        data.gv_t3 = atof(split_line[8].c_str());
      }

      // GOVERNOR_UO
      if (nstr > 9) {
        data.gv_uo = atof(split_line[9].c_str());
      }

      // GOVERNOR_UC
      if (nstr > 10) {
        data.gv_uc = atof(split_line[10].c_str());
      }

      // GOVERNOR_PMAX
      if (nstr > 11) {
        data.pmax = atof(split_line[11].c_str());
      }

      // GOVERNOR_PMIN
      if (nstr > 12) {
        data.pmin = atof(split_line[12].c_str());
      }

      // GOVERNOR_T4
      if (nstr > 13) {
        data.gv_t4 = atof(split_line[13].c_str());
      }

      // GOVERNOR_K1
      if (nstr > 14) {
        data.gv_k1 = atof(split_line[14].c_str());
      }

      // GOVERNOR_K2
      if (nstr > 15) {
        data.gv_k2 = atof(split_line[15].c_str());
      }

      // GOVERNOR_T5
      if (nstr > 16) {
        data.gv_t5 = atof(split_line[16].c_str());
      }

      // GOVERNOR_K3
      if (nstr > 17) {
        data.gv_k3 = atof(split_line[17].c_str());
      }

      // GOVERNOR_K4
      if (nstr > 18) {
        data.gv_k4 = atof(split_line[18].c_str());
      }

      // GOVERNOR_T6
      if (nstr > 19) {
        data.gv_t6 = atof(split_line[19].c_str());
      }

      // GOVERNOR_K5
      if (nstr > 20) {
        data.gv_k5 = atof(split_line[20].c_str());
      }

      // GOVERNOR_K6
      if (nstr > 21) {
        data.gv_k6 = atof(split_line[21].c_str());
      }

      // GOVERNOR_T7
      if (nstr > 22) {
        data.gv_t7 = atof(split_line[22].c_str());
      }

      // GOVERNOR_K7
      if (nstr > 23) {
        data.gv_k7 = atof(split_line[23].c_str());
      }

      // GOVERNOR_K8
      if (nstr > 24) {
        data.gv_k8 = atof(split_line[24].c_str());
      }

      // GOVERNOR_DB1
      if (nstr > 25) {
        data.db1 = atof(split_line[25].c_str());
      }

      // GOVERNOR_ERR
      if (nstr > 26) {
        data.err = atof(split_line[26].c_str());
      }

      // GOVERNOR_DB2
      if (nstr > 27) {
        data.db2 = atof(split_line[27].c_str());
      }

      // GOVERNOR_GV1
      if (nstr > 28) {
        data.gv1 = atof(split_line[28].c_str());
      }

      // GOVERNOR_PGV1
      if (nstr > 29) {
        data.pgv1 = atof(split_line[29].c_str());
      }

      // GOVERNOR_GV2
      if (nstr > 30) {
        data.gv2 = atof(split_line[30].c_str());
      }

      // GOVERNOR_PGV2
      if (nstr > 31) {
        data.pgv2 = atof(split_line[31].c_str());
      }

      // GOVERNOR_GV3
      if (nstr > 32) {
        data.gv3 = atof(split_line[32].c_str());
      }

      // GOVERNOR_PGV3
      if (nstr > 33) {
        data.pgv3 = atof(split_line[33].c_str());
      }

      // GOVERNOR_GV4
      if (nstr > 34) {
        data.gv4 = atof(split_line[34].c_str());
      }

      // GOVERNOR_PGV4
      if (nstr > 35) {
        data.pgv4 = atof(split_line[35].c_str());
      }

      // GOVERNOR_GV5
      if (nstr > 36) {
        data.gv5 = atof(split_line[36].c_str());
      }

      // GOVERNOR_PGV5
      if (nstr > 37) {
        data.pgv5 = atof(split_line[37].c_str());
      }

      // GOVERNOR_IBLOCK
      if (nstr > 38) {
        data.iblock = atoi(split_line[38].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
