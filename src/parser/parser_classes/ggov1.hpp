/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: June 20, 2016
 *      Author: Bruce Palmer
 */
#ifndef GGOV1_HPP
#define GGOV1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Ggov1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Ggov1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Ggov1Parser()
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

      // GOVERNOR_RSELECT
      if (!data->getValue(GOVERNOR_RSELECT,&rval,g_id)) {
        data->addValue(GOVERNOR_RSELECT, data_struct.rselect, g_id);
      } else {
        data->setValue(GOVERNOR_RSELECT, data_struct.rselect, g_id);
      }

      // GOVERNOR_FLAGSWITCH
      if (!data->getValue(GOVERNOR_FLAGSWITCH,&rval,g_id)) {
        data->addValue(GOVERNOR_FLAGSWITCH, data_struct.flagswitch, g_id);
      } else {
        data->setValue(GOVERNOR_FLAGSWITCH, data_struct.flagswitch, g_id);
      }

      // GOVERNOR_R
      if (!data->getValue(GOVERNOR_R,&rval,g_id)) {
        data->addValue(GOVERNOR_R, data_struct.gv_r, g_id);
      } else {
        data->setValue(GOVERNOR_R, data_struct.gv_r, g_id);
      }

      // GOVERNOR_TPELEC
      if (!data->getValue(GOVERNOR_TPELEC,&rval,g_id)) {
        data->addValue(GOVERNOR_TPELEC, data_struct.tpelec, g_id);
      } else {
        data->setValue(GOVERNOR_TPELEC, data_struct.tpelec, g_id);
      }

      // GOVERNOR_MAXERR
      if (!data->getValue(GOVERNOR_MAXERR,&rval,g_id)) {
        data->addValue(GOVERNOR_MAXERR, data_struct.maxerr, g_id);
      } else {
        data->setValue(GOVERNOR_MAXERR, data_struct.maxerr, g_id);
      }

      // GOVERNOR_MINERR
      if (!data->getValue(GOVERNOR_MINERR,&rval,g_id)) {
        data->addValue(GOVERNOR_MINERR, data_struct.minerr, g_id);
      } else {
        data->setValue(GOVERNOR_MINERR, data_struct.minerr, g_id);
      }

      // GOVERNOR_KPGOV
      if (!data->getValue(GOVERNOR_KPGOV,&rval,g_id)) {
        data->addValue(GOVERNOR_KPGOV, data_struct.kpgov, g_id);
      } else {
        data->setValue(GOVERNOR_KPGOV, data_struct.kpgov, g_id);
      }

      // GOVERNOR_KIGOV
      if (!data->getValue(GOVERNOR_KIGOV,&rval,g_id)) {
        data->addValue(GOVERNOR_KIGOV, data_struct.kigov, g_id);
      } else {
        data->setValue(GOVERNOR_KIGOV, data_struct.kigov, g_id);
      }

      // GOVERNOR_KDGOV
      if (!data->getValue(GOVERNOR_KDGOV,&rval,g_id)) {
        data->addValue(GOVERNOR_KDGOV, data_struct.kdgov, g_id);
      } else {
        data->setValue(GOVERNOR_KDGOV, data_struct.kdgov, g_id);
      }

      // GOVERNOR_TDGOV
      if (!data->getValue(GOVERNOR_TDGOV,&rval,g_id)) {
        data->addValue(GOVERNOR_TDGOV, data_struct.tdgov, g_id);
      } else {
        data->setValue(GOVERNOR_TDGOV, data_struct.tdgov, g_id);
      }

      // GOVERNOR_VMAX
      if (!data->getValue(GOVERNOR_VMAX,&rval,g_id)) {
        data->addValue(GOVERNOR_VMAX, data_struct.vmax, g_id);
      } else {
        data->setValue(GOVERNOR_VMAX, data_struct.vmax, g_id);
      }

      // GOVERNOR_VMIN
      if (!data->getValue(GOVERNOR_VMIN,&rval,g_id)) {
        data->addValue(GOVERNOR_VMIN, data_struct.vmin, g_id);
      } else {
        data->setValue(GOVERNOR_VMIN, data_struct.vmin, g_id);
      }

      // GOVERNOR_TACT
      if (!data->getValue(GOVERNOR_TACT,&rval,g_id)) {
        data->addValue(GOVERNOR_TACT, data_struct.tact, g_id);
      } else {
        data->setValue(GOVERNOR_TACT, data_struct.tact, g_id);
      }

      // GOVERNOR_KTURB
      if (!data->getValue(GOVERNOR_KTURB,&rval,g_id)) {
        data->addValue(GOVERNOR_KTURB, data_struct.kturb, g_id);
      } else {
        data->setValue(GOVERNOR_KTURB, data_struct.kturb, g_id);
      }

      // GOVERNOR_WFNL
      if (!data->getValue(GOVERNOR_WFNL,&rval,g_id)) {
        data->addValue(GOVERNOR_WFNL, data_struct.wfnl, g_id);
      } else {
        data->setValue(GOVERNOR_WFNL, data_struct.wfnl, g_id);
      }

      // GOVERNOR_TB
      if (!data->getValue(GOVERNOR_TB,&rval,g_id)) {
        data->addValue(GOVERNOR_TB, data_struct.gv_tb, g_id);
      } else {
        data->setValue(GOVERNOR_TB, data_struct.gv_tb, g_id);
      }

      // GOVERNOR_TC
      if (!data->getValue(GOVERNOR_TC,&rval,g_id)) {
        data->addValue(GOVERNOR_TC, data_struct.gv_tc, g_id);
      } else {
        data->setValue(GOVERNOR_TC, data_struct.gv_tc, g_id);
      }

      // GOVERNOR_TENG
      if (!data->getValue(GOVERNOR_TENG,&rval,g_id)) {
        data->addValue(GOVERNOR_TENG, data_struct.teng, g_id);
      } else {
        data->setValue(GOVERNOR_TENG, data_struct.teng, g_id);
      }

      // GOVERNOR_TFLOAD
      if (!data->getValue(GOVERNOR_TFLOAD,&rval,g_id)) {
        data->addValue(GOVERNOR_TFLOAD, data_struct.tfload, g_id);
      } else {
        data->setValue(GOVERNOR_TFLOAD, data_struct.tfload, g_id);
      }

      // GOVERNOR_KPLOAD
      if (!data->getValue(GOVERNOR_KPLOAD,&rval,g_id)) {
        data->addValue(GOVERNOR_KPLOAD, data_struct.kpload, g_id);
      } else {
        data->setValue(GOVERNOR_KPLOAD, data_struct.kpload, g_id);
      }

      // GOVERNOR_KILOAD
      if (!data->getValue(GOVERNOR_KILOAD,&rval,g_id)) {
        data->addValue(GOVERNOR_KILOAD, data_struct.kiload, g_id);
      } else {
        data->setValue(GOVERNOR_KILOAD, data_struct.kiload, g_id);
      }

      // GOVERNOR_LDREF
      if (!data->getValue(GOVERNOR_LDREF,&rval,g_id)) {
        data->addValue(GOVERNOR_LDREF, data_struct.ldref, g_id);
      } else {
        data->setValue(GOVERNOR_LDREF, data_struct.ldref, g_id);
      }

      // GOVERNOR_DM
      if (!data->getValue(GOVERNOR_DM,&rval,g_id)) {
        data->addValue(GOVERNOR_DM, data_struct.gv_dm, g_id);
      } else {
        data->setValue(GOVERNOR_DM, data_struct.gv_dm, g_id);
      }

      // GOVERNOR_ROPEN
      if (!data->getValue(GOVERNOR_ROPEN,&rval,g_id)) {
        data->addValue(GOVERNOR_ROPEN, data_struct.ropen, g_id);
      } else {
        data->setValue(GOVERNOR_ROPEN, data_struct.ropen, g_id);
      }

      // GOVERNOR_RCLOSE
      if (!data->getValue(GOVERNOR_RCLOSE,&rval,g_id)) {
        data->addValue(GOVERNOR_RCLOSE, data_struct.rclose, g_id);
      } else {
        data->setValue(GOVERNOR_RCLOSE, data_struct.rclose, g_id);
      }

      // GOVERNOR_KIMW
      if (!data->getValue(GOVERNOR_KIMW,&rval,g_id)) {
        data->addValue(GOVERNOR_KIMW, data_struct.kimw, g_id);
      } else {
        data->setValue(GOVERNOR_KIMW, data_struct.kimw, g_id);
      }

      // GOVERNOR_ASET
      if (!data->getValue(GOVERNOR_ASET,&rval,g_id)) {
        data->addValue(GOVERNOR_ASET, data_struct.aset, g_id);
      } else {
        data->setValue(GOVERNOR_ASET, data_struct.aset, g_id);
      }

      // GOVERNOR_KA
      if (!data->getValue(GOVERNOR_KA,&rval,g_id)) {
        data->addValue(GOVERNOR_KA, data_struct.gv_ka, g_id);
      } else {
        data->setValue(GOVERNOR_KA, data_struct.gv_ka, g_id);
      }

      // GOVERNOR_TA
      if (!data->getValue(GOVERNOR_TA,&rval,g_id)) {
        data->addValue(GOVERNOR_TA, data_struct.gv_ta, g_id);
      } else {
        data->setValue(GOVERNOR_TA, data_struct.gv_ta, g_id);
      }

      // GOVERNOR_TRATE
      if (!data->getValue(GOVERNOR_TRATE,&rval,g_id)) {
        data->addValue(GOVERNOR_TRATE, data_struct.trate, g_id);
      } else {
        data->setValue(GOVERNOR_TRATE, data_struct.trate, g_id);
      }

      // GOVERNOR_DB
      if (!data->getValue(GOVERNOR_DB,&rval,g_id)) {
        data->addValue(GOVERNOR_DB, data_struct.gv_db, g_id);
      } else {
        data->setValue(GOVERNOR_DB, data_struct.gv_db, g_id);
      }

      // GOVERNOR_TSA
      if (!data->getValue(GOVERNOR_TSA,&rval,g_id)) {
        data->addValue(GOVERNOR_TSA, data_struct.tsa, g_id);
      } else {
        data->setValue(GOVERNOR_TSA, data_struct.tsa, g_id);
      }

      // GOVERNOR_TSB
      if (!data->getValue(GOVERNOR_TSB,&rval,g_id)) {
        data->addValue(GOVERNOR_TSB, data_struct.tsb, g_id);
      } else {
        data->setValue(GOVERNOR_TSB, data_struct.tsb, g_id);
      }

      // GOVERNOR_RUP
      if (!data->getValue(GOVERNOR_RUP,&rval,g_id)) {
        data->addValue(GOVERNOR_RUP, data_struct.rup, g_id);
      } else {
        data->setValue(GOVERNOR_RUP, data_struct.rup, g_id);
      }

      // GOVERNOR_RDOWN
      if (!data->getValue(GOVERNOR_RDOWN,&rval,g_id)) {
        data->addValue(GOVERNOR_RDOWN, data_struct.rdown, g_id);
      } else {
        data->setValue(GOVERNOR_RDOWN, data_struct.rdown, g_id);
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

      // GOVERNOR_RSELECT
      if (nstr > 3) {
        if (!data->getValue(GOVERNOR_RSELECT,&rval,g_id)) {
          data->addValue(GOVERNOR_RSELECT,
              atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_RSELECT,
              atof(split_line[3].c_str()), g_id);
        }
      } 

      // GOVERNOR_FLAGSWITCH
      if (nstr > 4) {
        if (!data->getValue(GOVERNOR_FLAGSWITCH,&rval,g_id)) {
          data->addValue(GOVERNOR_FLAGSWITCH,
              atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_FLAGSWITCH,
              atof(split_line[4].c_str()), g_id);
        }
      } 

      // GOVERNOR_R
      if (nstr > 5) {
        if (!data->getValue(GOVERNOR_R,&rval,g_id)) {
          data->addValue(GOVERNOR_R,
              atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_R,
              atof(split_line[5].c_str()), g_id);
        }
      } 

      // GOVERNOR_TPELEC
      if (nstr > 6) {
        if (!data->getValue(GOVERNOR_TPELEC,&rval,g_id)) {
          data->addValue(GOVERNOR_TPELEC,
              atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TPELEC,
              atof(split_line[6].c_str()), g_id);
        }
      } 

      // GOVERNOR_MAXERR
      if (nstr > 7) {
        if (!data->getValue(GOVERNOR_MAXERR,&rval,g_id)) {
          data->addValue(GOVERNOR_MAXERR,
              atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_MAXERR,
              atof(split_line[7].c_str()), g_id);
        }
      } 

      // GOVERNOR_MINERR
      if (nstr > 8) {
        if (!data->getValue(GOVERNOR_MINERR,&rval,g_id)) {
          data->addValue(GOVERNOR_MINERR,
              atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_MINERR,
              atof(split_line[8].c_str()), g_id);
        }
      } 

      // GOVERNOR_KPGOV
      if (nstr > 9) {
        if (!data->getValue(GOVERNOR_KPGOV,&rval,g_id)) {
          data->addValue(GOVERNOR_KPGOV,
              atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_KPGOV,
              atof(split_line[9].c_str()), g_id);
        }
      }

      // GOVERNOR_KIGOV
      if (nstr > 10) {
        if (!data->getValue(GOVERNOR_KIGOV,&rval,g_id)) {
          data->addValue(GOVERNOR_KIGOV,
              atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_KIGOV,
              atof(split_line[10].c_str()), g_id);
        }
      } 

      // GOVERNOR_KDGOV
      if (nstr > 11) {
        if (!data->getValue(GOVERNOR_KDGOV,&rval,g_id)) {
          data->addValue(GOVERNOR_KDGOV,
              atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_KDGOV,
              atof(split_line[11].c_str()), g_id);
        }
      } 

      // GOVERNOR_TDGOV
      if (nstr > 12) {
        if (!data->getValue(GOVERNOR_TDGOV,&rval,g_id)) {
          data->addValue(GOVERNOR_TDGOV,
              atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TDGOV,
              atof(split_line[12].c_str()), g_id);
        }
      } 

      // GOVERNOR_VMAX
      if (nstr > 13) {
        if (!data->getValue(GOVERNOR_VMAX,&rval,g_id)) {
          data->addValue(GOVERNOR_VMAX,
              atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_VMAX,
              atof(split_line[13].c_str()), g_id);
        }
      } 

      // GOVERNOR_VMIN
      if (nstr > 14) {
        if (!data->getValue(GOVERNOR_VMIN,&rval,g_id)) {
          data->addValue(GOVERNOR_VMIN,
              atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_VMIN,
              atof(split_line[14].c_str()), g_id);
        }
      } 

      // GOVERNOR_TACT
      if (nstr > 15) {
        if (!data->getValue(GOVERNOR_TACT,&rval,g_id)) {
          data->addValue(GOVERNOR_TACT,
              atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TACT,
              atof(split_line[15].c_str()), g_id);
        }
      } 

      // GOVERNOR_KTURB
      if (nstr > 16) {
        if (!data->getValue(GOVERNOR_KTURB,&rval,g_id)) {
          data->addValue(GOVERNOR_KTURB,
              atof(split_line[16].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_KTURB,
              atof(split_line[16].c_str()), g_id);
        }
      } 

      // GOVERNOR_WFNL
      if (nstr > 17) {
        if (!data->getValue(GOVERNOR_WFNL,&rval,g_id)) {
          data->addValue(GOVERNOR_WFNL,
              atof(split_line[17].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_WFNL,
              atof(split_line[17].c_str()), g_id);
        }
      } 

      // GOVERNOR_TB
      if (nstr > 18) {
        if (!data->getValue(GOVERNOR_TB,&rval,g_id)) {
          data->addValue(GOVERNOR_TB,
              atof(split_line[18].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TB,
              atof(split_line[18].c_str()), g_id);
        }
      } 

      // GOVERNOR_TC
      if (nstr > 19) {
        if (!data->getValue(GOVERNOR_TC,&rval,g_id)) {
          data->addValue(GOVERNOR_TC,
              atof(split_line[19].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TC,
              atof(split_line[19].c_str()), g_id);
        }
      } 

      // GOVERNOR_TENG
      if (nstr > 20) {
        if (!data->getValue(GOVERNOR_TENG,&rval,g_id)) {
          data->addValue(GOVERNOR_TENG,
              atof(split_line[20].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TENG,
              atof(split_line[20].c_str()), g_id);
        }
      } 

      // GOVERNOR_TFLOAD
      if (nstr > 21) {
        if (!data->getValue(GOVERNOR_TFLOAD,&rval,g_id)) {
          data->addValue(GOVERNOR_TFLOAD,
              atof(split_line[21].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TFLOAD,
              atof(split_line[21].c_str()), g_id);
        }
      } 

      // GOVERNOR_KPLOAD
      if (nstr > 22) {
        if (!data->getValue(GOVERNOR_KPLOAD,&rval,g_id)) {
          data->addValue(GOVERNOR_KPLOAD,
              atof(split_line[22].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_KPLOAD,
              atof(split_line[22].c_str()), g_id);
        }
      } 

      // GOVERNOR_KILOAD
      if (nstr > 23) {
        if (!data->getValue(GOVERNOR_KILOAD,&rval,g_id)) {
          data->addValue(GOVERNOR_KILOAD,
              atof(split_line[23].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_KILOAD,
              atof(split_line[23].c_str()), g_id);
        }
      } 

      // GOVERNOR_LDREF
      if (nstr > 24) {
        if (!data->getValue(GOVERNOR_LDREF,&rval,g_id)) {
          data->addValue(GOVERNOR_LDREF,
              atof(split_line[24].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_LDREF,
              atof(split_line[24].c_str()), g_id);
        }
      } 

      // GOVERNOR_DM
      if (nstr > 25) {
        if (!data->getValue(GOVERNOR_DM,&rval,g_id)) {
          data->addValue(GOVERNOR_DM,
              atof(split_line[25].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_DM,
              atof(split_line[25].c_str()), g_id);
        }
      } 

      // GOVERNOR_ROPEN
      if (nstr > 26) {
        if (!data->getValue(GOVERNOR_ROPEN,&rval,g_id)) {
          data->addValue(GOVERNOR_ROPEN,
              atof(split_line[26].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_ROPEN,
              atof(split_line[26].c_str()), g_id);
        }
      } 

      // GOVERNOR_RCLOSE
      if (nstr > 27) {
        if (!data->getValue(GOVERNOR_RCLOSE,&rval,g_id)) {
          data->addValue(GOVERNOR_RCLOSE,
              atof(split_line[27].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_RCLOSE,
              atof(split_line[27].c_str()), g_id);
        }
      } 

      // GOVERNOR_KIMW
      if (nstr > 28) {
        if (!data->getValue(GOVERNOR_KIMW,&rval,g_id)) {
          data->addValue(GOVERNOR_KIMW,
              atof(split_line[28].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_KIMW,
              atof(split_line[28].c_str()), g_id);
        }
      } 

      // GOVERNOR_ASET
      if (nstr > 29) {
        if (!data->getValue(GOVERNOR_ASET,&rval,g_id)) {
          data->addValue(GOVERNOR_ASET,
              atof(split_line[29].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_ASET,
              atof(split_line[29].c_str()), g_id);
        }
      } 

      // GOVERNOR_KA
      if (nstr > 30) {
        if (!data->getValue(GOVERNOR_KA,&rval,g_id)) {
          data->addValue(GOVERNOR_KA,
              atof(split_line[30].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_KA,
              atof(split_line[30].c_str()), g_id);
        }
      } 

      // GOVERNOR_TA
      if (nstr > 31) {
        if (!data->getValue(GOVERNOR_TA,&rval,g_id)) {
          data->addValue(GOVERNOR_TA,
              atof(split_line[31].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TA,
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

      // GOVERNOR_DB
      if (nstr > 33) {
        if (!data->getValue(GOVERNOR_DB,&rval,g_id)) {
          data->addValue(GOVERNOR_DB,
              atof(split_line[33].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_DB,
              atof(split_line[33].c_str()), g_id);
        }
      } 

      // GOVERNOR_TSA
      if (nstr > 34) {
        if (!data->getValue(GOVERNOR_TSA,&rval,g_id)) {
          data->addValue(GOVERNOR_TSA,
              atof(split_line[34].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TSA,
              atof(split_line[34].c_str()), g_id);
        }
      } 

      // GOVERNOR_TSB
      if (nstr > 35) {
        if (!data->getValue(GOVERNOR_TSB,&rval,g_id)) {
          data->addValue(GOVERNOR_TSB,
              atof(split_line[35].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_TSB,
              atof(split_line[35].c_str()), g_id);
        }
      } 

      // GOVERNOR_RUP
      if (nstr > 36) {
        if (!data->getValue(GOVERNOR_RUP,&rval,g_id)) {
          data->addValue(GOVERNOR_RUP,
              atof(split_line[36].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_RUP,
              atof(split_line[36].c_str()), g_id);
        }
      } 

      // GOVERNOR_RDOWN
      if (nstr > 37) {
        if (!data->getValue(GOVERNOR_RDOWN,&rval,g_id)) {
          data->addValue(GOVERNOR_RDOWN,
              atof(split_line[37].c_str()), g_id);
        } else {
          data->setValue(GOVERNOR_RDOWN,
              atof(split_line[37].c_str()), g_id);
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
      // GOVERNOR_RSELECT
      if (nstr > 3) {
        data.rselect = atof(split_line[3].c_str());
      }

      // GOVERNOR_FLAGSWITCH
      if (nstr > 4) {
        data.flagswitch = atof(split_line[4].c_str());
      }

      // GOVERNOR_R
      if (nstr > 5) {
        data.gv_r = atof(split_line[5].c_str());
      }

      // GOVERNOR_TPELEC
      if (nstr > 6) {
        data.tpelec = atof(split_line[6].c_str());
      }

      // GOVERNOR_MAXERR
      if (nstr > 7) {
        data.maxerr = atof(split_line[7].c_str());
      }

      // GOVERNOR_MINERR
      if (nstr > 8) {
        data.minerr = atof(split_line[8].c_str());
      }

      // GOVERNOR_KPGOV
      if (nstr > 9) {
        data.kpgov = atof(split_line[9].c_str());
      }

      // GOVERNOR_KIGOV
      if (nstr > 10) {
        data.kigov = atof(split_line[10].c_str());
      }

      // GOVERNOR_KDGOV
      if (nstr > 11) {
        data.kdgov = atof(split_line[11].c_str());
      }

      // GOVERNOR_TDGOV
      if (nstr > 12) {
        data.tdgov = atof(split_line[12].c_str());
      }

      // GOVERNOR_VMAX
      if (nstr > 13) {
        data.vmax = atof(split_line[13].c_str());
      }

      // GOVERNOR_VMIN
      if (nstr > 14) {
        data.vmin = atof(split_line[14].c_str());
      }

      // GOVERNOR_TACT
      if (nstr > 15) {
        data.tact = atof(split_line[15].c_str());
      }

      // GOVERNOR_KTURB
      if (nstr > 16) {
        data.kturb = atof(split_line[16].c_str());
      }

      // GOVERNOR_WFNL
      if (nstr > 17) {
        data.wfnl = atof(split_line[17].c_str());
      }

      // GOVERNOR_TB
      if (nstr > 18) {
        data.gv_tb = atof(split_line[18].c_str());
      }

      // GOVERNOR_TC
      if (nstr > 19) {
        data.gv_tc = atof(split_line[19].c_str());
      }

      // GOVERNOR_TENG
      if (nstr > 20) {
        data.teng = atof(split_line[20].c_str());
      }

      // GOVERNOR_TFLOAD
      if (nstr > 21) {
        data.tfload = atof(split_line[21].c_str());
      }

      // GOVERNOR_KPLOAD
      if (nstr > 22) {
        data.kpload = atof(split_line[22].c_str());
      }

      // GOVERNOR_KILOAD
      if (nstr > 23) {
        data.kiload = atof(split_line[23].c_str());
      }

      // GOVERNOR_LDREF
      if (nstr > 24) {
        data.ldref = atof(split_line[24].c_str());
      }

      // GOVERNOR_DM
      if (nstr > 25) {
        data.gv_dm = atof(split_line[25].c_str());
      }

      // GOVERNOR_ROPEN
      if (nstr > 26) {
        data.ropen = atof(split_line[26].c_str());
      }

      // GOVERNOR_RCLOSE
      if (nstr > 27) {
        data.rclose = atof(split_line[27].c_str());
      }

      // GOVERNOR_KIMW
      if (nstr > 28) {
        data.kimw = atof(split_line[28].c_str());
      }

      // GOVERNOR_ASET
      if (nstr > 29) {
        data.aset = atof(split_line[29].c_str());
      }

      // GOVERNOR_KA
      if (nstr > 30) {
        data.gv_ka = atof(split_line[30].c_str());
      }

      // GOVERNOR_TA
      if (nstr > 31) {
        data.gv_ta = atof(split_line[31].c_str());
      }

      // GOVERNOR_TRATE
      if (nstr > 32) {
        data.trate = atof(split_line[32].c_str());
      }

      // GOVERNOR_DB
      if (nstr > 33) {
        data.gv_db = atof(split_line[33].c_str());
      }

      // GOVERNOR_TSA
      if (nstr > 34) {
        data.tsa = atof(split_line[34].c_str());
      }

      // GOVERNOR_TSB
      if (nstr > 35) {
        data.tsb = atof(split_line[35].c_str());
      }

      // GOVERNOR_RUP
      if (nstr > 36) {
        data.rup = atof(split_line[36].c_str());
      }

      // GOVERNOR_RDOWN
      if (nstr > 37) {
        data.rdown = atof(split_line[37].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
