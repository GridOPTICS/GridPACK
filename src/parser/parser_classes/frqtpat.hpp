/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: July 13, 2016
 *      Author: Bruce Palmer
 */
#ifndef FRQTPAT_HPP
#define FRQTPAT_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class FrqtpatParser
{
  public:
    /**
     * Constructor
     */
    explicit FrqtpatParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~FrqtpatParser()
    {
    }

    /**
     * Extract data from _data_struct and store it in data collection object
     * @param data_struct data struct object
     * @param data data collection object
     * @param gen_id index of generator
     */
    void extract(_data_struct &data_struct,
        gridpack::component::DataCollection *data)
    {
      double rval;
      int ival, r_id;
      bool bval;
      std::string stmp;
      // RELAY_NUMBER
      if (!data->getValue(RELAY_NUMBER,&ival)) {
        ival = 0;
        data->addValue(RELAY_NUMBER, ival+1);
      } else {
        data->setValue(RELAY_NUMBER, ival+1);
      }
      r_id = ival;

      // RELAY_MODEL
      if (!data->getValue(RELAY_MODEL,&stmp,r_id)) {
        data->addValue(RELAY_MODEL, data_struct.model, r_id);
      } else {
        data->setValue(RELAY_MODEL, data_struct.model, r_id);
      }

      // RELAY_GENID
      if (!data->getValue(RELAY_GENID,&stmp,r_id)) {
        data->addValue(RELAY_GENID, data_struct.tag, r_id);
      } else {
        data->setValue(RELAY_GENID, data_struct.tag, r_id);
      }

      // RELAY_MINS
      if (!data->getValue(RELAY_MINS,&ival,r_id)) {
        data->addValue(RELAY_MINS, data_struct.mins, r_id);
      } else {
        data->setValue(RELAY_MINS, data_struct.mins, r_id);
      }

      // RELAY_FREBUS
      if (!data->getValue(RELAY_FREBUS,&ival,r_id)) {
        data->addValue(RELAY_FREBUS, data_struct.frebus, r_id);
      } else {
        data->setValue(RELAY_FREBUS, data_struct.frebus, r_id);
      }

      // RELAY_FL
      if (!data->getValue(RELAY_FL,&rval,r_id)) {
        data->addValue(RELAY_FL, data_struct.fl, r_id);
      } else {
        data->setValue(RELAY_FL, data_struct.fl, r_id);
      }

      // RELAY_FU
      if (!data->getValue(RELAY_FU,&rval,r_id)) {
        data->addValue(RELAY_FU, data_struct.fu, r_id);
      } else {
        data->setValue(RELAY_FU, data_struct.fu, r_id);
      }

      // RELAY_TP
      if (!data->getValue(RELAY_TP,&rval,r_id)) {
        data->addValue(RELAY_TP, data_struct.tp, r_id);
      } else {
        data->setValue(RELAY_TP, data_struct.tp, r_id);
      }

      // RELAY_TB
      if (!data->getValue(RELAY_TB,&rval,r_id)) {
        data->addValue(RELAY_TB, data_struct.tb, r_id);
      } else {
        data->setValue(RELAY_TB, data_struct.tb, r_id);
      }
    }

    /**
     * Parser list of strings and store results in data collection object
     * @param split_line list of tokens from .dyr file
     * @param data data collection object
     * @param model name of generator model
     * @param gen_id index of generator
     */
    void parse(std::vector<std::string> &split_line,
        gridpack::component::DataCollection *data)
    {
      double rval;
      int nstr = split_line.size();
      int ival, r_id;

      // RELAY_NUMBER
      if (!data->getValue(RELAY_NUMBER,&ival)) {
        ival = 0;
        data->addValue(RELAY_NUMBER, ival+1);
      } else {
        data->setValue(RELAY_NUMBER, ival+1);
      }
      r_id = ival;

      // RELAY_MODEL
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(RELAY_MODEL,&stmp,r_id)) {
        data->addValue(RELAY_MODEL, model.c_str(), r_id);
      } else {
        data->setValue(RELAY_MODEL, model.c_str(), r_id);
      }

      // RELAY_GENID
      if (nstr > 4) {
        model = util.clean2Char(split_line[4]);
        if (!data->getValue(RELAY_GENID,&stmp,r_id)) {
          data->addValue(RELAY_GENID, model.c_str(), r_id);
        } else {
          data->setValue(RELAY_GENID, model.c_str(), r_id);
        }
      } 

      // RELAY_MINS
      if (nstr > 1) {
        if (!data->getValue(RELAY_MINS,&ival,r_id)) {
          data->addValue(RELAY_MINS, atoi(split_line[0].c_str()), r_id);
        } else {
          data->setValue(RELAY_MINS, atoi(split_line[0].c_str()), r_id);
        }
      }

      // RELAY_FREBUS
      if (nstr > 2) {
        if (!data->getValue(RELAY_FREBUS,&ival,r_id)) {
          data->addValue(RELAY_FREBUS, atoi(split_line[2].c_str()), r_id);
        } else {
          data->setValue(RELAY_FREBUS, atoi(split_line[2].c_str()), r_id);
        }
      }

      // RELAY_FL
      if (nstr > 5) {
        if (!data->getValue(RELAY_FL,&rval,r_id)) {
          data->addValue(RELAY_FL, atof(split_line[5].c_str()), r_id);
        } else {
          data->setValue(RELAY_FL, atof(split_line[5].c_str()), r_id);
        }
      }

      // RELAY_FU
      if (nstr > 6) {
        if (!data->getValue(RELAY_FU,&rval,r_id)) {
          data->addValue(RELAY_FU, atof(split_line[6].c_str()), r_id);
        } else {
          data->setValue(RELAY_FU, atof(split_line[6].c_str()), r_id);
        }
      }

      // RELAY_TP
      if (nstr > 7) {
        if (!data->getValue(RELAY_TP,&rval,r_id)) {
          data->addValue(RELAY_TP, atof(split_line[7].c_str()), r_id);
        } else {
          data->setValue(RELAY_TP, atof(split_line[7].c_str()), r_id);
        }
      }

      // RELAY_TB
      if (nstr > 8) {
        if (!data->getValue(RELAY_TB,&rval,r_id)) {
          data->addValue(RELAY_TB, atof(split_line[8].c_str()), r_id);
        } else {
          data->setValue(RELAY_TB, atof(split_line[8].c_str()), r_id);
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
      // RELAY_BUSNUMBER               "I"                   integer
      int o_idx;
      o_idx = atoi(split_line[3].c_str());
      data.bus_id = o_idx;

      std::string sval;
      gridpack::utility::StringUtils util;
      sval = util.trimQuotes(split_line[1]);
      util.toUpper(sval);

      // RELAY_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());
      int nstr = split_line.size();

      // RELAY_GENID
      if (nstr > 4) {
        sval = util.clean2Char(split_line[4]);
        strcpy(data.tag, sval.c_str());
      }

      // RELAY_MINS
      if (nstr > 1) {
        data.mins = atoi(split_line[0].c_str());
      }

      // RELAY_FREBUS
      if (nstr > 2) {
        data.frebus = atoi(split_line[2].c_str());
      }

      // RELAY_FL
      if (nstr > 5) {
        data.fl = atof(split_line[5].c_str());
      }

      // RELAY_FU
      if (nstr > 6) {
        data.fu = atof(split_line[6].c_str());
      }

      // RELAY_TP
      if (nstr > 7) {
        data.tp = atof(split_line[7].c_str());
      }

      // RELAY_TB
      if (nstr > 8) {
        data.tb = atof(split_line[8].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
