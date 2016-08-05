/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: July 13, 2016
 *      Author: Bruce Palmer
 */
#ifndef LVSHBL_HPP
#define LVSHBL_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class LvshblParser
{
  public:
    /**
     * Constructor
     */
    explicit LvshblParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~LvshblParser()
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

      // RELAY_LID
      if (!data->getValue(RELAY_LID,&stmp,r_id)) {
        data->addValue(RELAY_LID, data_struct.tag, r_id);
      } else {
        data->setValue(RELAY_LID, data_struct.tag, r_id);
      }

      // RELAY_JBUS
      if (!data->getValue(RELAY_JBUS,&ival,r_id)) {
        data->addValue(RELAY_JBUS, data_struct.jbus, r_id);
      } else {
        data->setValue(RELAY_JBUS, data_struct.jbus, r_id);
      }

      // RELAY_V1
      if (!data->getValue(RELAY_V1,&rval,r_id)) {
        data->addValue(RELAY_V1, data_struct.v1, r_id);
      } else {
        data->setValue(RELAY_V1, data_struct.v1, r_id);
      }

      // RELAY_T1
      if (!data->getValue(RELAY_T1,&rval,r_id)) {
        data->addValue(RELAY_T1, data_struct.t1, r_id);
      } else {
        data->setValue(RELAY_T1, data_struct.t1, r_id);
      }

      // RELAY_F1
      if (!data->getValue(RELAY_F1,&rval,r_id)) {
        data->addValue(RELAY_F1, data_struct.f1, r_id);
      } else {
        data->setValue(RELAY_F1, data_struct.f1, r_id);
      }

      // RELAY_V2
      if (!data->getValue(RELAY_V2,&rval,r_id)) {
        data->addValue(RELAY_V2, data_struct.v2, r_id);
      } else {
        data->setValue(RELAY_V2, data_struct.v2, r_id);
      }

      // RELAY_T2
      if (!data->getValue(RELAY_T2,&rval,r_id)) {
        data->addValue(RELAY_T2, data_struct.t2, r_id);
      } else {
        data->setValue(RELAY_T2, data_struct.t2, r_id);
      }

      // RELAY_F2
      if (!data->getValue(RELAY_F2,&rval,r_id)) {
        data->addValue(RELAY_F2, data_struct.f2, r_id);
      } else {
        data->setValue(RELAY_F2, data_struct.f2, r_id);
      }

      // RELAY_V3
      if (!data->getValue(RELAY_V3,&rval,r_id)) {
        data->addValue(RELAY_V3, data_struct.v3, r_id);
      } else {
        data->setValue(RELAY_V3, data_struct.v3, r_id);
      }

      // RELAY_T3
      if (!data->getValue(RELAY_T3,&rval,r_id)) {
        data->addValue(RELAY_T3, data_struct.t3, r_id);
      } else {
        data->setValue(RELAY_T3, data_struct.t3, r_id);
      }

      // RELAY_F3
      if (!data->getValue(RELAY_F3,&rval,r_id)) {
        data->addValue(RELAY_F3, data_struct.f3, r_id);
      } else {
        data->setValue(RELAY_F3, data_struct.f3, r_id);
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

      // RELAY_LID
      if (nstr > 2) {
        model = util.clean2Char(split_line[2]);
        if (!data->getValue(RELAY_LID,&stmp,r_id)) {
          data->addValue(RELAY_LID, model.c_str(), r_id);
        } else {
          data->setValue(RELAY_LID, model.c_str(), r_id);
        }
      } 

      // RELAY_JBUS
      if (nstr > 3) {
        if (!data->getValue(RELAY_JBUS,&ival,r_id)) {
          data->addValue(RELAY_JBUS, atoi(split_line[3].c_str()), r_id);
        } else {
          data->setValue(RELAY_JBUS, atoi(split_line[3].c_str()), r_id);
        }
      }

      // RELAY_V1
      if (nstr > 4) {
        if (!data->getValue(RELAY_V1,&rval,r_id)) {
          data->addValue(RELAY_V1, atof(split_line[4].c_str()), r_id);
        } else {
          data->setValue(RELAY_V1, atof(split_line[4].c_str()), r_id);
        }
      }

      // RELAY_T1
      if (nstr > 5) {
        if (!data->getValue(RELAY_T1,&rval,r_id)) {
          data->addValue(RELAY_T1, atof(split_line[5].c_str()), r_id);
        } else {
          data->setValue(RELAY_T1, atof(split_line[5].c_str()), r_id);
        }
      }

      // RELAY_F1
      if (nstr > 6) {
        if (!data->getValue(RELAY_F1,&rval,r_id)) {
          data->addValue(RELAY_F1, atof(split_line[6].c_str()), r_id);
        } else {
          data->setValue(RELAY_F1, atof(split_line[6].c_str()), r_id);
        }
      }

      // RELAY_V2
      if (nstr > 7) {
        if (!data->getValue(RELAY_V2,&rval,r_id)) {
          data->addValue(RELAY_V2, atof(split_line[7].c_str()), r_id);
        } else {
          data->setValue(RELAY_V2, atof(split_line[7].c_str()), r_id);
        }
      }

      // RELAY_T2
      if (nstr > 8) {
        if (!data->getValue(RELAY_T2,&rval,r_id)) {
          data->addValue(RELAY_T2, atof(split_line[8].c_str()), r_id);
        } else {
          data->setValue(RELAY_T2, atof(split_line[8].c_str()), r_id);
        }
      }

      // RELAY_F2
      if (nstr > 9) {
        if (!data->getValue(RELAY_F2,&rval,r_id)) {
          data->addValue(RELAY_F2, atof(split_line[9].c_str()), r_id);
        } else {
          data->setValue(RELAY_F2, atof(split_line[9].c_str()), r_id);
        }
      }

      // RELAY_V3
      if (nstr > 10) {
        if (!data->getValue(RELAY_V3,&rval,r_id)) {
          data->addValue(RELAY_V3, atof(split_line[10].c_str()), r_id);
        } else {
          data->setValue(RELAY_V3, atof(split_line[10].c_str()), r_id);
        }
      }

      // RELAY_T3
      if (nstr > 11) {
        if (!data->getValue(RELAY_T3,&rval,r_id)) {
          data->addValue(RELAY_T3, atof(split_line[11].c_str()), r_id);
        } else {
          data->setValue(RELAY_T3, atof(split_line[11].c_str()), r_id);
        }
      }

      // RELAY_F3
      if (nstr > 12) {
        if (!data->getValue(RELAY_F3,&rval,r_id)) {
          data->addValue(RELAY_F3, atof(split_line[12].c_str()), r_id);
        } else {
          data->setValue(RELAY_F3, atof(split_line[12].c_str()), r_id);
        }
      }

      // RELAY_TB
      if (nstr > 13) {
        if (!data->getValue(RELAY_TB,&rval,r_id)) {
          data->addValue(RELAY_TB, atof(split_line[13].c_str()), r_id);
        } else {
          data->setValue(RELAY_TB, atof(split_line[13].c_str()), r_id);
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
      o_idx = atoi(split_line[0].c_str());
      data.bus_id = o_idx;

      std::string sval;
      gridpack::utility::StringUtils util;
      sval = util.trimQuotes(split_line[1]);
      util.toUpper(sval);

      // RELAY_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());
      int nstr = split_line.size();

      // RELAY_LID
      if (nstr > 2) {
        sval = util.clean2Char(split_line[2]);
          strcpy(data.tag, sval.c_str());
      }

      // RELAY_JBUS
      if (nstr > 3) {
        data.jbus = atoi(split_line[3].c_str());
      }

      // RELAY_V1
      if (nstr > 4) {
        data.v1 = atof(split_line[4].c_str());
      }

      // RELAY_T1
      if (nstr > 5) {
        data.t1 = atof(split_line[5].c_str());
      }

      // RELAY_F1
      if (nstr > 6) {
        data.f1 = atof(split_line[6].c_str());
      }

      // RELAY_V2
      if (nstr > 7) {
        data.v2 = atof(split_line[7].c_str());
      }

      // RELAY_T2
      if (nstr > 8) {
        data.t2 = atof(split_line[8].c_str());
      }

      // RELAY_F2
      if (nstr > 9) {
        data.f2 = atof(split_line[9].c_str());
      }

      // RELAY_V3
      if (nstr > 10) {
        data.v3 = atof(split_line[10].c_str());
      }

      // RELAY_T3
      if (nstr > 11) {
        data.t3 = atof(split_line[11].c_str());
      }

      // RELAY_F3
      if (nstr > 12) {
        data.f3 = atof(split_line[12].c_str());
      }

      // RELAY_TB
      if (nstr > 13) {
        data.tb = atof(split_line[13].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
