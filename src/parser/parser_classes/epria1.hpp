/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: February, 7
 *      Author: Shrirang Abhyankar
 */
#ifndef EPRIA1_HPP
#define EPRIA1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Epria1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Epria1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Epria1Parser()
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
      // GENERATOR_MODEL              "MODEL"        string
      std::string stmp;
      if (!data->getValue(GENERATOR_MODEL,&stmp,g_id)) {
        data->addValue(GENERATOR_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(GENERATOR_MODEL, data_struct.model, g_id);
      }
	  
      int ival;

      if (!data->getValue(EPRIA1_PARAM1,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM1, data_struct.epria1_param1, g_id);
      } else {
        data->setValue(EPRIA1_PARAM1, data_struct.epria1_param1, g_id);
      }

      if (!data->getValue(EPRIA1_PARAM2,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM2, data_struct.epria1_param2, g_id);
      } else {
        data->setValue(EPRIA1_PARAM2, data_struct.epria1_param2, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM3,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM3, data_struct.epria1_param3, g_id);
      } else {
        data->setValue(EPRIA1_PARAM3, data_struct.epria1_param3, g_id);
      }

      if (!data->getValue(EPRIA1_PARAM4,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM4, data_struct.epria1_param4, g_id);
      } else {
        data->setValue(EPRIA1_PARAM4, data_struct.epria1_param4, g_id);
      }

      if (!data->getValue(EPRIA1_PARAM5,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM5, data_struct.epria1_param5, g_id);
      } else {
        data->setValue(EPRIA1_PARAM5, data_struct.epria1_param5, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM6,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM6, data_struct.epria1_param6, g_id);
      } else {
        data->setValue(EPRIA1_PARAM6, data_struct.epria1_param6, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM7,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM7, data_struct.epria1_param7, g_id);
      } else {
        data->setValue(EPRIA1_PARAM7, data_struct.epria1_param7, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM8,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM8, data_struct.epria1_param8, g_id);
      } else {
        data->setValue(EPRIA1_PARAM8, data_struct.epria1_param8, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM9,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM9, data_struct.epria1_param9, g_id);
      } else {
        data->setValue(EPRIA1_PARAM9, data_struct.epria1_param9, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM10,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM10, data_struct.epria1_param10, g_id);
      } else {
        data->setValue(EPRIA1_PARAM10, data_struct.epria1_param10, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM11,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM11, data_struct.epria1_param11, g_id);
      } else {
        data->setValue(EPRIA1_PARAM11, data_struct.epria1_param11, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM12,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM12, data_struct.epria1_param12, g_id);
      } else {
        data->setValue(EPRIA1_PARAM12, data_struct.epria1_param12, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM13,&ival,g_id)) {
        data->addValue(EPRIA1_PARAM13, data_struct.epria1_param13, g_id);
      } else {
        data->setValue(EPRIA1_PARAM13, data_struct.epria1_param13, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM14,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM14, data_struct.epria1_param14, g_id);
      } else {
        data->setValue(EPRIA1_PARAM14, data_struct.epria1_param14, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM15,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM15, data_struct.epria1_param15, g_id);
      } else {
        data->setValue(EPRIA1_PARAM15, data_struct.epria1_param15, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM16,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM16, data_struct.epria1_param16, g_id);
      } else {
        data->setValue(EPRIA1_PARAM16, data_struct.epria1_param16, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM17,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM17, data_struct.epria1_param17, g_id);
      } else {
        data->setValue(EPRIA1_PARAM17, data_struct.epria1_param17, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM18,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM18, data_struct.epria1_param18, g_id);
      } else {
        data->setValue(EPRIA1_PARAM18, data_struct.epria1_param18, g_id);
      }
      
      if (!data->getValue(EPRIA1_PARAM19,&rval,g_id)) {
        data->addValue(EPRIA1_PARAM19, data_struct.epria1_param19, g_id);
      } else {
        data->setValue(EPRIA1_PARAM19, data_struct.epria1_param19, g_id);
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
      // GENERATOR_MODEL              "MODEL"                  string
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(GENERATOR_MODEL, &stmp, g_id)) {
        data->addValue(GENERATOR_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(GENERATOR_MODEL, model.c_str(), g_id);
      }

      int ival;
      int ctr=3;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM1,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM1, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM1, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM2,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM2, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM2, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM3,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM3, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM3, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM4,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM4, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM4, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM5,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM5, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM5, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM6,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM6, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM6, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM7,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM7, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM7, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM8,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM8, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM8, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM9,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM9, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM9, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM10,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM10, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM10, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM11,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM11, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM11, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM12,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM12, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM12, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM13,&ival,g_id)) {
          data->addValue(EPRIA1_PARAM13, atoi(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM13, atoi(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM14,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM14, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM14, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM15,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM15, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM15, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM16,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM16, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM16, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM17,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM17, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM17, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM18,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM18, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM18, atof(split_line[ctr].c_str()), g_id);
        }
      }
      ctr += 1;
      if (nstr > ctr) {
        if (!data->getValue(EPRIA1_PARAM19,&rval,g_id)) {
          data->addValue(EPRIA1_PARAM19, atof(split_line[ctr].c_str()), g_id);
        } else {
          data->setValue(EPRIA1_PARAM19, atof(split_line[ctr].c_str()), g_id);
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
      // GENERATOR_BUSNUMBER               "I"                   integer
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

      // GENERATOR_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size();
      int ctr = 3;
      if (nstr > ctr) {
        data.epria1_param1 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param2 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param3 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param4 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param5 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param6 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param7 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param8 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param9 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param10 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param11 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param12 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param13 = atoi(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param14 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param15 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param16 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param17 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param18 = atof(split_line[ctr].c_str());
      }
      ctr += 1;
      if (nstr > ctr) {
        data.epria1_param19 = atof(split_line[ctr].c_str());
      }
      
    }	
};
}  // parser
}  // gridpack
#endif
