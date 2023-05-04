/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: Apr 03, 2023
 *      Author: Shuangshuang Jin
 */
#ifndef STAB1_HPP
#define STAB1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Stab1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Stab1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Stab1Parser()
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

      // HAS_PSS
      if (!data->getValue(HAS_PSS,&bval,g_id)) {
        data->addValue(HAS_PSS, true, g_id);
      } else {
        data->setValue(HAS_PSS, true, g_id);
      }

      // STAB1_MODEL
      std::string stmp;
      if (!data->getValue(PSS_MODEL, &stmp, g_id)) {
        data->addValue(PSS_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(PSS_MODEL, data_struct.model, g_id);
      }
        
        // STAB1_J
        if (!data->getValue(STAB1_J,&rval,g_id)) {
          data->addValue(STAB1_J, data_struct.stab1_j, g_id); 
        } else {
          data->setValue(STAB1_J, data_struct.stab1_j, g_id);
        }
      
      // STAB1_J1
      if (!data->getValue(STAB1_J1,&rval,g_id)) {
        data->addValue(STAB1_J1, data_struct.stab1_j1, g_id);
      } else {
        data->setValue(STAB1_J1, data_struct.stab1_j1, g_id);
      }

      // STAB1_J2
      if (!data->getValue(STAB1_J2,&rval,g_id)) {
        data->addValue(STAB1_J2, data_struct.stab1_j2, g_id);
      } else {
        data->setValue(STAB1_J2, data_struct.stab1_j2, g_id);
      }

      // STAB1_J3
      if (!data->getValue(STAB1_J3,&rval,g_id)) {
        data->addValue(STAB1_J3, data_struct.stab1_j3, g_id);
      } else {
        data->setValue(STAB1_J3, data_struct.stab1_j3, g_id);
      }

      // STAB1_J4
      if (!data->getValue(STAB1_J4,&rval,g_id)) {
        data->addValue(STAB1_J4, data_struct.stab1_j4, g_id);
      } else {
        data->setValue(STAB1_J4, data_struct.stab1_j4, g_id);
      }
      
      // STAB1_J5
      if (!data->getValue(STAB1_J5,&rval,g_id)) {
        data->addValue(STAB1_J5, data_struct.stab1_j5, g_id);
      } else {
        data->setValue(STAB1_J5, data_struct.stab1_j5, g_id);
      }
      
      // STAB1_J6
      if (!data->getValue(STAB1_J6,&rval,g_id)) {
        data->addValue(STAB1_J6, data_struct.stab1_j6, g_id);
      } else {
        data->setValue(STAB1_J6, data_struct.stab1_j6, g_id);
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
      // HAS_PSS
      if (!data->getValue(HAS_PSS,&bval,g_id)) {
        data->addValue(HAS_PSS, true, g_id);
      } else {
        data->setValue(HAS_PSS, true, g_id);
      }

      // STAB1_MODEL
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(PSS_MODEL,&stmp,g_id)) {
        data->addValue(PSS_MODEL, model.c_str(), g_id);
      } else {
        data->setValue(PSS_MODEL, model.c_str(), g_id);
      }
        
        // STAB1_J
        if (nstr > 3) {
          if (!data->getValue(STAB1_J,&rval,g_id)) {
            data->addValue(STAB1_J,
                atof(split_line[3].c_str()), g_id); 
            data->setValue(STAB1_J,
                atof(split_line[3].c_str()), g_id);
          }
        }

      // STAB1_J1
      if (nstr > 4) {
        if (!data->getValue(STAB1_J1,&rval,g_id)) {
          data->addValue(STAB1_J1,
              atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(STAB1_J1,
              atof(split_line[4].c_str()), g_id);
        }
      } 

      // STAB1_J2
      if (nstr > 5) {
        if (!data->getValue(STAB1_J2,&rval,g_id)) {
          data->addValue(STAB1_J2,
              atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(STAB1_J2,
              atof(split_line[5].c_str()), g_id);
        }
      } 

      // STAB1_J3
      if (nstr > 6) {
        if (!data->getValue(STAB1_J3,&rval,g_id)) {
          data->addValue(STAB1_J3,
              atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(STAB1_J3,
              atof(split_line[6].c_str()), g_id);
        }
      } 

      // STAB1_J4
      if (nstr > 7) {
        if (!data->getValue(STAB1_J4,&rval,g_id)) {
          data->addValue(STAB1_J4,
              atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(STAB1_J4,
              atof(split_line[7].c_str()), g_id);
        }
      } 

      // STAB1_J5
      if (nstr > 8) {
        if (!data->getValue(STAB1_J5,&rval,g_id)) {
          data->addValue(STAB1_J5,
              atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(STAB1_J5,
              atof(split_line[8].c_str()), g_id);
        }
      } 

      // STAB1_J6
      if (nstr > 9) {
        if (!data->getValue(STAB1_J6,&rval,g_id)) {
          data->addValue(STAB1_J6,
              atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(STAB1_J6,
              atof(split_line[9].c_str()), g_id);
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
      // STAB1_JNUMBER               "I"                   integer
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

      // STAB1_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());

      int nstr = split_line.size(); 
    
        // STAB1_J
        if (nstr > 3) {
          data.stab1_j= atof(split_line[3].c_str()); 
        }

      // STAB1_J1
      if (nstr > 4) {
        data.stab1_j1= atof(split_line[4].c_str()); 
      }

      // STAB1_J2
      if (nstr > 5) {
        data.stab1_j2 = atof(split_line[5].c_str()); 
      }

      // STAB1_J3
      if (nstr > 6) {
        data.stab1_j3 = atof(split_line[6].c_str()); 
      }

      // STAB1_J4
      if (nstr > 7) {
        data.stab1_j4 = atof(split_line[7].c_str()); 
      }

      // STAB1_J5
      if (nstr > 8) {
        data.stab1_j5 = atof(split_line[8].c_str()); 
      }

      // STAB1_J6
      if (nstr > 9) {
        data.stab1_j6 = atof(split_line[9].c_str()); 
      }

    }
};
}  // parser
}  // gridpack
#endif
