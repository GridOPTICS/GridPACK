/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: June 16, 2016
 *      Author: Bruce Palmer
 */
#ifndef GENROU_HPP
#define GENROU_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class GenrouParser
{
  public:
    /**
     * Constructor
     */
    explicit GenrouParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~GenrouParser()
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
      if (!data->getValue(GENERATOR_MODEL, &stmp, g_id)) {
        data->addValue(GENERATOR_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(GENERATOR_MODEL, data_struct.model, g_id);
      }

      // GENERATOR_TDOP
      if (!data->getValue(GENERATOR_TDOP,&rval,g_id)) {
        data->addValue(GENERATOR_TDOP,data_struct.tdop, g_id);
      } else {
        data->setValue(GENERATOR_TDOP, data_struct.tdop, g_id);
      }

      // GENERATOR_TDOPP
      if (!data->getValue(GENERATOR_TDOPP,&rval,g_id)) {
        data->addValue(GENERATOR_TDOPP, data_struct.tdopp, g_id);
      } else {
        data->setValue(GENERATOR_TDOPP, data_struct.tdopp, g_id);
      }

      // GENERATOR_TQOPP
      if (!data->getValue(GENERATOR_TQOPP,&rval,g_id)) {
        data->addValue(GENERATOR_TQOPP,
            data_struct.tqopp, g_id);
      } else {
        data->setValue(GENERATOR_TQOPP, data_struct.tqopp, g_id);
      }

      // GENERATOR_INERTIA_CONSTANT_H                           float
      if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H,&rval,g_id)) {
        data->addValue(GENERATOR_INERTIA_CONSTANT_H,
            data_struct.inertia, g_id);
      } else {
        data->setValue(GENERATOR_INERTIA_CONSTANT_H,
            data_struct.inertia, g_id);
      }

      // GENERATOR_DAMPING_COEFFICIENT_0                           float
      if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&rval,g_id)) {
        data->addValue(GENERATOR_DAMPING_COEFFICIENT_0,
            data_struct.damping, g_id);
      } else {
        data->setValue(GENERATOR_DAMPING_COEFFICIENT_0,
            data_struct.damping, g_id);
      }

      // GENERATOR_XD
      if (!data->getValue(GENERATOR_XD,&rval,g_id)) {
        data->addValue(GENERATOR_XD, data_struct.gn_xd, g_id);
      } else {
        data->setValue(GENERATOR_XD, data_struct.gn_xd, g_id);
      }

      // GENERATOR_XQ
      if (!data->getValue(GENERATOR_XQ,&rval,g_id)) {
        data->addValue(GENERATOR_XQ, data_struct.gn_xq, g_id);
      } else {
        data->setValue(GENERATOR_XQ, data_struct.gn_xq, g_id);
      }

      // GENERATOR_XDP
      if (!data->getValue(GENERATOR_XDP,&rval,g_id)) {
        data->addValue(GENERATOR_XDP, data_struct.xdp, g_id);
      } else {
        data->setValue(GENERATOR_XDP, data_struct.xdp, g_id);
      }

      // GENERATOR_XDPP
      if (!data->getValue(GENERATOR_XDPP,&rval,g_id)) {
        data->addValue(GENERATOR_XDPP, data_struct.xdpp, g_id);
      } else {
        data->setValue(GENERATOR_XDPP, data_struct.xdpp, g_id);
      }

      // GENERATOR_XL
      if (!data->getValue(GENERATOR_XL,&rval,g_id)) {
        data->addValue(GENERATOR_XL,
            data_struct.gn_xl, g_id);
      } else {
        data->setValue(GENERATOR_XL,
            data_struct.gn_xl, g_id);
      }

      // GENERATOR_S1
      if (!data->getValue(GENERATOR_S1,&rval,g_id)) {
        data->addValue(GENERATOR_XL, data_struct.gn_s1, g_id);
      } else {
        data->setValue(GENERATOR_S1, data_struct.gn_s1, g_id);
      }

      // GENERATOR_S12
      if (!data->getValue(GENERATOR_S12,&rval,g_id)) {
        data->addValue(GENERATOR_XL, data_struct.s12, g_id);
      } else {
        data->setValue(GENERATOR_S12, data_struct.s12, g_id);
      }

      // GENERATOR_TQOP
      if (!data->getValue(GENERATOR_TQOP,&rval,g_id)) {
        data->addValue(GENERATOR_TQOP,
            data_struct.tqop, g_id);
      } else {
        data->setValue(GENERATOR_TQOP, data_struct.tqop, g_id);
      }

      // GENERATOR_XQP
      if (!data->getValue(GENERATOR_XQP,&rval,g_id)) {
        data->addValue(GENERATOR_XQP,
            data_struct.xqp, g_id);
      } else {
        data->setValue(GENERATOR_XQP, data_struct.xqp, g_id);
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

      // GENERATOR_TDOP
      if (nstr > 3) {
        if (!data->getValue(GENERATOR_TDOP,&rval,g_id)) {
          data->addValue(GENERATOR_TDOP,
              atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_TDOP,
              atof(split_line[3].c_str()), g_id);
        }
      } 

      // GENERATOR_TDOPP
      if (nstr > 4) {
        if (!data->getValue(GENERATOR_TDOPP,&rval,g_id)) {
          data->addValue(GENERATOR_TDOPP,
              atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_TDOPP,
              atof(split_line[4].c_str()), g_id);
        }
      } 

      // GENERATOR_TQOP
      if (nstr > 5) {
        if (!data->getValue(GENERATOR_TQOP,&rval,g_id)) {
          data->addValue(GENERATOR_TQOP,
              atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_TQOP,
              atof(split_line[5].c_str()), g_id);
        }
      } 

      // GENERATOR_TQOPP
      if (nstr > 6) {
        if (!data->getValue(GENERATOR_TQOPP,&rval,g_id)) {
          data->addValue(GENERATOR_TQOPP,
              atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_TQOPP,
              atof(split_line[6].c_str()), g_id);
        }
      } 

      // GENERATOR_INERTIA_CONSTANT_H                           float
      if (nstr > 7) {
        if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H,&rval,g_id)) {
          data->addValue(GENERATOR_INERTIA_CONSTANT_H,
              atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_INERTIA_CONSTANT_H,
              atof(split_line[7].c_str()), g_id);
        }
      } 

      // GENERATOR_DAMPING_COEFFICIENT_0                           float
      if (nstr > 8) {
        if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&rval,g_id)) {
          data->addValue(GENERATOR_DAMPING_COEFFICIENT_0,
              atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_DAMPING_COEFFICIENT_0,
              atof(split_line[8].c_str()), g_id);
        }
      }

      // GENERATOR_XD
      if (nstr > 9) {
        if (!data->getValue(GENERATOR_XD,&rval,g_id)) {
          data->addValue(GENERATOR_XD,
              atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_XD,
              atof(split_line[9].c_str()), g_id);
        }
      } 

      // GENERATOR_XQ
      if (nstr > 10) {
        if (!data->getValue(GENERATOR_XQ,&rval,g_id)) {
          data->addValue(GENERATOR_XQ,
              atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_XQ,
              atof(split_line[9].c_str()), g_id);
        }
      } 

      // GENERATOR_XDP
      if (nstr > 11) {
        if (!data->getValue(GENERATOR_XDP,&rval,g_id)) {
          data->addValue(GENERATOR_XDP,
              atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_XDP,
              atof(split_line[11].c_str()), g_id);
        }
      } 

      // GENERATOR_XQP
      if (nstr > 12) {
        if (!data->getValue(GENERATOR_XQP,&rval,g_id)) {
          data->addValue(GENERATOR_XQP,
              atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_XQP,
              atof(split_line[12].c_str()), g_id);
        }
      } 

      // GENERATOR_XDPP
      if (nstr > 13) {
        if (!data->getValue(GENERATOR_XDPP,&rval,g_id)) {
          data->addValue(GENERATOR_XDPP,
              atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_XDPP,
              atof(split_line[13].c_str()), g_id);
        }
      } 

      // GENERATOR_XL
      if (nstr > 14) {
        if (!data->getValue(GENERATOR_XL,&rval,g_id)) {
          data->addValue(GENERATOR_XL,
              atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_XL,
              atof(split_line[14].c_str()), g_id);
        }
      } 

      // GENERATOR_S1
      if (nstr > 15) {
        if (!data->getValue(GENERATOR_S1,&rval,g_id)) {
          data->addValue(GENERATOR_XL,
              atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_S1,
              atof(split_line[15].c_str()), g_id);
        }
      } 

      // GENERATOR_S12
      if (nstr > 16) {
        if (!data->getValue(GENERATOR_S12,&rval,g_id)) {
          data->addValue(GENERATOR_XL,
              atof(split_line[16].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_S12,
              atof(split_line[16].c_str()), g_id);
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
      // GENERATOR_TDOP
      if (nstr > 3) {
        data.tdop = atof(split_line[3].c_str());
      } 

      // GENERATOR_TDOPP
      if (nstr > 4) {
        data.tdopp = atof(split_line[4].c_str());
      } 

      // GENERATOR_TQOP
      if (nstr > 5) {
        data.tqop = atof(split_line[5].c_str());
      } 

      // GENERATOR_TQOPP
      if (nstr > 6) {
        data.tqopp = atof(split_line[6].c_str());
      } 

      // GENERATOR_INERTIA_CONSTANT_H                           float
      if (nstr > 7) {
        data.inertia = atof(split_line[7].c_str());
      } 

      // GENERATOR_DAMPING_COEFFICIENT_0                           float
      if (nstr > 8) {
        data.damping = atof(split_line[8].c_str());
      }

      // GENERATOR_XD
      if (nstr > 9) {
        data.gn_xd = atof(split_line[9].c_str());
      } 

      // GENERATOR_XQ
      if (nstr > 10) {
        data.gn_xq = atof(split_line[10].c_str());
      } 

      // GENERATOR_XDP
      if (nstr > 11) {
        data.xdp = atof(split_line[11].c_str());
      } 

      // GENERATOR_XQP
      if (nstr > 12) {
        data.xqp = atof(split_line[12].c_str());
      } 

      // GENERATOR_XDPP
      if (nstr > 13) {
        data.xdpp = atof(split_line[13].c_str());
      } 

      // GENERATOR_XL
      if (nstr > 14) {
        data.gn_xl = atof(split_line[14].c_str());
      } 

      // GENERATOR_S1
      if (nstr > 15) {
        data.gn_s1 = atof(split_line[15].c_str());
      } 

      // GENERATOR_S12
      if (nstr > 16) {
        data.s12 = atof(split_line[16].c_str());
      } 
    }
};
}  // parser
}  // gridpack
#endif
