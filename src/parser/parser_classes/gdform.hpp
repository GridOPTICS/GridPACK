/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: June 16, 2016
 *      Author: Bruce Palmer
 */
#ifndef GDFORM_HPP
#define GDFORM_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class GdformParser
{
  public:
    /**
     * Constructor
     */
    explicit GdformParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~GdformParser()
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

      // GENERATOR_XL
      if (!data->getValue(GENERATOR_XL,&rval,g_id)) {
        data->addValue(GENERATOR_XL,data_struct.gn_xl, g_id);
      } else {
        data->setValue(GENERATOR_XL, data_struct.gn_xl, g_id);
      }

      // GENERATOR_VSET
      if (!data->getValue(GENERATOR_VSET,&rval,g_id)) {
        data->addValue(GENERATOR_VSET, data_struct.vset, g_id);
      } else {
        data->setValue(GENERATOR_VSET, data_struct.vset, g_id);
      }

      // GENERATOR_MQ
      if (!data->getValue(GENERATOR_MQ,&rval,g_id)) {
        data->addValue(GENERATOR_MQ, data_struct.mq, g_id);
      } else {
        data->setValue(GENERATOR_MQ, data_struct.mq, g_id);
      }

      // GENERATOR_KPV                           float
      if (!data->getValue(GENERATOR_KPV,&rval,g_id)) {
        data->addValue(GENERATOR_KPV, data_struct.kpv, g_id);
      } else {
        data->setValue(GENERATOR_KPV, data_struct.kpv, g_id);
      }

      // GENERATOR_KIV                           float
      if (!data->getValue(GENERATOR_KIV,&rval,g_id)) {
        data->addValue(GENERATOR_KIV, data_struct.kiv, g_id);
      } else {
        data->setValue(GENERATOR_KIV, data_struct.kiv, g_id);
      }

      // GENERATOR_EMAX
      if (!data->getValue(GENERATOR_EMAX,&rval,g_id)) {
        data->addValue(GENERATOR_EMAX, data_struct.emax, g_id);
      } else {
        data->setValue(GENERATOR_EMAX, data_struct.emax, g_id);
      }

      // GENERATOR_EMIN
      if (!data->getValue(GENERATOR_EMIN,&rval,g_id)) {
        data->addValue(GENERATOR_EMIN, data_struct.emin, g_id);
      } else {
        data->setValue(GENERATOR_EMIN, data_struct.emin, g_id);
      }

      // GENERATOR_MP
      if (!data->getValue(GENERATOR_MP,&rval,g_id)) {
        data->addValue(GENERATOR_MP, data_struct.mp, g_id);
      } else {
        data->setValue(GENERATOR_MP, data_struct.mp, g_id);
      }

      // GENERATOR_KPPMAX
      if (!data->getValue(GENERATOR_KPPMAX,&rval,g_id)) {
        data->addValue(GENERATOR_KPPMAX, data_struct.kppmax, g_id);
      } else {
        data->setValue(GENERATOR_KPPMAX, data_struct.kppmax, g_id);
      }

      // GENERATOR_KIPMAX
      if (!data->getValue(GENERATOR_KIPMAX,&rval,g_id)) {
        data->addValue(GENERATOR_KIPMAX, data_struct.kipmax, g_id);
      } else {
        data->setValue(GENERATOR_KIPMAX, data_struct.kipmax, g_id);
      }

      // GENERATOR_PSET
      if (!data->getValue(GENERATOR_PSET,&rval,g_id)) {
        data->addValue(GENERATOR_PSET, data_struct.pset, g_id);
      } else {
        data->setValue(GENERATOR_PSET, data_struct.pset, g_id);
      }

      // GENERATOR_S1
      if (!data->getValue(GENERATOR_PMAX,&rval,g_id)) {
        data->addValue(GENERATOR_PMAX, data_struct.pmax, g_id);
      } else {
        data->setValue(GENERATOR_PMAX, data_struct.pmax, g_id);
      }

      // GENERATOR_PMIN
      if (!data->getValue(GENERATOR_PMIN,&rval,g_id)) {
        data->addValue(GENERATOR_PMIN, data_struct.pmin, g_id);
      } else {
        data->setValue(GENERATOR_PMIN, data_struct.pmin, g_id);
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

      // GENERATOR_XL
      if (nstr > 3) {
        if (!data->getValue(GENERATOR_XL,&rval,g_id)) {
          data->addValue(GENERATOR_XL, atof(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_XL, atof(split_line[3].c_str()), g_id);
        }
      } 

      // GENERATOR_VSET
      if (nstr > 4) {
        if (!data->getValue(GENERATOR_VSET,&rval,g_id)) {
          data->addValue(GENERATOR_VSET, atof(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_VSET, atof(split_line[4].c_str()), g_id);
        }
      } 

      // GENERATOR_MQ
      if (nstr > 5) {
        if (!data->getValue(GENERATOR_MQ,&rval,g_id)) {
          data->addValue(GENERATOR_MQ, atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_MQ, atof(split_line[5].c_str()), g_id);
        }
      } 

      // GENERATOR_KPV                           float
      if (nstr > 6) {
        if (!data->getValue(GENERATOR_KPV,&rval,g_id)) {
          data->addValue(GENERATOR_KPV, atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_KPV, atof(split_line[6].c_str()), g_id);
        }
      } 

      // GENERATOR_KIV                           float
      if (nstr > 7) {
        if (!data->getValue(GENERATOR_KIV,&rval,g_id)) {
          data->addValue(GENERATOR_KIV, atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_KIV, atof(split_line[7].c_str()), g_id);
        }
      }

      // GENERATOR_EMAX
      if (nstr > 8) {
        if (!data->getValue(GENERATOR_EMAX,&rval,g_id)) {
          data->addValue(GENERATOR_EMAX, atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_EMAX, atof(split_line[8].c_str()), g_id);
        }
      } 

      // GENERATOR_EMIN
      if (nstr > 9) {
        if (!data->getValue(GENERATOR_EMIN,&rval,g_id)) {
          data->addValue(GENERATOR_EMIN, atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_EMIN, atof(split_line[9].c_str()), g_id);
        }
      } 

      // GENERATOR_MP
      if (nstr > 10) {
        if (!data->getValue(GENERATOR_MP,&rval,g_id)) {
          data->addValue(GENERATOR_MP, atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_MP, atof(split_line[10].c_str()), g_id);
        }
      } 

      // GENERATOR_KPPMAX
      if (nstr > 11) {
        if (!data->getValue(GENERATOR_KPPMAX,&rval,g_id)) {
          data->addValue(GENERATOR_KPPMAX, atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_KPPMAX, atof(split_line[11].c_str()), g_id);
        }
      } 

      // GENERATOR_KIPMAX
      if (nstr > 12) {
        if (!data->getValue(GENERATOR_KIPMAX,&rval,g_id)) {
          data->addValue(GENERATOR_XL, atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_XL, atof(split_line[12].c_str()), g_id);
        }
      } 

      // GENERATOR_PSET
      if (nstr > 13) {
        if (!data->getValue(GENERATOR_PSET,&rval,g_id)) {
          data->addValue(GENERATOR_PSET,
              atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_PSET,
              atof(split_line[13].c_str()), g_id);
        }
      } 

      // GENERATOR_PMAX
      if (nstr > 14) {
        if (!data->getValue(GENERATOR_PMAX,&rval,g_id)) {
          data->addValue(GENERATOR_PMAX, atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_PMAX, atof(split_line[14].c_str()), g_id);
        }
      } 

      // GENERATOR_PMIN
      if (nstr > 15) {
        if (!data->getValue(GENERATOR_PMIN,&rval,g_id)) {
          data->addValue(GENERATOR_PMIN, atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_PMIN, atof(split_line[15].c_str()), g_id);
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
      // GENERATOR_XL
      if (nstr > 3) {
        data.gn_xl = atof(split_line[3].c_str());
      } 

      // GENERATOR_VSET
      if (nstr > 4) {
        data.vset = atof(split_line[4].c_str());
      } 

      // GENERATOR_MQ
      if (nstr > 5) {
        data.mq = atof(split_line[5].c_str());
      } 

      // GENERATOR_KPV                           float
      if (nstr > 6) {
        data.kpv = atof(split_line[6].c_str());
      } 

      // GENERATOR_KIV                           float
      if (nstr > 7) {
        data.kiv = atof(split_line[7].c_str());
      }

      // GENERATOR_EMAX
      if (nstr > 8) {
        data.emax = atof(split_line[8].c_str());
      } 

      // GENERATOR_EMIN
      if (nstr > 9) {
        data.emin = atof(split_line[9].c_str());
      } 

      // GENERATOR_MP
      if (nstr > 10) {
        data.mp = atof(split_line[10].c_str());
      } 

      // GENERATOR_KPPMAX
      if (nstr > 11) {
        data.kppmax = atof(split_line[11].c_str());
      } 

      // GENERATOR_KIPMAX
      if (nstr > 12) {
        data.kipmax = atof(split_line[12].c_str());
      } 

      // GENERATOR_PSET
      if (nstr > 13) {
        data.pset = atof(split_line[13].c_str());
      } 

      // GENERATOR_PMAX
      if (nstr > 14) {
        data.pmax = atof(split_line[14].c_str());
      } 

      // GENERATOR_PMIN
      if (nstr > 15) {
        data.pmin = atof(split_line[15].c_str());
      } 
    }
};
}  // parser
}  // gridpack
#endif
