/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: June 16, 2016
 *      Author: Bruce Palmer
*
*  Created on: March 24, 2026
*       Author: Yuan Liu
 */
#ifndef REGFMA1_HPP
#define REGFMA1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Regfma1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Regfma1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Regfma1Parser()
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
      int ival;
      // GENERATOR_MODEL              "MODEL"        string
      std::string stmp;
      if (!data->getValue(GENERATOR_MODEL,&stmp,g_id)) {
        data->addValue(GENERATOR_MODEL, data_struct.model, g_id);
      } else {
        data->setValue(GENERATOR_MODEL, data_struct.model, g_id);
      }

      // GENERATOR_REGFMA1_VFLAG, integer
      if (!data->getValue(GENERATOR_REGFMA1_VFLAG,&ival,g_id)) {
        data->addValue(GENERATOR_REGFMA1_VFLAG, data_struct.vflag, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_VFLAG, data_struct.vflag, g_id);
      }

      // GENERATOR_REGFMA1_QVFLAG, integer
      if (!data->getValue(GENERATOR_REGFMA1_QVFLAG,&ival,g_id)) {
        data->addValue(GENERATOR_REGFMA1_QVFLAG, data_struct.qvflag, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_QVFLAG, data_struct.qvflag, g_id);
      }

      // GENERATOR_REGFMA1_TPF, float
      if (!data->getValue(GENERATOR_REGFMA1_TPF,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_TPF, data_struct.tpf, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_TPF, data_struct.tpf, g_id);
      }

      // GENERATOR_REGFMA1_TQF, float
      if (!data->getValue(GENERATOR_REGFMA1_TQF,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_TQF, data_struct.tqf, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_TQF, data_struct.tqf, g_id);
      }

      // GENERATOR_REGFMA1_TVF, float
      if (!data->getValue(GENERATOR_REGFMA1_TVF,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_TVF, data_struct.tvf, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_TVF, data_struct.tvf, g_id);
      }

      // GENERATOR_REGFMA1_IMAX, float
      if (!data->getValue(GENERATOR_REGFMA1_IMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_IMAX, data_struct.imax, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_IMAX, data_struct.imax, g_id);
      }

      // GENERATOR_REGFMA1_EMAX, float
      if (!data->getValue(GENERATOR_REGFMA1_EMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_EMAX, data_struct.emax, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_EMAX, data_struct.emax, g_id);
      }

      // GENERATOR_REGFMA1_EMIN, float
      if (!data->getValue(GENERATOR_REGFMA1_EMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_EMIN, data_struct.emin, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_EMIN, data_struct.emin, g_id);
      }

      // GENERATOR_REGFMA1_PMAX, float
      if (!data->getValue(GENERATOR_REGFMA1_PMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_PMAX, data_struct.pmax, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_PMAX, data_struct.pmax, g_id);
      }

      // GENERATOR_REGFMA1_PMIN, float
      if (!data->getValue(GENERATOR_REGFMA1_PMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_PMIN, data_struct.pmin, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_PMIN, data_struct.pmin, g_id);
      }

      // GENERATOR_REGFMA1_QMAX, float
      if (!data->getValue(GENERATOR_REGFMA1_QMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_QMAX, data_struct.qmax, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_QMAX, data_struct.qmax, g_id);
      }

      // GENERATOR_REGFMA1_QMIN, float
      if (!data->getValue(GENERATOR_REGFMA1_QMIN,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_QMIN, data_struct.qmin, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_QMIN, data_struct.qmin, g_id);
      }


      // GENERATOR_REGFMA1_MP, float
      if (!data->getValue(GENERATOR_REGFMA1_MP,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_MP, data_struct.mp, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_MP, data_struct.mp, g_id);
      }

      // GENERATOR_REGFMA1_MQ, float
      if (!data->getValue(GENERATOR_REGFMA1_MQ,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_MQ, data_struct.mq, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_MQ, data_struct.mq, g_id);
      }

      // GENERATOR_REGFMA1_KPPMAX, float
      if (!data->getValue(GENERATOR_REGFMA1_KPPMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_KPPMAX, data_struct.kppmax, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_KPPMAX, data_struct.kppmax, g_id);
      }

      // GENERATOR_REGFMA1_KIPMAX, float
      if (!data->getValue(GENERATOR_REGFMA1_KIPMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_KIPMAX, data_struct.kipmax, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_KIPMAX, data_struct.kipmax, g_id);
      }

      // GENERATOR_REGFMA1_KPQMAX, float
      if (!data->getValue(GENERATOR_REGFMA1_KPQMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_KPQMAX, data_struct.kpqmax, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_KPQMAX, data_struct.kpqmax, g_id);
      }

      // GENERATOR_REGFMA1_KIQMAX, float
      if (!data->getValue(GENERATOR_REGFMA1_KIQMAX,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_KIQMAX, data_struct.kiqmax, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_KIQMAX, data_struct.kiqmax, g_id);
      }

      // GENERATOR_REGFMA1_KPV, float
      if (!data->getValue(GENERATOR_REGFMA1_KPV,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_KPV, data_struct.kpv, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_KPV, data_struct.kpv, g_id);
      }

      // GENERATOR_REGFMA1_KIV, float
      if (!data->getValue(GENERATOR_REGFMA1_KIV,&rval,g_id)) {
        data->addValue(GENERATOR_REGFMA1_KIV, data_struct.kiv, g_id);
      } else {
        data->setValue(GENERATOR_REGFMA1_KIV, data_struct.kiv, g_id);
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
      int ival;
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

      // GENERATOR_REGFMA1_VFLAG, integer
      if (nstr > 3) {
        if (!data->getValue(GENERATOR_REGFMA1_VFLAG,&ival,g_id)) {
          data->addValue(GENERATOR_REGFMA1_VFLAG, atoi(split_line[3].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_VFLAG, atoi(split_line[3].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_QVFLAG, integer
      if (nstr > 4) {
        if (!data->getValue(GENERATOR_REGFMA1_QVFLAG,&ival,g_id)) {
          data->addValue(GENERATOR_REGFMA1_QVFLAG, atoi(split_line[4].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_QVFLAG, atoi(split_line[4].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_TPF, float
      if (nstr > 5) {
        if (!data->getValue(GENERATOR_REGFMA1_TPF,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_TPF, atof(split_line[5].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_TPF, atof(split_line[5].c_str()), g_id);
        }
      }

      // GENERATOR_REGFMA1_TQF, float
      if (nstr > 6) {
        if (!data->getValue(GENERATOR_REGFMA1_TQF,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_TQF, atof(split_line[6].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_TQF, atof(split_line[6].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_TVF, float
      if (nstr > 7) {
        if (!data->getValue(GENERATOR_REGFMA1_TVF,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_TVF, atof(split_line[7].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_TVF, atof(split_line[7].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_IMAX, float
      if (nstr > 8) {
        if (!data->getValue(GENERATOR_REGFMA1_IMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_IMAX, atof(split_line[8].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_IMAX, atof(split_line[8].c_str()), g_id);
        }
      }


      // GENERATOR_REGFMA1_EMAX, float
      if (nstr > 9) {
        if (!data->getValue(GENERATOR_REGFMA1_EMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_EMAX, atof(split_line[9].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_EMAX, atof(split_line[9].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_EMIN, float
      if (nstr > 10) {
        if (!data->getValue(GENERATOR_REGFMA1_EMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_EMIN, atof(split_line[10].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_EMIN, atof(split_line[10].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_PMAX, float
      if (nstr > 11) {
        if (!data->getValue(GENERATOR_REGFMA1_PMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_PMAX, atof(split_line[11].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_PMAX, atof(split_line[11].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_PMIN, float
      if (nstr > 12) {
        if (!data->getValue(GENERATOR_REGFMA1_PMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_PMIN, atof(split_line[12].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_PMIN, atof(split_line[12].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_QMAX, float
      if (nstr > 13) {
        if (!data->getValue(GENERATOR_REGFMA1_QMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_QMAX, atof(split_line[13].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_QMAX, atof(split_line[13].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_QMIN, float
      if (nstr > 14) {
        if (!data->getValue(GENERATOR_REGFMA1_QMIN,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_QMIN, atof(split_line[14].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_QMIN, atof(split_line[14].c_str()), g_id);
        }
      }

      // GENERATOR_REGFMA1_MP, float
      if (nstr > 15) {
        if (!data->getValue(GENERATOR_REGFMA1_MP,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_MP, atof(split_line[15].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_MP, atof(split_line[15].c_str()), g_id);
        }
      }

      // GENERATOR_REGFMA1_MQ, float
      if (nstr > 16) {
        if (!data->getValue(GENERATOR_REGFMA1_MQ,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_MQ, atof(split_line[16].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_MQ, atof(split_line[16].c_str()), g_id);
        }
      }

      // GENERATOR_REGFMA1_KPPMAX, float
      if (nstr > 17) {
        if (!data->getValue(GENERATOR_REGFMA1_KPPMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_KPPMAX, atof(split_line[17].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_KPPMAX, atof(split_line[17].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_KIPMAX, float
      if (nstr > 18) {
        if (!data->getValue(GENERATOR_REGFMA1_KIPMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_KIPMAX, atof(split_line[18].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_KIPMAX, atof(split_line[18].c_str()), g_id);
        }
      }

      // GENERATOR_REGFMA1_KPQMAX, float
      if (nstr > 19) {
        if (!data->getValue(GENERATOR_REGFMA1_KPQMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_KPQMAX, atof(split_line[19].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_KPQMAX, atof(split_line[19].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_KIQMAX, float
      if (nstr > 20) {
        if (!data->getValue(GENERATOR_REGFMA1_KIQMAX,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_KIQMAX, atof(split_line[20].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_KIQMAX, atof(split_line[20].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_KPV, float
      if (nstr > 21) {
        if (!data->getValue(GENERATOR_REGFMA1_KPV,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_KPV, atof(split_line[21].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_KPV, atof(split_line[21].c_str()), g_id);
        }
      } 

      // GENERATOR_REGFMA1_KIV, float
      if (nstr > 22) {
        if (!data->getValue(GENERATOR_REGFMA1_KIV,&rval,g_id)) {
          data->addValue(GENERATOR_REGFMA1_KIV, atof(split_line[22].c_str()), g_id);
        } else {
          data->setValue(GENERATOR_REGFMA1_KIV, atof(split_line[22].c_str()), g_id);
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

      // GENERATOR_REGFMA1_VFLAG, integer
      if (nstr > 3) {
        data.vflag = atoi(split_line[3].c_str());
      } 

      // GENERATOR_REGFMA1_QVFLAG, integer
      if (nstr > 4) {
        data.qvflag = atoi(split_line[4].c_str());
      } 

      // GENERATOR_REGFMA1_TPF, float
      if (nstr > 5) {
        data.tpf = atof(split_line[5].c_str());
      }

      // GENERATOR_REGFMA1_TQF, float
      if (nstr > 6) {
        data.tqf = atof(split_line[6].c_str());
      } 

      // GENERATOR_REGFMA1_TVF, float
      if (nstr > 7) {
        data.tvf = atof(split_line[7].c_str());
      } 

      // GENERATOR_REGFMA1_IMAX, float
      if (nstr > 8) {
        data.imax = atof(split_line[8].c_str());
      } 

      // GENERATOR_REGFMA1_EMAX, float
      if (nstr > 9) {
        data.emax = atof(split_line[9].c_str());
      } 

      // GENERATOR_REGFMA1_EMIN, float
      if (nstr > 10) {
        data.emin = atof(split_line[10].c_str());
      } 


      // GENERATOR_REGFMA1_PMAX, float
      if (nstr > 11) {
        data.pmax = atof(split_line[11].c_str());
      } 

      // GENERATOR_REGFMA1_PMIN, float
      if (nstr > 12) {
        data.pmin = atof(split_line[12].c_str());
      } 

      // GENERATOR_REGFMA1_QMAX, float
      if (nstr > 13) {
        data.qmax = atof(split_line[13].c_str());
      } 

      // GENERATOR_REGFMA1_QMIN, float
      if (nstr > 14) {
        data.qmin = atof(split_line[14].c_str());
      }

      // GENERATOR_REGFMA1_MP, float
      if (nstr > 15) {
        data.mp = atof(split_line[15].c_str());
      } 

      // GENERATOR_REGFMA1_MQ, float
      if (nstr > 16) {
        data.mq = atof(split_line[16].c_str());
      }

      // GENERATOR_REGFMA1_KPPMAX, float
      if (nstr > 17) {
        data.kppmax = atof(split_line[17].c_str());
      } 

      // GENERATOR_REGFMA1_KIPMAX, float
      if (nstr > 18) {
        data.kipmax = atof(split_line[18].c_str());
      }

      // GENERATOR_REGFMA1_KPQMAX, float
      if (nstr > 19) {
        data.kpqmax = atof(split_line[19].c_str());
      } 

      // GENERATOR_REGFMA1_KIQMAX, float
      if (nstr > 20) {
        data.kiqmax = atof(split_line[20].c_str());
      }

      // GENERATOR_REGFMA1_KPV, float
      if (nstr > 21) {
        data.kpv = atof(split_line[21].c_str());
      } 

      // GENERATOR_REGFMA1_KIV, float
      if (nstr > 22) {
        data.kiv = atof(split_line[22].c_str());
      } 
    }
};
}  // parser
}  // gridpack
#endif
