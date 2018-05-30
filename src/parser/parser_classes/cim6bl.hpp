/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: September 8, 2016
 *      Author: Bruce Palmer
 */
#ifndef CIM6BL_HPP
#define CIM6BL_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Cim6blParser
{
  public:
    /**
     * Constructor
     */
    explicit Cim6blParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Cim6blParser()
    {
    }

    /**
     * Extract data from _data_struct and store it in data collection object
     * @param data_struct data struct object
     * @param data data collection object
     * @param gen_id index of generator
     */
    void extract(_data_struct &data_struct,
        gridpack::component::DataCollection *data, int l_id)
    {
      double rval;
      int ival;
      std::string stmp;

      // LOAD_MODEL
      if (!data->getValue(LOAD_MODEL,&stmp,l_id)) {
        data->addValue(LOAD_MODEL, data_struct.model, l_id);
      } else {
        data->setValue(LOAD_MODEL, data_struct.model, l_id);
      }

      // LOAD_IT
      if (!data->getValue(LOAD_IT,&ival,l_id)) {
        data->addValue(LOAD_IT, data_struct.it, l_id);
      } else {
        data->setValue(LOAD_IT, data_struct.it, l_id);
      }

      // LOAD_RA
      if (!data->getValue(LOAD_RA,&rval,l_id)) {
        data->addValue(LOAD_RA, data_struct.ra, l_id);
      } else {
        data->setValue(LOAD_RA, data_struct.ra, l_id);
      }

      // LOAD_XA
      if (!data->getValue(LOAD_XA,&rval,l_id)) {
        data->addValue(LOAD_XA, data_struct.xa, l_id);
      } else {
        data->setValue(LOAD_XA, data_struct.xa, l_id);
      }

      // LOAD_XM
      if (!data->getValue(LOAD_XM,&rval,l_id)) {
        data->addValue(LOAD_XM, data_struct.xm, l_id);
      } else {
        data->setValue(LOAD_XM, data_struct.xm, l_id);
      }

      // LOAD_R1
      if (!data->getValue(LOAD_R1,&rval,l_id)) {
        data->addValue(LOAD_R1, data_struct.r1, l_id);
      } else {
        data->setValue(LOAD_R1, data_struct.r1, l_id);
      }

      // LOAD_X1
      if (!data->getValue(LOAD_X1,&rval,l_id)) {
        data->addValue(LOAD_X1, data_struct.r1, l_id);
      } else {
        data->setValue(LOAD_X1, data_struct.r1, l_id);
      }

      // LOAD_R2
      if (!data->getValue(LOAD_R2,&rval,l_id)) {
        data->addValue(LOAD_R2, data_struct.r2, l_id);
      } else {
        data->setValue(LOAD_R2, data_struct.r2, l_id);
      }

      // LOAD_X2
      if (!data->getValue(LOAD_X2,&rval,l_id)) {
        data->addValue(LOAD_X2, data_struct.x2, l_id);
      } else {
        data->setValue(LOAD_X2, data_struct.x2, l_id);
      }

      // LOAD_E1
      if (!data->getValue(LOAD_E1,&rval,l_id)) {
        data->addValue(LOAD_E1, data_struct.e1, l_id);
      } else {
        data->setValue(LOAD_E1, data_struct.e1, l_id);
      }

      // LOAD_SE1
      if (!data->getValue(LOAD_SE1,&rval,l_id)) {
        data->addValue(LOAD_SE1, data_struct.se1, l_id);
      } else {
        data->setValue(LOAD_SE1, data_struct.se1, l_id);
      }

      // LOAD_E2
      if (!data->getValue(LOAD_E2,&rval,l_id)) {
        data->addValue(LOAD_E2, data_struct.se2, l_id);
      } else {
        data->setValue(LOAD_E2, data_struct.se2, l_id);
      }

      // LOAD_SE2
      if (!data->getValue(LOAD_SE2,&rval,l_id)) {
        data->addValue(LOAD_SE2, data_struct.se2, l_id);
      } else {
        data->setValue(LOAD_SE2, data_struct.se2, l_id);
      }

      // LOAD_MBASE
      if (!data->getValue(LOAD_MBASE,&rval,l_id)) {
        data->addValue(LOAD_MBASE, data_struct.mbase, l_id);
      } else {
        data->setValue(LOAD_MBASE, data_struct.mbase, l_id);
      }

      // LOAD_PMULT
      if (!data->getValue(LOAD_PMULT,&rval,l_id)) {
        data->addValue(LOAD_PMULT, data_struct.pmult, l_id);
      } else {
        data->setValue(LOAD_PMULT, data_struct.pmult, l_id);
      }

      // LOAD_H
      if (!data->getValue(LOAD_H,&rval,l_id)) {
        data->addValue(LOAD_H, data_struct.h, l_id);
      } else {
        data->setValue(LOAD_H, data_struct.h, l_id);
      }

      // LOAD_VI
      if (!data->getValue(LOAD_VI,&rval,l_id)) {
        data->addValue(LOAD_VI, data_struct.vi, l_id);
      } else {
        data->setValue(LOAD_VI, data_struct.vi, l_id);
      }

      // LOAD_TI
      if (!data->getValue(LOAD_TI,&rval,l_id)) {
        data->addValue(LOAD_TI, data_struct.ti, l_id);
      } else {
        data->setValue(LOAD_TI, data_struct.ti, l_id);
      }

      // LOAD_TB
      if (!data->getValue(LOAD_TB,&rval,l_id)) {
        data->addValue(LOAD_TB, data_struct.tb, l_id);
      } else {
        data->setValue(LOAD_TB, data_struct.tb, l_id);
      }

      // LOAD_A
      if (!data->getValue(LOAD_A,&rval,l_id)) {
        data->addValue(LOAD_A, data_struct.a, l_id);
      } else {
        data->setValue(LOAD_A, data_struct.a, l_id);
      }

      // LOAD_B
      if (!data->getValue(LOAD_B,&rval,l_id)) {
        data->addValue(LOAD_B, data_struct.b, l_id);
      } else {
        data->setValue(LOAD_B, data_struct.b, l_id);
      }

      // LOAD_D
      if (!data->getValue(LOAD_D,&rval,l_id)) {
        data->addValue(LOAD_D, data_struct.d, l_id);
      } else {
        data->setValue(LOAD_D, data_struct.d, l_id);
      }

      // LOAD_E
      if (!data->getValue(LOAD_E,&rval,l_id)) {
        data->addValue(LOAD_E, data_struct.e, l_id);
      } else {
        data->setValue(LOAD_E, data_struct.e, l_id);
      }

      // LOAD_C0
      if (!data->getValue(LOAD_C0,&rval,l_id)) {
        data->addValue(LOAD_C0, data_struct.c0, l_id);
      } else {
        data->setValue(LOAD_C0, data_struct.c0, l_id);
      }

      // LOAD_TNOM
      if (!data->getValue(LOAD_TNOM,&rval,l_id)) {
        data->addValue(LOAD_TNOM, data_struct.tnom, l_id);
      } else {
        data->setValue(LOAD_TNOM, data_struct.tnom, l_id);
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
        gridpack::component::DataCollection *data, int l_id)
    {
      double rval;
      int nstr = split_line.size();
      int ival;

      // LOAD_MODEL
      std::string stmp, model;
      gridpack::utility::StringUtils util;
      model = util.trimQuotes(split_line[1]);
      util.toUpper(model);
      if (!data->getValue(LOAD_MODEL,&stmp,l_id)) {
        data->addValue(LOAD_MODEL, model.c_str(), l_id);
      } else {
        data->setValue(LOAD_MODEL, model.c_str(), l_id);
      }

      // LOAD_IT
      if (nstr > 3) {
        if (!data->getValue(LOAD_IT,&ival,l_id)) {
          data->addValue(LOAD_IT, atoi(split_line[3].c_str()), l_id);
        } else {
          data->setValue(LOAD_IT, atoi(split_line[3].c_str()), l_id);
        }
      }

      // LOAD_RA
      if (nstr > 4) {
        if (!data->getValue(LOAD_RA,&rval,l_id)) {
          data->addValue(LOAD_RA, atof(split_line[4].c_str()), l_id);
        } else {
          data->setValue(LOAD_RA, atof(split_line[4].c_str()), l_id);
        }
      } 

      // LOAD_XA
      if (nstr > 5) {
        if (!data->getValue(LOAD_XA,&rval,l_id)) {
          data->addValue(LOAD_XA, atof(split_line[5].c_str()), l_id);
        } else {
          data->setValue(LOAD_XA, atof(split_line[5].c_str()), l_id);
        }
      } 

      // LOAD_XM
      if (nstr > 6) {
        if (!data->getValue(LOAD_XM,&rval,l_id)) {
          data->addValue(LOAD_XM, atof(split_line[6].c_str()), l_id);
        } else {
          data->setValue(LOAD_XM, atof(split_line[6].c_str()), l_id);
        }
      } 

      // LOAD_R1
      if (nstr > 7) {
        if (!data->getValue(LOAD_R1,&rval,l_id)) {
          data->addValue(LOAD_R1, atof(split_line[7].c_str()), l_id);
        } else {
          data->setValue(LOAD_R1, atof(split_line[7].c_str()), l_id);
        }
      } 

      // LOAD_X1
      if (nstr > 8) {
        if (!data->getValue(LOAD_X1,&rval,l_id)) {
          data->addValue(LOAD_X1, atof(split_line[8].c_str()), l_id);
        } else {
          data->setValue(LOAD_X1, atof(split_line[8].c_str()), l_id);
        }
      } 

      // LOAD_R2
      if (nstr > 9) {
        if (!data->getValue(LOAD_R2,&rval,l_id)) {
          data->addValue(LOAD_R2, atof(split_line[9].c_str()), l_id);
        } else {
          data->setValue(LOAD_R2, atof(split_line[9].c_str()), l_id);
        }
      } 

      // LOAD_X2
      if (nstr > 10) {
        if (!data->getValue(LOAD_X2,&rval,l_id)) {
          data->addValue(LOAD_X2, atof(split_line[10].c_str()), l_id);
        } else {
          data->setValue(LOAD_X2, atof(split_line[10].c_str()), l_id);
        }
      } 

      // LOAD_E1
      if (nstr > 11) {
        if (!data->getValue(LOAD_E1,&rval,l_id)) {
          data->addValue(LOAD_E1, atof(split_line[11].c_str()), l_id);
        } else {
          data->setValue(LOAD_E1, atof(split_line[11].c_str()), l_id);
        }
      } 

      // LOAD_SE1
      if (nstr > 12) {
        if (!data->getValue(LOAD_SE1,&rval,l_id)) {
          data->addValue(LOAD_SE1, atof(split_line[12].c_str()), l_id);
        } else {
          data->setValue(LOAD_SE1, atof(split_line[12].c_str()), l_id);
        }
      } 

      // LOAD_E2
      if (nstr > 13) {
        if (!data->getValue(LOAD_E2,&rval,l_id)) {
          data->addValue(LOAD_E2, atof(split_line[13].c_str()), l_id);
        } else {
          data->setValue(LOAD_E2, atof(split_line[13].c_str()), l_id);
        }
      } 

      // LOAD_SE2
      if (nstr > 14) {
        if (!data->getValue(LOAD_SE2,&rval,l_id)) {
          data->addValue(LOAD_SE2, atof(split_line[14].c_str()), l_id);
        } else {
          data->setValue(LOAD_SE2, atof(split_line[14].c_str()), l_id);
        }
      } 

      // LOAD_MBASE
      if (nstr > 15) {
        if (!data->getValue(LOAD_MBASE,&rval,l_id)) {
          data->addValue(LOAD_MBASE, atof(split_line[15].c_str()), l_id);
        } else {
          data->setValue(LOAD_MBASE, atof(split_line[15].c_str()), l_id);
        }
      } 

      // LOAD_PMULT
      if (nstr > 16) {
        if (!data->getValue(LOAD_PMULT,&rval,l_id)) {
          data->addValue(LOAD_PMULT, atof(split_line[16].c_str()), l_id);
        } else {
          data->setValue(LOAD_PMULT, atof(split_line[16].c_str()), l_id);
        }
      } 

      // LOAD_H
      if (nstr > 17) {
        if (!data->getValue(LOAD_H,&rval,l_id)) {
          data->addValue(LOAD_H, atof(split_line[17].c_str()), l_id);
        } else {
          data->setValue(LOAD_H, atof(split_line[17].c_str()), l_id);
        }
      }

      // LOAD_VI
      if (nstr > 18) {
        if (!data->getValue(LOAD_VI,&rval,l_id)) {
          data->addValue(LOAD_VI, atof(split_line[18].c_str()), l_id);
        } else {
          data->setValue(LOAD_VI, atof(split_line[18].c_str()), l_id);
        }
      }

      // LOAD_TI
      if (nstr > 19) {
        if (!data->getValue(LOAD_TI,&rval,l_id)) {
          data->addValue(LOAD_TI, atof(split_line[19].c_str()), l_id);
        } else {
          data->setValue(LOAD_TI, atof(split_line[19].c_str()), l_id);
        }
      }

      // LOAD_TB
      if (nstr > 20) {
        if (!data->getValue(LOAD_TB,&rval,l_id)) {
          data->addValue(LOAD_TB, atof(split_line[20].c_str()), l_id);
        } else {
          data->setValue(LOAD_TB, atof(split_line[20].c_str()), l_id);
        }
      }

      // LOAD_A
      if (nstr > 21) {
        if (!data->getValue(LOAD_A,&rval,l_id)) {
          data->addValue(LOAD_A, atof(split_line[21].c_str()), l_id);
        } else {
          data->setValue(LOAD_A, atof(split_line[21].c_str()), l_id);
        }
      }

      // LOAD_B
      if (nstr > 22) {
        if (!data->getValue(LOAD_B,&rval,l_id)) {
          data->addValue(LOAD_B, atof(split_line[22].c_str()), l_id);
        } else {
          data->setValue(LOAD_B, atof(split_line[22].c_str()), l_id);
        }
      }

      // LOAD_D
      if (nstr > 23) {
        if (!data->getValue(LOAD_D,&rval,l_id)) {
          data->addValue(LOAD_D, atof(split_line[23].c_str()), l_id);
        } else {
          data->setValue(LOAD_D, atof(split_line[23].c_str()), l_id);
        }
      }

      // LOAD_E
      if (nstr > 24) {
        if (!data->getValue(LOAD_E,&rval,l_id)) {
          data->addValue(LOAD_E, atof(split_line[24].c_str()), l_id);
        } else {
          data->setValue(LOAD_E, atof(split_line[24].c_str()), l_id);
        }
      } 

      // LOAD_C0
      if (nstr > 25) {
        if (!data->getValue(LOAD_C0,&rval,l_id)) {
          data->addValue(LOAD_C0, atof(split_line[25].c_str()), l_id);
        } else {
          data->setValue(LOAD_C0, atof(split_line[25].c_str()), l_id);
        }
      } 

      // LOAD_TNOM
      if (nstr > 26) {
        if (!data->getValue(LOAD_TNOM,&rval,l_id)) {
          data->addValue(LOAD_TNOM, atof(split_line[26].c_str()), l_id);
        } else {
          data->setValue(LOAD_TNOM, atof(split_line[26].c_str()), l_id);
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
      std::string sval;
      gridpack::utility::StringUtils util;
      sval = util.trimQuotes(split_line[1]);
      util.toUpper(sval);

      // LOAD_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());
      int nstr = split_line.size();

      // LOAD_ID
      if (nstr > 2) {
        sval = util.clean2Char(split_line[2]);
        strcpy(data.id, sval.c_str());
      }

      // LOAD_IT
      if (nstr > 3) {
        data.it = atoi(split_line[3].c_str());
      }

      // LOAD_RA
      if (nstr > 4) {
        data.ra = atof(split_line[4].c_str());
      }
      // LOAD_XA
      if (nstr > 5) {
        data.xa = atof(split_line[5].c_str());
      }

      // LOAD_XM
      if (nstr > 6) {
        data.xm = atof(split_line[6].c_str());
      }

      // LOAD_R1
      if (nstr > 7) {
        data.r1 = atof(split_line[7].c_str());
      }

      // LOAD_X1
      if (nstr > 8) {
        data.x1 = atof(split_line[8].c_str());
      }

      // LOAD_R2
      if (nstr > 9) {
        data.r2 = atof(split_line[9].c_str());
      }

      // LOAD_X2
      if (nstr > 10) {
        data.x2 = atof(split_line[10].c_str());
      }

      // LOAD_E1
      if (nstr > 11) {
        data.e1 = atof(split_line[11].c_str());
      }

      // LOAD_SE1
      if (nstr > 12) {
        data.se1 = atof(split_line[12].c_str());
      }

      // LOAD_E2
      if (nstr > 13) {
        data.e2 = atof(split_line[13].c_str());
      }

      // LOAD_SE2
      if (nstr > 14) {
        data.se2 = atof(split_line[14].c_str());
      }

      // LOAD_MBASE
      if (nstr > 15) {
        data.mbase = atof(split_line[15].c_str());
      }

      // LOAD_PMULT
      if (nstr > 16) {
        data.pmult = atof(split_line[16].c_str());
      }

      // LOAD_H
      if (nstr > 17) {
        data.h = atof(split_line[17].c_str());
      }

      // LOAD_VI
      if (nstr > 18) {
        data.vi = atof(split_line[18].c_str());
      }

      // LOAD_TI
      if (nstr > 19) {
        data.ti = atof(split_line[19].c_str());
      }

      // LOAD_TB
      if (nstr > 20) {
        data.tb = atof(split_line[20].c_str());
      }

      // LOAD_A
      if (nstr > 21) {
        data.a = atof(split_line[21].c_str());
      }

      // LOAD_B
      if (nstr > 22) {
        data.b = atof(split_line[22].c_str());
      }

      // LOAD_D
      if (nstr > 23) {
        data.d = atof(split_line[23].c_str());
      }

      // LOAD_E
      if (nstr > 24) {
        data.e = atof(split_line[24].c_str());
      }

      // LOAD_C0
      if (nstr > 25) {
        data.c0 = atof(split_line[25].c_str());
      }

      // LOAD_TNOM
      if (nstr > 26) {
        data.tnom = atof(split_line[26].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
