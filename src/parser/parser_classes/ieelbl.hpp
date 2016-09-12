/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: September 8, 2016
 *      Author: Bruce Palmer
 */
#ifndef IEELBL_HPP
#define IEELBL_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class IeelblParser
{
  public:
    /**
     * Constructor
     */
    explicit IeelblParser()
    {
    }

    /**
     * Destructor
     */
    virtual ~IeelblParser()
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

      // LOAD_A1
      if (!data->getValue(LOAD_A1,&rval,l_id)) {
        data->addValue(LOAD_A1, data_struct.a1, l_id);
      } else {
        data->setValue(LOAD_A1, data_struct.a1, l_id);
      }

      // LOAD_A2
      if (!data->getValue(LOAD_A2,&rval,l_id)) {
        data->addValue(LOAD_A2, data_struct.a2, l_id);
      } else {
        data->setValue(LOAD_A2, data_struct.a2, l_id);
      }

      // LOAD_A3
      if (!data->getValue(LOAD_A3,&rval,l_id)) {
        data->addValue(LOAD_A3, data_struct.a3, l_id);
      } else {
        data->setValue(LOAD_A3, data_struct.a3, l_id);
      }

      // LOAD_A4
      if (!data->getValue(LOAD_A4,&rval,l_id)) {
        data->addValue(LOAD_A4, data_struct.a4, l_id);
      } else {
        data->setValue(LOAD_A4, data_struct.a4, l_id);
      }

      // LOAD_A5
      if (!data->getValue(LOAD_A5,&rval,l_id)) {
        data->addValue(LOAD_A5, data_struct.a5, l_id);
      } else {
        data->setValue(LOAD_A5, data_struct.a5, l_id);
      }

      // LOAD_A6
      if (!data->getValue(LOAD_A6,&rval,l_id)) {
        data->addValue(LOAD_A6, data_struct.a6, l_id);
      } else {
        data->setValue(LOAD_A6, data_struct.a6, l_id);
      }

      // LOAD_A7
      if (!data->getValue(LOAD_A7,&rval,l_id)) {
        data->addValue(LOAD_A7, data_struct.a7, l_id);
      } else {
        data->setValue(LOAD_A7, data_struct.a7, l_id);
      }

      // LOAD_A8
      if (!data->getValue(LOAD_A8,&rval,l_id)) {
        data->addValue(LOAD_A8, data_struct.a8, l_id);
      } else {
        data->setValue(LOAD_A8, data_struct.a8, l_id);
      }

      // LOAD_N1
      if (!data->getValue(LOAD_N1,&rval,l_id)) {
        data->addValue(LOAD_N1, data_struct.n1, l_id);
      } else {
        data->setValue(LOAD_N1, data_struct.n1, l_id);
      }

      // LOAD_N2
      if (!data->getValue(LOAD_N2,&rval,l_id)) {
        data->addValue(LOAD_N2, data_struct.n2, l_id);
      } else {
        data->setValue(LOAD_N2, data_struct.n2, l_id);
      }

      // LOAD_N3
      if (!data->getValue(LOAD_N3,&rval,l_id)) {
        data->addValue(LOAD_N3, data_struct.n3, l_id);
      } else {
        data->setValue(LOAD_N3, data_struct.n3, l_id);
      }

      // LOAD_N4
      if (!data->getValue(LOAD_N4,&rval,l_id)) {
        data->addValue(LOAD_N4, data_struct.n4, l_id);
      } else {
        data->setValue(LOAD_N4, data_struct.n4, l_id);
      }

      // LOAD_N5
      if (!data->getValue(LOAD_N5,&rval,l_id)) {
        data->addValue(LOAD_N5, data_struct.n5, l_id);
      } else {
        data->setValue(LOAD_N5, data_struct.n5, l_id);
      }

      // LOAD_N6
      if (!data->getValue(LOAD_N6,&rval,l_id)) {
        data->addValue(LOAD_N6, data_struct.n6, l_id);
      } else {
        data->setValue(LOAD_N6, data_struct.n6, l_id);
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

      // LOAD_A1
      if (nstr > 3) {
        if (!data->getValue(LOAD_A1,&ival,l_id)) {
          data->addValue(LOAD_A1, atof(split_line[3].c_str()), l_id);
        } else {
          data->setValue(LOAD_A1, atof(split_line[3].c_str()), l_id);
        }
      }

      // LOAD_A2
      if (nstr > 4) {
        if (!data->getValue(LOAD_A2,&rval,l_id)) {
          data->addValue(LOAD_A2, atof(split_line[4].c_str()), l_id);
        } else {
          data->setValue(LOAD_A2, atof(split_line[4].c_str()), l_id);
        }
      } 

      // LOAD_A3
      if (nstr > 5) {
        if (!data->getValue(LOAD_A3,&rval,l_id)) {
          data->addValue(LOAD_A3, atof(split_line[5].c_str()), l_id);
        } else {
          data->setValue(LOAD_A3, atof(split_line[5].c_str()), l_id);
        }
      } 

      // LOAD_A4
      if (nstr > 6) {
        if (!data->getValue(LOAD_A4,&rval,l_id)) {
          data->addValue(LOAD_A4, atof(split_line[6].c_str()), l_id);
        } else {
          data->setValue(LOAD_A4, atof(split_line[6].c_str()), l_id);
        }
      } 

      // LOAD_A5
      if (nstr > 7) {
        if (!data->getValue(LOAD_A5,&rval,l_id)) {
          data->addValue(LOAD_A5, atof(split_line[7].c_str()), l_id);
        } else {
          data->setValue(LOAD_A5, atof(split_line[7].c_str()), l_id);
        }
      } 

      // LOAD_A6
      if (nstr > 8) {
        if (!data->getValue(LOAD_A6,&rval,l_id)) {
          data->addValue(LOAD_A6, atof(split_line[8].c_str()), l_id);
        } else {
          data->setValue(LOAD_A6, atof(split_line[8].c_str()), l_id);
        }
      } 

      // LOAD_A7
      if (nstr > 9) {
        if (!data->getValue(LOAD_A7,&rval,l_id)) {
          data->addValue(LOAD_A7, atof(split_line[9].c_str()), l_id);
        } else {
          data->setValue(LOAD_A7, atof(split_line[9].c_str()), l_id);
        }
      } 

      // LOAD_A8
      if (nstr > 10) {
        if (!data->getValue(LOAD_A8,&rval,l_id)) {
          data->addValue(LOAD_A8, atof(split_line[10].c_str()), l_id);
        } else {
          data->setValue(LOAD_A8, atof(split_line[10].c_str()), l_id);
        }
      } 

      // LOAD_N1
      if (nstr > 11) {
        if (!data->getValue(LOAD_N1,&rval,l_id)) {
          data->addValue(LOAD_N1, atof(split_line[11].c_str()), l_id);
        } else {
          data->setValue(LOAD_N1, atof(split_line[11].c_str()), l_id);
        }
      } 

      // LOAD_N2
      if (nstr > 12) {
        if (!data->getValue(LOAD_N2,&rval,l_id)) {
          data->addValue(LOAD_N2, atof(split_line[12].c_str()), l_id);
        } else {
          data->setValue(LOAD_N2, atof(split_line[12].c_str()), l_id);
        }
      } 

      // LOAD_N3
      if (nstr > 13) {
        if (!data->getValue(LOAD_N3,&rval,l_id)) {
          data->addValue(LOAD_N3, atof(split_line[13].c_str()), l_id);
        } else {
          data->setValue(LOAD_N3, atof(split_line[13].c_str()), l_id);
        }
      } 

      // LOAD_N4
      if (nstr > 14) {
        if (!data->getValue(LOAD_N4,&rval,l_id)) {
          data->addValue(LOAD_N4, atof(split_line[14].c_str()), l_id);
        } else {
          data->setValue(LOAD_N4, atof(split_line[14].c_str()), l_id);
        }
      } 

      // LOAD_N5
      if (nstr > 15) {
        if (!data->getValue(LOAD_N5,&rval,l_id)) {
          data->addValue(LOAD_N5, atof(split_line[15].c_str()), l_id);
        } else {
          data->setValue(LOAD_N5, atof(split_line[15].c_str()), l_id);
        }
      } 

      // LOAD_N6
      if (nstr > 16) {
        if (!data->getValue(LOAD_N6,&rval,l_id)) {
          data->addValue(LOAD_N6, atof(split_line[16].c_str()), l_id);
        } else {
          data->setValue(LOAD_N6, atof(split_line[16].c_str()), l_id);
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

      // LOAD_A1
      if (nstr > 3) {
        data.a1 = atof(split_line[3].c_str());
      }

      // LOAD_A2
      if (nstr > 4) {
        data.a2 = atof(split_line[4].c_str());
      }
      // LOAD_A3
      if (nstr > 5) {
        data.a3 = atof(split_line[5].c_str());
      }

      // LOAD_A4
      if (nstr > 6) {
        data.a4 = atof(split_line[6].c_str());
      }

      // LOAD_A5
      if (nstr > 7) {
        data.a5 = atof(split_line[7].c_str());
      }

      // LOAD_A6
      if (nstr > 8) {
        data.a6 = atof(split_line[8].c_str());
      }

      // LOAD_A7
      if (nstr > 9) {
        data.a7 = atof(split_line[9].c_str());
      }

      // LOAD_A8
      if (nstr > 10) {
        data.a8 = atof(split_line[10].c_str());
      }

      // LOAD_N1
      if (nstr > 11) {
        data.n1 = atof(split_line[11].c_str());
      }

      // LOAD_N2
      if (nstr > 12) {
        data.n2 = atof(split_line[12].c_str());
      }

      // LOAD_N3
      if (nstr > 13) {
        data.n3 = atof(split_line[13].c_str());
      }

      // LOAD_N4
      if (nstr > 14) {
        data.n4 = atof(split_line[14].c_str());
      }

      // LOAD_N5
      if (nstr > 15) {
        data.n5 = atof(split_line[15].c_str());
      }

      // LOAD_N6
      if (nstr > 16) {
        data.n6 = atof(split_line[16].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
