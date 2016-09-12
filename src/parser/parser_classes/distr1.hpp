/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: July 13, 2016
 *      Author: Bruce Palmer
 */
#ifndef DISTR1_HPP
#define DISTR1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Distr1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Distr1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Distr1Parser()
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

      // RELAY_ID
      if (!data->getValue(RELAY_ID,&stmp,r_id)) {
        data->addValue(RELAY_ID, data_struct.id, r_id);
      } else {
        data->setValue(RELAY_ID, data_struct.id, r_id);
      }

      // RELAY_RS
      if (!data->getValue(RELAY_RS,&ival,r_id)) {
        data->addValue(RELAY_RS, data_struct.rs, r_id);
      } else {
        data->setValue(RELAY_RS, data_struct.rs, r_id);
      }

      // RELAY_MTYPE
      if (!data->getValue(RELAY_MTYPE,&ival,r_id)) {
        data->addValue(RELAY_MTYPE, data_struct.mtype, r_id);
      } else {
        data->setValue(RELAY_MTYPE, data_struct.mtype, r_id);
      }

      // RELAY_BMON
      if (!data->getValue(RELAY_BMON,&ival,r_id)) {
        data->addValue(RELAY_BMON, data_struct.bmon, r_id);
      } else {
        data->setValue(RELAY_BMON, data_struct.bmon, r_id);
      }

      // RELAY_IBUS1
      if (!data->getValue(RELAY_IBUS1,&ival,r_id)) {
        data->addValue(RELAY_IBUS1, data_struct.ibus1, r_id);
      } else {
        data->setValue(RELAY_IBUS1, data_struct.ibus1, r_id);
      }

      // RELAY_JBUS1
      if (!data->getValue(RELAY_JBUS1,&ival,r_id)) {
        data->addValue(RELAY_JBUS1, data_struct.jbus1, r_id);
      } else {
        data->setValue(RELAY_JBUS1, data_struct.jbus1, r_id);
      }

      // RELAY_ID1
      if (!data->getValue(RELAY_ID1,&stmp,r_id)) {
        data->addValue(RELAY_ID1, data_struct.id1, r_id);
      } else {
        data->setValue(RELAY_ID1, data_struct.id1, r_id);
      }

      // RELAY_IBUS2
      if (!data->getValue(RELAY_IBUS2,&ival,r_id)) {
        data->addValue(RELAY_IBUS2, data_struct.ibus2, r_id);
      } else {
        data->setValue(RELAY_IBUS2, data_struct.ibus2, r_id);
      }

      // RELAY_JBUS2
      if (!data->getValue(RELAY_JBUS2,&ival,r_id)) {
        data->addValue(RELAY_JBUS2, data_struct.jbus2, r_id);
      } else {
        data->setValue(RELAY_JBUS2, data_struct.jbus2, r_id);
      }

      // RELAY_ID2
      if (!data->getValue(RELAY_ID2,&stmp,r_id)) {
        data->addValue(RELAY_ID2, data_struct.id2, r_id);
      } else {
        data->setValue(RELAY_ID2, data_struct.id2, r_id);
      }

      // RELAY_IBUS3
      if (!data->getValue(RELAY_IBUS3,&ival,r_id)) {
        data->addValue(RELAY_IBUS3, data_struct.ibus3, r_id);
      } else {
        data->setValue(RELAY_IBUS3, data_struct.ibus3, r_id);
      }

      // RELAY_JBUS3
      if (!data->getValue(RELAY_JBUS3,&ival,r_id)) {
        data->addValue(RELAY_JBUS3, data_struct.jbus3, r_id);
      } else {
        data->setValue(RELAY_JBUS3, data_struct.jbus3, r_id);
      }

      // RELAY_ID3
      if (!data->getValue(RELAY_ID3,&stmp,r_id)) {
        data->addValue(RELAY_ID3, data_struct.id3, r_id);
      } else {
        data->setValue(RELAY_ID3, data_struct.id3, r_id);
      }

      // RELAY_ZONE1_TIME
      if (!data->getValue(RELAY_ZONE1_TIME,&rval,r_id)) {
        data->addValue(RELAY_ZONE1_TIME, data_struct.zone1_time, r_id);
      } else {
        data->setValue(RELAY_ZONE1_TIME, data_struct.zone1_time, r_id);
      }

      // RELAY_ZONE1_REACH
      if (!data->getValue(RELAY_ZONE1_REACH,&rval,r_id)) {
        data->addValue(RELAY_ZONE1_REACH, data_struct.zone1_reach, r_id);
      } else {
        data->setValue(RELAY_ZONE1_REACH, data_struct.zone1_reach, r_id);
      }

      // RELAY_ZONE1_CENANG
      if (!data->getValue(RELAY_ZONE1_CENANG,&rval,r_id)) {
        data->addValue(RELAY_ZONE1_CENANG, data_struct.zone1_cenang, r_id);
      } else {
        data->setValue(RELAY_ZONE1_CENANG, data_struct.zone1_cenang, r_id);
      }

      // RELAY_ZONE1_CENDIS
      if (!data->getValue(RELAY_ZONE1_CENDIS,&rval,r_id)) {
        data->addValue(RELAY_ZONE1_CENDIS, data_struct.zone1_cendis, r_id);
      } else {
        data->setValue(RELAY_ZONE1_CENDIS, data_struct.zone1_cendis, r_id);
      }

      // RELAY_ZONE2_TIME
      if (!data->getValue(RELAY_ZONE2_TIME,&rval,r_id)) {
        data->addValue(RELAY_ZONE2_TIME, data_struct.zone2_time, r_id);
      } else {
        data->setValue(RELAY_ZONE2_TIME, data_struct.zone2_time, r_id);
      }

      // RELAY_ZONE2_REACH
      if (!data->getValue(RELAY_ZONE2_REACH,&rval,r_id)) {
        data->addValue(RELAY_ZONE2_REACH, data_struct.zone2_reach, r_id);
      } else {
        data->setValue(RELAY_ZONE2_REACH, data_struct.zone2_reach, r_id);
      }

      // RELAY_ZONE2_CENANG
      if (!data->getValue(RELAY_ZONE2_CENANG,&rval,r_id)) {
        data->addValue(RELAY_ZONE2_CENANG, data_struct.zone2_cenang, r_id);
      } else {
        data->setValue(RELAY_ZONE2_CENANG, data_struct.zone2_cenang, r_id);
      }

      // RELAY_ZONE2_CENDIS
      if (!data->getValue(RELAY_ZONE2_CENDIS,&rval,r_id)) {
        data->addValue(RELAY_ZONE2_CENDIS, data_struct.zone2_cendis, r_id);
      } else {
        data->setValue(RELAY_ZONE2_CENDIS, data_struct.zone2_cendis, r_id);
      }

      // RELAY_ZONE3_TIME
      if (!data->getValue(RELAY_ZONE3_TIME,&rval,r_id)) {
        data->addValue(RELAY_ZONE3_TIME, data_struct.zone2_time, r_id);
      } else {
        data->setValue(RELAY_ZONE3_TIME, data_struct.zone2_time, r_id);
      }

      // RELAY_ZONE3_REACH
      if (!data->getValue(RELAY_ZONE3_REACH,&rval,r_id)) {
        data->addValue(RELAY_ZONE3_REACH, data_struct.zone3_reach, r_id);
      } else {
        data->setValue(RELAY_ZONE3_REACH, data_struct.zone3_reach, r_id);
      }

      // RELAY_ZONE3_CENANG
      if (!data->getValue(RELAY_ZONE3_CENANG,&rval,r_id)) {
        data->addValue(RELAY_ZONE3_CENANG, data_struct.zone3_cenang, r_id);
      } else {
        data->setValue(RELAY_ZONE3_CENANG, data_struct.zone3_cenang, r_id);
      }

      // RELAY_ZONE3_CENDIS
      if (!data->getValue(RELAY_ZONE3_CENDIS,&rval,r_id)) {
        data->addValue(RELAY_ZONE3_CENDIS, data_struct.zone3_cendis, r_id);
      } else {
        data->setValue(RELAY_ZONE3_CENDIS, data_struct.zone3_cendis, r_id);
      }

      // RELAY_DIRANG
      if (!data->getValue(RELAY_DIRANG,&rval,r_id)) {
        data->addValue(RELAY_DIRANG, data_struct.dirang, r_id);
      } else {
        data->setValue(RELAY_DIRANG, data_struct.dirang, r_id);
      }

      // RELAY_THCUR
      if (!data->getValue(RELAY_THCUR,&rval,r_id)) {
        data->addValue(RELAY_THCUR, data_struct.thcur, r_id);
      } else {
        data->setValue(RELAY_THCUR, data_struct.thcur, r_id);
      }

      // RELAY_SEBTIME
      if (!data->getValue(RELAY_SEBTIME,&rval,r_id)) {
        data->addValue(RELAY_SEBTIME, data_struct.sebtime, r_id);
      } else {
        data->setValue(RELAY_SEBTIME, data_struct.sebtime, r_id);
      }

      // RELAY_SERCTIME
      if (!data->getValue(RELAY_SERCTIME,&rval,r_id)) {
        data->addValue(RELAY_SERCTIME, data_struct.serctime, r_id);
      } else {
        data->setValue(RELAY_SERCTIME, data_struct.serctime, r_id);
      }

      // RELAY_TRBTIME
      if (!data->getValue(RELAY_TRBTIME,&rval,r_id)) {
        data->addValue(RELAY_TRBTIME, data_struct.trbtime, r_id);
      } else {
        data->setValue(RELAY_TRBTIME, data_struct.trbtime, r_id);
      }

      // RELAY_TRRCTIME
      if (!data->getValue(RELAY_TRRCTIME,&rval,r_id)) {
        data->addValue(RELAY_TRRCTIME, data_struct.trrctime, r_id);
      } else {
        data->setValue(RELAY_TRRCTIME, data_struct.trrctime, r_id);
      }

      // RELAY_BLTYPE1
      if (!data->getValue(RELAY_BLTYPE1,&ival,r_id)) {
        data->addValue(RELAY_BLTYPE1, data_struct.bltype1, r_id);
      } else {
        data->setValue(RELAY_BLTYPE1, data_struct.bltype1, r_id);
      }

      // RELAY_BLINT1
      if (!data->getValue(RELAY_BLINT1,&rval,r_id)) {
        data->addValue(RELAY_BLINT1, data_struct.blint1, r_id);
      } else {
        data->setValue(RELAY_BLINT1, data_struct.blint1, r_id);
      }

      // RELAY_BLRO1
      if (!data->getValue(RELAY_BLRO1,&rval,r_id)) {
        data->addValue(RELAY_BLRO1, data_struct.blro1, r_id);
      } else {
        data->setValue(RELAY_BLRO1, data_struct.blro1, r_id);
      }

      // RELAY_BLTYPE2
      if (!data->getValue(RELAY_BLTYPE2,&ival,r_id)) {
        data->addValue(RELAY_BLTYPE2, data_struct.bltype2, r_id);
      } else {
        data->setValue(RELAY_BLTYPE2, data_struct.bltype2, r_id);
      }

      // RELAY_BLINT2
      if (!data->getValue(RELAY_BLINT2,&rval,r_id)) {
        data->addValue(RELAY_BLINT2, data_struct.blint2, r_id);
      } else {
        data->setValue(RELAY_BLINT2, data_struct.blint2, r_id);
      }

      // RELAY_BLRO2
      if (!data->getValue(RELAY_BLRO2,&rval,r_id)) {
        data->addValue(RELAY_BLRO2, data_struct.blro2, r_id);
      } else {
        data->setValue(RELAY_BLRO2, data_struct.blro2, r_id);
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

      // RELAY_ID
      if (nstr > 3) {
        model = util.clean2Char(split_line[3]);
        if (!data->getValue(RELAY_ID,&stmp,r_id)) {
          data->addValue(RELAY_ID, model.c_str(), r_id);
        } else {
          data->setValue(RELAY_ID, model.c_str(), r_id);
        }
      } 

      // RELAY_RS
      if (nstr > 4) {
        if (!data->getValue(RELAY_RS,&ival,r_id)) {
          data->addValue(RELAY_RS, atoi(split_line[4].c_str()), r_id);
        } else {
          data->setValue(RELAY_RS, atoi(split_line[4].c_str()), r_id);
        }
      }

      // RELAY_MTYPE
      if (nstr > 5) {
        if (!data->getValue(RELAY_MTYPE,&ival,r_id)) {
          data->addValue(RELAY_MTYPE, atoi(split_line[5].c_str()), r_id);
        } else {
          data->setValue(RELAY_MTYPE, atoi(split_line[5].c_str()), r_id);
        }
      }

      // RELAY_BMON
      if (nstr > 6) {
        if (!data->getValue(RELAY_BMON,&ival,r_id)) {
          data->addValue(RELAY_BMON, atoi(split_line[6].c_str()), r_id);
        } else {
          data->setValue(RELAY_BMON, atoi(split_line[6].c_str()), r_id);
        }
      }

      // RELAY_IBUS1
      if (nstr > 7) {
        if (!data->getValue(RELAY_IBUS1,&ival,r_id)) {
          data->addValue(RELAY_IBUS1, atoi(split_line[7].c_str()), r_id);
        } else {
          data->setValue(RELAY_IBUS1, atoi(split_line[7].c_str()), r_id);
        }
      }

      // RELAY_JBUS1
      if (nstr > 8) {
        if (!data->getValue(RELAY_JBUS1,&ival,r_id)) {
          data->addValue(RELAY_JBUS1, atoi(split_line[8].c_str()), r_id);
        } else {
          data->setValue(RELAY_JBUS1, atoi(split_line[8].c_str()), r_id);
        }
      }

      // RELAY_ID1
      if (nstr > 9) {
        model = util.clean2Char(split_line[9]);
        if (!data->getValue(RELAY_ID1,&stmp,r_id)) {
          data->addValue(RELAY_ID1, model.c_str(), r_id);
        } else {
          data->setValue(RELAY_ID1, model.c_str(), r_id);
        }
      } 

      // RELAY_IBUS2
      if (nstr > 10) {
        if (!data->getValue(RELAY_IBUS2,&ival,r_id)) {
          data->addValue(RELAY_IBUS2, atoi(split_line[10].c_str()), r_id);
        } else {
          data->setValue(RELAY_IBUS2, atoi(split_line[10].c_str()), r_id);
        }
      }

      // RELAY_JBUS2
      if (nstr > 11) {
        if (!data->getValue(RELAY_JBUS2,&ival,r_id)) {
          data->addValue(RELAY_JBUS2, atoi(split_line[11].c_str()), r_id);
        } else {
          data->setValue(RELAY_JBUS2, atoi(split_line[11].c_str()), r_id);
        }
      }

      // RELAY_ID2
      if (nstr > 12) {
        model = util.clean2Char(split_line[12]);
        if (!data->getValue(RELAY_ID2,&stmp,r_id)) {
          data->addValue(RELAY_ID2, model.c_str(), r_id);
        } else {
          data->setValue(RELAY_ID2, model.c_str(), r_id);
        }
      } 

      // RELAY_IBUS3
      if (nstr > 13) {
        if (!data->getValue(RELAY_IBUS3,&ival,r_id)) {
          data->addValue(RELAY_IBUS3, atoi(split_line[13].c_str()), r_id);
        } else {
          data->setValue(RELAY_IBUS3, atoi(split_line[13].c_str()), r_id);
        }
      }

      // RELAY_JBUS3
      if (nstr > 14) {
        if (!data->getValue(RELAY_JBUS3,&ival,r_id)) {
          data->addValue(RELAY_JBUS3, atoi(split_line[14].c_str()), r_id);
        } else {
          data->setValue(RELAY_JBUS3, atoi(split_line[14].c_str()), r_id);
        }
      }

      // RELAY_ID3
      if (nstr > 15) {
        model = util.clean2Char(split_line[15]);
        if (!data->getValue(RELAY_ID3,&stmp,r_id)) {
          data->addValue(RELAY_ID3, model.c_str(), r_id);
        } else {
          data->setValue(RELAY_ID3, model.c_str(), r_id);
        }
      } 

      // RELAY_ZONE1_TIME
      if (nstr > 16) {
        if (!data->getValue(RELAY_ZONE1_TIME,&rval,r_id)) {
          data->addValue(RELAY_ZONE1_TIME, atof(split_line[16].c_str()), r_id);
        } else {
          data->setValue(RELAY_ZONE1_TIME, atof(split_line[16].c_str()), r_id);
        }
      } 

      // RELAY_ZONE1_REACH
      if (nstr > 17) {
        if (!data->getValue(RELAY_ZONE1_REACH,&rval,r_id)) {
          data->addValue(RELAY_ZONE1_REACH, atof(split_line[17].c_str()), r_id);
        } else {
          data->setValue(RELAY_ZONE1_REACH, atof(split_line[17].c_str()), r_id);
        }
      } 

      // RELAY_ZONE1_CENANG
      if (nstr > 18) {
        if (!data->getValue(RELAY_ZONE1_CENANG,&rval,r_id)) {
          data->addValue(RELAY_ZONE1_CENANG, atof(split_line[18].c_str()), r_id);
        } else {
          data->setValue(RELAY_ZONE1_CENANG, atof(split_line[18].c_str()), r_id);
        }
      } 

      // RELAY_ZONE1_CENDIS
      if (nstr > 19) {
        if (!data->getValue(RELAY_ZONE1_CENDIS,&rval,r_id)) {
          data->addValue(RELAY_ZONE1_CENDIS, atof(split_line[19].c_str()), r_id);
        } else {
          data->setValue(RELAY_ZONE1_CENDIS, atof(split_line[19].c_str()), r_id);
        }
      } 

      // RELAY_ZONE2_TIME
      if (nstr > 20) {
        if (!data->getValue(RELAY_ZONE2_TIME,&rval,r_id)) {
          data->addValue(RELAY_ZONE2_TIME, atof(split_line[20].c_str()), r_id);
        } else {
          data->setValue(RELAY_ZONE2_TIME, atof(split_line[20].c_str()), r_id);
        }
      } 

      // RELAY_ZONE2_REACH
      if (nstr > 21) {
        if (!data->getValue(RELAY_ZONE2_REACH,&rval,r_id)) {
          data->addValue(RELAY_ZONE2_REACH, atof(split_line[21].c_str()), r_id);
        } else {
          data->setValue(RELAY_ZONE2_REACH, atof(split_line[21].c_str()), r_id);
        }
      } 

      // RELAY_ZONE2_CENANG
      if (nstr > 22) {
        if (!data->getValue(RELAY_ZONE2_CENANG,&rval,r_id)) {
          data->addValue(RELAY_ZONE2_CENANG, atof(split_line[22].c_str()), r_id);
        } else {
          data->setValue(RELAY_ZONE2_CENANG, atof(split_line[22].c_str()), r_id);
        }
      } 

      // RELAY_ZONE2_CENDIS
      if (nstr > 23) {
        if (!data->getValue(RELAY_ZONE2_CENDIS,&rval,r_id)) {
          data->addValue(RELAY_ZONE2_CENDIS, atof(split_line[23].c_str()), r_id);
        } else {
          data->setValue(RELAY_ZONE2_CENDIS, atof(split_line[23].c_str()), r_id);
        }
      } 

      // RELAY_ZONE3_TIME
      if (nstr > 24) {
        if (!data->getValue(RELAY_ZONE3_TIME,&rval,r_id)) {
          data->addValue(RELAY_ZONE3_TIME, atof(split_line[24].c_str()), r_id);
        } else {
          data->setValue(RELAY_ZONE3_TIME, atof(split_line[24].c_str()), r_id);
        }
      } 

      // RELAY_ZONE3_REACH
      if (nstr > 25) {
        if (!data->getValue(RELAY_ZONE3_REACH,&rval,r_id)) {
          data->addValue(RELAY_ZONE3_REACH, atof(split_line[25].c_str()), r_id);
        } else {
          data->setValue(RELAY_ZONE3_REACH, atof(split_line[25].c_str()), r_id);
        }
      } 

      // RELAY_ZONE3_CENANG
      if (nstr > 26) {
        if (!data->getValue(RELAY_ZONE3_CENANG,&rval,r_id)) {
          data->addValue(RELAY_ZONE3_CENANG, atof(split_line[26].c_str()), r_id);
        } else {
          data->setValue(RELAY_ZONE3_CENANG, atof(split_line[26].c_str()), r_id);
        }
      } 

      // RELAY_ZONE3_CENDIS
      if (nstr > 27) {
        if (!data->getValue(RELAY_ZONE3_CENDIS,&rval,r_id)) {
          data->addValue(RELAY_ZONE3_CENDIS, atof(split_line[27].c_str()), r_id);
        } else {
          data->setValue(RELAY_ZONE3_CENDIS, atof(split_line[27].c_str()), r_id);
        }
      } 

      // RELAY_DIRANG
      if (nstr > 28) {
        if (!data->getValue(RELAY_DIRANG,&rval,r_id)) {
          data->addValue(RELAY_DIRANG, atof(split_line[28].c_str()), r_id);
        } else {
          data->setValue(RELAY_DIRANG, atof(split_line[28].c_str()), r_id);
        }
      }

      // RELAY_THCUR
      if (nstr > 29) {
        if (!data->getValue(RELAY_THCUR,&rval,r_id)) {
          data->addValue(RELAY_THCUR, atof(split_line[29].c_str()), r_id);
        } else {
          data->setValue(RELAY_THCUR, atof(split_line[29].c_str()), r_id);
        }
      }

      // RELAY_SEBTIME
      if (nstr > 30) {
        if (!data->getValue(RELAY_SEBTIME,&rval,r_id)) {
          data->addValue(RELAY_SEBTIME, atof(split_line[30].c_str()), r_id);
        } else {
          data->setValue(RELAY_SEBTIME, atof(split_line[30].c_str()), r_id);
        }
      }

      // RELAY_SERCTIME
      if (nstr > 31) {
        if (!data->getValue(RELAY_SERCTIME,&rval,r_id)) {
          data->addValue(RELAY_SERCTIME, atof(split_line[31].c_str()), r_id);
        } else {
          data->setValue(RELAY_SERCTIME, atof(split_line[31].c_str()), r_id);
        }
      }

      // RELAY_TRBTIME
      if (nstr > 32) {
        if (!data->getValue(RELAY_TRBTIME,&rval,r_id)) {
          data->addValue(RELAY_TRBTIME, atof(split_line[32].c_str()), r_id);
        } else {
          data->setValue(RELAY_TRBTIME, atof(split_line[32].c_str()), r_id);
        }
      }

      // RELAY_TRRCTIME
      if (nstr > 33) {
        if (!data->getValue(RELAY_TRRCTIME,&rval,r_id)) {
          data->addValue(RELAY_TRRCTIME, atof(split_line[33].c_str()), r_id);
        } else {
          data->setValue(RELAY_TRRCTIME, atof(split_line[33].c_str()), r_id);
        }
      }

      // RELAY_BLTYPE1
      if (nstr > 34) {
        if (!data->getValue(RELAY_BLTYPE1,&ival,r_id)) {
          data->addValue(RELAY_BLTYPE1, atoi(split_line[34].c_str()), r_id);
        } else {
          data->setValue(RELAY_BLTYPE1, atoi(split_line[34].c_str()), r_id);
        }
      }

      // RELAY_BLINT1
      if (nstr > 35) {
        if (!data->getValue(RELAY_BLINT1,&rval,r_id)) {
          data->addValue(RELAY_BLINT1, atof(split_line[35].c_str()), r_id);
        } else {
          data->setValue(RELAY_BLINT1, atof(split_line[35].c_str()), r_id);
        }
      } 

      // RELAY_BLRO1
      if (nstr > 36) {
        if (!data->getValue(RELAY_BLRO1,&rval,r_id)) {
          data->addValue(RELAY_BLRO1, atof(split_line[36].c_str()), r_id);
        } else {
          data->setValue(RELAY_BLRO1, atof(split_line[36].c_str()), r_id);
        }
      } 

      // RELAY_BLTYPE2
      if (nstr > 37) {
        if (!data->getValue(RELAY_BLTYPE2,&ival,r_id)) {
          data->addValue(RELAY_BLTYPE2, atoi(split_line[37].c_str()), r_id);
        } else {
          data->setValue(RELAY_BLTYPE2, atoi(split_line[37].c_str()), r_id);
        }
      }

      // RELAY_BLINT2
      if (nstr > 38) {
        if (!data->getValue(RELAY_BLINT2,&rval,r_id)) {
          data->addValue(RELAY_BLINT2, atof(split_line[38].c_str()), r_id);
        } else {
          data->setValue(RELAY_BLINT2, atof(split_line[38].c_str()), r_id);
        }
      } 

      // RELAY_BLRO2
      if (nstr > 39) {
        if (!data->getValue(RELAY_BLRO2,&rval,r_id)) {
          data->addValue(RELAY_BLRO2, atof(split_line[39].c_str()), r_id);
        } else {
          data->setValue(RELAY_BLRO2, atof(split_line[39].c_str()), r_id);
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
      // RELAY_FROM_BUS
      int from_idx, to_idx;
      from_idx = atoi(split_line[0].c_str());
      data.from_bus = from_idx;
      to_idx = atoi(split_line[2].c_str());
      data.to_bus = to_idx;

      std::string sval;
      gridpack::utility::StringUtils util;
      sval = util.trimQuotes(split_line[1]);
      util.toUpper(sval);

      // RELAY_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());
      int nstr = split_line.size();

      // RELAY_ID
      if (nstr > 3) {
        sval = util.clean2Char(split_line[3]);
        strcpy(data.id, sval.c_str());
      }

      // RELAY_RS
      if (nstr > 4) {
        data.rs = atoi(split_line[4].c_str());
      }

      // RELAY_MTYPE
      if (nstr > 5) {
        data.mtype = atoi(split_line[5].c_str());
      }

      // RELAY_BMON
      if (nstr > 6) {
        data.bmon = atoi(split_line[6].c_str());
      }

      // RELAY_IBUS1
      if (nstr > 7) {
        data.ibus1 = atoi(split_line[7].c_str());
      }

      // RELAY_JBUS1
      if (nstr > 8) {
        data.jbus1 = atoi(split_line[8].c_str());
      }

      // RELAY_ID1
      if (nstr > 9) {
        sval = util.clean2Char(split_line[9]);
        strcpy(data.id1, sval.c_str());
      }

      // RELAY_IBUS2
      if (nstr > 10) {
        data.ibus2 = atoi(split_line[10].c_str());
      }

      // RELAY_JBUS2
      if (nstr > 11) {
        data.jbus2 = atoi(split_line[11].c_str());
      }

      // RELAY_ID2
      if (nstr > 12) {
        sval = util.clean2Char(split_line[12]);
        strcpy(data.id2, sval.c_str());
      }

      // RELAY_IBUS3
      if (nstr > 13) {
        data.ibus3 = atoi(split_line[13].c_str());
      }

      // RELAY_JBUS3
      if (nstr > 14) {
        data.jbus3 = atoi(split_line[14].c_str());
      }

      // RELAY_ID3
      if (nstr > 15) {
        sval = util.clean2Char(split_line[15]);
        strcpy(data.id3, sval.c_str());
      }

      // RELAY_ZONE1_TIME
      if (nstr > 16) {
        data.zone1_time = atof(split_line[16].c_str());
      }

      // RELAY_ZONE1_REACH
      if (nstr > 17) {
        data.zone1_reach = atof(split_line[17].c_str());
      }

      // RELAY_ZONE1_CENANG
      if (nstr > 18) {
        data.zone1_cenang = atof(split_line[18].c_str());
      }

      // RELAY_ZONE1_CENDIS
      if (nstr > 19) {
        data.zone1_cendis = atof(split_line[19].c_str());
      }

      // RELAY_ZONE2_TIME
      if (nstr > 20) {
        data.zone2_time = atof(split_line[20].c_str());
      }

      // RELAY_ZONE2_REACH
      if (nstr > 21) {
        data.zone2_reach = atof(split_line[21].c_str());
      }

      // RELAY_ZONE2_CENANG
      if (nstr > 22) {
        data.zone2_cenang = atof(split_line[22].c_str());
      }

      // RELAY_ZONE2_CENDIS
      if (nstr > 23) {
        data.zone2_cendis = atof(split_line[23].c_str());
      }

      // RELAY_ZONE3_TIME
      if (nstr > 24) {
        data.zone3_time = atof(split_line[24].c_str());
      }

      // RELAY_ZONE3_REACH
      if (nstr > 25) {
        data.zone3_reach = atof(split_line[25].c_str());
      }

      // RELAY_ZONE3_CENANG
      if (nstr > 26) {
        data.zone3_cenang = atof(split_line[26].c_str());
      }

      // RELAY_ZONE3_CENDIS
      if (nstr > 27) {
        data.zone3_cendis = atof(split_line[27].c_str());
      }

      // RELAY_DIRANG
      if (nstr > 28) {
        data.dirang = atof(split_line[28].c_str());
      }

      // RELAY_THCUR
      if (nstr > 29) {
        data.thcur = atof(split_line[29].c_str());
      }

      // RELAY_SEBTIME
      if (nstr > 30) {
        data.sebtime = atof(split_line[30].c_str());
      }

      // RELAY_SERCTIME
      if (nstr > 31) {
        data.serctime = atof(split_line[31].c_str());
      }

      // RELAY_TRBTIME
      if (nstr > 32) {
        data.trbtime = atof(split_line[32].c_str());
      }

      // RELAY_TRRCTIME
      if (nstr > 33) {
        data.trrctime = atof(split_line[33].c_str());
      }

      // RELAY_BLTYPE1
      if (nstr > 34) {
        data.bltype1 = atoi(split_line[34].c_str());
      }

      // RELAY_BLINT1
      if (nstr > 35) {
        data.blint1 = atof(split_line[35].c_str());
      }

      // RELAY_BLRO1
      if (nstr > 36) {
        data.blro1 = atof(split_line[36].c_str());
      }

      // RELAY_BLTYPE2
      if (nstr > 37) {
        data.bltype2 = atoi(split_line[37].c_str());
      }

      // RELAY_BLINT2
      if (nstr > 38) {
        data.blint2 = atof(split_line[38].c_str());
      }

      // RELAY_BLRO2
      if (nstr > 39) {
        data.blro2 = atof(split_line[39].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
