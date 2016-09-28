/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: September 8, 2016
 *      Author: Bruce Palmer
 */
#ifndef ACMTBLU1_HPP
#define ACMTBLU1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Acmtblu1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Acmtblu1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Acmtblu1Parser()
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

      // LOAD_TSTALL
      if (!data->getValue(LOAD_TSTALL,&rval,l_id)) {
        data->addValue(LOAD_TSTALL, data_struct.tstall, l_id);
      } else {
        data->setValue(LOAD_TSTALL, data_struct.tstall, l_id);
      }

      // LOAD_TRESTART
      if (!data->getValue(LOAD_TRESTART,&rval,l_id)) {
        data->addValue(LOAD_TRESTART, data_struct.trestart, l_id);
      } else {
        data->setValue(LOAD_TRESTART, data_struct.trestart, l_id);
      }

      // LOAD_TV
      if (!data->getValue(LOAD_TV,&rval,l_id)) {
        data->addValue(LOAD_TV, data_struct.tv, l_id);
      } else {
        data->setValue(LOAD_TV, data_struct.tv, l_id);
      }

      // LOAD_TF
      if (!data->getValue(LOAD_TF,&rval,l_id)) {
        data->addValue(LOAD_TF, data_struct.tf, l_id);
      } else {
        data->setValue(LOAD_TF, data_struct.tf, l_id);
      }

      // LOAD_COMPLF
      if (!data->getValue(LOAD_COMPLF,&rval,l_id)) {
        data->addValue(LOAD_COMPLF, data_struct.complf, l_id);
      } else {
        data->setValue(LOAD_COMPLF, data_struct.complf, l_id);
      }

      // LOAD_COMPPF
      if (!data->getValue(LOAD_COMPPF,&rval,l_id)) {
        data->addValue(LOAD_COMPPF, data_struct.comppf, l_id);
      } else {
        data->setValue(LOAD_COMPPF, data_struct.comppf, l_id);
      }

      // LOAD_VSTALL
      if (!data->getValue(LOAD_VSTALL,&rval,l_id)) {
        data->addValue(LOAD_VSTALL, data_struct.vstall, l_id);
      } else {
        data->setValue(LOAD_VSTALL, data_struct.vstall, l_id);
      }

      // LOAD_RSTALL
      if (!data->getValue(LOAD_RSTALL,&rval,l_id)) {
        data->addValue(LOAD_RSTALL, data_struct.rstall, l_id);
      } else {
        data->setValue(LOAD_RSTALL, data_struct.rstall, l_id);
      }

      // LOAD_XSTALL
      if (!data->getValue(LOAD_XSTALL,&rval,l_id)) {
        data->addValue(LOAD_XSTALL, data_struct.xstall, l_id);
      } else {
        data->setValue(LOAD_XSTALL, data_struct.xstall, l_id);
      }

      // LOAD_LFADJ
      if (!data->getValue(LOAD_LFADJ,&rval,l_id)) {
        data->addValue(LOAD_LFADJ, data_struct.lfadj, l_id);
      } else {
        data->setValue(LOAD_LFADJ, data_struct.lfadj, l_id);
      }

      // LOAD_KP1
      if (!data->getValue(LOAD_KP1,&rval,l_id)) {
        data->addValue(LOAD_KP1, data_struct.kp1, l_id);
      } else {
        data->setValue(LOAD_KP1, data_struct.kp1, l_id);
      }

      // LOAD_NP1
      if (!data->getValue(LOAD_NP1,&rval,l_id)) {
        data->addValue(LOAD_NP1, data_struct.np1, l_id);
      } else {
        data->setValue(LOAD_NP1, data_struct.np1, l_id);
      }

      // LOAD_KQ1
      if (!data->getValue(LOAD_KQ1,&rval,l_id)) {
        data->addValue(LOAD_KQ1, data_struct.kq1, l_id);
      } else {
        data->setValue(LOAD_KQ1, data_struct.kq1, l_id);
      }

      // LOAD_NQ1
      if (!data->getValue(LOAD_NQ1,&rval,l_id)) {
        data->addValue(LOAD_NQ1, data_struct.kq1, l_id);
      } else {
        data->setValue(LOAD_NQ1, data_struct.kq1, l_id);
      }

      // LOAD_KP2
      if (!data->getValue(LOAD_KP2,&rval,l_id)) {
        data->addValue(LOAD_KP2, data_struct.kp2, l_id);
      } else {
        data->setValue(LOAD_KP2, data_struct.kp2, l_id);
      }

      // LOAD_NP2
      if (!data->getValue(LOAD_NP2,&rval,l_id)) {
        data->addValue(LOAD_NP2, data_struct.np2, l_id);
      } else {
        data->setValue(LOAD_NP2, data_struct.np2, l_id);
      }

      // LOAD_KQ2
      if (!data->getValue(LOAD_KQ2,&rval,l_id)) {
        data->addValue(LOAD_KQ2, data_struct.kq2, l_id);
      } else {
        data->setValue(LOAD_KQ2, data_struct.kq2, l_id);
      }

      // LOAD_NQ2
      if (!data->getValue(LOAD_NQ2,&rval,l_id)) {
        data->addValue(LOAD_NQ2, data_struct.kq2, l_id);
      } else {
        data->setValue(LOAD_NQ2, data_struct.kq2, l_id);
      }

      // LOAD_VBRK
      if (!data->getValue(LOAD_VBRK,&rval,l_id)) {
        data->addValue(LOAD_VBRK, data_struct.vbrk, l_id);
      } else {
        data->setValue(LOAD_VBRK, data_struct.vbrk, l_id);
      }

      // LOAD_FRST
      if (!data->getValue(LOAD_FRST,&rval,l_id)) {
        data->addValue(LOAD_FRST, data_struct.frst, l_id);
      } else {
        data->setValue(LOAD_FRST, data_struct.frst, l_id);
      }

      // LOAD_VRST
      if (!data->getValue(LOAD_VRST,&rval,l_id)) {
        data->addValue(LOAD_VRST, data_struct.vrst, l_id);
      } else {
        data->setValue(LOAD_VRST, data_struct.vrst, l_id);
      }

      // LOAD_CMPKPF
      if (!data->getValue(LOAD_CMPKPF,&rval,l_id)) {
        data->addValue(LOAD_CMPKPF, data_struct.cmpkpf, l_id);
      } else {
        data->setValue(LOAD_CMPKPF, data_struct.cmpkpf, l_id);
      }

      // LOAD_CMPKQF
      if (!data->getValue(LOAD_CMPKQF,&rval,l_id)) {
        data->addValue(LOAD_CMPKQF, data_struct.cmpkqf, l_id);
      } else {
        data->setValue(LOAD_CMPKQF, data_struct.cmpkqf, l_id);
      }

      // LOAD_VC1OFF
      if (!data->getValue(LOAD_VC1OFF,&rval,l_id)) {
        data->addValue(LOAD_VC1OFF, data_struct.vc1off, l_id);
      } else {
        data->setValue(LOAD_VC1OFF, data_struct.vc1off, l_id);
      }

      // LOAD_VC2OFF
      if (!data->getValue(LOAD_VC2OFF,&rval,l_id)) {
        data->addValue(LOAD_VC2OFF, data_struct.vc2off, l_id);
      } else {
        data->setValue(LOAD_VC2OFF, data_struct.vc2off, l_id);
      }

      // LOAD_VC1ON
      if (!data->getValue(LOAD_VC1ON,&rval,l_id)) {
        data->addValue(LOAD_VC1ON, data_struct.vc1on, l_id);
      } else {
        data->setValue(LOAD_VC1ON, data_struct.vc1on, l_id);
      }

      // LOAD_VC2ON
      if (!data->getValue(LOAD_VC2ON,&rval,l_id)) {
        data->addValue(LOAD_VC2ON, data_struct.vc2on, l_id);
      } else {
        data->setValue(LOAD_VC2ON, data_struct.vc2on, l_id);
      }

      // LOAD_TTH
      if (!data->getValue(LOAD_TTH,&rval,l_id)) {
        data->addValue(LOAD_TTH, data_struct.tth, l_id);
      } else {
        data->setValue(LOAD_TTH, data_struct.tth, l_id);
      }

      // LOAD_TH1T
      if (!data->getValue(LOAD_TH1T,&rval,l_id)) {
        data->addValue(LOAD_TH1T, data_struct.th1t, l_id);
      } else {
        data->setValue(LOAD_TH1T, data_struct.th1t, l_id);
      }

      // LOAD_TH2T
      if (!data->getValue(LOAD_TH2T,&rval,l_id)) {
        data->addValue(LOAD_TH2T, data_struct.th2t, l_id);
      } else {
        data->setValue(LOAD_TH2T, data_struct.th2t, l_id);
      }

      // LOAD_FUVR
      if (!data->getValue(LOAD_FUVR,&rval,l_id)) {
        data->addValue(LOAD_FUVR, data_struct.fuvr, l_id);
      } else {
        data->setValue(LOAD_FUVR, data_struct.fuvr, l_id);
      }

      // LOAD_UVTR1
      if (!data->getValue(LOAD_UVTR1,&rval,l_id)) {
        data->addValue(LOAD_UVTR1, data_struct.uvtr1, l_id);
      } else {
        data->setValue(LOAD_UVTR1, data_struct.uvtr1, l_id);
      }

      // LOAD_TTR1
      if (!data->getValue(LOAD_TTR1,&rval,l_id)) {
        data->addValue(LOAD_TTR1, data_struct.ttr1, l_id);
      } else {
        data->setValue(LOAD_TTR1, data_struct.ttr1, l_id);
      }

      // LOAD_UVTR2
      if (!data->getValue(LOAD_UVTR2,&rval,l_id)) {
        data->addValue(LOAD_UVTR2, data_struct.uvtr2, l_id);
      } else {
        data->setValue(LOAD_UVTR2, data_struct.uvtr2, l_id);
      }

      // LOAD_TTR2
      if (!data->getValue(LOAD_TTR2,&rval,l_id)) {
        data->addValue(LOAD_TTR2, data_struct.ttr2, l_id);
      } else {
        data->setValue(LOAD_TTR2, data_struct.ttr2, l_id);
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
      model = util.trimQuotes(split_line[3]);
      util.toUpper(model);
      if (!data->getValue(LOAD_MODEL,&stmp,l_id)) {
        data->addValue(LOAD_MODEL, model.c_str(), l_id);
      } else {
        data->setValue(LOAD_MODEL, model.c_str(), l_id);
      }

      // LOAD_TSTALL
      if (nstr > 11) {
        if (!data->getValue(LOAD_TSTALL,&rval,l_id)) {
          data->addValue(LOAD_TSTALL, atof(split_line[11].c_str()), l_id);
        } else {
          data->setValue(LOAD_TSTALL, atof(split_line[11].c_str()), l_id);
        }
      } 

      // LOAD_TRESTART
      if (nstr > 12) {
        if (!data->getValue(LOAD_TRESTART,&rval,l_id)) {
          data->addValue(LOAD_TRESTART, atof(split_line[12].c_str()), l_id);
        } else {
          data->setValue(LOAD_TRESTART, atof(split_line[12].c_str()), l_id);
        }
      } 

      // LOAD_TV
      if (nstr > 13) {
        if (!data->getValue(LOAD_TV,&rval,l_id)) {
          data->addValue(LOAD_TV, atof(split_line[13].c_str()), l_id);
        } else {
          data->setValue(LOAD_TV, atof(split_line[13].c_str()), l_id);
        }
      } 

      // LOAD_TF
      if (nstr > 14) {
        if (!data->getValue(LOAD_TF,&rval,l_id)) {
          data->addValue(LOAD_TF, atof(split_line[14].c_str()), l_id);
        } else {
          data->setValue(LOAD_TF, atof(split_line[14].c_str()), l_id);
        }
      } 

      // LOAD_COMPLF
      if (nstr > 15) {
        if (!data->getValue(LOAD_COMPLF,&rval,l_id)) {
          data->addValue(LOAD_COMPLF, atof(split_line[15].c_str()), l_id);
        } else {
          data->setValue(LOAD_COMPLF, atof(split_line[15].c_str()), l_id);
        }
      } 

      // LOAD_COMPPF
      if (nstr > 16) {
        if (!data->getValue(LOAD_COMPPF,&rval,l_id)) {
          data->addValue(LOAD_COMPPF, atof(split_line[16].c_str()), l_id);
        } else {
          data->setValue(LOAD_COMPPF, atof(split_line[16].c_str()), l_id);
        }
      } 

      // LOAD_VSTALL
      if (nstr > 17) {
        if (!data->getValue(LOAD_VSTALL,&rval,l_id)) {
          data->addValue(LOAD_VSTALL, atof(split_line[17].c_str()), l_id);
        } else {
          data->setValue(LOAD_VSTALL, atof(split_line[17].c_str()), l_id);
        }
      }

      // LOAD_RSTALL
      if (nstr > 18) {
        if (!data->getValue(LOAD_RSTALL,&rval,l_id)) {
          data->addValue(LOAD_RSTALL, atof(split_line[18].c_str()), l_id);
        } else {
          data->setValue(LOAD_RSTALL, atof(split_line[18].c_str()), l_id);
        }
      }

      // LOAD_XSTALL
      if (nstr > 19) {
        if (!data->getValue(LOAD_XSTALL,&rval,l_id)) {
          data->addValue(LOAD_XSTALL, atof(split_line[19].c_str()), l_id);
        } else {
          data->setValue(LOAD_XSTALL, atof(split_line[19].c_str()), l_id);
        }
      }

      // LOAD_LFADJ
      if (nstr > 20) {
        if (!data->getValue(LOAD_LFADJ,&rval,l_id)) {
          data->addValue(LOAD_LFADJ, atof(split_line[20].c_str()), l_id);
        } else {
          data->setValue(LOAD_LFADJ, atof(split_line[20].c_str()), l_id);
        }
      }

      // LOAD_KP1
      if (nstr > 21) {
        if (!data->getValue(LOAD_KP1,&rval,l_id)) {
          data->addValue(LOAD_KP1, atof(split_line[21].c_str()), l_id);
        } else {
          data->setValue(LOAD_KP1, atof(split_line[21].c_str()), l_id);
        }
      }

      // LOAD_NP1
      if (nstr > 22) {
        if (!data->getValue(LOAD_NP1,&rval,l_id)) {
          data->addValue(LOAD_NP1, atof(split_line[22].c_str()), l_id);
        } else {
          data->setValue(LOAD_NP1, atof(split_line[22].c_str()), l_id);
        }
      }

      // LOAD_KQ1
      if (nstr > 23) {
        if (!data->getValue(LOAD_KQ1,&rval,l_id)) {
          data->addValue(LOAD_KQ1, atof(split_line[23].c_str()), l_id);
        } else {
          data->setValue(LOAD_KQ1, atof(split_line[23].c_str()), l_id);
        }
      }

      // LOAD_NQ1
      if (nstr > 24) {
        if (!data->getValue(LOAD_NQ1,&rval,l_id)) {
          data->addValue(LOAD_NQ1, atof(split_line[24].c_str()), l_id);
        } else {
          data->setValue(LOAD_NQ1, atof(split_line[24].c_str()), l_id);
        }
      } 

      // LOAD_KP2
      if (nstr > 25) {
        if (!data->getValue(LOAD_KP2,&rval,l_id)) {
          data->addValue(LOAD_KP2, atof(split_line[25].c_str()), l_id);
        } else {
          data->setValue(LOAD_KP2, atof(split_line[25].c_str()), l_id);
        }
      } 

      // LOAD_NP2
      if (nstr > 26) {
        if (!data->getValue(LOAD_NP2,&rval,l_id)) {
          data->addValue(LOAD_NP2, atof(split_line[26].c_str()), l_id);
        } else {
          data->setValue(LOAD_NP2, atof(split_line[26].c_str()), l_id);
        }
      }

      // LOAD_KQ2
      if (nstr > 27) {
        if (!data->getValue(LOAD_KQ2,&rval,l_id)) {
          data->addValue(LOAD_KQ2, atof(split_line[27].c_str()), l_id);
        } else {
          data->setValue(LOAD_KQ2, atof(split_line[27].c_str()), l_id);
        }
      }

      // LOAD_NQ2
      if (nstr > 28) {
        if (!data->getValue(LOAD_NQ2,&rval,l_id)) {
          data->addValue(LOAD_NQ2, atof(split_line[28].c_str()), l_id);
        } else {
          data->setValue(LOAD_NQ2, atof(split_line[28].c_str()), l_id);
        }
      }

      // LOAD_VBRK
      if (nstr > 29) {
        if (!data->getValue(LOAD_VBRK,&rval,l_id)) {
          data->addValue(LOAD_VBRK, atof(split_line[29].c_str()), l_id);
        } else {
          data->setValue(LOAD_VBRK, atof(split_line[29].c_str()), l_id);
        }
      }

      // LOAD_FRST
      if (nstr > 30) {
        if (!data->getValue(LOAD_FRST,&rval,l_id)) {
          data->addValue(LOAD_FRST, atof(split_line[30].c_str()), l_id);
        } else {
          data->setValue(LOAD_FRST, atof(split_line[30].c_str()), l_id);
        }
      }

      // LOAD_VRST
      if (nstr > 31) {
        if (!data->getValue(LOAD_VRST,&rval,l_id)) {
          data->addValue(LOAD_VRST, atof(split_line[31].c_str()), l_id);
        } else {
          data->setValue(LOAD_VRST, atof(split_line[31].c_str()), l_id);
        }
      }

      // LOAD_CMPKPF
      if (nstr > 32) {
        if (!data->getValue(LOAD_CMPKPF,&rval,l_id)) {
          data->addValue(LOAD_CMPKPF, atof(split_line[32].c_str()), l_id);
        } else {
          data->setValue(LOAD_CMPKPF, atof(split_line[32].c_str()), l_id);
        }
      }

      // LOAD_CMPKQF
      if (nstr > 33) {
        if (!data->getValue(LOAD_CMPKQF,&rval,l_id)) {
          data->addValue(LOAD_CMPKQF, atof(split_line[33].c_str()), l_id);
        } else {
          data->setValue(LOAD_CMPKQF, atof(split_line[33].c_str()), l_id);
        }
      }

      // LOAD_VC1OFF
      if (nstr > 34) {
        if (!data->getValue(LOAD_VC1OFF,&rval,l_id)) {
          data->addValue(LOAD_VC1OFF, atof(split_line[34].c_str()), l_id);
        } else {
          data->setValue(LOAD_VC1OFF, atof(split_line[34].c_str()), l_id);
        }
      }

      // LOAD_VC2OFF
      if (nstr > 35) {
        if (!data->getValue(LOAD_VC2OFF,&rval,l_id)) {
          data->addValue(LOAD_VC2OFF, atof(split_line[35].c_str()), l_id);
        } else {
          data->setValue(LOAD_VC2OFF, atof(split_line[35].c_str()), l_id);
        }
      }

      // LOAD_VC1ON
      if (nstr > 36) {
        if (!data->getValue(LOAD_VC1ON,&rval,l_id)) {
          data->addValue(LOAD_VC1ON, atof(split_line[36].c_str()), l_id);
        } else {
          data->setValue(LOAD_VC1ON, atof(split_line[36].c_str()), l_id);
        }
      }

      // LOAD_VC2ON
      if (nstr > 37) {
        if (!data->getValue(LOAD_VC2ON,&rval,l_id)) {
          data->addValue(LOAD_VC2ON, atof(split_line[37].c_str()), l_id);
        } else {
          data->setValue(LOAD_VC2ON, atof(split_line[37].c_str()), l_id);
        }
      }

      // LOAD_TTH
      if (nstr > 38) {
        if (!data->getValue(LOAD_TTH,&rval,l_id)) {
          data->addValue(LOAD_TTH, atof(split_line[38].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTH, atof(split_line[38].c_str()), l_id);
        }
      }

      // LOAD_TH1T
      if (nstr > 39) {
        if (!data->getValue(LOAD_TH1T,&rval,l_id)) {
          data->addValue(LOAD_TH1T, atof(split_line[39].c_str()), l_id);
        } else {
          data->setValue(LOAD_TH1T, atof(split_line[39].c_str()), l_id);
        }
      }

      // LOAD_TH2T
      if (nstr > 40) {
        if (!data->getValue(LOAD_TH2T,&rval,l_id)) {
          data->addValue(LOAD_TH2T, atof(split_line[40].c_str()), l_id);
        } else {
          data->setValue(LOAD_TH2T, atof(split_line[40].c_str()), l_id);
        }
      }

      // LOAD_FUVR
      if (nstr > 41) {
        if (!data->getValue(LOAD_FUVR,&rval,l_id)) {
          data->addValue(LOAD_FUVR, atof(split_line[41].c_str()), l_id);
        } else {
          data->setValue(LOAD_FUVR, atof(split_line[41].c_str()), l_id);
        }
      }

      // LOAD_UVTR1
      if (nstr > 42) {
        if (!data->getValue(LOAD_UVTR1,&rval,l_id)) {
          data->addValue(LOAD_UVTR1, atof(split_line[42].c_str()), l_id);
        } else {
          data->setValue(LOAD_UVTR1, atof(split_line[42].c_str()), l_id);
        }
      }

      // LOAD_TTR1
      if (nstr > 43) {
        if (!data->getValue(LOAD_TTR1,&rval,l_id)) {
          data->addValue(LOAD_TTR1, atof(split_line[43].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTR1, atof(split_line[43].c_str()), l_id);
        }
      }

      // LOAD_UVTR1
      if (nstr > 44) {
        if (!data->getValue(LOAD_UVTR1,&rval,l_id)) {
          data->addValue(LOAD_UVTR1, atof(split_line[44].c_str()), l_id);
        } else {
          data->setValue(LOAD_UVTR1, atof(split_line[44].c_str()), l_id);
        }
      }

      // LOAD_TTR2
      if (nstr > 45) {
        if (!data->getValue(LOAD_TTR2,&rval,l_id)) {
          data->addValue(LOAD_TTR2, atof(split_line[45].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTR2, atof(split_line[45].c_str()), l_id);
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
      sval = util.trimQuotes(split_line[3]);
      util.toUpper(sval);

      // LOAD_MODEL              "MODEL"                  integer
      strcpy(data.model, sval.c_str());
      int nstr = split_line.size();

      // LOAD_ID
      if (nstr > 2) {
        sval = util.clean2Char(split_line[2]);
        strcpy(data.id, sval.c_str());
      }

      // LOAD_TSTALL
      if (nstr > 11) {
        data.tstall = atof(split_line[11].c_str());
      }

      // LOAD_TRESTART
      if (nstr > 12) {
        data.trestart = atof(split_line[12].c_str());
      }

      // LOAD_TV
      if (nstr > 13) {
        data.tv = atof(split_line[13].c_str());
      }

      // LOAD_TF
      if (nstr > 14) {
        data.tf = atof(split_line[14].c_str());
      }

      // LOAD_COMPLF
      if (nstr > 15) {
        data.complf = atof(split_line[15].c_str());
      }

      // LOAD_COMPPF
      if (nstr > 16) {
        data.comppf = atof(split_line[16].c_str());
      }

      // LOAD_VSTALL
      if (nstr > 17) {
        data.vstall = atof(split_line[17].c_str());
      }

      // LOAD_RSTALL
      if (nstr > 18) {
        data.rstall = atof(split_line[18].c_str());
      }

      // LOAD_XSTALL
      if (nstr > 19) {
        data.xstall = atof(split_line[19].c_str());
      }

      // LOAD_LFADJ
      if (nstr > 20) {
        data.lfadj = atof(split_line[20].c_str());
      }

      // LOAD_KP1
      if (nstr > 21) {
        data.kp1 = atof(split_line[21].c_str());
      }

      // LOAD_NP1
      if (nstr > 22) {
        data.np1 = atof(split_line[22].c_str());
      }

      // LOAD_KQ1
      if (nstr > 23) {
        data.kq1 = atof(split_line[23].c_str());
      }

      // LOAD_NQ1
      if (nstr > 24) {
        data.nq1 = atof(split_line[24].c_str());
      }

      // LOAD_KP2
      if (nstr > 25) {
        data.kp2 = atof(split_line[25].c_str());
      }

      // LOAD_NP2
      if (nstr > 26) {
        data.np2 = atof(split_line[26].c_str());
      }

      // LOAD_KQ2
      if (nstr > 27) {
        data.kq2 = atof(split_line[27].c_str());
      }

      // LOAD_NQ2
      if (nstr > 28) {
        data.nq2 = atof(split_line[28].c_str());
      }

      // LOAD_VBRK
      if (nstr > 29) {
        data.vbrk = atof(split_line[29].c_str());
      }

      // LOAD_FRST
      if (nstr > 30) {
        data.frst = atof(split_line[30].c_str());
      }

      // LOAD_VRST
      if (nstr > 31) {
        data.vrst = atof(split_line[31].c_str());
      }

      // LOAD_CMPKPF
      if (nstr > 32) {
        data.cmpkpf = atof(split_line[32].c_str());
      }

      // LOAD_CMPKQF
      if (nstr > 33) {
        data.cmpkqf = atof(split_line[33].c_str());
      }

      // LOAD_VC1OFF
      if (nstr > 34) {
        data.vc1off = atof(split_line[34].c_str());
      }

      // LOAD_VC2OFF
      if (nstr > 35) {
        data.vc2off = atof(split_line[35].c_str());
      }

      // LOAD_VC1ON
      if (nstr > 36) {
        data.vc1on = atof(split_line[36].c_str());
      }

      // LOAD_VC2ON
      if (nstr > 37) {
        data.vc2on = atof(split_line[37].c_str());
      }

      // LOAD_TTH
      if (nstr > 38) {
        data.tth = atof(split_line[38].c_str());
      }

      // LOAD_TH1T
      if (nstr > 39) {
        data.th1t = atof(split_line[39].c_str());
      }

      // LOAD_TH2T
      if (nstr > 40) {
        data.th2t = atof(split_line[40].c_str());
      }

      // LOAD_FUVR
      if (nstr > 41) {
        data.fuvr = atof(split_line[41].c_str());
      }

      // LOAD_UVTR1
      if (nstr > 42) {
        data.uvtr1 = atof(split_line[42].c_str());
      }

      // LOAD_TTR1
      if (nstr > 43) {
        data.ttr1 = atof(split_line[43].c_str());
      }

      // LOAD_UVTR2
      if (nstr > 44) {
        data.uvtr2 = atof(split_line[44].c_str());
      }

      // LOAD_TTR2
      if (nstr > 45) {
        data.ttr2 = atof(split_line[45].c_str());
      }
    }
};
}  // parser
}  // gridpack
#endif
