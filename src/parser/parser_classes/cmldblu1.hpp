/*
 *    Copyright (c) 2013 Battelle Memorial Institute
 *    Licensed under modified BSD License. A copy of this license can be found
 *    in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: September 8, 2016
 *      Author: Bruce Palmer
 */
#ifndef CMLDBLU1_HPP
#define CMLDBLU1_HPP
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/string_utils.hpp"
namespace gridpack {
namespace parser {
template <class _data_struct> class Cmldblu1Parser
{
  public:
    /**
     * Constructor
     */
    explicit Cmldblu1Parser()
    {
    }

    /**
     * Destructor
     */
    virtual ~Cmldblu1Parser()
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

      // LOAD_MVA
      if (!data->getValue(LOAD_MVA,&rval,l_id)) {
        data->addValue(LOAD_MVA, data_struct.mva, l_id);
      } else {
        data->setValue(LOAD_MVA, data_struct.mva, l_id);
      }

      // LOAD_BSS
      if (!data->getValue(LOAD_BSS,&rval,l_id)) {
        data->addValue(LOAD_BSS, data_struct.bss, l_id);
      } else {
        data->setValue(LOAD_BSS, data_struct.bss, l_id);
      }

      // LOAD_RFDR
      if (!data->getValue(LOAD_RFDR,&rval,l_id)) {
        data->addValue(LOAD_RFDR, data_struct.rfdr, l_id);
      } else {
        data->setValue(LOAD_RFDR, data_struct.rfdr, l_id);
      }

      // LOAD_XFDR
      if (!data->getValue(LOAD_XFDR,&rval,l_id)) {
        data->addValue(LOAD_XFDR, data_struct.xfdr, l_id);
      } else {
        data->setValue(LOAD_XFDR, data_struct.xfdr, l_id);
      }

      // LOAD_FB
      if (!data->getValue(LOAD_FB,&rval,l_id)) {
        data->addValue(LOAD_FB, data_struct.fb, l_id);
      } else {
        data->setValue(LOAD_FB, data_struct.fb, l_id);
      }

      // LOAD_XXF
      if (!data->getValue(LOAD_XXF,&rval,l_id)) {
        data->addValue(LOAD_XXF, data_struct.xxf, l_id);
      } else {
        data->setValue(LOAD_XXF, data_struct.xxf, l_id);
      }

      // LOAD_TFIXHS
      if (!data->getValue(LOAD_TFIXHS,&rval,l_id)) {
        data->addValue(LOAD_TFIXHS, data_struct.tfixhs, l_id);
      } else {
        data->setValue(LOAD_TFIXHS, data_struct.tfixhs, l_id);
      }

      // LOAD_TFIXLS
      if (!data->getValue(LOAD_TFIXLS,&rval,l_id)) {
        data->addValue(LOAD_TFIXLS, data_struct.tfixls, l_id);
      } else {
        data->setValue(LOAD_TFIXLS, data_struct.tfixls, l_id);
      }

      // LOAD_LTC
      if (!data->getValue(LOAD_LTC,&rval,l_id)) {
        data->addValue(LOAD_LTC, data_struct.ltc, l_id);
      } else {
        data->setValue(LOAD_LTC, data_struct.ltc, l_id);
      }

      // LOAD_TMIN
      if (!data->getValue(LOAD_TMIN,&rval,l_id)) {
        data->addValue(LOAD_TMIN, data_struct.tmin, l_id);
      } else {
        data->setValue(LOAD_TMIN, data_struct.tmin, l_id);
      }

      // LOAD_TMAX
      if (!data->getValue(LOAD_TMAX,&rval,l_id)) {
        data->addValue(LOAD_TMAX, data_struct.tmax, l_id);
      } else {
        data->setValue(LOAD_TMAX, data_struct.tmax, l_id);
      }

      // LOAD_STEP
      if (!data->getValue(LOAD_STEP,&rval,l_id)) {
        data->addValue(LOAD_STEP, data_struct.step, l_id);
      } else {
        data->setValue(LOAD_STEP, data_struct.step, l_id);
      }

      // LOAD_VMIN
      if (!data->getValue(LOAD_VMIN,&rval,l_id)) {
        data->addValue(LOAD_VMIN, data_struct.vmin, l_id);
      } else {
        data->setValue(LOAD_VMIN, data_struct.vmin, l_id);
      }

      // LOAD_VMAX
      if (!data->getValue(LOAD_VMAX,&rval,l_id)) {
        data->addValue(LOAD_VMAX, data_struct.vmax, l_id);
      } else {
        data->setValue(LOAD_VMAX, data_struct.vmax, l_id);
      }

      // LOAD_TDEL
      if (!data->getValue(LOAD_TDEL,&rval,l_id)) {
        data->addValue(LOAD_TDEL, data_struct.tdel, l_id);
      } else {
        data->setValue(LOAD_TDEL, data_struct.tdel, l_id);
      }

      // LOAD_TTAP
      if (!data->getValue(LOAD_TTAP,&rval,l_id)) {
        data->addValue(LOAD_TTAP, data_struct.ttap, l_id);
      } else {
        data->setValue(LOAD_TTAP, data_struct.ttap, l_id);
      }

      // LOAD_RCOMP
      if (!data->getValue(LOAD_RCOMP,&rval,l_id)) {
        data->addValue(LOAD_RCOMP, data_struct.rcomp, l_id);
      } else {
        data->setValue(LOAD_RCOMP, data_struct.rcomp, l_id);
      }

      // LOAD_XCOMP
      if (!data->getValue(LOAD_XCOMP,&rval,l_id)) {
        data->addValue(LOAD_XCOMP, data_struct.xcomp, l_id);
      } else {
        data->setValue(LOAD_XCOMP, data_struct.xcomp, l_id);
      }

      // LOAD_FMA
      if (!data->getValue(LOAD_FMA,&rval,l_id)) {
        data->addValue(LOAD_FMA, data_struct.fma, l_id);
      } else {
        data->setValue(LOAD_FMA, data_struct.fma, l_id);
      }

      // LOAD_FMB
      if (!data->getValue(LOAD_FMB,&rval,l_id)) {
        data->addValue(LOAD_FMB, data_struct.fmb, l_id);
      } else {
        data->setValue(LOAD_FMB, data_struct.fmb, l_id);
      }

      // LOAD_FMC
      if (!data->getValue(LOAD_FMC,&rval,l_id)) {
        data->addValue(LOAD_FMC, data_struct.fmc, l_id);
      } else {
        data->setValue(LOAD_FMC, data_struct.fmc, l_id);
      }

      // LOAD_FMD
      if (!data->getValue(LOAD_FMD,&rval,l_id)) {
        data->addValue(LOAD_FMD, data_struct.fmd, l_id);
      } else {
        data->setValue(LOAD_FMD, data_struct.fmd, l_id);
      }

      // LOAD_FEL
      if (!data->getValue(LOAD_FEL,&rval,l_id)) {
        data->addValue(LOAD_FEL, data_struct.fel, l_id);
      } else {
        data->setValue(LOAD_FEL, data_struct.fel, l_id);
      }

      // LOAD_PFEL
      if (!data->getValue(LOAD_PFEL,&rval,l_id)) {
        data->addValue(LOAD_PFEL, data_struct.pfel, l_id);
      } else {
        data->setValue(LOAD_PFEL, data_struct.pfel, l_id);
      }

      // LOAD_VD1
      if (!data->getValue(LOAD_VD1,&rval,l_id)) {
        data->addValue(LOAD_VD1, data_struct.vd1, l_id);
      } else {
        data->setValue(LOAD_VD1, data_struct.vd1, l_id);
      }

      // LOAD_VD2
      if (!data->getValue(LOAD_VD2,&rval,l_id)) {
        data->addValue(LOAD_VD2, data_struct.vd2, l_id);
      } else {
        data->setValue(LOAD_VD2, data_struct.vd2, l_id);
      }

      // LOAD_PFS
      if (!data->getValue(LOAD_PFS,&rval,l_id)) {
        data->addValue(LOAD_PFS, data_struct.pfs, l_id);
      } else {
        data->setValue(LOAD_PFS, data_struct.pfs, l_id);
      }

      // LOAD_P1E
      if (!data->getValue(LOAD_P1E,&rval,l_id)) {
        data->addValue(LOAD_P1E, data_struct.p1e, l_id);
      } else {
        data->setValue(LOAD_P1E, data_struct.p1e, l_id);
      }

      // LOAD_P1C
      if (!data->getValue(LOAD_P1C,&rval,l_id)) {
        data->addValue(LOAD_P1C, data_struct.p1c, l_id);
      } else {
        data->setValue(LOAD_P1C, data_struct.p1c, l_id);
      }

      // LOAD_P2E
      if (!data->getValue(LOAD_P2E,&rval,l_id)) {
        data->addValue(LOAD_P2E, data_struct.p2e, l_id);
      } else {
        data->setValue(LOAD_P2E, data_struct.p2e, l_id);
      }

      // LOAD_P2C
      if (!data->getValue(LOAD_P2C,&rval,l_id)) {
        data->addValue(LOAD_P2C, data_struct.p2c, l_id);
      } else {
        data->setValue(LOAD_P2C, data_struct.p2c, l_id);
      }

      // LOAD_PFREQ
      if (!data->getValue(LOAD_PFREQ,&rval,l_id)) {
        data->addValue(LOAD_PFREQ, data_struct.pfreq, l_id);
      } else {
        data->setValue(LOAD_PFREQ, data_struct.pfreq, l_id);
      }

      // LOAD_Q1E
      if (!data->getValue(LOAD_Q1E,&rval,l_id)) {
        data->addValue(LOAD_Q1E, data_struct.q1e, l_id);
      } else {
        data->setValue(LOAD_Q1E, data_struct.q1e, l_id);
      }

      // LOAD_Q1C
      if (!data->getValue(LOAD_Q1C,&rval,l_id)) {
        data->addValue(LOAD_Q1C, data_struct.q1c, l_id);
      } else {
        data->setValue(LOAD_Q1C, data_struct.q1c, l_id);
      }

      // LOAD_Q2E
      if (!data->getValue(LOAD_Q2E,&rval,l_id)) {
        data->addValue(LOAD_Q2E, data_struct.q2e, l_id);
      } else {
        data->setValue(LOAD_Q2E, data_struct.q2e, l_id);
      }

      // LOAD_Q2C
      if (!data->getValue(LOAD_Q2C,&rval,l_id)) {
        data->addValue(LOAD_Q2C, data_struct.q2c, l_id);
      } else {
        data->setValue(LOAD_Q2C, data_struct.q2c, l_id);
      }

      // LOAD_QFREQ
      if (!data->getValue(LOAD_QFREQ,&rval,l_id)) {
        data->addValue(LOAD_QFREQ, data_struct.qfreq, l_id);
      } else {
        data->setValue(LOAD_QFREQ, data_struct.qfreq, l_id);
      }

      // LOAD_MTPA
      if (!data->getValue(LOAD_MTPA,&ival,l_id)) {
        data->addValue(LOAD_MTPA, data_struct.mtpa, l_id);
      } else {
        data->setValue(LOAD_MTPA, data_struct.mtpa, l_id);
      }

      // LOAD_LFMA
      if (!data->getValue(LOAD_LFMA,&rval,l_id)) {
        data->addValue(LOAD_LFMA, data_struct.lfma, l_id);
      } else {
        data->setValue(LOAD_LFMA, data_struct.lfma, l_id);
      }

      // LOAD_RSA
      if (!data->getValue(LOAD_RSA,&rval,l_id)) {
        data->addValue(LOAD_RSA, data_struct.rsa, l_id);
      } else {
        data->setValue(LOAD_RSA, data_struct.rsa, l_id);
      }

      // LOAD_LSA
      if (!data->getValue(LOAD_LSA,&rval,l_id)) {
        data->addValue(LOAD_LSA, data_struct.lsa, l_id);
      } else {
        data->setValue(LOAD_LSA, data_struct.lsa, l_id);
      }

      // LOAD_LPA
      if (!data->getValue(LOAD_LPA,&rval,l_id)) {
        data->addValue(LOAD_LPA, data_struct.lpa, l_id);
      } else {
        data->setValue(LOAD_LPA, data_struct.lpa, l_id);
      }

      // LOAD_LPPA
      if (!data->getValue(LOAD_LPPA,&rval,l_id)) {
        data->addValue(LOAD_LPPA, data_struct.lppa, l_id);
      } else {
        data->setValue(LOAD_LPPA, data_struct.lppa, l_id);
      }

      // LOAD_TPOA
      if (!data->getValue(LOAD_TPOA,&rval,l_id)) {
        data->addValue(LOAD_TPOA, data_struct.tpoa, l_id);
      } else {
        data->setValue(LOAD_TPOA, data_struct.tpoa, l_id);
      }

      // LOAD_TPPOA
      if (!data->getValue(LOAD_TPPOA,&rval,l_id)) {
        data->addValue(LOAD_TPPOA, data_struct.tppoa, l_id);
      } else {
        data->setValue(LOAD_TPPOA, data_struct.tppoa, l_id);
      }

      // LOAD_HA
      if (!data->getValue(LOAD_HA,&rval,l_id)) {
        data->addValue(LOAD_HA, data_struct.ha, l_id);
      } else {
        data->setValue(LOAD_HA, data_struct.ha, l_id);
      }

      // LOAD_ETRQA
      if (!data->getValue(LOAD_ETRQA,&rval,l_id)) {
        data->addValue(LOAD_ETRQA, data_struct.etrqa, l_id);
      } else {
        data->setValue(LOAD_ETRQA, data_struct.etrqa, l_id);
      }

      // LOAD_VTR1A
      if (!data->getValue(LOAD_VTR1A,&rval,l_id)) {
        data->addValue(LOAD_VTR1A, data_struct.vtr1a, l_id);
      } else {
        data->setValue(LOAD_VTR1A, data_struct.vtr1a, l_id);
      }

      // LOAD_TTR1A
      if (!data->getValue(LOAD_TTR1A,&rval,l_id)) {
        data->addValue(LOAD_TTR1A, data_struct.ttr1a, l_id);
      } else {
        data->setValue(LOAD_TTR1A, data_struct.ttr1a, l_id);
      }

      // LOAD_FTR1A
      if (!data->getValue(LOAD_FTR1A,&rval,l_id)) {
        data->addValue(LOAD_FTR1A, data_struct.ftr1a, l_id);
      } else {
        data->setValue(LOAD_FTR1A, data_struct.ftr1a, l_id);
      }

      // LOAD_VRC1A
      if (!data->getValue(LOAD_VRC1A,&rval,l_id)) {
        data->addValue(LOAD_VRC1A, data_struct.vrc1a, l_id);
      } else {
        data->setValue(LOAD_VRC1A, data_struct.vrc1a, l_id);
      }

      // LOAD_TRC1A
      if (!data->getValue(LOAD_TRC1A,&rval,l_id)) {
        data->addValue(LOAD_TRC1A, data_struct.trc1a, l_id);
      } else {
        data->setValue(LOAD_TRC1A, data_struct.trc1a, l_id);
      }

      // LOAD_VTR2A
      if (!data->getValue(LOAD_VTR2A,&rval,l_id)) {
        data->addValue(LOAD_VTR2A, data_struct.vtr2a, l_id);
      } else {
        data->setValue(LOAD_VTR2A, data_struct.vtr2a, l_id);
      }

      // LOAD_TTR2A
      if (!data->getValue(LOAD_TTR2A,&rval,l_id)) {
        data->addValue(LOAD_TTR2A, data_struct.ttr2a, l_id);
      } else {
        data->setValue(LOAD_TTR2A, data_struct.ttr2a, l_id);
      }

      // LOAD_FTR2A
      if (!data->getValue(LOAD_FTR2A,&rval,l_id)) {
        data->addValue(LOAD_FTR2A, data_struct.ftr2a, l_id);
      } else {
        data->setValue(LOAD_FTR2A, data_struct.ftr2a, l_id);
      }

      // LOAD_VRC2A
      if (!data->getValue(LOAD_VRC2A,&rval,l_id)) {
        data->addValue(LOAD_VRC2A, data_struct.vrc2a, l_id);
      } else {
        data->setValue(LOAD_VRC2A, data_struct.vrc2a, l_id);
      }

      // LOAD_TRC2A
      if (!data->getValue(LOAD_TRC2A,&rval,l_id)) {
        data->addValue(LOAD_TRC2A, data_struct.trc2a, l_id);
      } else {
        data->setValue(LOAD_TRC2A, data_struct.trc2a, l_id);
      }

      // LOAD_MTPB
      if (!data->getValue(LOAD_MTPB,&ival,l_id)) {
        data->addValue(LOAD_MTPB, data_struct.mtpb, l_id);
      } else {
        data->setValue(LOAD_MTPB, data_struct.mtpb, l_id);
      }

      // LOAD_LFMB
      if (!data->getValue(LOAD_LFMB,&rval,l_id)) {
        data->addValue(LOAD_LFMB, data_struct.lfmb, l_id);
      } else {
        data->setValue(LOAD_LFMB, data_struct.lfmb, l_id);
      }

      // LOAD_RSB
      if (!data->getValue(LOAD_RSB,&rval,l_id)) {
        data->addValue(LOAD_RSB, data_struct.rsb, l_id);
      } else {
        data->setValue(LOAD_RSB, data_struct.rsb, l_id);
      }

      // LOAD_LSB
      if (!data->getValue(LOAD_LSB,&rval,l_id)) {
        data->addValue(LOAD_LSB, data_struct.lsb, l_id);
      } else {
        data->setValue(LOAD_LSB, data_struct.lsb, l_id);
      }

      // LOAD_LPB
      if (!data->getValue(LOAD_LPB,&rval,l_id)) {
        data->addValue(LOAD_LPB, data_struct.lpb, l_id);
      } else {
        data->setValue(LOAD_LPB, data_struct.lpb, l_id);
      }

      // LOAD_LPPB
      if (!data->getValue(LOAD_LPPB,&rval,l_id)) {
        data->addValue(LOAD_LPPB, data_struct.lppb, l_id);
      } else {
        data->setValue(LOAD_LPPB, data_struct.lppb, l_id);
      }

      // LOAD_TPOB
      if (!data->getValue(LOAD_TPOB,&rval,l_id)) {
        data->addValue(LOAD_TPOB, data_struct.tpob, l_id);
      } else {
        data->setValue(LOAD_TPOB, data_struct.tpob, l_id);
      }

      // LOAD_TPPOB
      if (!data->getValue(LOAD_TPPOB,&rval,l_id)) {
        data->addValue(LOAD_TPPOB, data_struct.tppob, l_id);
      } else {
        data->setValue(LOAD_TPPOB, data_struct.tppob, l_id);
      }

      // LOAD_HB
      if (!data->getValue(LOAD_HB,&rval,l_id)) {
        data->addValue(LOAD_HB, data_struct.hb, l_id);
      } else {
        data->setValue(LOAD_HB, data_struct.hb, l_id);
      }

      // LOAD_ETRQB
      if (!data->getValue(LOAD_ETRQB,&rval,l_id)) {
        data->addValue(LOAD_ETRQB, data_struct.etrqb, l_id);
      } else {
        data->setValue(LOAD_ETRQB, data_struct.etrqb, l_id);
      }

      // LOAD_VTR1B
      if (!data->getValue(LOAD_VTR1B,&rval,l_id)) {
        data->addValue(LOAD_VTR1B, data_struct.vtr1b, l_id);
      } else {
        data->setValue(LOAD_VTR1B, data_struct.vtr1b, l_id);
      }

      // LOAD_TTR1B
      if (!data->getValue(LOAD_TTR1B,&rval,l_id)) {
        data->addValue(LOAD_TTR1B, data_struct.ttr1b, l_id);
      } else {
        data->setValue(LOAD_TTR1B, data_struct.ttr1b, l_id);
      }

      // LOAD_FTR1B
      if (!data->getValue(LOAD_FTR1B,&rval,l_id)) {
        data->addValue(LOAD_FTR1B, data_struct.ftr1b, l_id);
      } else {
        data->setValue(LOAD_FTR1B, data_struct.ftr1b, l_id);
      }

      // LOAD_VRC1B
      if (!data->getValue(LOAD_VRC1B,&rval,l_id)) {
        data->addValue(LOAD_VRC1B, data_struct.vrc1b, l_id);
      } else {
        data->setValue(LOAD_VRC1B, data_struct.vrc1b, l_id);
      }

      // LOAD_TRC1B
      if (!data->getValue(LOAD_TRC1B,&rval,l_id)) {
        data->addValue(LOAD_TRC1B, data_struct.trc1b, l_id);
      } else {
        data->setValue(LOAD_TRC1B, data_struct.trc1b, l_id);
      }

      // LOAD_VTR2B
      if (!data->getValue(LOAD_VTR2B,&rval,l_id)) {
        data->addValue(LOAD_VTR2B, data_struct.vtr2b, l_id);
      } else {
        data->setValue(LOAD_VTR2B, data_struct.vtr2b, l_id);
      }

      // LOAD_TTR2B
      if (!data->getValue(LOAD_TTR2B,&rval,l_id)) {
        data->addValue(LOAD_TTR2B, data_struct.ttr2b, l_id);
      } else {
        data->setValue(LOAD_TTR2B, data_struct.ttr2b, l_id);
      }

      // LOAD_FTR2B
      if (!data->getValue(LOAD_FTR2B,&rval,l_id)) {
        data->addValue(LOAD_FTR2B, data_struct.ftr2b, l_id);
      } else {
        data->setValue(LOAD_FTR2B, data_struct.ftr2b, l_id);
      }

      // LOAD_VRC2B
      if (!data->getValue(LOAD_VRC2B,&rval,l_id)) {
        data->addValue(LOAD_VRC2B, data_struct.vrc2b, l_id);
      } else {
        data->setValue(LOAD_VRC2B, data_struct.vrc2b, l_id);
      }

      // LOAD_TRC2B
      if (!data->getValue(LOAD_TRC2B,&rval,l_id)) {
        data->addValue(LOAD_TRC2B, data_struct.trc2b, l_id);
      } else {
        data->setValue(LOAD_TRC2B, data_struct.trc2b, l_id);
      }

      // LOAD_MTPC
      if (!data->getValue(LOAD_MTPC,&ival,l_id)) {
        data->addValue(LOAD_MTPC, data_struct.mtpc, l_id);
      } else {
        data->setValue(LOAD_MTPC, data_struct.mtpc, l_id);
      }

      // LOAD_LFMC
      if (!data->getValue(LOAD_LFMC,&rval,l_id)) {
        data->addValue(LOAD_LFMC, data_struct.lfmc, l_id);
      } else {
        data->setValue(LOAD_LFMC, data_struct.lfmc, l_id);
      }

      // LOAD_RSC
      if (!data->getValue(LOAD_RSC,&rval,l_id)) {
        data->addValue(LOAD_RSC, data_struct.rsc, l_id);
      } else {
        data->setValue(LOAD_RSC, data_struct.rsc, l_id);
      }

      // LOAD_LSC
      if (!data->getValue(LOAD_LSC,&rval,l_id)) {
        data->addValue(LOAD_LSC, data_struct.lsc, l_id);
      } else {
        data->setValue(LOAD_LSC, data_struct.lsc, l_id);
      }

      // LOAD_LPC
      if (!data->getValue(LOAD_LPC,&rval,l_id)) {
        data->addValue(LOAD_LPC, data_struct.lpc, l_id);
      } else {
        data->setValue(LOAD_LPC, data_struct.lpc, l_id);
      }

      // LOAD_LPPC
      if (!data->getValue(LOAD_LPPC,&rval,l_id)) {
        data->addValue(LOAD_LPPC, data_struct.lppc, l_id);
      } else {
        data->setValue(LOAD_LPPC, data_struct.lppc, l_id);
      }

      // LOAD_TPOC
      if (!data->getValue(LOAD_TPOC,&rval,l_id)) {
        data->addValue(LOAD_TPOC, data_struct.tpoc, l_id);
      } else {
        data->setValue(LOAD_TPOC, data_struct.tpoc, l_id);
      }

      // LOAD_TPPOC
      if (!data->getValue(LOAD_TPPOC,&rval,l_id)) {
        data->addValue(LOAD_TPPOC, data_struct.tppoc, l_id);
      } else {
        data->setValue(LOAD_TPPOC, data_struct.tppoc, l_id);
      }

      // LOAD_HC
      if (!data->getValue(LOAD_HC,&rval,l_id)) {
        data->addValue(LOAD_HC, data_struct.hc, l_id);
      } else {
        data->setValue(LOAD_HC, data_struct.hc, l_id);
      }

      // LOAD_ETRQC
      if (!data->getValue(LOAD_ETRQC,&rval,l_id)) {
        data->addValue(LOAD_ETRQC, data_struct.etrqc, l_id);
      } else {
        data->setValue(LOAD_ETRQC, data_struct.etrqc, l_id);
      }

      // LOAD_VTR1C
      if (!data->getValue(LOAD_VTR1C,&rval,l_id)) {
        data->addValue(LOAD_VTR1C, data_struct.vtr1c, l_id);
      } else {
        data->setValue(LOAD_VTR1C, data_struct.vtr1c, l_id);
      }

      // LOAD_TTR1C
      if (!data->getValue(LOAD_TTR1C,&rval,l_id)) {
        data->addValue(LOAD_TTR1C, data_struct.ttr1c, l_id);
      } else {
        data->setValue(LOAD_TTR1C, data_struct.ttr1c, l_id);
      }

      // LOAD_FTR1C
      if (!data->getValue(LOAD_FTR1C,&rval,l_id)) {
        data->addValue(LOAD_FTR1C, data_struct.ftr1c, l_id);
      } else {
        data->setValue(LOAD_FTR1C, data_struct.ftr1c, l_id);
      }

      // LOAD_VRC1C
      if (!data->getValue(LOAD_VRC1C,&rval,l_id)) {
        data->addValue(LOAD_VRC1C, data_struct.vrc1c, l_id);
      } else {
        data->setValue(LOAD_VRC1C, data_struct.vrc1c, l_id);
      }

      // LOAD_TRC1C
      if (!data->getValue(LOAD_TRC1C,&rval,l_id)) {
        data->addValue(LOAD_TRC1C, data_struct.trc1c, l_id);
      } else {
        data->setValue(LOAD_TRC1C, data_struct.trc1c, l_id);
      }

      // LOAD_VTR2C
      if (!data->getValue(LOAD_VTR2C,&rval,l_id)) {
        data->addValue(LOAD_VTR2C, data_struct.vtr2c, l_id);
      } else {
        data->setValue(LOAD_VTR2C, data_struct.vtr2c, l_id);
      }

      // LOAD_TTR2C
      if (!data->getValue(LOAD_TTR2C,&rval,l_id)) {
        data->addValue(LOAD_TTR2C, data_struct.ttr2c, l_id);
      } else {
        data->setValue(LOAD_TTR2C, data_struct.ttr2c, l_id);
      }

      // LOAD_FTR2C
      if (!data->getValue(LOAD_FTR2C,&rval,l_id)) {
        data->addValue(LOAD_FTR2C, data_struct.ftr2c, l_id);
      } else {
        data->setValue(LOAD_FTR2C, data_struct.ftr2c, l_id);
      }

      // LOAD_VRC2C
      if (!data->getValue(LOAD_VRC2C,&rval,l_id)) {
        data->addValue(LOAD_VRC2C, data_struct.vrc2c, l_id);
      } else {
        data->setValue(LOAD_VRC2C, data_struct.vrc2c, l_id);
      }

      // LOAD_TRC2C
      if (!data->getValue(LOAD_TRC2C,&rval,l_id)) {
        data->addValue(LOAD_TRC2C, data_struct.trc2c, l_id);
      } else {
        data->setValue(LOAD_TRC2C, data_struct.trc2c, l_id);
      }

      // LOAD_TSTALL
      if (!data->getValue(LOAD_TSTALL,&rval,l_id)) {
        data->addValue(LOAD_TSTALL, data_struct.tstall, l_id);
      } else {
        data->setValue(LOAD_TSTALL, data_struct.tstall, l_id);
      }

      // LOAD_TRST
      if (!data->getValue(LOAD_TRST,&rval,l_id)) {
        data->addValue(LOAD_TRST, data_struct.trst, l_id);
      } else {
        data->setValue(LOAD_TRST, data_struct.trst, l_id);
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

      // LOAD_LFMD
      if (!data->getValue(LOAD_LFMD,&rval,l_id)) {
        data->addValue(LOAD_LFMD, data_struct.lfmd, l_id);
      } else {
        data->setValue(LOAD_LFMD, data_struct.lfmd, l_id);
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
        data->addValue(LOAD_NQ1, data_struct.nq1, l_id);
      } else {
        data->setValue(LOAD_NQ1, data_struct.nq1, l_id);
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
        data->addValue(LOAD_NQ2, data_struct.nq2, l_id);
      } else {
        data->setValue(LOAD_NQ2, data_struct.nq2, l_id);
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
        data->addValue(LOAD_TTR1, data_struct.ttr2, l_id);
      } else {
        data->setValue(LOAD_TTR1, data_struct.ttr2, l_id);
      }

      // LOAD_UVTR2
      if (!data->getValue(LOAD_UVTR2,&rval,l_id)) {
        data->addValue(LOAD_UVTR2, data_struct.uvtr2, l_id);
      } else {
        data->setValue(LOAD_VC2ON, data_struct.uvtr2, l_id);
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

      // LOAD_IT
      if (nstr > 5) {
        if (!data->getValue(LOAD_IT,&ival,l_id)) {
          data->addValue(LOAD_IT, atoi(split_line[5].c_str()), l_id);
        } else {
          data->setValue(LOAD_IT, atoi(split_line[5].c_str()), l_id);
        }
      } 

      // LOAD_MVA
      if (nstr > 11) {
        if (!data->getValue(LOAD_MVA,&rval,l_id)) {
          data->addValue(LOAD_MVA, atof(split_line[11].c_str()), l_id);
        } else {
          data->setValue(LOAD_MVA, atof(split_line[11].c_str()), l_id);
        }
      } 

      // LOAD_BSS
      if (nstr > 12) {
        if (!data->getValue(LOAD_BSS,&rval,l_id)) {
          data->addValue(LOAD_BSS, atof(split_line[12].c_str()), l_id);
        } else {
          data->setValue(LOAD_BSS, atof(split_line[12].c_str()), l_id);
        }
      } 

      // LOAD_RFDR
      if (nstr > 13) {
        if (!data->getValue(LOAD_RFDR,&rval,l_id)) {
          data->addValue(LOAD_RFDR, atof(split_line[13].c_str()), l_id);
        } else {
          data->setValue(LOAD_RFDR, atof(split_line[13].c_str()), l_id);
        }
      } 

      // LOAD_XFDR
      if (nstr > 14) {
        if (!data->getValue(LOAD_XFDR,&rval,l_id)) {
          data->addValue(LOAD_XFDR, atof(split_line[14].c_str()), l_id);
        } else {
          data->setValue(LOAD_XFDR, atof(split_line[14].c_str()), l_id);
        }
      } 

      // LOAD_FB
      if (nstr > 15) {
        if (!data->getValue(LOAD_FB,&rval,l_id)) {
          data->addValue(LOAD_FB, atof(split_line[15].c_str()), l_id);
        } else {
          data->setValue(LOAD_FB, atof(split_line[15].c_str()), l_id);
        }
      } 

      // LOAD_XXF
      if (nstr > 16) {
        if (!data->getValue(LOAD_XXF,&rval,l_id)) {
          data->addValue(LOAD_XXF, atof(split_line[16].c_str()), l_id);
        } else {
          data->setValue(LOAD_XXF, atof(split_line[16].c_str()), l_id);
        }
      } 

      // LOAD_TFIXHS
      if (nstr > 17) {
        if (!data->getValue(LOAD_TFIXHS,&rval,l_id)) {
          data->addValue(LOAD_TFIXHS, atof(split_line[17].c_str()), l_id);
        } else {
          data->setValue(LOAD_TFIXHS, atof(split_line[17].c_str()), l_id);
        }
      }

      // LOAD_TFIXLS
      if (nstr > 18) {
        if (!data->getValue(LOAD_TFIXLS,&rval,l_id)) {
          data->addValue(LOAD_TFIXLS, atof(split_line[18].c_str()), l_id);
        } else {
          data->setValue(LOAD_TFIXLS, atof(split_line[18].c_str()), l_id);
        }
      }

      // LOAD_LTC
      if (nstr > 19) {
        if (!data->getValue(LOAD_LTC,&rval,l_id)) {
          data->addValue(LOAD_LTC, atof(split_line[19].c_str()), l_id);
        } else {
          data->setValue(LOAD_LTC, atof(split_line[19].c_str()), l_id);
        }
      }

      // LOAD_TMIN
      if (nstr > 20) {
        if (!data->getValue(LOAD_TMIN,&rval,l_id)) {
          data->addValue(LOAD_TMIN, atof(split_line[20].c_str()), l_id);
        } else {
          data->setValue(LOAD_TMIN, atof(split_line[20].c_str()), l_id);
        }
      }

      // LOAD_TMAX
      if (nstr > 21) {
        if (!data->getValue(LOAD_TMAX,&rval,l_id)) {
          data->addValue(LOAD_TMAX, atof(split_line[21].c_str()), l_id);
        } else {
          data->setValue(LOAD_TMAX, atof(split_line[21].c_str()), l_id);
        }
      }

      // LOAD_STEP
      if (nstr > 22) {
        if (!data->getValue(LOAD_STEP,&rval,l_id)) {
          data->addValue(LOAD_STEP, atof(split_line[22].c_str()), l_id);
        } else {
          data->setValue(LOAD_STEP, atof(split_line[22].c_str()), l_id);
        }
      }

      // LOAD_VMIN
      if (nstr > 23) {
        if (!data->getValue(LOAD_VMIN,&rval,l_id)) {
          data->addValue(LOAD_VMIN, atof(split_line[23].c_str()), l_id);
        } else {
          data->setValue(LOAD_VMIN, atof(split_line[23].c_str()), l_id);
        }
      }

      // LOAD_VMAX
      if (nstr > 24) {
        if (!data->getValue(LOAD_VMAX,&rval,l_id)) {
          data->addValue(LOAD_VMAX, atof(split_line[24].c_str()), l_id);
        } else {
          data->setValue(LOAD_VMAX, atof(split_line[24].c_str()), l_id);
        }
      } 

      // LOAD_TDEL
      if (nstr > 25) {
        if (!data->getValue(LOAD_TDEL,&rval,l_id)) {
          data->addValue(LOAD_TDEL, atof(split_line[25].c_str()), l_id);
        } else {
          data->setValue(LOAD_TDEL, atof(split_line[25].c_str()), l_id);
        }
      } 

      // LOAD_TTAP
      if (nstr > 26) {
        if (!data->getValue(LOAD_TTAP,&rval,l_id)) {
          data->addValue(LOAD_TTAP, atof(split_line[26].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTAP, atof(split_line[26].c_str()), l_id);
        }
      }

      // LOAD_RCOMP
      if (nstr > 27) {
        if (!data->getValue(LOAD_RCOMP,&rval,l_id)) {
          data->addValue(LOAD_RCOMP, atof(split_line[27].c_str()), l_id);
        } else {
          data->setValue(LOAD_RCOMP, atof(split_line[27].c_str()), l_id);
        }
      }

      // LOAD_XCOMP
      if (nstr > 28) {
        if (!data->getValue(LOAD_XCOMP,&rval,l_id)) {
          data->addValue(LOAD_XCOMP, atof(split_line[28].c_str()), l_id);
        } else {
          data->setValue(LOAD_XCOMP, atof(split_line[28].c_str()), l_id);
        }
      }

      // LOAD_FMA
      if (nstr > 29) {
        if (!data->getValue(LOAD_FMA,&rval,l_id)) {
          data->addValue(LOAD_FMA, atof(split_line[29].c_str()), l_id);
        } else {
          data->setValue(LOAD_FMA, atof(split_line[29].c_str()), l_id);
        }
      }

      // LOAD_FMB
      if (nstr > 30) {
        if (!data->getValue(LOAD_FMB,&rval,l_id)) {
          data->addValue(LOAD_FMB, atof(split_line[30].c_str()), l_id);
        } else {
          data->setValue(LOAD_FMB, atof(split_line[30].c_str()), l_id);
        }
      }

      // LOAD_FMC
      if (nstr > 31) {
        if (!data->getValue(LOAD_FMC,&rval,l_id)) {
          data->addValue(LOAD_FMC, atof(split_line[31].c_str()), l_id);
        } else {
          data->setValue(LOAD_FMC, atof(split_line[31].c_str()), l_id);
        }
      }

      // LOAD_FMD
      if (nstr > 32) {
        if (!data->getValue(LOAD_FMD,&rval,l_id)) {
          data->addValue(LOAD_FMD, atof(split_line[32].c_str()), l_id);
        } else {
          data->setValue(LOAD_FMD, atof(split_line[32].c_str()), l_id);
        }
      }

      // LOAD_FEL
      if (nstr > 33) {
        if (!data->getValue(LOAD_FEL,&rval,l_id)) {
          data->addValue(LOAD_FEL, atof(split_line[33].c_str()), l_id);
        } else {
          data->setValue(LOAD_FEL, atof(split_line[33].c_str()), l_id);
        }
      }

      // LOAD_PFEL
      if (nstr > 34) {
        if (!data->getValue(LOAD_PFEL,&rval,l_id)) {
          data->addValue(LOAD_PFEL, atof(split_line[34].c_str()), l_id);
        } else {
          data->setValue(LOAD_PFEL, atof(split_line[34].c_str()), l_id);
        }
      }

      // LOAD_VD1
      if (nstr > 35) {
        if (!data->getValue(LOAD_VD1,&rval,l_id)) {
          data->addValue(LOAD_VD1, atof(split_line[35].c_str()), l_id);
        } else {
          data->setValue(LOAD_VD1, atof(split_line[35].c_str()), l_id);
        }
      }

      // LOAD_VD2
      if (nstr > 36) {
        if (!data->getValue(LOAD_VD2,&rval,l_id)) {
          data->addValue(LOAD_VD2, atof(split_line[36].c_str()), l_id);
        } else {
          data->setValue(LOAD_VD2, atof(split_line[36].c_str()), l_id);
        }
      }

      // LOAD_PFS
      if (nstr > 37) {
        if (!data->getValue(LOAD_PFS,&rval,l_id)) {
          data->addValue(LOAD_PFS, atof(split_line[37].c_str()), l_id);
        } else {
          data->setValue(LOAD_PFS, atof(split_line[37].c_str()), l_id);
        }
      }

      // LOAD_P1E
      if (nstr > 38) {
        if (!data->getValue(LOAD_P1E,&rval,l_id)) {
          data->addValue(LOAD_P1E, atof(split_line[38].c_str()), l_id);
        } else {
          data->setValue(LOAD_P1E, atof(split_line[38].c_str()), l_id);
        }
      }

      // LOAD_P1C
      if (nstr > 39) {
        if (!data->getValue(LOAD_P1C,&rval,l_id)) {
          data->addValue(LOAD_P1C, atof(split_line[39].c_str()), l_id);
        } else {
          data->setValue(LOAD_P1C, atof(split_line[39].c_str()), l_id);
        }
      }

      // LOAD_P2E
      if (nstr > 40) {
        if (!data->getValue(LOAD_P2E,&rval,l_id)) {
          data->addValue(LOAD_P2E, atof(split_line[40].c_str()), l_id);
        } else {
          data->setValue(LOAD_P2E, atof(split_line[40].c_str()), l_id);
        }
      }

      // LOAD_P2C
      if (nstr > 41) {
        if (!data->getValue(LOAD_P2C,&rval,l_id)) {
          data->addValue(LOAD_P2C, atof(split_line[41].c_str()), l_id);
        } else {
          data->setValue(LOAD_P2C, atof(split_line[41].c_str()), l_id);
        }
      }

      // LOAD_PFREQ
      if (nstr > 42) {
        if (!data->getValue(LOAD_PFREQ,&rval,l_id)) {
          data->addValue(LOAD_PFREQ, atof(split_line[42].c_str()), l_id);
        } else {
          data->setValue(LOAD_PFREQ, atof(split_line[42].c_str()), l_id);
        }
      }

      // LOAD_Q1E
      if (nstr > 43) {
        if (!data->getValue(LOAD_Q1E,&rval,l_id)) {
          data->addValue(LOAD_Q1E, atof(split_line[43].c_str()), l_id);
        } else {
          data->setValue(LOAD_Q1E, atof(split_line[43].c_str()), l_id);
        }
      }

      // LOAD_Q1C
      if (nstr > 44) {
        if (!data->getValue(LOAD_Q1C,&rval,l_id)) {
          data->addValue(LOAD_Q1C, atof(split_line[44].c_str()), l_id);
        } else {
          data->setValue(LOAD_Q1C, atof(split_line[44].c_str()), l_id);
        }
      }

      // LOAD_Q2E
      if (nstr > 45) {
        if (!data->getValue(LOAD_Q2E,&rval,l_id)) {
          data->addValue(LOAD_Q2E, atof(split_line[45].c_str()), l_id);
        } else {
          data->setValue(LOAD_Q2E, atof(split_line[45].c_str()), l_id);
        }
      }

      // LOAD_Q2C
      if (nstr > 46) {
        if (!data->getValue(LOAD_Q2C,&rval,l_id)) {
          data->addValue(LOAD_Q2C, atof(split_line[46].c_str()), l_id);
        } else {
          data->setValue(LOAD_Q2C, atof(split_line[46].c_str()), l_id);
        }
      }

      // LOAD_QFREQ
      if (nstr > 47) {
        if (!data->getValue(LOAD_QFREQ,&rval,l_id)) {
          data->addValue(LOAD_QFREQ, atof(split_line[47].c_str()), l_id);
        } else {
          data->setValue(LOAD_QFREQ, atof(split_line[47].c_str()), l_id);
        }
      }

      // LOAD_MTPA
      if (nstr > 48) {
        if (!data->getValue(LOAD_MTPA,&ival,l_id)) {
          data->addValue(LOAD_MTPA, atoi(split_line[48].c_str()), l_id);
        } else {
          data->setValue(LOAD_MTPA, atoi(split_line[48].c_str()), l_id);
        }
      }

      // LOAD_LFMA
      if (nstr > 49) {
        if (!data->getValue(LOAD_LFMA,&rval,l_id)) {
          data->addValue(LOAD_LFMA, atof(split_line[49].c_str()), l_id);
        } else {
          data->setValue(LOAD_LFMA, atof(split_line[49].c_str()), l_id);
        }
      }

      // LOAD_RSA
      if (nstr > 50) {
        if (!data->getValue(LOAD_RSA,&rval,l_id)) {
          data->addValue(LOAD_RSA, atof(split_line[50].c_str()), l_id);
        } else {
          data->setValue(LOAD_RSA, atof(split_line[50].c_str()), l_id);
        }
      }

      // LOAD_LSA
      if (nstr > 51) {
        if (!data->getValue(LOAD_LSA,&rval,l_id)) {
          data->addValue(LOAD_LSA, atof(split_line[51].c_str()), l_id);
        } else {
          data->setValue(LOAD_LSA, atof(split_line[51].c_str()), l_id);
        }
      }

      // LOAD_LPA
      if (nstr > 52) {
        if (!data->getValue(LOAD_LPA,&rval,l_id)) {
          data->addValue(LOAD_LPA, atof(split_line[52].c_str()), l_id);
        } else {
          data->setValue(LOAD_LPA, atof(split_line[52].c_str()), l_id);
        }
      }

      // LOAD_LPPA
      if (nstr > 53) {
        if (!data->getValue(LOAD_LPPA,&rval,l_id)) {
          data->addValue(LOAD_LPPA, atof(split_line[53].c_str()), l_id);
        } else {
          data->setValue(LOAD_LPPA, atof(split_line[53].c_str()), l_id);
        }
      }

      // LOAD_TPOA
      if (nstr > 54) {
        if (!data->getValue(LOAD_TPOA,&rval,l_id)) {
          data->addValue(LOAD_TPOA, atof(split_line[54].c_str()), l_id);
        } else {
          data->setValue(LOAD_TPOA, atof(split_line[54].c_str()), l_id);
        }
      }

      // LOAD_TPPOA
      if (nstr > 55) {
        if (!data->getValue(LOAD_TPPOA,&rval,l_id)) {
          data->addValue(LOAD_TPPOA, atof(split_line[55].c_str()), l_id);
        } else {
          data->setValue(LOAD_TPPOA, atof(split_line[55].c_str()), l_id);
        }
      }

      // LOAD_HA
      if (nstr > 56) {
        if (!data->getValue(LOAD_HA,&rval,l_id)) {
          data->addValue(LOAD_HA, atof(split_line[56].c_str()), l_id);
        } else {
          data->setValue(LOAD_HA, atof(split_line[56].c_str()), l_id);
        }
      }

      // LOAD_ETRQA
      if (nstr > 57) {
        if (!data->getValue(LOAD_ETRQA,&rval,l_id)) {
          data->addValue(LOAD_ETRQA, atof(split_line[57].c_str()), l_id);
        } else {
          data->setValue(LOAD_ETRQA, atof(split_line[57].c_str()), l_id);
        }
      }

      // LOAD_VTR1A
      if (nstr > 58) {
        if (!data->getValue(LOAD_VTR1A,&rval,l_id)) {
          data->addValue(LOAD_VTR1A, atof(split_line[58].c_str()), l_id);
        } else {
          data->setValue(LOAD_VTR1A, atof(split_line[58].c_str()), l_id);
        }
      }

      // LOAD_TTR1A
      if (nstr > 59) {
        if (!data->getValue(LOAD_TTR1A,&rval,l_id)) {
          data->addValue(LOAD_TTR1A, atof(split_line[59].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTR1A, atof(split_line[59].c_str()), l_id);
        }
      }

      // LOAD_FTR1A
      if (nstr > 60) {
        if (!data->getValue(LOAD_FTR1A,&rval,l_id)) {
          data->addValue(LOAD_FTR1A, atof(split_line[60].c_str()), l_id);
        } else {
          data->setValue(LOAD_FTR1A, atof(split_line[60].c_str()), l_id);
        }
      }

      // LOAD_VRC1A
      if (nstr > 61) {
        if (!data->getValue(LOAD_VRC1A,&rval,l_id)) {
          data->addValue(LOAD_VRC1A, atof(split_line[61].c_str()), l_id);
        } else {
          data->setValue(LOAD_VRC1A, atof(split_line[61].c_str()), l_id);
        }
      }

      // LOAD_TRC1A
      if (nstr > 62) {
        if (!data->getValue(LOAD_TRC1A,&rval,l_id)) {
          data->addValue(LOAD_TRC1A, atof(split_line[62].c_str()), l_id);
        } else {
          data->setValue(LOAD_TRC1A, atof(split_line[62].c_str()), l_id);
        }
      }

      // LOAD_VTR2A
      if (nstr > 63) {
        if (!data->getValue(LOAD_VTR2A,&rval,l_id)) {
          data->addValue(LOAD_VTR2A, atof(split_line[63].c_str()), l_id);
        } else {
          data->setValue(LOAD_VTR2A, atof(split_line[63].c_str()), l_id);
        }
      }

      // LOAD_TTR2A
      if (nstr > 64) {
        if (!data->getValue(LOAD_TTR2A,&rval,l_id)) {
          data->addValue(LOAD_TTR2A, atof(split_line[64].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTR2A, atof(split_line[64].c_str()), l_id);
        }
      }

      // LOAD_FTR2A
      if (nstr > 65) {
        if (!data->getValue(LOAD_FTR2A,&rval,l_id)) {
          data->addValue(LOAD_FTR2A, atof(split_line[65].c_str()), l_id);
        } else {
          data->setValue(LOAD_FTR2A, atof(split_line[65].c_str()), l_id);
        }
      }

      // LOAD_VRC2A
      if (nstr > 66) {
        if (!data->getValue(LOAD_VRC2A,&rval,l_id)) {
          data->addValue(LOAD_VRC2A, atof(split_line[66].c_str()), l_id);
        } else {
          data->setValue(LOAD_VRC2A, atof(split_line[66].c_str()), l_id);
        }
      }

      // LOAD_TRC2A
      if (nstr > 67) {
        if (!data->getValue(LOAD_TRC2A,&rval,l_id)) {
          data->addValue(LOAD_TRC2A, atof(split_line[67].c_str()), l_id);
        } else {
          data->setValue(LOAD_TRC2A, atof(split_line[67].c_str()), l_id);
        }
      }

      // LOAD_MTPB
      if (nstr > 68) {
        if (!data->getValue(LOAD_MTPB,&ival,l_id)) {
          data->addValue(LOAD_MTPB, atoi(split_line[68].c_str()), l_id);
        } else {
          data->setValue(LOAD_MTPB, atoi(split_line[68].c_str()), l_id);
        }
      }

      // LOAD_LFMB
      if (nstr > 69) {
        if (!data->getValue(LOAD_LFMB,&rval,l_id)) {
          data->addValue(LOAD_LFMB, atof(split_line[69].c_str()), l_id);
        } else {
          data->setValue(LOAD_LFMB, atof(split_line[69].c_str()), l_id);
        }
      }

      // LOAD_RSB
      if (nstr > 71) {
        if (!data->getValue(LOAD_RSB,&rval,l_id)) {
          data->addValue(LOAD_RSB, atof(split_line[70].c_str()), l_id);
        } else {
          data->setValue(LOAD_RSB, atof(split_line[70].c_str()), l_id);
        }
      }

      // LOAD_LSB
      if (nstr > 72) {
        if (!data->getValue(LOAD_LSB,&rval,l_id)) {
          data->addValue(LOAD_LSB, atof(split_line[71].c_str()), l_id);
        } else {
          data->setValue(LOAD_LSB, atof(split_line[71].c_str()), l_id);
        }
      }

      // LOAD_LPB
      if (nstr > 72) {
        if (!data->getValue(LOAD_LPB,&rval,l_id)) {
          data->addValue(LOAD_LPB, atof(split_line[72].c_str()), l_id);
        } else {
          data->setValue(LOAD_LPB, atof(split_line[72].c_str()), l_id);
        }
      }

      // LOAD_LPPB
      if (nstr > 73) {
        if (!data->getValue(LOAD_LPPB,&rval,l_id)) {
          data->addValue(LOAD_LPPB, atof(split_line[73].c_str()), l_id);
        } else {
          data->setValue(LOAD_LPPB, atof(split_line[73].c_str()), l_id);
        }
      }

      // LOAD_TPOB
      if (nstr > 74) {
        if (!data->getValue(LOAD_TPOB,&rval,l_id)) {
          data->addValue(LOAD_TPOB, atof(split_line[74].c_str()), l_id);
        } else {
          data->setValue(LOAD_TPOB, atof(split_line[74].c_str()), l_id);
        }
      }

      // LOAD_TPPOB
      if (nstr > 75) {
        if (!data->getValue(LOAD_TPPOB,&rval,l_id)) {
          data->addValue(LOAD_TPPOB, atof(split_line[75].c_str()), l_id);
        } else {
          data->setValue(LOAD_TPPOB, atof(split_line[75].c_str()), l_id);
        }
      }

      // LOAD_HB
      if (nstr > 76) {
        if (!data->getValue(LOAD_HB,&rval,l_id)) {
          data->addValue(LOAD_HB, atof(split_line[76].c_str()), l_id);
        } else {
          data->setValue(LOAD_HB, atof(split_line[76].c_str()), l_id);
        }
      }

      // LOAD_ETRQB
      if (nstr > 77) {
        if (!data->getValue(LOAD_ETRQB,&rval,l_id)) {
          data->addValue(LOAD_ETRQB, atof(split_line[77].c_str()), l_id);
        } else {
          data->setValue(LOAD_ETRQB, atof(split_line[77].c_str()), l_id);
        }
      }

      // LOAD_VTR1B
      if (nstr > 78) {
        if (!data->getValue(LOAD_VTR1B,&rval,l_id)) {
          data->addValue(LOAD_VTR1B, atof(split_line[78].c_str()), l_id);
        } else {
          data->setValue(LOAD_VTR1B, atof(split_line[78].c_str()), l_id);
        }
      }

      // LOAD_TTR1B
      if (nstr > 79) {
        if (!data->getValue(LOAD_TTR1B,&rval,l_id)) {
          data->addValue(LOAD_TTR1B, atof(split_line[79].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTR1B, atof(split_line[79].c_str()), l_id);
        }
      }

      // LOAD_FTR1B
      if (nstr > 80) {
        if (!data->getValue(LOAD_FTR1B,&rval,l_id)) {
          data->addValue(LOAD_FTR1B, atof(split_line[80].c_str()), l_id);
        } else {
          data->setValue(LOAD_FTR1B, atof(split_line[80].c_str()), l_id);
        }
      }

      // LOAD_VRC1B
      if (nstr > 81) {
        if (!data->getValue(LOAD_VRC1B,&rval,l_id)) {
          data->addValue(LOAD_VRC1B, atof(split_line[81].c_str()), l_id);
        } else {
          data->setValue(LOAD_VRC1B, atof(split_line[81].c_str()), l_id);
        }
      }

      // LOAD_TRC1B
      if (nstr > 82) {
        if (!data->getValue(LOAD_TRC1B,&rval,l_id)) {
          data->addValue(LOAD_TRC1B, atof(split_line[82].c_str()), l_id);
        } else {
          data->setValue(LOAD_TRC1B, atof(split_line[82].c_str()), l_id);
        }
      }

      // LOAD_VTR2B
      if (nstr > 83) {
        if (!data->getValue(LOAD_VTR2B,&rval,l_id)) {
          data->addValue(LOAD_VTR2B, atof(split_line[83].c_str()), l_id);
        } else {
          data->setValue(LOAD_VTR2B, atof(split_line[83].c_str()), l_id);
        }
      }

      // LOAD_TTR2B
      if (nstr > 84) {
        if (!data->getValue(LOAD_TTR2B,&rval,l_id)) {
          data->addValue(LOAD_TTR2B, atof(split_line[84].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTR2B, atof(split_line[84].c_str()), l_id);
        }
      }

      // LOAD_FTR2B
      if (nstr > 85) {
        if (!data->getValue(LOAD_FTR2B,&rval,l_id)) {
          data->addValue(LOAD_FTR2B, atof(split_line[85].c_str()), l_id);
        } else {
          data->setValue(LOAD_FTR2B, atof(split_line[85].c_str()), l_id);
        }
      }

      // LOAD_VRC2B
      if (nstr > 86) {
        if (!data->getValue(LOAD_VRC2B,&rval,l_id)) {
          data->addValue(LOAD_VRC2B, atof(split_line[86].c_str()), l_id);
        } else {
          data->setValue(LOAD_VRC2B, atof(split_line[86].c_str()), l_id);
        }
      }

      // LOAD_TRC2B
      if (nstr > 87) {
        if (!data->getValue(LOAD_TRC2B,&rval,l_id)) {
          data->addValue(LOAD_TRC2B, atof(split_line[87].c_str()), l_id);
        } else {
          data->setValue(LOAD_TRC2B, atof(split_line[87].c_str()), l_id);
        }
      }

      // LOAD_MTPC
      if (nstr > 88) {
        if (!data->getValue(LOAD_MTPC,&ival,l_id)) {
          data->addValue(LOAD_MTPC, atoi(split_line[88].c_str()), l_id);
        } else {
          data->setValue(LOAD_MTPC, atoi(split_line[88].c_str()), l_id);
        }
      }

      // LOAD_LFMC
      if (nstr > 89) {
        if (!data->getValue(LOAD_LFMC,&rval,l_id)) {
          data->addValue(LOAD_LFMC, atof(split_line[89].c_str()), l_id);
        } else {
          data->setValue(LOAD_LFMC, atof(split_line[89].c_str()), l_id);
        }
      }

      // LOAD_RSC
      if (nstr > 90) {
        if (!data->getValue(LOAD_RSC,&rval,l_id)) {
          data->addValue(LOAD_RSC, atof(split_line[90].c_str()), l_id);
        } else {
          data->setValue(LOAD_RSC, atof(split_line[90].c_str()), l_id);
        }
      }

      // LOAD_LSC
      if (nstr > 91) {
        if (!data->getValue(LOAD_LSC,&rval,l_id)) {
          data->addValue(LOAD_LSC, atof(split_line[91].c_str()), l_id);
        } else {
          data->setValue(LOAD_LSC, atof(split_line[91].c_str()), l_id);
        }
      }

      // LOAD_LPC
      if (nstr > 92) {
        if (!data->getValue(LOAD_LPC,&rval,l_id)) {
          data->addValue(LOAD_LPC, atof(split_line[92].c_str()), l_id);
        } else {
          data->setValue(LOAD_LPC, atof(split_line[92].c_str()), l_id);
        }
      }

      // LOAD_LPPC
      if (nstr > 93) {
        if (!data->getValue(LOAD_LPPC,&rval,l_id)) {
          data->addValue(LOAD_LPPC, atof(split_line[93].c_str()), l_id);
        } else {
          data->setValue(LOAD_LPPC, atof(split_line[93].c_str()), l_id);
        }
      }

      // LOAD_TPOC
      if (nstr > 94) {
        if (!data->getValue(LOAD_TPOC,&rval,l_id)) {
          data->addValue(LOAD_TPOC, atof(split_line[94].c_str()), l_id);
        } else {
          data->setValue(LOAD_TPOC, atof(split_line[94].c_str()), l_id);
        }
      }

      // LOAD_TPPOC
      if (nstr > 95) {
        if (!data->getValue(LOAD_TPPOC,&rval,l_id)) {
          data->addValue(LOAD_TPPOC, atof(split_line[95].c_str()), l_id);
        } else {
          data->setValue(LOAD_TPPOC, atof(split_line[95].c_str()), l_id);
        }
      }

      // LOAD_HC
      if (nstr > 96) {
        if (!data->getValue(LOAD_HC,&rval,l_id)) {
          data->addValue(LOAD_HC, atof(split_line[96].c_str()), l_id);
        } else {
          data->setValue(LOAD_HC, atof(split_line[96].c_str()), l_id);
        }
      }

      // LOAD_ETRQC
      if (nstr > 97) {
        if (!data->getValue(LOAD_ETRQC,&rval,l_id)) {
          data->addValue(LOAD_ETRQC, atof(split_line[97].c_str()), l_id);
        } else {
          data->setValue(LOAD_ETRQC, atof(split_line[97].c_str()), l_id);
        }
      }

      // LOAD_VTR1C
      if (nstr > 98) {
        if (!data->getValue(LOAD_VTR1C,&rval,l_id)) {
          data->addValue(LOAD_VTR1C, atof(split_line[98].c_str()), l_id);
        } else {
          data->setValue(LOAD_VTR1C, atof(split_line[98].c_str()), l_id);
        }
      }

      // LOAD_TTR1C
      if (nstr > 99) {
        if (!data->getValue(LOAD_TTR1C,&rval,l_id)) {
          data->addValue(LOAD_TTR1C, atof(split_line[99].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTR1C, atof(split_line[99].c_str()), l_id);
        }
      }

      // LOAD_FTR1C
      if (nstr > 100) {
        if (!data->getValue(LOAD_FTR1C,&rval,l_id)) {
          data->addValue(LOAD_FTR1C, atof(split_line[100].c_str()), l_id);
        } else {
          data->setValue(LOAD_FTR1C, atof(split_line[100].c_str()), l_id);
        }
      }

      // LOAD_VRC1C
      if (nstr > 101) {
        if (!data->getValue(LOAD_VRC1C,&rval,l_id)) {
          data->addValue(LOAD_VRC1C, atof(split_line[101].c_str()), l_id);
        } else {
          data->setValue(LOAD_VRC1C, atof(split_line[101].c_str()), l_id);
        }
      }

      // LOAD_TRC1C
      if (nstr > 102) {
        if (!data->getValue(LOAD_TRC1C,&rval,l_id)) {
          data->addValue(LOAD_TRC1C, atof(split_line[102].c_str()), l_id);
        } else {
          data->setValue(LOAD_TRC1C, atof(split_line[102].c_str()), l_id);
        }
      }

      // LOAD_VTR2C
      if (nstr > 103) {
        if (!data->getValue(LOAD_VTR2C,&rval,l_id)) {
          data->addValue(LOAD_VTR2C, atof(split_line[103].c_str()), l_id);
        } else {
          data->setValue(LOAD_VTR2C, atof(split_line[103].c_str()), l_id);
        }
      }

      // LOAD_TTR2C
      if (nstr > 104) {
        if (!data->getValue(LOAD_TTR2C,&rval,l_id)) {
          data->addValue(LOAD_TTR2C, atof(split_line[104].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTR2C, atof(split_line[104].c_str()), l_id);
        }
      }

      // LOAD_FTR2C
      if (nstr > 105) {
        if (!data->getValue(LOAD_FTR2C,&rval,l_id)) {
          data->addValue(LOAD_FTR2C, atof(split_line[105].c_str()), l_id);
        } else {
          data->setValue(LOAD_FTR2C, atof(split_line[105].c_str()), l_id);
        }
      }

      // LOAD_VRC2C
      if (nstr > 106) {
        if (!data->getValue(LOAD_VRC2C,&rval,l_id)) {
          data->addValue(LOAD_VRC2C, atof(split_line[106].c_str()), l_id);
        } else {
          data->setValue(LOAD_VRC2C, atof(split_line[106].c_str()), l_id);
        }
      }

      // LOAD_TRC2C
      if (nstr > 107) {
        if (!data->getValue(LOAD_TRC2C,&rval,l_id)) {
          data->addValue(LOAD_TRC2C, atof(split_line[107].c_str()), l_id);
        } else {
          data->setValue(LOAD_TRC2C, atof(split_line[107].c_str()), l_id);
        }
      }

      // LOAD_TSTALL
      if (nstr > 108) {
        if (!data->getValue(LOAD_TSTALL,&rval,l_id)) {
          data->addValue(LOAD_TSTALL, atof(split_line[108].c_str()), l_id);
        } else {
          data->setValue(LOAD_TSTALL, atof(split_line[108].c_str()), l_id);
        }
      }

      // LOAD_TRST
      if (nstr > 109) {
        if (!data->getValue(LOAD_TRST,&rval,l_id)) {
          data->addValue(LOAD_TRST, atof(split_line[109].c_str()), l_id);
        } else {
          data->setValue(LOAD_TRST, atof(split_line[109].c_str()), l_id);
        }
      }

      // LOAD_TV
      if (nstr > 110) {
        if (!data->getValue(LOAD_TV,&rval,l_id)) {
          data->addValue(LOAD_TV, atof(split_line[110].c_str()), l_id);
        } else {
          data->setValue(LOAD_TV, atof(split_line[110].c_str()), l_id);
        }
      }

      // LOAD_TF
      if (nstr > 111) {
        if (!data->getValue(LOAD_TF,&rval,l_id)) {
          data->addValue(LOAD_TF, atof(split_line[111].c_str()), l_id);
        } else {
          data->setValue(LOAD_TF, atof(split_line[111].c_str()), l_id);
        }
      }

      // LOAD_LFMD
      if (nstr > 112) {
        if (!data->getValue(LOAD_LFMD,&rval,l_id)) {
          data->addValue(LOAD_LFMD, atof(split_line[112].c_str()), l_id);
        } else {
          data->setValue(LOAD_LFMD, atof(split_line[112].c_str()), l_id);
        }
      }

      // LOAD_COMPPF
      if (nstr > 113) {
        if (!data->getValue(LOAD_COMPPF,&rval,l_id)) {
          data->addValue(LOAD_COMPPF, atof(split_line[113].c_str()), l_id);
        } else {
          data->setValue(LOAD_COMPPF, atof(split_line[113].c_str()), l_id);
        }
      }

      // LOAD_VSTALL
      if (nstr > 114) {
        if (!data->getValue(LOAD_VSTALL,&rval,l_id)) {
          data->addValue(LOAD_VSTALL, atof(split_line[114].c_str()), l_id);
        } else {
          data->setValue(LOAD_VSTALL, atof(split_line[114].c_str()), l_id);
        }
      }

      // LOAD_RSTALL
      if (nstr > 115) {
        if (!data->getValue(LOAD_RSTALL,&rval,l_id)) {
          data->addValue(LOAD_RSTALL, atof(split_line[115].c_str()), l_id);
        } else {
          data->setValue(LOAD_RSTALL, atof(split_line[115].c_str()), l_id);
        }
      }

      // LOAD_XSTALL
      if (nstr > 116) {
        if (!data->getValue(LOAD_XSTALL,&rval,l_id)) {
          data->addValue(LOAD_XSTALL, atof(split_line[116].c_str()), l_id);
        } else {
          data->setValue(LOAD_XSTALL, atof(split_line[116].c_str()), l_id);
        }
      }

      // LOAD_LFADJ
      if (nstr > 117) {
        if (!data->getValue(LOAD_LFADJ,&rval,l_id)) {
          data->addValue(LOAD_LFADJ, atof(split_line[117].c_str()), l_id);
        } else {
          data->setValue(LOAD_LFADJ, atof(split_line[117].c_str()), l_id);
        }
      }

      // LOAD_KP1
      if (nstr > 118) {
        if (!data->getValue(LOAD_KP1,&rval,l_id)) {
          data->addValue(LOAD_KP1, atof(split_line[118].c_str()), l_id);
        } else {
          data->setValue(LOAD_KP1, atof(split_line[118].c_str()), l_id);
        }
      }

      // LOAD_NP1
      if (nstr > 119) {
        if (!data->getValue(LOAD_NP1,&rval,l_id)) {
          data->addValue(LOAD_NP1, atof(split_line[119].c_str()), l_id);
        } else {
          data->setValue(LOAD_NP1, atof(split_line[119].c_str()), l_id);
        }
      }

      // LOAD_KQ1
      if (nstr > 120) {
        if (!data->getValue(LOAD_KQ1,&rval,l_id)) {
          data->addValue(LOAD_KQ1, atof(split_line[120].c_str()), l_id);
        } else {
          data->setValue(LOAD_KQ1, atof(split_line[120].c_str()), l_id);
        }
      }

      // LOAD_NQ1
      if (nstr > 121) {
        if (!data->getValue(LOAD_NQ1,&rval,l_id)) {
          data->addValue(LOAD_NQ1, atof(split_line[121].c_str()), l_id);
        } else {
          data->setValue(LOAD_NQ1, atof(split_line[121].c_str()), l_id);
        }
      }

      // LOAD_KP2
      if (nstr > 122) {
        if (!data->getValue(LOAD_KP2,&rval,l_id)) {
          data->addValue(LOAD_KP2, atof(split_line[122].c_str()), l_id);
        } else {
          data->setValue(LOAD_KP2, atof(split_line[122].c_str()), l_id);
        }
      }

      // LOAD_NP2
      if (nstr > 123) {
        if (!data->getValue(LOAD_NP2,&rval,l_id)) {
          data->addValue(LOAD_NP2, atof(split_line[123].c_str()), l_id);
        } else {
          data->setValue(LOAD_NP2, atof(split_line[123].c_str()), l_id);
        }
      }

      // LOAD_KQ2
      if (nstr > 124) {
        if (!data->getValue(LOAD_KQ2,&rval,l_id)) {
          data->addValue(LOAD_KQ2, atof(split_line[124].c_str()), l_id);
        } else {
          data->setValue(LOAD_KQ2, atof(split_line[124].c_str()), l_id);
        }
      }

      // LOAD_NQ2
      if (nstr > 125) {
        if (!data->getValue(LOAD_NQ2,&rval,l_id)) {
          data->addValue(LOAD_NQ2, atof(split_line[125].c_str()), l_id);
        } else {
          data->setValue(LOAD_NQ2, atof(split_line[125].c_str()), l_id);
        }
      }

      // LOAD_VBRK
      if (nstr > 126) {
        if (!data->getValue(LOAD_VBRK,&rval,l_id)) {
          data->addValue(LOAD_VBRK, atof(split_line[126].c_str()), l_id);
        } else {
          data->setValue(LOAD_VBRK, atof(split_line[126].c_str()), l_id);
        }
      }

      // LOAD_FRST
      if (nstr > 127) {
        if (!data->getValue(LOAD_FRST,&rval,l_id)) {
          data->addValue(LOAD_FRST, atof(split_line[127].c_str()), l_id);
        } else {
          data->setValue(LOAD_FRST, atof(split_line[127].c_str()), l_id);
        }
      }

      // LOAD_VRST
      if (nstr > 128) {
        if (!data->getValue(LOAD_VRST,&rval,l_id)) {
          data->addValue(LOAD_VRST, atof(split_line[128].c_str()), l_id);
        } else {
          data->setValue(LOAD_VRST, atof(split_line[128].c_str()), l_id);
        }
      }

      // LOAD_CMPKPF
      if (nstr > 129) {
        if (!data->getValue(LOAD_CMPKPF,&rval,l_id)) {
          data->addValue(LOAD_CMPKPF, atof(split_line[129].c_str()), l_id);
        } else {
          data->setValue(LOAD_CMPKPF, atof(split_line[129].c_str()), l_id);
        }
      }

      // LOAD_CMPKQF
      if (nstr > 130) {
        if (!data->getValue(LOAD_CMPKQF,&rval,l_id)) {
          data->addValue(LOAD_CMPKQF, atof(split_line[130].c_str()), l_id);
        } else {
          data->setValue(LOAD_CMPKQF, atof(split_line[130].c_str()), l_id);
        }
      }

      // LOAD_VC1OFF
      if (nstr > 131) {
        if (!data->getValue(LOAD_VC1OFF,&rval,l_id)) {
          data->addValue(LOAD_VC1OFF, atof(split_line[131].c_str()), l_id);
        } else {
          data->setValue(LOAD_VC1OFF, atof(split_line[131].c_str()), l_id);
        }
      }

      // LOAD_VC2OFF
      if (nstr > 132) {
        if (!data->getValue(LOAD_VC2OFF,&rval,l_id)) {
          data->addValue(LOAD_VC2OFF, atof(split_line[132].c_str()), l_id);
        } else {
          data->setValue(LOAD_VC2OFF, atof(split_line[132].c_str()), l_id);
        }
      }

      // LOAD_VC1ON
      if (nstr > 133) {
        if (!data->getValue(LOAD_VC1ON,&rval,l_id)) {
          data->addValue(LOAD_VC1ON, atof(split_line[133].c_str()), l_id);
        } else {
          data->setValue(LOAD_VC1ON, atof(split_line[133].c_str()), l_id);
        }
      }

      // LOAD_VC2ON
      if (nstr > 134) {
        if (!data->getValue(LOAD_VC2ON,&rval,l_id)) {
          data->addValue(LOAD_VC2ON, atof(split_line[134].c_str()), l_id);
        } else {
          data->setValue(LOAD_VC2ON, atof(split_line[134].c_str()), l_id);
        }
      }

      // LOAD_TTH
      if (nstr > 135) {
        if (!data->getValue(LOAD_TTH,&rval,l_id)) {
          data->addValue(LOAD_TTH, atof(split_line[135].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTH, atof(split_line[135].c_str()), l_id);
        }
      }

      // LOAD_TH1T
      if (nstr > 136) {
        if (!data->getValue(LOAD_TH1T,&rval,l_id)) {
          data->addValue(LOAD_TH1T, atof(split_line[136].c_str()), l_id);
        } else {
          data->setValue(LOAD_TH1T, atof(split_line[136].c_str()), l_id);
        }
      }

      // LOAD_TH2T
      if (nstr > 137) {
        if (!data->getValue(LOAD_TH2T,&rval,l_id)) {
          data->addValue(LOAD_TH2T, atof(split_line[137].c_str()), l_id);
        } else {
          data->setValue(LOAD_TH2T, atof(split_line[137].c_str()), l_id);
        }
      }

      // LOAD_FUVR
      if (nstr > 138) {
        if (!data->getValue(LOAD_FUVR,&rval,l_id)) {
          data->addValue(LOAD_FUVR, atof(split_line[138].c_str()), l_id);
        } else {
          data->setValue(LOAD_FUVR, atof(split_line[138].c_str()), l_id);
        }
      }

      // LOAD_UVTR1
      if (nstr > 139) {
        if (!data->getValue(LOAD_UVTR1,&rval,l_id)) {
          data->addValue(LOAD_UVTR1, atof(split_line[139].c_str()), l_id);
        } else {
          data->setValue(LOAD_UVTR1, atof(split_line[139].c_str()), l_id);
        }
      }

      // LOAD_TTR1
      if (nstr > 140) {
        if (!data->getValue(LOAD_TTR1,&rval,l_id)) {
          data->addValue(LOAD_TTR1, atof(split_line[140].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTR1, atof(split_line[140].c_str()), l_id);
        }
      }

      // LOAD_UVTR2
      if (nstr > 141) {
        if (!data->getValue(LOAD_UVTR2,&rval,l_id)) {
          data->addValue(LOAD_UVTR2, atof(split_line[141].c_str()), l_id);
        } else {
          data->setValue(LOAD_UVTR2, atof(split_line[141].c_str()), l_id);
        }
      }

      // LOAD_TTR2
      if (nstr > 142) {
        if (!data->getValue(LOAD_TTR2,&rval,l_id)) {
          data->addValue(LOAD_TTR2, atof(split_line[142].c_str()), l_id);
        } else {
          data->setValue(LOAD_TTR2, atof(split_line[142].c_str()), l_id);
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

      // LOAD_IT
      if (nstr > 5) {
        data.it = atoi(split_line[5].c_str());
      }

      // LOAD_MVA
      if (nstr > 11) {
        data.mva = atof(split_line[11].c_str());
      }

      // LOAD_BSS
      if (nstr > 12) {
        data.bss = atof(split_line[12].c_str());
      }

      // LOAD_RFDR
      if (nstr > 13) {
        data.rfdr = atof(split_line[13].c_str());
      }

      // LOAD_XFDR
      if (nstr > 14) {
        data.xfdr = atof(split_line[14].c_str());
      }

      // LOAD_FB
      if (nstr > 15) {
        data.fb = atof(split_line[15].c_str());
      }

      // LOAD_XXF
      if (nstr > 16) {
        data.xxf = atof(split_line[16].c_str());
      }

      // LOAD_TFIXHS
      if (nstr > 17) {
        data.tfixhs = atof(split_line[17].c_str());
      }

      // LOAD_TFIXLS
      if (nstr > 18) {
        data.tfixls = atof(split_line[18].c_str());
      }

      // LOAD_LTC
      if (nstr > 19) {
        data.ltc = atof(split_line[19].c_str());
      }

      // LOAD_TMIN
      if (nstr > 20) {
        data.tmin = atof(split_line[20].c_str());
      }

      // LOAD_TMAX
      if (nstr > 21) {
        data.tmax = atof(split_line[21].c_str());
      }

      // LOAD_STEP
      if (nstr > 22) {
        data.step = atof(split_line[22].c_str());
      }

      // LOAD_VMIN
      if (nstr > 23) {
        data.vmin = atof(split_line[23].c_str());
      }

      // LOAD_VMAX
      if (nstr > 24) {
        data.vmax = atof(split_line[24].c_str());
      }

      // LOAD_TDEL
      if (nstr > 25) {
        data.tdel = atof(split_line[25].c_str());
      }

      // LOAD_TTAP
      if (nstr > 26) {
        data.ttap = atof(split_line[26].c_str());
      }

      // LOAD_RCOMP
      if (nstr > 27) {
        data.rcomp = atof(split_line[27].c_str());
      }

      // LOAD_XCOMP
      if (nstr > 28) {
        data.xcomp = atof(split_line[28].c_str());
      }

      // LOAD_FMA
      if (nstr > 29) {
        data.fma = atof(split_line[29].c_str());
      }

      // LOAD_FMB
      if (nstr > 30) {
        data.fmb = atof(split_line[30].c_str());
      }

      // LOAD_FMC
      if (nstr > 31) {
        data.fmc = atof(split_line[31].c_str());
      }

      // LOAD_FMD
      if (nstr > 32) {
        data.fmd = atof(split_line[32].c_str());
      }

      // LOAD_FEL
      if (nstr > 33) {
        data.fel = atof(split_line[33].c_str());
      }

      // LOAD_PFEL
      if (nstr > 34) {
        data.pfel = atof(split_line[34].c_str());
      }

      // LOAD_VD1
      if (nstr > 35) {
        data.vd1 = atof(split_line[35].c_str());
      }

      // LOAD_VD2
      if (nstr > 36) {
        data.vd2 = atof(split_line[36].c_str());
      }

      // LOAD_PFS
      if (nstr > 37) {
        data.pfs = atof(split_line[37].c_str());
      }

      // LOAD_P1E
      if (nstr > 38) {
        data.p1e = atof(split_line[38].c_str());
      }

      // LOAD_P1C
      if (nstr > 39) {
        data.p1c = atof(split_line[39].c_str());
      }

      // LOAD_P2E
      if (nstr > 40) {
        data.p2e = atof(split_line[40].c_str());
      }

      // LOAD_P2C
      if (nstr > 41) {
        data.p2c = atof(split_line[41].c_str());
      }

      // LOAD_PFREQ
      if (nstr > 42) {
        data.pfreq = atof(split_line[42].c_str());
      }

      // LOAD_Q1E
      if (nstr > 43) {
        data.q1e = atof(split_line[43].c_str());
      }

      // LOAD_Q1C
      if (nstr > 44) {
        data.q1c = atof(split_line[44].c_str());
      }

      // LOAD_Q2E
      if (nstr > 45) {
        data.q2e = atof(split_line[45].c_str());
      }

      // LOAD_Q2C
      if (nstr > 46) {
        data.q2c = atof(split_line[46].c_str());
      }

      // LOAD_QFREQ
      if (nstr > 47) {
        data.qfreq = atof(split_line[47].c_str());
      }

      // LOAD_MTPA
      if (nstr > 48) {
        data.mtpa = atoi(split_line[48].c_str());
      }

      // LOAD_LFMA
      if (nstr > 49) {
        data.lfma = atof(split_line[49].c_str());
      }

      // LOAD_RSA
      if (nstr > 50) {
        data.rsa = atof(split_line[50].c_str());
      }

      // LOAD_LSA
      if (nstr > 51) {
        data.lsa = atof(split_line[51].c_str());
      }

      // LOAD_LPA
      if (nstr > 52) {
        data.lpa = atof(split_line[52].c_str());
      }

      // LOAD_LPPA
      if (nstr > 53) {
        data.lppa = atof(split_line[53].c_str());
      }

      // LOAD_TPOA
      if (nstr > 54) {
        data.tpoa = atof(split_line[54].c_str());
      }

      // LOAD_TPPOA
      if (nstr > 55) {
        data.tppoa = atof(split_line[55].c_str());
      }

      // LOAD_HA
      if (nstr > 56) {
        data.ha = atof(split_line[56].c_str());
      }

      // LOAD_ETRQA
      if (nstr > 57) {
        data.etrqa = atof(split_line[57].c_str());
      }

      // LOAD_VTR1A
      if (nstr > 58) {
        data.vtr1a = atof(split_line[58].c_str());
      }

      // LOAD_TTR1A
      if (nstr > 59) {
        data.ttr1a = atof(split_line[59].c_str());
      }

      // LOAD_FTR1A
      if (nstr > 60) {
        data.ftr1a = atof(split_line[60].c_str());
      }

      // LOAD_VRC1A
      if (nstr > 61) {
        data.vrc1a = atof(split_line[61].c_str());
      }

      // LOAD_TRC1A
      if (nstr > 62) {
        data.trc1a = atof(split_line[62].c_str());
      }

      // LOAD_VTR2A
      if (nstr > 63) {
        data.vtr2a = atof(split_line[63].c_str());
      }

      // LOAD_TTR2A
      if (nstr > 64) {
        data.ttr2a = atof(split_line[64].c_str());
      }

      // LOAD_FTR2A
      if (nstr > 65) {
        data.ftr2a = atof(split_line[65].c_str());
      }

      // LOAD_VRC2A
      if (nstr > 66) {
        data.vrc2a = atof(split_line[66].c_str());
      }

      // LOAD_TRC2A
      if (nstr > 67) {
        data.trc2a = atof(split_line[67].c_str());
      }

      // LOAD_MTPB
      if (nstr > 68) {
        data.mtpb = atoi(split_line[68].c_str());
      }

      // LOAD_LFMB
      if (nstr > 69) {
        data.lfmb = atof(split_line[69].c_str());
      }

      // LOAD_RSB
      if (nstr > 70) {
        data.rsb = atof(split_line[70].c_str());
      }

      // LOAD_LSB
      if (nstr > 71) {
        data.lsb = atof(split_line[71].c_str());
      }

      // LOAD_LPB
      if (nstr > 72) {
        data.lpb = atof(split_line[72].c_str());
      }

      // LOAD_LPPB
      if (nstr > 73) {
        data.lppb = atof(split_line[73].c_str());
      }

      // LOAD_TPOB
      if (nstr > 74) {
        data.tpob = atof(split_line[74].c_str());
      }

      // LOAD_TPPOB
      if (nstr > 75) {
        data.tppob = atof(split_line[75].c_str());
      }

      // LOAD_HB
      if (nstr > 76) {
        data.hb = atof(split_line[76].c_str());
      }

      // LOAD_ETRQB
      if (nstr > 77) {
        data.etrqb = atof(split_line[77].c_str());
      }

      // LOAD_VTR1B
      if (nstr > 78) {
        data.vtr1b = atof(split_line[78].c_str());
      }

      // LOAD_TTR1B
      if (nstr > 79) {
        data.ttr1b = atof(split_line[79].c_str());
      }

      // LOAD_FTR1B
      if (nstr > 80) {
        data.ftr1b = atof(split_line[80].c_str());
      }

      // LOAD_VRC1B
      if (nstr > 81) {
        data.vrc1b = atof(split_line[81].c_str());
      }

      // LOAD_TRC1B
      if (nstr > 82) {
        data.trc1b = atof(split_line[82].c_str());
      }

      // LOAD_VTR2B
      if (nstr > 83) {
        data.vtr2b = atof(split_line[83].c_str());
      }

      // LOAD_TTR2B
      if (nstr > 84) {
        data.ttr2b = atof(split_line[84].c_str());
      }

      // LOAD_FTR2B
      if (nstr > 85) {
        data.ftr2b = atof(split_line[85].c_str());
      }

      // LOAD_VRC2B
      if (nstr > 86) {
        data.vrc2b = atof(split_line[86].c_str());
      }

      // LOAD_TRC2B
      if (nstr > 87) {
        data.trc2b = atof(split_line[87].c_str());
      }

      // LOAD_MTPC
      if (nstr > 88) {
        data.mtpc = atoi(split_line[88].c_str());
      }

      // LOAD_LFMC
      if (nstr > 89) {
        data.lfmc = atof(split_line[89].c_str());
      }

      // LOAD_RSC
      if (nstr > 90) {
        data.rsc = atof(split_line[90].c_str());
      }

      // LOAD_LSC
      if (nstr > 91) {
        data.lsc = atof(split_line[91].c_str());
      }

      // LOAD_LPC
      if (nstr > 92) {
        data.lpc = atof(split_line[92].c_str());
      }

      // LOAD_LPPC
      if (nstr > 93) {
        data.lppc = atof(split_line[93].c_str());
      }

      // LOAD_TPOC
      if (nstr > 94) {
        data.tpoc = atof(split_line[94].c_str());
      }

      // LOAD_TPPOC
      if (nstr > 95) {
        data.tppoc = atof(split_line[95].c_str());
      }

      // LOAD_HC
      if (nstr > 96) {
        data.hc = atof(split_line[96].c_str());
      }

      // LOAD_ETRQC
      if (nstr > 97) {
        data.etrqc = atof(split_line[97].c_str());
      }

      // LOAD_VTR1C
      if (nstr > 98) {
        data.vtr1c = atof(split_line[98].c_str());
      }

      // LOAD_TTR1C
      if (nstr > 99) {
        data.ttr1c = atof(split_line[99].c_str());
      }

      // LOAD_FTR1C
      if (nstr > 100) {
        data.ftr1c = atof(split_line[100].c_str());
      }

      // LOAD_VRC1C
      if (nstr > 101) {
        data.vrc1c = atof(split_line[101].c_str());
      }

      // LOAD_TRC1C
      if (nstr > 102) {
        data.trc1c = atof(split_line[102].c_str());
      }

      // LOAD_VTR2C
      if (nstr > 103) {
        data.vtr2c = atof(split_line[103].c_str());
      }

      // LOAD_TTR2C
      if (nstr > 104) {
        data.ttr2c = atof(split_line[104].c_str());
      }

      // LOAD_FTR2C
      if (nstr > 105) {
        data.ftr2c = atof(split_line[105].c_str());
      }

      // LOAD_VRC2C
      if (nstr > 106) {
        data.vrc2c = atof(split_line[106].c_str());
      }

      // LOAD_TRC2C
      if (nstr > 107) {
        data.trc2c = atof(split_line[107].c_str());
      }

      // LOAD_TSTALL
      if (nstr > 108) {
        data.tstall = atof(split_line[108].c_str());
      }

      // LOAD_TRST
      if (nstr > 109) {
        data.trst = atof(split_line[109].c_str());
      }

      // LOAD_TV
      if (nstr > 110) {
        data.tv = atof(split_line[110].c_str());
      }

      // LOAD_TF
      if (nstr > 111) {
        data.tf = atof(split_line[111].c_str());
      }

      // LOAD_LFMD
      if (nstr > 112) {
        data.lfmd = atof(split_line[112].c_str());
      }

      // LOAD_COMPPF
      if (nstr > 113) {
        data.comppf = atof(split_line[113].c_str());
      }

      // LOAD_VSTALL
      if (nstr > 114) {
        data.vstall = atof(split_line[114].c_str());
      }

      // LOAD_RSTALL
      if (nstr > 115) {
        data.rstall = atof(split_line[115].c_str());
      }

      // LOAD_XSTALL
      if (nstr > 116) {
        data.xstall = atof(split_line[116].c_str());
      }

      // LOAD_LFADJ
      if (nstr > 117) {
        data.lfadj = atof(split_line[117].c_str());
      }

      // LOAD_KP1
      if (nstr > 118) {
        data.kp1 = atof(split_line[118].c_str());
      }

      // LOAD_NP1
      if (nstr > 119) {
        data.np1 = atof(split_line[119].c_str());
      }

      // LOAD_KQ1
      if (nstr > 120) {
        data.kq1 = atof(split_line[120].c_str());
      }

      // LOAD_NQ1
      if (nstr > 121) {
        data.nq1 = atof(split_line[121].c_str());
      }

      // LOAD_KP2
      if (nstr > 122) {
        data.kp2 = atof(split_line[122].c_str());
      }

      // LOAD_NP2
      if (nstr > 123) {
        data.np2 = atof(split_line[123].c_str());
      }

      // LOAD_KQ2
      if (nstr > 124) {
        data.kq2 = atof(split_line[124].c_str());
      }

      // LOAD_NQ2
      if (nstr > 125) {
        data.nq2 = atof(split_line[125].c_str());
      }

      // LOAD_VBRK
      if (nstr > 126) {
        data.vbrk = atof(split_line[126].c_str());
      }

      // LOAD_FRST
      if (nstr > 127) {
        data.frst = atof(split_line[127].c_str());
      }

      // LOAD_VRST
      if (nstr > 128) {
        data.vrst = atof(split_line[128].c_str());
      }

      // LOAD_CMPKPF
      if (nstr > 129) {
        data.cmpkpf = atof(split_line[129].c_str());
      }

      // LOAD_CMPKQF
      if (nstr > 130) {
        data.cmpkqf = atof(split_line[130].c_str());
      }

      // LOAD_VC1OFF
      if (nstr > 131) {
        data.vc1off = atof(split_line[131].c_str());
      }

      // LOAD_VC2OFF
      if (nstr > 132) {
        data.vc2off = atof(split_line[132].c_str());
      }

      // LOAD_VC1ON
      if (nstr > 133) {
        data.vc1on = atof(split_line[133].c_str());
      }

      // LOAD_VC2ON
      if (nstr > 134) {
        data.vc2on = atof(split_line[134].c_str());
      }

      // LOAD_TTH
      if (nstr > 135) {
        data.tth = atof(split_line[135].c_str());
      }

      // LOAD_TH1T
      if (nstr > 136) {
        data.th1t = atof(split_line[136].c_str());
      }

      // LOAD_TH2T
      if (nstr > 137) {
        data.th2t = atof(split_line[137].c_str());
      }

      // LOAD_FUVR
      if (nstr > 138) {
        data.fuvr = atof(split_line[138].c_str());
      }

      // LOAD_UVTR1
      if (nstr > 139) {
        data.uvtr1 = atof(split_line[139].c_str());
      }

      // LOAD_TTR1
      if (nstr > 140) {
        data.ttr1= atof(split_line[140].c_str());
      }

      // LOAD_UVTR2
      if (nstr > 141) {
        data.uvtr2 = atof(split_line[141].c_str());
      }

      // LOAD_TTR2
      if (nstr > 142) {
        data.ttr2 = atof(split_line[142].c_str());
      }
    }

    // Transfer data from composite load collection to new transformer
    // (branch) collection
    // @param comp_data data collection object from bus containing composite
    //     load object
    // @param t_data data collection object from new branch containing
    //     transformer
    // @param l_idx index of original composite load
    void setTransformer(gridpack::component::DataCollection *comp_data,
        gridpack::component::DataCollection *t_data, int l_idx) {
      int ival;
      double rval;
      ival = 1;
      t_data->addValue(BRANCH_NUM_ELEMENTS, ival);
      t_data->addValue(BRANCH_STATUS, ival, 0);
      t_data->addValue(BRANCH_SWITCHED, false, 0);
      t_data->addValue(BRANCH_CKT, " 1", 0);
      if (comp_data->getValue(LOAD_XXF, &rval, l_idx)) {
        t_data->addValue(BRANCH_X, rval, 0);
      }
      rval = 0.0;
      t_data->addValue(BRANCH_SHIFT, rval, 0);
      t_data->addValue(BRANCH_R, rval, 0);
      t_data->addValue(BRANCH_B, rval, 0);
      t_data->addValue(BRANCH_SHUNT_ADMTTNC_G1, rval, 0);
      t_data->addValue(BRANCH_SHUNT_ADMTTNC_B1, rval, 0);
      t_data->addValue(BRANCH_SHUNT_ADMTTNC_G2, rval, 0);
      t_data->addValue(BRANCH_SHUNT_ADMTTNC_B2, rval, 0);
      rval = 1.0;
      t_data->addValue(BRANCH_TAP, rval, 0);

      if (comp_data->getValue(LOAD_MVA, &rval, l_idx)) {
        t_data->addValue(LOAD_MVA, rval, 0);
      }
      if (comp_data->getValue(LOAD_XXF, &rval, l_idx)) {
        t_data->addValue(LOAD_XXF, rval, 0);
      }
      if (comp_data->getValue(LOAD_TFIXHS, &rval, l_idx)) {
        t_data->addValue(LOAD_TFIXHS, rval, 0);
      }
      if (comp_data->getValue(LOAD_TFIXLS, &rval, l_idx)) {
        t_data->addValue(LOAD_TFIXLS, rval, 0);
      }
      if (comp_data->getValue(LOAD_LTC, &rval, l_idx)) {
        t_data->addValue(LOAD_LTC, rval, 0);
      }
      if (comp_data->getValue(LOAD_TMIN, &rval, l_idx)) {
        t_data->addValue(LOAD_TMIN, rval, 0);
      }
      if (comp_data->getValue(LOAD_TMAX, &rval, l_idx)) {
        t_data->addValue(LOAD_TMAX, rval, 0);
      }
      if (comp_data->getValue(LOAD_STEP, &rval, l_idx)) {
        t_data->addValue(LOAD_STEP, rval, 0);
      }
      if (comp_data->getValue(LOAD_VMIN, &rval, l_idx)) {
        t_data->addValue(LOAD_VMIN, rval, 0);
      }
      if (comp_data->getValue(LOAD_VMAX, &rval, l_idx)) {
        t_data->addValue(LOAD_VMAX, rval, 0);
      }
      if (comp_data->getValue(LOAD_TDEL, &rval, l_idx)) {
        t_data->addValue(LOAD_TDEL, rval, 0);
      }
      if (comp_data->getValue(LOAD_TTAP, &rval, l_idx)) {
        t_data->addValue(LOAD_TTAP, rval, 0);
      }
      if (comp_data->getValue(LOAD_RCOMP, &rval, l_idx)) {
        t_data->addValue(LOAD_RCOMP, rval, 0);
      }
      if (comp_data->getValue(LOAD_XCOMP, &rval, l_idx)) {
        t_data->addValue(LOAD_XCOMP, rval, 0);
      }
      if (comp_data->getValue(LOAD_FMA, &rval, l_idx)) {
        t_data->addValue(LOAD_FMA, rval, 0);
      }
      if (comp_data->getValue(LOAD_FMB, &rval, l_idx)) {
        t_data->addValue(LOAD_FMB, rval, 0);
      }
      if (comp_data->getValue(LOAD_FMC, &rval, l_idx)) {
        t_data->addValue(LOAD_FMC, rval, 0);
      }
      if (comp_data->getValue(LOAD_FMD, &rval, l_idx)) {
        t_data->addValue(LOAD_FMD, rval, 0);
      }
    }

    // Transfer data from composite load collection to new feeder
    // (branch) collection
    // @param comp_data data collection object from bus containing composite
    //     load object
    // @param t_data data collection object from new branch representing
    //     feeder
    // @param l_idx index of original composite load
    void setFeeder( gridpack::component::DataCollection *comp_data,
        gridpack::component::DataCollection *t_data, int l_idx) {
      int ival;
      double rval;
      ival = 1;
      t_data->addValue(BRANCH_NUM_ELEMENTS, ival);
      t_data->addValue(BRANCH_STATUS,ival,0);
      t_data->addValue(BRANCH_SWITCHED, false, 0);
      t_data->addValue(BRANCH_CKT, " 1", 0);
      if (comp_data->getValue(LOAD_RFDR, &rval, l_idx)) {
        t_data->addValue(BRANCH_R, rval, 0);
      }
      if (comp_data->getValue(LOAD_XFDR, &rval, l_idx)) {
        t_data->addValue(BRANCH_X, rval, 0);
      }
      rval = 0.0;
      t_data->addValue(BRANCH_B,rval,0);
      t_data->addValue(BRANCH_SHUNT_ADMTTNC_G1,rval,0);
      t_data->addValue(BRANCH_SHUNT_ADMTTNC_B1,rval,0);
      t_data->addValue(BRANCH_SHUNT_ADMTTNC_G2,rval,0);
      t_data->addValue(BRANCH_SHUNT_ADMTTNC_B2,rval,0);
      t_data->addValue(BRANCH_SHIFT,rval,0);
      rval = 1.0;
      t_data->addValue(BRANCH_TAP,rval,0);

      if (comp_data->getValue(LOAD_MVA, &rval, l_idx)) {
        t_data->addValue(LOAD_MVA, rval, 0);
      }
      if (comp_data->getValue(LOAD_RFDR, &rval, l_idx)) {
        t_data->addValue(LOAD_RFDR, rval, 0);
      }
      if (comp_data->getValue(LOAD_XFDR, &rval, l_idx)) {
        t_data->addValue(LOAD_XFDR, rval, 0);
      }
      if (comp_data->getValue(LOAD_FB, &rval, l_idx)) {
        t_data->addValue(LOAD_FB, rval, 0);
      }
      if (comp_data->getValue(LOAD_FMA, &rval, l_idx)) {
        t_data->addValue(LOAD_FMA, rval, 0);
      }
      if (comp_data->getValue(LOAD_FMB, &rval, l_idx)) {
        t_data->addValue(LOAD_FMB, rval, 0);
      }
      if (comp_data->getValue(LOAD_FMC, &rval, l_idx)) {
        t_data->addValue(LOAD_FMC, rval, 0);
      }
      if (comp_data->getValue(LOAD_FMD, &rval, l_idx)) {
        t_data->addValue(LOAD_FMD, rval, 0);
      }
    }

    // Transfer data from composite load collection to load bus
    // collection
    // @param comp_data data collection object from bus containing composite
    //     load object
    // @param t_data data collection object from new bus representing
    //     load bus
    // @param l_idx index of original composite load
    void setLoadBus( gridpack::component::DataCollection *comp_data,
        gridpack::component::DataCollection *t_data, int l_idx) {
      int ival;
      double rval;
      // Parameters needed by the bus as a whole
      if (comp_data->getValue(CASE_SBASE, &rval)) {
        t_data->addValue(CASE_SBASE, rval);
      }
      rval = 0.0;
      t_data->addValue(BUS_SHUNT_GL, rval, 0);
      t_data->addValue(BUS_SHUNT_BL, rval, 0);
      t_data->addValue(SHUNT_BINIT, rval);
      ival = 1;
      t_data->addValue(BUS_TYPE, ival);

      if (comp_data->getValue(LOAD_P1C, &rval, l_idx)) {
        t_data->addValue(LOAD_P1C, rval);
      }
      if (comp_data->getValue(LOAD_P1E, &rval, l_idx)) {
        t_data->addValue(LOAD_P1E, rval);
      }
      if (comp_data->getValue(LOAD_P2C, &rval, l_idx)) {
        t_data->addValue(LOAD_P2C, rval);
      }
      if (comp_data->getValue(LOAD_P2E, &rval, l_idx)) {
        t_data->addValue(LOAD_P2E, rval);
      }
      if (comp_data->getValue(LOAD_PFREQ, &rval, l_idx)) {
        t_data->addValue(LOAD_PFREQ, rval);
      }
      if (comp_data->getValue(LOAD_PFS, &rval, l_idx)) {
        t_data->addValue(LOAD_PFS, rval);
      }
      if (comp_data->getValue(LOAD_PL, &rval, l_idx)) {
        t_data->addValue(LOAD_PL, rval);
      }
      if (comp_data->getValue(LOAD_Q1C, &rval, l_idx)) {
        t_data->addValue(LOAD_Q1C, rval);
      }
      if (comp_data->getValue(LOAD_Q1E, &rval, l_idx)) {
        t_data->addValue(LOAD_Q1E, rval);
      }
      if (comp_data->getValue(LOAD_Q2C, &rval, l_idx)) {
        t_data->addValue(LOAD_Q2C, rval);
      }
      if (comp_data->getValue(LOAD_Q2E, &rval, l_idx)) {
        t_data->addValue(LOAD_Q2E, rval);
      }
      if (comp_data->getValue(LOAD_QFREQ, &rval, l_idx)) {
        t_data->addValue(LOAD_QFREQ, rval);
      }
      if (comp_data->getValue(LOAD_QL, &rval, l_idx)) {
        t_data->addValue(LOAD_QL, rval);
      }

      if (comp_data->getValue(LOAD_MVA, &rval, l_idx)) {
        t_data->addValue(LOAD_MVA, rval);
      }
      if (comp_data->getValue(LOAD_BSS, &rval, l_idx)) {
        t_data->addValue(LOAD_BSS, rval);
      }
      if (comp_data->getValue(LOAD_XXF, &rval, l_idx)) {
        t_data->addValue(LOAD_XXF, rval);
      }
      if (comp_data->getValue(LOAD_RFDR, &rval, l_idx)) {
        t_data->addValue(LOAD_RFDR, rval);
      }
      if (comp_data->getValue(LOAD_XFDR, &rval, l_idx)) {
        t_data->addValue(LOAD_XFDR, rval);
      }
      if (comp_data->getValue(LOAD_TFIXHS, &rval, l_idx)) {
        t_data->addValue(LOAD_TFIXHS, rval);
      }
      if (comp_data->getValue(LOAD_TFIXLS, &rval, l_idx)) {
        t_data->addValue(LOAD_TFIXLS, rval);
      }
      if (comp_data->getValue(LOAD_VMIN, &rval, l_idx)) {
        t_data->addValue(LOAD_VMIN, rval);
      }
      if (comp_data->getValue(LOAD_VMAX, &rval, l_idx)) {
        t_data->addValue(LOAD_VMAX, rval);
      }
      if (comp_data->getValue(LOAD_TMIN, &rval, l_idx)) {
        t_data->addValue(LOAD_TMIN, rval);
      }
      if (comp_data->getValue(LOAD_TMAX, &rval, l_idx)) {
        t_data->addValue(LOAD_TMAX, rval);
      }
      if (comp_data->getValue(LOAD_STEP, &rval, l_idx)) {
        t_data->addValue(LOAD_STEP, rval);
      }
      if (comp_data->getValue(LOAD_FMA, &rval, l_idx)) {
        t_data->addValue(LOAD_FMA, rval);
      }
      if (comp_data->getValue(LOAD_FMB, &rval, l_idx)) {
        t_data->addValue(LOAD_FMB, rval);
      }
      if (comp_data->getValue(LOAD_FMC, &rval, l_idx)) {
        t_data->addValue(LOAD_FMC, rval);
      }
      if (comp_data->getValue(LOAD_FMD, &rval, l_idx)) {
        t_data->addValue(LOAD_FMD, rval);
      }

      // Parameters for different loads
      t_data->addValue(LOAD_NUMBER, 4);
      t_data->addValue(LOAD_ID, "M1", 0);
      if (comp_data->getValue(LOAD_MTPA, &ival, l_idx)) {
        t_data->addValue(LOAD_MTP, ival, 0);
      }
      if (comp_data->getValue(LOAD_LFMA, &rval, l_idx)) {
        t_data->addValue(LOAD_LFM, rval, 0);
      }
      if (comp_data->getValue(LOAD_RSA, &rval, l_idx)) {
        t_data->addValue(LOAD_RS, rval, 0);
      }
      if (comp_data->getValue(LOAD_LSA, &rval, l_idx)) {
        t_data->addValue(LOAD_LS, rval, 0);
      }
      if (comp_data->getValue(LOAD_LPA, &rval, l_idx)) {
        t_data->addValue(LOAD_LP, rval, 0);
      }
      if (comp_data->getValue(LOAD_LPPA, &rval, l_idx)) {
        t_data->addValue(LOAD_LPP, rval, 0);
      }
      if (comp_data->getValue(LOAD_TPOA, &rval, l_idx)) {
        t_data->addValue(LOAD_TPO, rval, 0);
      }
      if (comp_data->getValue(LOAD_TPPOA, &rval, l_idx)) {
        t_data->addValue(LOAD_TPPO, rval, 0);
      }
      if (comp_data->getValue(LOAD_HA, &rval, l_idx)) {
        t_data->addValue(LOAD_H, rval, 0);
      }
      if (comp_data->getValue(LOAD_ETRQA, &rval, l_idx)) {
        t_data->addValue(LOAD_ETRQ, rval, 0);
      }
      if (comp_data->getValue(LOAD_VTR1A, &rval, l_idx)) {
        t_data->addValue(LOAD_VTR1, rval, 0);
      }
      if (comp_data->getValue(LOAD_TTR1A, &rval, l_idx)) {
        t_data->addValue(LOAD_TTR1, rval, 0);
      }
      if (comp_data->getValue(LOAD_FTR1A, &rval, l_idx)) {
        t_data->addValue(LOAD_FTR1, rval, 0);
      }
      if (comp_data->getValue(LOAD_VRC1A, &rval, l_idx)) {
        t_data->addValue(LOAD_VRC1, rval, 0);
      }
      if (comp_data->getValue(LOAD_TRC1A, &rval, l_idx)) {
        t_data->addValue(LOAD_TRC1, rval, 0);
      }
      if (comp_data->getValue(LOAD_VTR2A, &rval, l_idx)) {
        t_data->addValue(LOAD_VTR2, rval, 0);
      }
      if (comp_data->getValue(LOAD_TTR2A, &rval, l_idx)) {
        t_data->addValue(LOAD_TTR2, rval, 0);
      }
      if (comp_data->getValue(LOAD_FTR2A, &rval, l_idx)) {
        t_data->addValue(LOAD_FTR2, rval, 0);
      }
      if (comp_data->getValue(LOAD_VRC2A, &rval, l_idx)) {
        t_data->addValue(LOAD_VRC2, rval, 0);
      }
      if (comp_data->getValue(LOAD_TRC2A, &rval, l_idx)) {
        t_data->addValue(LOAD_TRC2, rval, 0);
      }

      t_data->addValue(LOAD_ID, "M2", 1);
      if (comp_data->getValue(LOAD_MTPB, &ival, l_idx)) {
        t_data->addValue(LOAD_MTP, ival, 1);
      }
      if (comp_data->getValue(LOAD_LFMB, &rval, l_idx)) {
        t_data->addValue(LOAD_LFM, rval, 1);
      }
      if (comp_data->getValue(LOAD_RSB, &rval, l_idx)) {
        t_data->addValue(LOAD_RS, rval, 1);
      }
      if (comp_data->getValue(LOAD_LSB, &rval, l_idx)) {
        t_data->addValue(LOAD_LS, rval, 1);
      }
      if (comp_data->getValue(LOAD_LPB, &rval, l_idx)) {
        t_data->addValue(LOAD_LP, rval, 1);
      }
      if (comp_data->getValue(LOAD_LPPB, &rval, l_idx)) {
        t_data->addValue(LOAD_LPP, rval, 1);
      }
      if (comp_data->getValue(LOAD_TPOB, &rval, l_idx)) {
        t_data->addValue(LOAD_TPO, rval, 1);
      }
      if (comp_data->getValue(LOAD_TPPOB, &rval, l_idx)) {
        t_data->addValue(LOAD_TPPO, rval, 1);
      }
      if (comp_data->getValue(LOAD_HB, &rval, l_idx)) {
        t_data->addValue(LOAD_H, rval, 1);
      }
      if (comp_data->getValue(LOAD_ETRQB, &rval, l_idx)) {
        t_data->addValue(LOAD_ETRQ, rval, 1);
      }
      if (comp_data->getValue(LOAD_VTR1B, &rval, l_idx)) {
        t_data->addValue(LOAD_VTR1, rval, 1);
      }
      if (comp_data->getValue(LOAD_TTR1B, &rval, l_idx)) {
        t_data->addValue(LOAD_TTR1, rval, 1);
      }
      if (comp_data->getValue(LOAD_FTR1B, &rval, l_idx)) {
        t_data->addValue(LOAD_FTR1, rval, 1);
      }
      if (comp_data->getValue(LOAD_VRC1B, &rval, l_idx)) {
        t_data->addValue(LOAD_VRC1, rval, 1);
      }
      if (comp_data->getValue(LOAD_TRC1B, &rval, l_idx)) {
        t_data->addValue(LOAD_TRC1, rval, 1);
      }
      if (comp_data->getValue(LOAD_VTR2B, &rval, l_idx)) {
        t_data->addValue(LOAD_VTR2, rval, 1);
      }
      if (comp_data->getValue(LOAD_TTR2B, &rval, l_idx)) {
        t_data->addValue(LOAD_TTR2, rval, 1);
      }
      if (comp_data->getValue(LOAD_FTR2B, &rval, l_idx)) {
        t_data->addValue(LOAD_FTR2, rval, 1);
      }
      if (comp_data->getValue(LOAD_VRC2B, &rval, l_idx)) {
        t_data->addValue(LOAD_VRC2, rval, 1);
      }
      if (comp_data->getValue(LOAD_TRC2B, &rval, l_idx)) {
        t_data->addValue(LOAD_TRC2, rval, 1);
      }

      t_data->addValue(LOAD_ID, "M3", 2);
      if (comp_data->getValue(LOAD_MTPC, &ival, l_idx)) {
        t_data->addValue(LOAD_MTP, ival, 2);
      }
      if (comp_data->getValue(LOAD_LFMC, &rval, l_idx)) {
        t_data->addValue(LOAD_LFM, rval, 2);
      }
      if (comp_data->getValue(LOAD_RSC, &rval, l_idx)) {
        t_data->addValue(LOAD_RS, rval, 2);
      }
      if (comp_data->getValue(LOAD_LSC, &rval, l_idx)) {
        t_data->addValue(LOAD_LS, rval, 2);
      }
      if (comp_data->getValue(LOAD_LPC, &rval, l_idx)) {
        t_data->addValue(LOAD_LP, rval, 2);
      }
      if (comp_data->getValue(LOAD_LPPC, &rval, l_idx)) {
        t_data->addValue(LOAD_LPP, rval, 2);
      }
      if (comp_data->getValue(LOAD_TPOC, &rval, l_idx)) {
        t_data->addValue(LOAD_TPO, rval, 2);
      }
      if (comp_data->getValue(LOAD_TPPOC, &rval, l_idx)) {
        t_data->addValue(LOAD_TPPO, rval, 2);
      }
      if (comp_data->getValue(LOAD_HC, &rval, l_idx)) {
        t_data->addValue(LOAD_H, rval, 2);
      }
      if (comp_data->getValue(LOAD_ETRQC, &rval, l_idx)) {
        t_data->addValue(LOAD_ETRQ, rval, 2);
      }
      if (comp_data->getValue(LOAD_VTR1C, &rval, l_idx)) {
        t_data->addValue(LOAD_VTR1, rval, 2);
      }
      if (comp_data->getValue(LOAD_TTR1C, &rval, l_idx)) {
        t_data->addValue(LOAD_TTR1, rval, 2);
      }
      if (comp_data->getValue(LOAD_FTR1C, &rval, l_idx)) {
        t_data->addValue(LOAD_FTR1, rval, 2);
      }
      if (comp_data->getValue(LOAD_VRC1C, &rval, l_idx)) {
        t_data->addValue(LOAD_VRC1, rval, 2);
      }
      if (comp_data->getValue(LOAD_TRC1C, &rval, l_idx)) {
        t_data->addValue(LOAD_TRC1, rval, 2);
      }
      if (comp_data->getValue(LOAD_VTR2C, &rval, l_idx)) {
        t_data->addValue(LOAD_VTR2, rval, 2);
      }
      if (comp_data->getValue(LOAD_TTR2C, &rval, l_idx)) {
        t_data->addValue(LOAD_TTR2, rval, 2);
      }
      if (comp_data->getValue(LOAD_FTR2C, &rval, l_idx)) {
        t_data->addValue(LOAD_FTR2, rval, 2);
      }
      if (comp_data->getValue(LOAD_VRC2C, &rval, l_idx)) {
        t_data->addValue(LOAD_VRC2, rval, 2);
      }
      if (comp_data->getValue(LOAD_TRC2C, &rval, l_idx)) {
        t_data->addValue(LOAD_TRC2, rval, 2);
      }

      t_data->addValue(LOAD_ID, "M4", 3);
      if (comp_data->getValue(LOAD_MVA, &rval, l_idx)) {
        t_data->addValue(LOAD_MVA, rval, 3);
      }
      t_data->addValue(LOAD_MTP, 1, 3);
      if (comp_data->getValue(LOAD_TSTALL, &rval, l_idx)) {
        t_data->addValue(LOAD_TSTALL, rval, 3);
      }
      if (comp_data->getValue(LOAD_TRST, &rval, l_idx)) {
        t_data->addValue(LOAD_TRST, rval, 3);
      }
      if (comp_data->getValue(LOAD_TV, &rval, l_idx)) {
        t_data->addValue(LOAD_TV, rval, 3);
      }
      if (comp_data->getValue(LOAD_TF, &rval, l_idx)) {
        t_data->addValue(LOAD_TF, rval, 3);
      }
      if (comp_data->getValue(LOAD_LFMD, &rval, l_idx)) {
        t_data->addValue(LOAD_LFM, rval, 3);
      }
      if (comp_data->getValue(LOAD_COMPPF, &rval, l_idx)) {
        t_data->addValue(LOAD_COMPPF, rval, 3);
      }
      if (comp_data->getValue(LOAD_VSTALL, &rval, l_idx)) {
        t_data->addValue(LOAD_VSTALL, rval, 3);
      }
      if (comp_data->getValue(LOAD_RSTALL, &rval, l_idx)) {
        t_data->addValue(LOAD_RSTALL, rval, 3);
      }
      if (comp_data->getValue(LOAD_XSTALL, &rval, l_idx)) {
        t_data->addValue(LOAD_XSTALL, rval, 3);
      }
      if (comp_data->getValue(LOAD_LFADJ, &rval, l_idx)) {
        t_data->addValue(LOAD_LFADJ, rval, 3);
      }
      if (comp_data->getValue(LOAD_KP1, &rval, l_idx)) {
        t_data->addValue(LOAD_KP1, rval, 3);
      }
      if (comp_data->getValue(LOAD_NP1, &rval, l_idx)) {
        t_data->addValue(LOAD_NP1, rval, 3);
      }
      if (comp_data->getValue(LOAD_KQ1, &rval, l_idx)) {
        t_data->addValue(LOAD_KQ1, rval, 3);
      }
      if (comp_data->getValue(LOAD_NQ1, &rval, l_idx)) {
        t_data->addValue(LOAD_NQ1, rval, 3);
      }
      if (comp_data->getValue(LOAD_KP2, &rval, l_idx)) {
        t_data->addValue(LOAD_KP2, rval, 3);
      }
      if (comp_data->getValue(LOAD_NP2, &rval, l_idx)) {
        t_data->addValue(LOAD_NP2, rval, 3);
      }
      if (comp_data->getValue(LOAD_KQ2, &rval, l_idx)) {
        t_data->addValue(LOAD_KQ2, rval, 3);
      }
      if (comp_data->getValue(LOAD_NQ2, &rval, l_idx)) {
        t_data->addValue(LOAD_NQ2, rval, 3);
      }
      if (comp_data->getValue(LOAD_VBRK, &rval, l_idx)) {
        t_data->addValue(LOAD_VBRK, rval, 3);
      }
      if (comp_data->getValue(LOAD_FRST, &rval, l_idx)) {
        t_data->addValue(LOAD_FRST, rval, 3);
      }
      if (comp_data->getValue(LOAD_VRST, &rval, l_idx)) {
        t_data->addValue(LOAD_VRST, rval, 3);
      }
      if (comp_data->getValue(LOAD_CMPKPF, &rval, l_idx)) {
        t_data->addValue(LOAD_CMPKPF, rval, 3);
      }
      if (comp_data->getValue(LOAD_CMPKQF, &rval, l_idx)) {
        t_data->addValue(LOAD_CMPKQF, rval, 3);
      }
      if (comp_data->getValue(LOAD_VC1OFF, &rval, l_idx)) {
        t_data->addValue(LOAD_VC1OFF, rval, 3);
      }
      if (comp_data->getValue(LOAD_VC2OFF, &rval, l_idx)) {
        t_data->addValue(LOAD_VC2OFF, rval, 3);
      }
      if (comp_data->getValue(LOAD_VC1ON, &rval, l_idx)) {
        t_data->addValue(LOAD_VC1ON, rval, 3);
      }
      if (comp_data->getValue(LOAD_VC2ON, &rval, l_idx)) {
        t_data->addValue(LOAD_VC2ON, rval, 3);
      }
      if (comp_data->getValue(LOAD_TTH, &rval, l_idx)) {
        t_data->addValue(LOAD_TTH, rval, 3);
      }
      if (comp_data->getValue(LOAD_TH1T, &rval, l_idx)) {
        t_data->addValue(LOAD_TH1T, rval, 3);
      }
      if (comp_data->getValue(LOAD_TH2T, &rval, l_idx)) {
        t_data->addValue(LOAD_TH2T, rval, 3);
      }
      if (comp_data->getValue(LOAD_FUVR, &rval, l_idx)) {
        t_data->addValue(LOAD_FUVR, rval, 3);
      }
      if (comp_data->getValue(LOAD_UVTR1, &rval, l_idx)) {
        t_data->addValue(LOAD_UVTR1, rval, 3);
      }
      if (comp_data->getValue(LOAD_TTR1, &rval, l_idx)) {
        t_data->addValue(LOAD_TTR1, rval, 3);
      }
      if (comp_data->getValue(LOAD_UVTR2, &rval, l_idx)) {
        t_data->addValue(LOAD_UVTR2, rval, 3);
      }
      if (comp_data->getValue(LOAD_TTR2, &rval, l_idx)) {
        t_data->addValue(LOAD_TTR2, rval, 3);
      }
    }

    // Transfer data from composite load collection to low side bus
    // collection
    // @param comp_data data collection object from bus containing composite
    //     load object
    // @param t_data data collection object from new bus representing
    //     low side bus
    // @param l_idx index of original composite load
    void setLowSideBus( gridpack::component::DataCollection *comp_data,
        gridpack::component::DataCollection *t_data, int l_idx) {
      double rval;
      int ival;
      if (comp_data->getValue(CASE_SBASE, &rval)) {
        t_data->addValue(CASE_SBASE, rval);
      }
      rval = 0.0;
      t_data->addValue(BUS_SHUNT_GL, rval, 0);
      t_data->addValue(SHUNT_BINIT, rval);
      if (comp_data->getValue(LOAD_BSS, &rval, l_idx)) {
        t_data->addValue(BUS_SHUNT_BL, rval, 0);
      }
      if (comp_data->getValue(LOAD_MVA, &rval, l_idx)) {
        t_data->addValue(LOAD_MVA, rval, 0);
      }
      if (comp_data->getValue(LOAD_BSS, &rval, l_idx)) {
        t_data->addValue(LOAD_BSS, rval, 0);
      }
      ival = 1;
      t_data->addValue(BUS_TYPE, ival);
    }

    // Expand composite model by creating new buses and branches
    // @param comp_data data collection object from bus containing composite
    //     load object
    // @param new_buses vector of data collection objects representing new buses
    // @param new_branches vector of data collection objects representing new
    //     branches
    // @param l_idx index of original composite load
    void expandModel(gridpack::component::DataCollection *comp_data,
        std::vector<gridpack::component::DataCollection*> &new_buses,
        std::vector<gridpack::component::DataCollection*> &new_branches,
        int l_idx)
    {
      new_buses.clear();
      new_branches.clear();
      gridpack::component::DataCollection *bus_ptr;
      gridpack::component::DataCollection *branch_ptr;
      // For convenience, add comp_data as first component in new_buses
      new_buses.push_back(comp_data);
      // create new data collection to represent low side bus
      bus_ptr = new gridpack::component::DataCollection;
      setLowSideBus(comp_data, bus_ptr, l_idx);
      bus_ptr->addValue("NEW_BUS_TYPE","LOW_SIDE_BUS");
      new_buses.push_back(bus_ptr);
      // create new data collection to represent load bus
      bus_ptr = new gridpack::component::DataCollection;
      setLoadBus(comp_data, bus_ptr, l_idx);
      bus_ptr->addValue("NEW_BUS_TYPE","LOAD_BUS");
      new_buses.push_back(bus_ptr);
      // create new data collection to represent transformer branch
      branch_ptr = new gridpack::component::DataCollection;
      setTransformer(comp_data, branch_ptr, l_idx);
      branch_ptr->addValue("NEW_BRANCH_TYPE","TRANSFORMER");
      new_branches.push_back(branch_ptr);
      // create new data collection to represent feeder branch
      branch_ptr = new gridpack::component::DataCollection;
      setFeeder(comp_data, branch_ptr, l_idx);
      branch_ptr->addValue("NEW_BRANCH_TYPE","FEEDER");
      new_branches.push_back(branch_ptr);
      // Assign endpoint indices (based on local ordering of new buses new_buses
      // vector) of branches so that we can update neighbor information in the
      // parser. Start with the transformer branch
      branch_ptr = new_branches[0];
      branch_ptr->addValue(BRANCH_FROMBUS,0);
      branch_ptr->addValue(BRANCH_TOBUS,1);
      branch_ptr = new_branches[1];
      branch_ptr->addValue(BRANCH_FROMBUS,1);
      branch_ptr->addValue(BRANCH_TOBUS,2);
    }

};
}  // parser
}  // gridpack
#endif
