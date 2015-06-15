/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 *  Created on: December 30, 2014
 *      Author: Kevin Glass, Bruce Palmer
 */

#ifndef BASEPTIPARSER_HPP_
#define BASEPTIPARSER_HPP_

#define OLD_MAP

#include <iostream>
#include <string>
#include <vector>
#include <map>
#ifndef OLD_MAP
#include <boost/unordered_map.hpp>
#endif


#include "gridpack/component/base_component.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/utilities/exception.hpp"
#include "gridpack/utilities/string_utils.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/parser/base_parser.hpp"
#include "gridpack/parser/hash_distr.hpp"

namespace gridpack {
namespace parser {

template <class _network>
class BasePTIParser : public BaseParser<_network>
{
  public:

    /**
     * Constructor
     */
    explicit BasePTIParser()
    {
      p_timer = gridpack::utility::CoarseTimer::instance();
    }


    /**
     * Destructor
     */
    virtual ~BasePTIParser(){}

    /**
     * Parse a second file after original network has been distributed. This
     * requires the data in the second file to be distributed to all network
     * objects that need the data
     * @param fileName name of file
     */
    void externalParse(const std::string &fileName)
    {
      std::string ext = getExtension(fileName);
      if (ext == "dyr") {
        getDSExternal(fileName);
      } else if (ext == "uc") {
        getUCExternal(fileName);
      }
    }

  protected:

    /* ************************************************************************
     **************************************************************************
     ***** PROTECTED SCOPE
     **************************************************************************
     *********************************************************************** */

    /**
     * Assign network to internal network pointer variable
     */
    void setNetwork(boost::shared_ptr<_network> network)
    {
      p_network = network;
      BaseParser<_network>::setNetwork(network);
    }

    /**
     * This routine opens up a .dyr file with parameters for dynamic
     * simulation. It assumes that a .raw file has already been parsed
     */
    void getDS(const std::string & fileName)
    {
      int t_ds = p_timer->createCategory("Parser:getDS");
      p_timer->start(t_ds);
      int me(p_network->communicator().rank());

      if (me == 0) {
        std::ifstream input;
        input.open(fileName.c_str());
        if (!input.is_open()) {
          p_timer->stop(t_ds);
          return;
        }
        find_ds_par(input);
        input.close();
      }
      p_timer->stop(t_ds);
#if 0
      int i;
      printf("BUS data size: %d\n",p_busData.size());
      for (i=0; i<p_network->numBuses(); i++) {
        printf("Dumping bus: %d\n",i);
        p_network->getBusData(i)->dump();
      }
#endif
    }

    struct ds_params{
      // Generator parameters
      int bus_id; // ID of bus that owns generator
      char gen_id[3]; // Generator ID
      char gen_model[8];  // Generator model
      double inertia;  // Inertia constant 0
      double damping;  // Damping coefficient
      double reactance; // Transient reactance
      double tdop;
      double tdopp;
      double tqop;
      double tqopp;
      double xd;
      double xq;
      double xdp;
      double xqp;
      double xdpp;
      double xl;
      double s1;
      double s12;
      // Exciter parameters
      char ex_model[8];  // Exciter model
      bool has_exciter;
      int jbus;
      int m;
      double k;
      double t1;
      double t2;
      double t3;
      double uo;
      double uc;
      double pmax;
      double pmin;
      double t4;
      double k1;
      double k2;
      double t5;
      double k3;
      double k4;
      double t6;
      double k5;
      double k6;
      double t7;
      double k7;
      double k8;
      double db1;
      double err;
      double db2;
      double gv1;
      double pgv1;
      double gv2;
      double pgv2;
      double gv3;
      double pgv3;
      double gv4;
      double pgv4;
      double gv5;
      double pgv5;
      int iblock;
      // Governor parameters
      char gov_model[8];  // Exciter model
      bool has_governor;
      double tr;
      double ka;
      double ta;
      double tb;
      double tc;
      double vrmax;
      double vrmin;
      double ke;
      double te;
      double kf;
      double tf1;
      double rswitch;
      double e1;
      double se1;
      double e2;
      double se2;
    };

    /**
     * This routine opens up a .dyr file with parameters for dynamic
     * simulation and distributes the parameters to whatever processor holds the
     * corresponding buses. It assumes that a .raw file has already been parsed
     */
    void getDSExternal(const std::string & fileName)
    {

      //      int t_ds = p_timer->createCategory("Parser:getDS");
      //      p_timer->start(t_ds);
      int me(p_network->communicator().rank());

      std::vector<ds_params> ds_data;
      if (me == 0) {
        std::ifstream            input;
        input.open(fileName.c_str());
        if (!input.is_open()) {
          // p_timer->stop(t_ds);
          return;
        }
        find_ds_vector(input, &ds_data);
        input.close();
      }
      int nsize = ds_data.size();
      std::vector<int> buses;
      int i;
      for (i=0; i<nsize; i++) {
        buses.push_back(ds_data[i].bus_id);
      }
      gridpack::hash_distr::HashDistribution<_network,ds_params,ds_params>
        distr(p_network);
      distr.distributeBusValues(buses,ds_data);
      // Now match data with corresponding data collection objects
      gridpack::component::DataCollection *data;
      nsize = buses.size();
      for (i=0; i<nsize; i++) {
        int l_idx = buses[i];
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());

        // Find out how many generators are already on bus
        int ngen;
        if (!data->getValue(GENERATOR_NUMBER, &ngen)) continue;
        // Identify index of generator to which this data applies
        int g_id = -1;
        // Clean up 2 character tag for generator ID
        std::string tag = ds_data[i].gen_id;
        int j;
        for (j=0; j<ngen; j++) {
          std::string t_id;
          data->getValue(GENERATOR_ID,&t_id,j);
          if (tag == t_id) {
            g_id = j;
            break;
          }
        }
        if (g_id == -1) continue;

        std::string sval;
        double rval;
        int ival;
        bool bval=true;
        // GENERATOR_MODEL              "MODEL"        string
        if (!data->getValue(GENERATOR_MODEL,&sval,g_id)) {
          data->addValue(GENERATOR_MODEL, ds_data[i].gen_model, g_id);
        } else {
          data->setValue(GENERATOR_MODEL, ds_data[i].gen_model, g_id);
        }

        if (!strcmp(ds_data[i].gen_model,"GENCLS")) {

          // GENERATOR_INERTIA_CONSTANT_H                float
          if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H,&rval,g_id)) {
            data->addValue(GENERATOR_INERTIA_CONSTANT_H,
                ds_data[i].inertia, g_id);
          } else {
            data->setValue(GENERATOR_INERTIA_CONSTANT_H,
                ds_data[i].inertia, g_id);
          }

          // GENERATOR_DAMPING_COEFFICIENT_0             float
          if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&rval,g_id)) {
            data->addValue(GENERATOR_DAMPING_COEFFICIENT_0,
                ds_data[i].damping, g_id);
          } else {
            data->setValue(GENERATOR_DAMPING_COEFFICIENT_0,
                ds_data[i].damping, g_id);
          }

          // GENERATOR_TRANSIENT_REACTANCE               float
          if (!data->getValue(GENERATOR_TRANSIENT_REACTANCE,&rval,g_id)) {
            data->addValue(GENERATOR_TRANSIENT_REACTANCE,
                ds_data[i].reactance, g_id);
          } else {
            data->setValue(GENERATOR_TRANSIENT_REACTANCE,
                ds_data[i].reactance, g_id);
          }
        } else if (!strcmp(ds_data[i].gen_model,"GENSAL") ||
            !strcmp(ds_data[i].gen_model,"GENROU")) {
          // GENERATOR_TDOP
          if (!data->getValue(GENERATOR_TDOP,&rval,g_id)) {
            data->addValue(GENERATOR_TDOP,ds_data[i].tdop, g_id);
          } else {
            data->setValue(GENERATOR_TDOP, ds_data[i].tdop, g_id);
          }

          // GENERATOR_TDOPP
          if (!data->getValue(GENERATOR_TDOPP,&rval,g_id)) {
            data->addValue(GENERATOR_TDOPP, ds_data[i].tdopp, g_id);
          } else {
            data->setValue(GENERATOR_TDOPP, ds_data[i].tdopp, g_id);
          }

          // GENERATOR_TQOPP
          if (!data->getValue(GENERATOR_TQOPP,&rval,g_id)) {
            data->addValue(GENERATOR_TQOPP,
                ds_data[i].tqopp, g_id);
          } else {
            data->setValue(GENERATOR_TQOPP, ds_data[i].tqopp, g_id);
          }

          // GENERATOR_INERTIA_CONSTANT_H                           float
          if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H,&rval,g_id)) {
            data->addValue(GENERATOR_INERTIA_CONSTANT_H,
                ds_data[i].inertia, g_id);
          } else {
            data->setValue(GENERATOR_INERTIA_CONSTANT_H,
                ds_data[i].inertia, g_id);
          }

          // GENERATOR_DAMPING_COEFFICIENT_0                           float
          if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&rval,g_id)) {
            data->addValue(GENERATOR_DAMPING_COEFFICIENT_0,
                ds_data[i].damping, g_id);
          } else {
            data->setValue(GENERATOR_DAMPING_COEFFICIENT_0,
                ds_data[i].damping, g_id);
          }

          // GENERATOR_XD
          if (!data->getValue(GENERATOR_XD,&rval,g_id)) {
            data->addValue(GENERATOR_XD, ds_data[i].xd, g_id);
          } else {
            data->setValue(GENERATOR_XD, ds_data[i].xd, g_id);
          }

          // GENERATOR_XQ
          if (!data->getValue(GENERATOR_XQ,&rval,g_id)) {
            data->addValue(GENERATOR_XQ, ds_data[i].xq, g_id);
          } else {
            data->setValue(GENERATOR_XQ, ds_data[i].xq, g_id);
          }

          // GENERATOR_XDP
          if (!data->getValue(GENERATOR_XDP,&rval,g_id)) {
            data->addValue(GENERATOR_XDP, ds_data[i].xdp, g_id);
          } else {
            data->setValue(GENERATOR_XDP, ds_data[i].xdp, g_id);
          }

          // GENERATOR_XDPP
          if (!data->getValue(GENERATOR_XDPP,&rval,g_id)) {
            data->addValue(GENERATOR_XDPP, ds_data[i].xdpp, g_id);
          } else {
            data->setValue(GENERATOR_XDPP, ds_data[i].xdpp, g_id);
          }

          // GENERATOR_XL
          if (!data->getValue(GENERATOR_XL,&rval,g_id)) {
            data->addValue(GENERATOR_XL,
                ds_data[i].xl, g_id);
          } else {
            data->setValue(GENERATOR_XL,
                ds_data[i].xl, g_id);
          }

          // GENERATOR_S1
          if (!data->getValue(GENERATOR_S1,&rval,g_id)) {
            data->addValue(GENERATOR_XL, ds_data[i].s1, g_id);
          } else {
            data->setValue(GENERATOR_S1, ds_data[i].s1, g_id);
          }

          // GENERATOR_S12
          if (!data->getValue(GENERATOR_S12,&rval,g_id)) {
            data->addValue(GENERATOR_XL, ds_data[i].s12, g_id);
          } else {
            data->setValue(GENERATOR_S12, ds_data[i].s12, g_id);
          }

          if (!strcmp(ds_data[i].gen_model,"GENROU")) {
            // GENERATOR_TQOP
            if (!data->getValue(GENERATOR_TQOP,&rval,g_id)) {
              data->addValue(GENERATOR_TQOP,
                  ds_data[i].tqop, g_id);
            } else {
              data->setValue(GENERATOR_TQOP, ds_data[i].tqop, g_id);
            }

            // GENERATOR_XQP
            if (!data->getValue(GENERATOR_XQP,&rval,g_id)) {
              data->addValue(GENERATOR_XQP,
                  ds_data[i].xqp, g_id);
            } else {
              data->setValue(GENERATOR_XQP, ds_data[i].xqp, g_id);
            }
          }
        } else if (!strcmp(ds_data[i].gen_model,"WSIEG1")) {
          // HAS_GOVERNOR
          if (!data->getValue(HAS_GOVERNOR,&bval,g_id)) {
            data->addValue(HAS_GOVERNOR, true, g_id);
          } else {
            data->setValue(HAS_GOVERNOR, true, g_id);
          }

          // GOVERNOR_NAME
          if (!data->getValue(GOVERNOR_MODEL,&sval,g_id)) {
            data->addValue(GOVERNOR_MODEL, ds_data[i].gen_model, g_id);
          } else {
            data->setValue(GOVERNOR_MODEL, ds_data[i].gen_model, g_id);
          }

          // GOVERNOR_JBUS
          if (!data->getValue(GOVERNOR_JBUS,&ival,g_id)) {
            data->addValue(GOVERNOR_JBUS, ds_data[i].jbus, g_id);
          } else {
            data->setValue(GOVERNOR_JBUS, ds_data[i].jbus, g_id);
          }

          // GOVERNOR_M
          if (!data->getValue(GOVERNOR_M,&ival,g_id)) {
            data->addValue(GOVERNOR_M, ds_data[i].m, g_id);
          } else {
            data->setValue(GOVERNOR_M, ds_data[i].m, g_id);
          }

          // GOVERNOR_K
          if (!data->getValue(GOVERNOR_K,&rval,g_id)) {
            data->addValue(GOVERNOR_K, ds_data[i].k, g_id);
          } else {
            data->setValue(GOVERNOR_K, ds_data[i].k, g_id);
          }

          // GOVERNOR_T1
          if (!data->getValue(GOVERNOR_T1,&rval,g_id)) {
            data->addValue(GOVERNOR_T1, ds_data[i].t1, g_id);
          } else {
            data->setValue(GOVERNOR_T1, ds_data[i].t1, g_id);
          }

          // GOVERNOR_T2
          if (!data->getValue(GOVERNOR_T2,&rval,g_id)) {
            data->addValue(GOVERNOR_T2, ds_data[i].t2, g_id);
          } else {
            data->setValue(GOVERNOR_T2, ds_data[i].t2, g_id);
          }

          // GOVERNOR_T3
          if (!data->getValue(GOVERNOR_T3,&rval,g_id)) {
            data->addValue(GOVERNOR_T3, ds_data[i].t2, g_id);
          } else {
            data->setValue(GOVERNOR_T3, ds_data[i].t3, g_id);
          }

          // GOVERNOR_UO
          if (!data->getValue(GOVERNOR_UO,&rval,g_id)) {
            data->addValue(GOVERNOR_UO, ds_data[i].uo, g_id);
          } else {
            data->setValue(GOVERNOR_UO, ds_data[i].uo, g_id);
          }

          // GOVERNOR_UC
          if (!data->getValue(GOVERNOR_UC,&rval,g_id)) {
            data->addValue(GOVERNOR_UC, ds_data[i].uc, g_id);
          } else {
            data->setValue(GOVERNOR_UC, ds_data[i].uc, g_id);
          }

          // GOVERNOR_PMAX
          if (!data->getValue(GOVERNOR_PMAX,&rval,g_id)) {
            data->addValue(GOVERNOR_PMAX, ds_data[i].pmax, g_id);
          } else {
            data->setValue(GOVERNOR_PMAX, ds_data[i].pmax, g_id);
          }

          // GOVERNOR_PMIN
          if (!data->getValue(GOVERNOR_PMIN,&rval,g_id)) {
            data->addValue(GOVERNOR_PMIN, ds_data[i].pmin, g_id);
          } else {
            data->setValue(GOVERNOR_PMIN, ds_data[i].pmin, g_id);
          }

          // GOVERNOR_T4
          if (!data->getValue(GOVERNOR_T4,&rval,g_id)) {
            data->addValue(GOVERNOR_T4, ds_data[i].t4, g_id);
          } else {
            data->setValue(GOVERNOR_T4, ds_data[i].t4, g_id);
          }

          // GOVERNOR_K1
          if (!data->getValue(GOVERNOR_K1,&rval,g_id)) {
            data->addValue(GOVERNOR_K1, ds_data[i].k1, g_id);
          } else {
            data->setValue(GOVERNOR_K1, ds_data[i].k1, g_id);
          }

          // GOVERNOR_K2
          if (!data->getValue(GOVERNOR_K2,&rval,g_id)) {
            data->addValue(GOVERNOR_K2, ds_data[i].k2, g_id);
          } else {
            data->setValue(GOVERNOR_K2, ds_data[i].k2, g_id);
          }

          // GOVERNOR_T5
          if (!data->getValue(GOVERNOR_T5,&rval,g_id)) {
            data->addValue(GOVERNOR_T5, ds_data[i].t5, g_id);
          } else {
            data->setValue(GOVERNOR_T5, ds_data[i].t5, g_id);
          }

          // GOVERNOR_K3
          if (!data->getValue(GOVERNOR_K3,&rval,g_id)) {
            data->addValue(GOVERNOR_K3, ds_data[i].k3, g_id);
          } else {
            data->setValue(GOVERNOR_K3, ds_data[i].k3, g_id);
          }

          // GOVERNOR_K4
          if (!data->getValue(GOVERNOR_K4,&rval,g_id)) {
            data->addValue(GOVERNOR_K4, ds_data[i].k4, g_id);
          } else {
            data->setValue(GOVERNOR_K4, ds_data[i].k4, g_id);
          }

          // GOVERNOR_T6
          if (!data->getValue(GOVERNOR_T6,&rval,g_id)) {
            data->addValue(GOVERNOR_T6, ds_data[i].t6, g_id);
          } else {
            data->setValue(GOVERNOR_T6, ds_data[i].t6, g_id);
          }

          // GOVERNOR_K5
          if (!data->getValue(GOVERNOR_K5,&rval,g_id)) {
            data->addValue(GOVERNOR_K5, ds_data[i].k5, g_id);
          } else {
            data->setValue(GOVERNOR_K5, ds_data[i].k5, g_id);
          }

          // GOVERNOR_K6
          if (!data->getValue(GOVERNOR_K6,&rval,g_id)) {
            data->addValue(GOVERNOR_K6, ds_data[i].k6, g_id);
          } else {
            data->setValue(GOVERNOR_K6, ds_data[i].k6, g_id);
          }

          // GOVERNOR_T7
          if (!data->getValue(GOVERNOR_T7,&rval,g_id)) {
            data->addValue(GOVERNOR_T7, ds_data[i].t7, g_id);
          } else {
            data->setValue(GOVERNOR_T7, ds_data[i].t7, g_id);
          }

          // GOVERNOR_K7
          if (!data->getValue(GOVERNOR_K7,&rval,g_id)) {
            data->addValue(GOVERNOR_K7, ds_data[i].k7, g_id);
          } else {
            data->setValue(GOVERNOR_K7, ds_data[i].k7, g_id);
          }

          // GOVERNOR_K8
          if (!data->getValue(GOVERNOR_K8,&rval,g_id)) {
            data->addValue(GOVERNOR_K8, ds_data[i].k8, g_id);
          } else {
            data->setValue(GOVERNOR_K8, ds_data[i].k8, g_id);
          }

          // GOVERNOR_DB1
          if (!data->getValue(GOVERNOR_DB1,&rval,g_id)) {
            data->addValue(GOVERNOR_DB1, ds_data[i].db1, g_id);
          } else {
            data->setValue(GOVERNOR_DB1, ds_data[i].db1, g_id);
          }

          // GOVERNOR_ERR
          if (!data->getValue(GOVERNOR_ERR,&rval,g_id)) {
            data->addValue(GOVERNOR_ERR, ds_data[i].err, g_id);
          } else {
            data->setValue(GOVERNOR_ERR, ds_data[i].err, g_id);
          }

          // GOVERNOR_DB2
          if (!data->getValue(GOVERNOR_DB2,&rval,g_id)) {
            data->addValue(GOVERNOR_DB2, ds_data[i].db2, g_id);
          } else {
            data->setValue(GOVERNOR_DB2, ds_data[i].db2, g_id);
          }

          // GOVERNOR_GV1
          if (!data->getValue(GOVERNOR_GV1,&rval,g_id)) {
            data->addValue(GOVERNOR_GV1, ds_data[i].gv1, g_id);
          } else {
            data->setValue(GOVERNOR_GV1, ds_data[i].gv1, g_id);
          }

          // GOVERNOR_PGV1
          if (!data->getValue(GOVERNOR_PGV1,&rval,g_id)) {
            data->addValue(GOVERNOR_PGV1, ds_data[i].pgv1, g_id);
          } else {
            data->setValue(GOVERNOR_PGV1, ds_data[i].pgv1, g_id);
          }

          // GOVERNOR_GV2
          if (!data->getValue(GOVERNOR_GV2,&rval,g_id)) {
            data->addValue(GOVERNOR_GV2, ds_data[i].gv2, g_id);
          } else {
            data->setValue(GOVERNOR_GV2, ds_data[i].gv2, g_id);
          }

          // GOVERNOR_PGV2
          if (!data->getValue(GOVERNOR_PGV2,&rval,g_id)) {
            data->addValue(GOVERNOR_PGV2, ds_data[i].pgv2, g_id);
          } else {
            data->setValue(GOVERNOR_PGV2, ds_data[i].pgv2, g_id);
          }

          // GOVERNOR_GV3
          if (!data->getValue(GOVERNOR_GV3,&rval,g_id)) {
            data->addValue(GOVERNOR_GV3, ds_data[i].gv3, g_id);
          } else {
            data->setValue(GOVERNOR_GV3, ds_data[i].gv3, g_id);
          }

          // GOVERNOR_PGV3
          if (!data->getValue(GOVERNOR_PGV3,&rval,g_id)) {
            data->addValue(GOVERNOR_PGV3, ds_data[i].pgv3, g_id);
          } else {
            data->setValue(GOVERNOR_PGV3, ds_data[i].pgv3, g_id);
          }

          // GOVERNOR_GV4
          if (!data->getValue(GOVERNOR_GV4,&rval,g_id)) {
            data->addValue(GOVERNOR_GV4, ds_data[i].gv4, g_id);
          } else {
            data->setValue(GOVERNOR_GV4, ds_data[i].gv4, g_id);
          }

          // GOVERNOR_PGV4
          if (!data->getValue(GOVERNOR_PGV4,&rval,g_id)) {
            data->addValue(GOVERNOR_PGV4, ds_data[i].pgv4, g_id);
          } else {
            data->setValue(GOVERNOR_PGV4, ds_data[i].pgv4, g_id);
          }

          // GOVERNOR_GV5
          if (!data->getValue(GOVERNOR_GV5,&rval,g_id)) {
            data->addValue(GOVERNOR_GV5, ds_data[i].gv5, g_id);
          } else {
            data->setValue(GOVERNOR_GV5, ds_data[i].gv5, g_id);
          }

          // GOVERNOR_PGV5
          if (!data->getValue(GOVERNOR_PGV5,&rval,g_id)) {
            data->addValue(GOVERNOR_PGV5, ds_data[i].pgv5, g_id);
          } else {
            data->setValue(GOVERNOR_PGV5, ds_data[i].pgv5, g_id);
          }

          // GOVERNOR_IBLOCK
          if (!data->getValue(GOVERNOR_IBLOCK,&ival,g_id)) {
            data->addValue(GOVERNOR_IBLOCK, ds_data[i].iblock, g_id);
          } else {
            data->setValue(GOVERNOR_IBLOCK, ds_data[i].iblock, g_id);
          }
        } else if (!strcmp(ds_data[i].gen_model,"EXDC1") ||
            !strcmp(ds_data[i].gen_model,"EXDC2")) {
          // HAS_EXCITER
          if (!data->getValue(HAS_EXCITER,&bval,g_id)) {
            data->addValue(HAS_EXCITER, true, g_id);
          } else {
            data->setValue(HAS_EXCITER, true, g_id);
          }

          // EXCITER_MODEL
          if (!data->getValue(EXCITER_MODEL,&sval,g_id)) {
            data->addValue(EXCITER_MODEL, ds_data[i].gen_model, g_id);
          } else {
            data->setValue(EXCITER_MODEL, ds_data[i].gen_model, g_id);
          }

          // EXCITER_TR
          if (!data->getValue(EXCITER_TR,&rval,g_id)) {
            data->addValue(EXCITER_TR, ds_data[i].tr, g_id);
          } else {
            data->setValue(EXCITER_TR, ds_data[i].tr, g_id);
          }

          // EXCITER_KA
          if (!data->getValue(EXCITER_KA,&rval,g_id)) {
            data->addValue(EXCITER_KA, ds_data[i].ka, g_id);
          } else {
            data->setValue(EXCITER_KA, ds_data[i].ka, g_id);
          }

          // EXCITER_TA
          if (!data->getValue(EXCITER_TA,&rval,g_id)) {
            data->addValue(EXCITER_TA, ds_data[i].ta, g_id);
          } else {
            data->setValue(EXCITER_TA, ds_data[i].ta, g_id);
          }

          // EXCITER_TB
          if (!data->getValue(EXCITER_TB,&rval,g_id)) {
            data->addValue(EXCITER_TB, ds_data[i].tb, g_id);
          } else {
            data->setValue(EXCITER_TB, ds_data[i].tb, g_id);
          }

          // EXCITER_TC
          if (!data->getValue(EXCITER_TC,&rval,g_id)) {
            data->addValue(EXCITER_TC, ds_data[i].tc, g_id);
          } else {
            data->setValue(EXCITER_TC, ds_data[i].tc, g_id);
          }

          // EXCITER_VRMAX
          if (!data->getValue(EXCITER_VRMAX,&rval,g_id)) {
            data->addValue(EXCITER_VRMAX, ds_data[i].vrmax, g_id);
          } else {
            data->setValue(EXCITER_VRMAX, ds_data[i].vrmax, g_id);
          }

          // EXCITER_VRMIN
          if (!data->getValue(EXCITER_VRMIN,&rval,g_id)) {
            data->addValue(EXCITER_VRMIN, ds_data[i].vrmin, g_id);
          } else {
            data->setValue(EXCITER_VRMIN, ds_data[i].vrmin, g_id);
          }

          // EXCITER_KE
          if (!data->getValue(EXCITER_KE,&rval,g_id)) {
            data->addValue(EXCITER_KE, ds_data[i].ke, g_id);
          } else {
            data->setValue(EXCITER_KE, ds_data[i].ke, g_id);
          }

          // EXCITER_TE
          if (!data->getValue(EXCITER_TE,&rval,g_id)) {
            data->addValue(EXCITER_TE, ds_data[i].te, g_id);
          } else {
            data->setValue(EXCITER_TE, ds_data[i].te, g_id);
          }

          // EXCITER_KF
          if (!data->getValue(EXCITER_KF,&rval,g_id)) {
            data->addValue(EXCITER_KF, ds_data[i].kf, g_id);
          } else {
            data->setValue(EXCITER_KF, ds_data[i].kf, g_id);
          }

          // EXCITER_TF1
          if (!data->getValue(EXCITER_TF1,&rval,g_id)) {
            data->addValue(EXCITER_TF1, ds_data[i].tf1, g_id);
          } else {
            data->setValue(EXCITER_TF1, ds_data[i].tf1, g_id);
          }

          // EXCITER_SWITCH
          if (!data->getValue(EXCITER_SWITCH,&rval,g_id)) {
            data->addValue(EXCITER_SWITCH, ds_data[i].rswitch, g_id);
          } else {
            data->setValue(EXCITER_SWITCH, ds_data[i].rswitch, g_id);
          }

          // EXCITER_E1
          if (!data->getValue(EXCITER_E1,&rval,g_id)) {
            data->addValue(EXCITER_E1, ds_data[i].e1, g_id);
          } else {
            data->setValue(EXCITER_E1, ds_data[i].e1, g_id);
          }

          // EXCITER_SE1
          if (!data->getValue(EXCITER_SE1,&rval,g_id)) {
            data->addValue(EXCITER_SE1, ds_data[i].se1, g_id);
          } else {
            data->setValue(EXCITER_SE1, ds_data[i].se1, g_id);
          }

          // EXCITER_E2
          if (!data->getValue(EXCITER_E2,&rval,g_id)) {
            data->addValue(EXCITER_E2, ds_data[i].e2, g_id);
          } else {
            data->setValue(EXCITER_E2, ds_data[i].e2, g_id);
          }

          // EXCITER_SE2
          if (!data->getValue(EXCITER_SE2,&rval,g_id)) {
            data->addValue(EXCITER_SE2, ds_data[i].se2, g_id);
          } else {
            data->setValue(EXCITER_SE2, ds_data[i].se2, g_id);
          }
        }
      }
      //      p_timer->stop(t_ds);
    }

    struct uc_params{
      int bus_id; // ID of bus that owns generator
      char gen_id[3]; // Generator ID
      int type;
      double init_level; // Initial production level
      double min_gen; // Minimum generation
      double max_gen; // Maximum generation
      double max_oper; // Maximum operating generation
      int min_up;
      int min_down;
      double ramp_up;
      double ramp_down;
      double start_up; // Start up cost
      double const_cost; // Constant cost
      double lin_cost; // Linear cost
      double co_2_cost;
      double init_prd; // Init periods
      double start_cap; // Startup cap
      double shut_cap; // Shutdown cap
    };

    /**
     * This routine opens up a .uc file with parameters for a unit commitment
     * calculationand distributes the parameters to whatever processor holds the
     * corresponding buses. It assumes that a .raw file has already been parsed
     */
    void getUCExternal(const std::string & fileName)
    {

      int me(p_network->communicator().rank());

      std::vector<uc_params> uc_data;
      if (me == 0) {
        std::ifstream            input;
        input.open(fileName.c_str());
        if (!input.is_open()) {
          return;
        }
        find_uc_vector(input, &uc_data);
        input.close();
      }
      int nsize = uc_data.size();
      std::vector<int> buses;
      int i;
      for (i=0; i<nsize; i++) {
        buses.push_back(uc_data[i].bus_id);
      }
      gridpack::hash_distr::HashDistribution<_network,uc_params,uc_params>
        distr(p_network);
      distr.distributeBusValues(buses,uc_data);
      // Now match data with corresponding data collection objects
      gridpack::component::DataCollection *data;
      nsize = buses.size();
      for (i=0; i<nsize; i++) {
        int l_idx = buses[i];
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());

        // Find out how many generators are already on bus
        int ngen;
        if (!data->getValue(GENERATOR_NUMBER, &ngen)) continue;
        // Identify index of generator to which this data applies
        int g_id = -1;
        // Clean up 2 character tag for generator ID
        std::string tag = uc_data[i].gen_id;
        int j;
        for (j=0; j<ngen; j++) {
          std::string t_id;
          data->getValue(GENERATOR_ID,&t_id,j);
          if (tag == t_id) {
            g_id = j;
            break;
          }
        }
        if (g_id == -1) continue;

        double rval;
        int ival;
        if (!data->getValue("GENERATOR_TYPE",&ival,g_id)) {
          data->addValue("GENERATOR_TYPE", uc_data[i].type, g_id);
        } else {
          data->setValue("GENERATOR_TYPE", uc_data[i].type, g_id);
        }

        if (!data->getValue("GENERATOR_INIT_LEVEL",&rval,g_id)) {
          data->addValue("GENERATOR_INIT_LEVEL", uc_data[i].init_level, g_id);
        } else {
          data->setValue("GENERATOR_INIT_LEVEL", uc_data[i].init_level, g_id);
        }

        if (!data->getValue("GENERATOR_MIN_GEN",&rval,g_id)) {
          data->addValue("GENERATOR_MIN_GEN", uc_data[i].min_gen, g_id);
        } else {
          data->setValue("GENERATOR_MIN_GEN", uc_data[i].min_gen, g_id);
        }

        if (!data->getValue("GENERATOR_MAX_GEN",&rval,g_id)) {
          data->addValue("GENERATOR_MAX_GEN", uc_data[i].max_gen, g_id);
        } else {
          data->setValue("GENERATOR_MAX_GEN", uc_data[i].max_gen, g_id);
        }

        if (!data->getValue("GENERATOR_MIN_UP",&rval,g_id)) {
          data->addValue("GENERATOR_MIN_UP", uc_data[i].min_up, g_id);
        } else {
          data->setValue("GENERATOR_MIN_UP", uc_data[i].min_up, g_id);
        }

        if (!data->getValue("GENERATOR_MIN_DOWN",&rval,g_id)) {
          data->addValue("GENERATOR_MIN_DOWN", uc_data[i].min_down, g_id);
        } else {
          data->setValue("GENERATOR_MIN_DOWN", uc_data[i].min_down, g_id);
        }

        if (!data->getValue("GENERATOR_RAMP_UP",&rval,g_id)) {
          data->addValue("GENERATOR_RAMP_UP", uc_data[i].ramp_up, g_id);
        } else {
          data->setValue("GENERATOR_RAMP_UP", uc_data[i].ramp_up, g_id);
        }

        if (!data->getValue("GENERATOR_RAMP_DOWN",&rval,g_id)) {
          data->addValue("GENERATOR_RAMP_DOWN", uc_data[i].ramp_down, g_id);
        } else {
          data->setValue("GENERATOR_RAMP_DOWN", uc_data[i].ramp_down, g_id);
        }

        if (!data->getValue("GENERATOR_START_UP",&rval,g_id)) {
          data->addValue("GENERATOR_START_UP", uc_data[i].start_up, g_id);
        } else {
          data->setValue("GENERATOR_START_UP", uc_data[i].start_up, g_id);
        }

        if (!data->getValue("GENERATOR_CONST_COST",&rval,g_id)) {
          data->addValue("GENERATOR_CONST_COST", uc_data[i].const_cost, g_id);
        } else {
          data->setValue("GENERATOR_CONST_COST", uc_data[i].const_cost, g_id);
        }

        if (!data->getValue("GENERATOR_CO_2_COST",&rval,g_id)) {
          data->addValue("GENERATOR_CO_2_COST", uc_data[i].co_2_cost, g_id);
        } else {
          data->setValue("GENERATOR_CO_2_COST", uc_data[i].co_2_cost, g_id);
        }

        if (!data->getValue("GENERATOR_INIT_PRD",&rval,g_id)) {
          data->addValue("GENERATOR_INIT_PRD", uc_data[i].init_prd, g_id);
        } else {
          data->setValue("GENERATOR_INIT_PRD", uc_data[i].init_prd, g_id);
        }

        if (!data->getValue("GENERATOR_START_CAP",&rval,g_id)) {
          data->addValue("GENERATOR_START_CAP", uc_data[i].start_cap, g_id);
        } else {
          data->setValue("GENERATOR_START_CAP", uc_data[i].start_cap, g_id);
        }

        if (!data->getValue("GENERATOR_SHUT_CAP",&rval,g_id)) {
          data->addValue("GENERATOR_SHUT_CAP", uc_data[i].shut_cap, g_id);
        } else {
          data->setValue("GENERATOR_SHUT_CAP", uc_data[i].shut_cap, g_id);
        }
      }
    }

    // Extract extension from file name and convert it to lower case
    std::string getExtension(const std::string file)
    {
      std::string ret;
      std::string line = file;
      int ntok1 = line.find('.',0);
      if (ntok1 == std::string::npos) return
        ret;
      ntok1++;
      int ntok2 = line.find(' ',ntok1);
      if (ntok2 == std::string::npos)
        ntok2 = line.size();
      // get extension
      ret = line.substr(ntok1,ntok2-ntok1);
      // convert all characters to lower case 
      int size = ret.size();
      int i;
      for (i=0; i<size; i++) {
        if (isalpha(ret[i])) {
          ret[i] = tolower(ret[i]);
        }
      }
      return ret;
    }

    void find_ds_par(std::ifstream & input)
    {
      std::string          line;
      gridpack::component::DataCollection *data;
      while(std::getline(input,line)) {
        std::string record = line;
        int idx = line.find('/');
        while (idx == std::string::npos) {
          std::getline(input,line);
          idx = line.find('/');
          record.append(line);
        }
        idx = record.find('/');
        if (idx != std::string::npos) record.erase(idx,record.length()-idx);
        std::vector<std::string>  split_line;
        boost::split(split_line, record, boost::algorithm::is_any_of(","),
            boost::token_compress_on);

        // GENERATOR_BUSNUMBER               "I"                   integer
        int l_idx, o_idx;
        o_idx = atoi(split_line[0].c_str());
#ifdef OLD_MAP
        std::map<int, int>::iterator it;
#else
        boost::unordered_map<int, int>::iterator it;
#endif
        int nstr = split_line.size();
        it = p_busMap->find(o_idx);
        if (it != p_busMap->end()) {
          l_idx = it->second;
        } else {
          continue;
        }
        data = dynamic_cast<gridpack::component::DataCollection*>
          (p_network->getBusData(l_idx).get());

        // Find out how many generators are already on bus
        int ngen;
        if (!data->getValue(GENERATOR_NUMBER, &ngen)) continue;
        // Identify index of generator to which this data applies
        int g_id = -1;
        // Clean up 2 character tag for generator ID
        gridpack::utility::StringUtils util;
        std::string tag = util.clean2Char(split_line[2]);
        int i;
        for (i=0; i<ngen; i++) {
          std::string t_id;
          data->getValue(GENERATOR_ID,&t_id,i);
          if (tag == t_id) {
            g_id = i;
            break;
          }
        }
        if (g_id == -1) continue;

        std::string sval;
        double rval;
        int ival;
        bool bval;

        // GENERATOR_MODEL              "MODEL"                  string
        sval = util.trimQuotes(split_line[1]);
        util.toUpper(sval);

        if (sval == "GENCLS") {
          // GENERATOR_MODEL              "MODEL"                  string
          if (!data->getValue(GENERATOR_MODEL,&sval,g_id)) {
            data->addValue(GENERATOR_MODEL, sval.c_str(), g_id);
          } else {
            data->setValue(GENERATOR_MODEL, sval.c_str(), g_id);
          }

          // GENERATOR_INERTIA_CONSTANT_H                           float
          if (nstr > 3) {
            if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H,&rval,g_id)) {
              data->addValue(GENERATOR_INERTIA_CONSTANT_H,
                  atof(split_line[3].c_str()), g_id);
            } else {
              data->setValue(GENERATOR_INERTIA_CONSTANT_H,
                  atof(split_line[3].c_str()), g_id);
            }
          } 

          // GENERATOR_DAMPING_COEFFICIENT_0                           float
          if (nstr > 4) {
            if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&rval,g_id)) {
              data->addValue(GENERATOR_DAMPING_COEFFICIENT_0,
                  atof(split_line[4].c_str()), g_id);
            } else {
              data->setValue(GENERATOR_DAMPING_COEFFICIENT_0,
                  atof(split_line[4].c_str()), g_id);
            }
          }
        } else if (sval == "GENSAL") {
          // GENERATOR_MODEL              "MODEL"                  string
          if (!data->getValue(GENERATOR_MODEL,&sval,g_id)) {
            data->addValue(GENERATOR_MODEL, sval.c_str(), g_id);
          } else {
            data->setValue(GENERATOR_MODEL, sval.c_str(), g_id);
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

          // GENERATOR_TQOPP
          if (nstr > 5) {
            if (!data->getValue(GENERATOR_TQOPP,&rval,g_id)) {
              data->addValue(GENERATOR_TQOPP,
                  atof(split_line[5].c_str()), g_id);
            } else {
              data->setValue(GENERATOR_TQOPP,
                  atof(split_line[5].c_str()), g_id);
            }
          } 

          // GENERATOR_INERTIA_CONSTANT_H                           float
          if (nstr > 6) {
            if (!data->getValue(GENERATOR_INERTIA_CONSTANT_H,&rval,g_id)) {
              data->addValue(GENERATOR_INERTIA_CONSTANT_H,
                  atof(split_line[6].c_str()), g_id);
            } else {
              data->setValue(GENERATOR_INERTIA_CONSTANT_H,
                  atof(split_line[6].c_str()), g_id);
            }
          } 

          // GENERATOR_DAMPING_COEFFICIENT_0                           float
          if (nstr > 7) {
            if (!data->getValue(GENERATOR_DAMPING_COEFFICIENT_0,&rval,g_id)) {
              data->addValue(GENERATOR_DAMPING_COEFFICIENT_0,
                  atof(split_line[7].c_str()), g_id);
            } else {
              data->setValue(GENERATOR_DAMPING_COEFFICIENT_0,
                  atof(split_line[7].c_str()), g_id);
            }
          }

          // GENERATOR_XD
          if (nstr > 8) {
            if (!data->getValue(GENERATOR_XD,&rval,g_id)) {
              data->addValue(GENERATOR_XD,
                  atof(split_line[8].c_str()), g_id);
            } else {
              data->setValue(GENERATOR_XD,
                  atof(split_line[8].c_str()), g_id);
            }
          } 

          // GENERATOR_XQ
          if (nstr > 9) {
            if (!data->getValue(GENERATOR_XQ,&rval,g_id)) {
              data->addValue(GENERATOR_XQ,
                  atof(split_line[9].c_str()), g_id);
            } else {
              data->setValue(GENERATOR_XQ,
                  atof(split_line[9].c_str()), g_id);
            }
          } 

          // GENERATOR_XDP
          if (nstr > 10) {
            if (!data->getValue(GENERATOR_XDP,&rval,g_id)) {
              data->addValue(GENERATOR_XDP,
                  atof(split_line[10].c_str()), g_id);
            } else {
              data->setValue(GENERATOR_XDP,
                  atof(split_line[10].c_str()), g_id);
            }
          } 

          // GENERATOR_XDPP
          if (nstr > 11) {
            if (!data->getValue(GENERATOR_XDPP,&rval,g_id)) {
              data->addValue(GENERATOR_XDPP,
                  atof(split_line[11].c_str()), g_id);
            } else {
              data->setValue(GENERATOR_XDPP,
                  atof(split_line[11].c_str()), g_id);
            }
          } 

          // GENERATOR_XL
          if (nstr > 12) {
            if (!data->getValue(GENERATOR_XL,&rval,g_id)) {
              data->addValue(GENERATOR_XL,
                  atof(split_line[12].c_str()), g_id);
            } else {
              data->setValue(GENERATOR_XL,
                  atof(split_line[12].c_str()), g_id);
            }
          } 

          // GENERATOR_S1
          if (nstr > 13) {
            if (!data->getValue(GENERATOR_S1,&rval,g_id)) {
              data->addValue(GENERATOR_XL,
                  atof(split_line[13].c_str()), g_id);
            } else {
              data->setValue(GENERATOR_S1,
                  atof(split_line[13].c_str()), g_id);
            }
          } 

          // GENERATOR_S12
          if (nstr > 14) {
            if (!data->getValue(GENERATOR_S12,&rval,g_id)) {
              data->addValue(GENERATOR_XL,
                  atof(split_line[14].c_str()), g_id);
            } else {
              data->setValue(GENERATOR_S12,
                  atof(split_line[14].c_str()), g_id);
            }
          } 
        } else if (sval == "GENROU") {
          // GENERATOR_MODEL              "MODEL"                  string
          if (!data->getValue(GENERATOR_MODEL,&sval,g_id)) {
            data->addValue(GENERATOR_MODEL, sval.c_str(), g_id);
          } else {
            data->setValue(GENERATOR_MODEL, sval.c_str(), g_id);
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
        } else if (sval == "WSIEG1") {
          // HAS_GOVERNOR
          if (!data->getValue(HAS_GOVERNOR,&bval,g_id)) {
            data->addValue(HAS_GOVERNOR, true, g_id);
          } else {
            data->setValue(HAS_GOVERNOR, true, g_id);
          }

          // GOVERNOR_MODEL
          std::string stmp;
          if (!data->getValue(GOVERNOR_MODEL,&stmp,g_id)) {
            data->addValue(GOVERNOR_MODEL, sval.c_str(), g_id);
          } else {
            data->setValue(GOVERNOR_MODEL, sval.c_str(), g_id);
          }

          // GOVERNOR_JBUS
          if (nstr > 3) {
            if (!data->getValue(GOVERNOR_JBUS,&ival,g_id)) {
              data->addValue(GOVERNOR_JBUS,
                  atoi(split_line[3].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_JBUS,
                  atoi(split_line[3].c_str()), g_id);
            }
          } 

          // GOVERNOR_M
          if (nstr > 4) {
            if (!data->getValue(GOVERNOR_M,&ival,g_id)) {
              data->addValue(GOVERNOR_M,
                  atoi(split_line[4].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_M,
                  atoi(split_line[4].c_str()), g_id);
            }
          } 

          // GOVERNOR_K
          if (nstr > 5) {
            if (!data->getValue(GOVERNOR_K,&rval,g_id)) {
              data->addValue(GOVERNOR_K,
                  atof(split_line[5].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_K,
                  atof(split_line[5].c_str()), g_id);
            }
          } 

          // GOVERNOR_T1
          if (nstr > 6) {
            if (!data->getValue(GOVERNOR_T1,&rval,g_id)) {
              data->addValue(GOVERNOR_T1,
                  atof(split_line[6].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_T1,
                  atof(split_line[6].c_str()), g_id);
            }
          } 

          // GOVERNOR_T2
          if (nstr > 7) {
            if (!data->getValue(GOVERNOR_T2,&rval,g_id)) {
              data->addValue(GOVERNOR_T2,
                  atof(split_line[7].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_T2,
                  atof(split_line[7].c_str()), g_id);
            }
          } 

          // GOVERNOR_T3
          if (nstr > 8) {
            if (!data->getValue(GOVERNOR_T3,&rval,g_id)) {
              data->addValue(GOVERNOR_T3,
                  atof(split_line[8].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_T3,
                  atof(split_line[8].c_str()), g_id);
            }
          } 

          // GOVERNOR_UO
          if (nstr > 9) {
            if (!data->getValue(GOVERNOR_UO,&rval,g_id)) {
              data->addValue(GOVERNOR_UO,
                  atof(split_line[9].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_UO,
                  atof(split_line[9].c_str()), g_id);
            }
          }

          // GOVERNOR_UC
          if (nstr > 10) {
            if (!data->getValue(GOVERNOR_UC,&rval,g_id)) {
              data->addValue(GOVERNOR_UC,
                  atof(split_line[10].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_UC,
                  atof(split_line[10].c_str()), g_id);
            }
          } 

          // GOVERNOR_PMAX
          if (nstr > 11) {
            if (!data->getValue(GOVERNOR_PMAX,&rval,g_id)) {
              data->addValue(GOVERNOR_PMAX,
                  atof(split_line[11].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_PMAX,
                  atof(split_line[11].c_str()), g_id);
            }
          } 

          // GOVERNOR_PMIN
          if (nstr > 12) {
            if (!data->getValue(GOVERNOR_PMIN,&rval,g_id)) {
              data->addValue(GOVERNOR_PMIN,
                  atof(split_line[12].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_PMIN,
                  atof(split_line[12].c_str()), g_id);
            }
          } 

          // GOVERNOR_T4
          if (nstr > 13) {
            if (!data->getValue(GOVERNOR_T4,&rval,g_id)) {
              data->addValue(GOVERNOR_T4,
                  atof(split_line[13].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_T4,
                  atof(split_line[13].c_str()), g_id);
            }
          } 

          // GOVERNOR_K1
          if (nstr > 14) {
            if (!data->getValue(GOVERNOR_K1,&rval,g_id)) {
              data->addValue(GOVERNOR_K1,
                  atof(split_line[14].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_K1,
                  atof(split_line[14].c_str()), g_id);
            }
          } 

          // GOVERNOR_K2
          if (nstr > 15) {
            if (!data->getValue(GOVERNOR_K2,&rval,g_id)) {
              data->addValue(GOVERNOR_K2,
                  atof(split_line[15].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_K2,
                  atof(split_line[15].c_str()), g_id);
            }
          } 

          // GOVERNOR_T5
          if (nstr > 16) {
            if (!data->getValue(GOVERNOR_T5,&rval,g_id)) {
              data->addValue(GOVERNOR_T5,
                  atof(split_line[16].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_T5,
                  atof(split_line[16].c_str()), g_id);
            }
          } 

          // GOVERNOR_K3
          if (nstr > 17) {
            if (!data->getValue(GOVERNOR_K3,&rval,g_id)) {
              data->addValue(GOVERNOR_K3,
                  atof(split_line[17].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_K3,
                  atof(split_line[17].c_str()), g_id);
            }
          } 

          // GOVERNOR_K4
          if (nstr > 18) {
            if (!data->getValue(GOVERNOR_K4,&rval,g_id)) {
              data->addValue(GOVERNOR_K4,
                  atof(split_line[18].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_K4,
                  atof(split_line[18].c_str()), g_id);
            }
          } 

          // GOVERNOR_T6
          if (nstr > 19) {
            if (!data->getValue(GOVERNOR_T6,&rval,g_id)) {
              data->addValue(GOVERNOR_T6,
                  atof(split_line[19].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_T6,
                  atof(split_line[19].c_str()), g_id);
            }
          } 

          // GOVERNOR_K5
          if (nstr > 20) {
            if (!data->getValue(GOVERNOR_K5,&rval,g_id)) {
              data->addValue(GOVERNOR_K5,
                  atof(split_line[20].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_K5,
                  atof(split_line[20].c_str()), g_id);
            }
          } 

          // GOVERNOR_K6
          if (nstr > 21) {
            if (!data->getValue(GOVERNOR_K6,&rval,g_id)) {
              data->addValue(GOVERNOR_K6,
                  atof(split_line[21].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_K6,
                  atof(split_line[21].c_str()), g_id);
            }
          } 

          // GOVERNOR_T7
          if (nstr > 22) {
            if (!data->getValue(GOVERNOR_T7,&rval,g_id)) {
              data->addValue(GOVERNOR_T7,
                  atof(split_line[22].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_T7,
                  atof(split_line[22].c_str()), g_id);
            }
          } 

          // GOVERNOR_K7
          if (nstr > 23) {
            if (!data->getValue(GOVERNOR_K7,&rval,g_id)) {
              data->addValue(GOVERNOR_K7,
                  atof(split_line[23].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_K7,
                  atof(split_line[23].c_str()), g_id);
            }
          } 

          // GOVERNOR_K8
          if (nstr > 24) {
            if (!data->getValue(GOVERNOR_K8,&rval,g_id)) {
              data->addValue(GOVERNOR_K8,
                  atof(split_line[24].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_K8,
                  atof(split_line[24].c_str()), g_id);
            }
          } 

          // GOVERNOR_DB1
          if (nstr > 25) {
            if (!data->getValue(GOVERNOR_DB1,&rval,g_id)) {
              data->addValue(GOVERNOR_DB1,
                  atof(split_line[25].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_DB1,
                  atof(split_line[25].c_str()), g_id);
            }
          } 

          // GOVERNOR_ERR
          if (nstr > 26) {
            if (!data->getValue(GOVERNOR_ERR,&rval,g_id)) {
              data->addValue(GOVERNOR_ERR,
                  atof(split_line[26].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_ERR,
                  atof(split_line[26].c_str()), g_id);
            }
          } 

          // GOVERNOR_DB2
          if (nstr > 27) {
            if (!data->getValue(GOVERNOR_DB2,&rval,g_id)) {
              data->addValue(GOVERNOR_DB2,
                  atof(split_line[27].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_DB2,
                  atof(split_line[27].c_str()), g_id);
            }
          } 

          // GOVERNOR_GV1
          if (nstr > 28) {
            if (!data->getValue(GOVERNOR_GV1,&rval,g_id)) {
              data->addValue(GOVERNOR_GV1,
                  atof(split_line[28].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_GV1,
                  atof(split_line[28].c_str()), g_id);
            }
          } 

          // GOVERNOR_PGV1
          if (nstr > 29) {
            if (!data->getValue(GOVERNOR_PGV1,&rval,g_id)) {
              data->addValue(GOVERNOR_PGV1,
                  atof(split_line[29].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_PGV1,
                  atof(split_line[29].c_str()), g_id);
            }
          } 

          // GOVERNOR_GV2
          if (nstr > 30) {
            if (!data->getValue(GOVERNOR_GV2,&rval,g_id)) {
              data->addValue(GOVERNOR_GV2,
                  atof(split_line[30].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_GV2,
                  atof(split_line[30].c_str()), g_id);
            }
          } 

          // GOVERNOR_PGV2
          if (nstr > 31) {
            if (!data->getValue(GOVERNOR_PGV2,&rval,g_id)) {
              data->addValue(GOVERNOR_PGV2,
                  atof(split_line[31].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_PGV2,
                  atof(split_line[31].c_str()), g_id);
            }
          } 

          // GOVERNOR_GV3
          if (nstr > 32) {
            if (!data->getValue(GOVERNOR_GV3,&rval,g_id)) {
              data->addValue(GOVERNOR_GV3,
                  atof(split_line[32].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_GV3,
                  atof(split_line[32].c_str()), g_id);
            }
          } 

          // GOVERNOR_PGV3
          if (nstr > 33) {
            if (!data->getValue(GOVERNOR_PGV3,&rval,g_id)) {
              data->addValue(GOVERNOR_PGV3,
                  atof(split_line[33].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_PGV3,
                  atof(split_line[33].c_str()), g_id);
            }
          } 

          // GOVERNOR_GV4
          if (nstr > 34) {
            if (!data->getValue(GOVERNOR_GV4,&rval,g_id)) {
              data->addValue(GOVERNOR_GV4,
                  atof(split_line[34].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_GV4,
                  atof(split_line[34].c_str()), g_id);
            }
          } 

          // GOVERNOR_PGV4
          if (nstr > 35) {
            if (!data->getValue(GOVERNOR_PGV4,&rval,g_id)) {
              data->addValue(GOVERNOR_PGV4,
                  atof(split_line[35].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_PGV4,
                  atof(split_line[35].c_str()), g_id);
            }
          } 

          // GOVERNOR_GV5
          if (nstr > 36) {
            if (!data->getValue(GOVERNOR_GV5,&rval,g_id)) {
              data->addValue(GOVERNOR_GV5,
                  atof(split_line[36].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_GV4,
                  atof(split_line[36].c_str()), g_id);
            }
          } 

          // GOVERNOR_PGV5
          if (nstr > 37) {
            if (!data->getValue(GOVERNOR_PGV5,&rval,g_id)) {
              data->addValue(GOVERNOR_PGV5,
                  atof(split_line[37].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_PGV5,
                  atof(split_line[37].c_str()), g_id);
            }
          } 

          // GOVERNOR_IBLOCK
          if (nstr > 38) {
            if (!data->getValue(GOVERNOR_IBLOCK,&ival,g_id)) {
              data->addValue(GOVERNOR_PGV5,
                  atoi(split_line[38].c_str()), g_id);
            } else {
              data->setValue(GOVERNOR_IBLOCK,
                  atoi(split_line[38].c_str()), g_id);
            }
          } 
        } else if (sval == "EXDC1" || sval == "EXDC2") {
          // HAS_EXCITER
          if (!data->getValue(HAS_EXCITER,&bval,g_id)) {
            data->addValue(HAS_EXCITER, true, g_id);
          } else {
            data->setValue(HAS_EXCITER, true, g_id);
          }

          // EXCITER_MODEL
          std::string stmp;
          if (!data->getValue(EXCITER_MODEL,&stmp,g_id)) {
            data->addValue(EXCITER_MODEL, sval.c_str(), g_id);
          } else {
            data->setValue(EXCITER_MODEL, sval.c_str(), g_id);
          }

          // EXCITER_TR
          if (nstr > 3) {
            if (!data->getValue(EXCITER_TR,&rval,g_id)) {
              data->addValue(EXCITER_TR,
                  atof(split_line[3].c_str()), g_id);
            } else {
              data->setValue(EXCITER_TR,
                  atof(split_line[3].c_str()), g_id);
            }
          } 

          // EXCITER_KA
          if (nstr > 4) {
            if (!data->getValue(EXCITER_KA,&rval,g_id)) {
              data->addValue(EXCITER_KA,
                  atof(split_line[4].c_str()), g_id);
            } else {
              data->setValue(EXCITER_KA,
                  atof(split_line[4].c_str()), g_id);
            }
          } 

          // EXCITER_TA
          if (nstr > 5) {
            if (!data->getValue(EXCITER_TA,&rval,g_id)) {
              data->addValue(EXCITER_TA,
                  atof(split_line[5].c_str()), g_id);
            } else {
              data->setValue(EXCITER_TA,
                  atof(split_line[5].c_str()), g_id);
            }
          }

          // EXCITER_TB
          if (nstr > 6) {
            if (!data->getValue(EXCITER_TB,&rval,g_id)) {
              data->addValue(EXCITER_TB,
                  atof(split_line[6].c_str()), g_id);
            } else {
              data->setValue(EXCITER_TB,
                  atof(split_line[6].c_str()), g_id);
            }
          }

          // EXCITER_TC
          if (nstr > 7) {
            if (!data->getValue(EXCITER_TC,&rval,g_id)) {
              data->addValue(EXCITER_TC,
                  atof(split_line[7].c_str()), g_id);
            } else {
              data->setValue(EXCITER_TC,
                  atof(split_line[7].c_str()), g_id);
            }
          }

          // EXCITER_VRMAX
          if (nstr > 8) {
            if (!data->getValue(EXCITER_VRMAX,&rval,g_id)) {
              data->addValue(EXCITER_VRMAX,
                  atof(split_line[8].c_str()), g_id);
            } else {
              data->setValue(EXCITER_VRMAX,
                  atof(split_line[8].c_str()), g_id);
            }
          } 

          // EXCITER_VRMIN
          if (nstr > 9) {
            if (!data->getValue(EXCITER_VRMIN,&rval,g_id)) {
              data->addValue(EXCITER_VRMIN,
                  atof(split_line[9].c_str()), g_id);
            } else {
              data->setValue(EXCITER_VRMIN,
                  atof(split_line[9].c_str()), g_id);
            }
          } 

          // EXCITER_KE
          if (nstr > 10) {
            if (!data->getValue(EXCITER_KE,&rval,g_id)) {
              data->addValue(EXCITER_KE,
                  atof(split_line[10].c_str()), g_id);
            } else {
              data->setValue(EXCITER_KE,
                  atof(split_line[10].c_str()), g_id);
            }
          } 

          // EXCITER_TE
          if (nstr > 11) {
            if (!data->getValue(EXCITER_TE,&rval,g_id)) {
              data->addValue(EXCITER_TE,
                  atof(split_line[11].c_str()), g_id);
            } else {
              data->setValue(EXCITER_TE,
                  atof(split_line[11].c_str()), g_id);
            }
          } 

          // EXCITER_KF
          if (nstr > 12) {
            if (!data->getValue(EXCITER_KF,&rval,g_id)) {
              data->addValue(EXCITER_KF,
                  atof(split_line[12].c_str()), g_id);
            } else {
              data->setValue(EXCITER_KF,
                  atof(split_line[12].c_str()), g_id);
            }
          } 

          // EXCITER_TF1
          if (nstr > 13) {
            if (!data->getValue(EXCITER_TF1,&rval,g_id)) {
              data->addValue(EXCITER_TF1,
                  atof(split_line[13].c_str()), g_id);
            } else {
              data->setValue(EXCITER_TF1,
                  atof(split_line[13].c_str()), g_id);
            }
          } 

          // EXCITER_SWITCH
          if (nstr > 14) {
            if (!data->getValue(EXCITER_SWITCH,&rval,g_id)) {
              data->addValue(EXCITER_SWITCH,
                  atof(split_line[14].c_str()), g_id);
            } else {
              data->setValue(EXCITER_SWITCH,
                  atof(split_line[14].c_str()), g_id);
            }
          } 

          // EXCITER_E1
          if (nstr > 15) {
            if (!data->getValue(EXCITER_E1,&rval,g_id)) {
              data->addValue(EXCITER_E1,
                  atof(split_line[15].c_str()), g_id);
            } else {
              data->setValue(EXCITER_E1,
                  atof(split_line[15].c_str()), g_id);
            }
          } 

          // EXCITER_SE1
          if (nstr > 16) {
            if (!data->getValue(EXCITER_SE1,&rval,g_id)) {
              data->addValue(EXCITER_SE1,
                  atof(split_line[16].c_str()), g_id);
            } else {
              data->setValue(EXCITER_SE1,
                  atof(split_line[16].c_str()), g_id);
            }
          } 

          // EXCITER_E2
          if (nstr > 17) {
            if (!data->getValue(EXCITER_E2,&rval,g_id)) {
              data->addValue(EXCITER_E2,
                  atof(split_line[17].c_str()), g_id);
            } else {
              data->setValue(EXCITER_E2,
                  atof(split_line[17].c_str()), g_id);
            }
          } 

          // EXCITER_SE2
          if (nstr > 18) {
            if (!data->getValue(EXCITER_SE2,&rval,g_id)) {
              data->addValue(EXCITER_SE2,
                  atof(split_line[18].c_str()), g_id);
            } else {
              data->setValue(EXCITER_SE2,
                  atof(split_line[18].c_str()), g_id);
            }
          } 
        }
      }
    }

    void find_ds_vector(std::ifstream & input, std::vector<ds_params> *ds_vector)
    {
      std::string          line;
      ds_vector->clear();
      while(std::getline(input,line)) {
        std::string record = line;
        int idx = line.find('/');
        while (idx == std::string::npos) {
          std::getline(input,line);
          idx = line.find('/');
          record.append(line);
        }
        idx = record.find('/');
        if (idx != std::string::npos) record.erase(idx,record.length()-idx);
        std::vector<std::string>  split_line;
        boost::split(split_line, record, boost::algorithm::is_any_of(","),
            boost::token_compress_on);

        ds_params data;

        // GENERATOR_BUSNUMBER               "I"                   integer
        int o_idx;
        o_idx = atoi(split_line[0].c_str());
        data.bus_id = o_idx;

        int nstr = split_line.size();

        // Clean up 2 character tag for generator ID
        gridpack::utility::StringUtils util;
        std::string tag = util.clean2Char(split_line[2]);
        strcpy(data.gen_id, tag.c_str());

        std::string sval;
        double rval;
        int ival;

        sval = util.trimQuotes(split_line[1]);
        util.toUpper(sval);

        // GENERATOR_MODEL              "MODEL"                  integer
        strcpy(data.gen_model, sval.c_str());

        if (sval == "GENCLS") {
          // GENERATOR_INERTIA_CONSTANT_H                           float
          if (nstr > 3) {
            data.inertia = atof(split_line[3].c_str());
          } 

          // GENERATOR_DAMPING_COEFFICIENT_0                           float
          if (nstr > 4) {
            data.damping = atof(split_line[4].c_str());
          }

          // GENERATOR_TRANSIENT_REACTANCE                             float
          if (nstr > 5) {
            data.reactance = atof(split_line[5].c_str());
          }
        } else if (sval == "GENSAL") {
          // GENERATOR_TDOP
          if (nstr > 3) {
            data.tdop = atof(split_line[3].c_str());
          } 

          // GENERATOR_TDOPP
          if (nstr > 4) {
            data.tdopp = atof(split_line[4].c_str());
          } 

          // GENERATOR_TQOPP
          if (nstr > 5) {
            data.tqopp = atof(split_line[5].c_str());
          } 

          // GENERATOR_INERTIA_CONSTANT_H                           float
          if (nstr > 6) {
            data.inertia = atof(split_line[6].c_str());
          } 

          // GENERATOR_DAMPING_COEFFICIENT_0                           float
          if (nstr > 7) {
            data.damping = atof(split_line[7].c_str());
          }

          // GENERATOR_XD
          if (nstr > 8) {
            data.xd = atof(split_line[8].c_str());
          } 

          // GENERATOR_XQ
          if (nstr > 9) {
            data.xq = atof(split_line[9].c_str());
          } 

          // GENERATOR_XDP
          if (nstr > 10) {
            data.xdp = atof(split_line[10].c_str());
          } 

          // GENERATOR_XDPP
          if (nstr > 11) {
            data.xdpp = atof(split_line[11].c_str());
          } 

          // GENERATOR_XL
          if (nstr > 12) {
            data.xl = atof(split_line[12].c_str());
          } 

          // GENERATOR_S1
          if (nstr > 13) {
            data.s1 = atof(split_line[13].c_str());
          } 

          // GENERATOR_S12
          if (nstr > 14) {
            data.s12 = atof(split_line[14].c_str());
          } 
        } else if (sval == "GENROU") {
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
            data.xd = atof(split_line[9].c_str());
          } 

          // GENERATOR_XQ
          if (nstr > 10) {
            data.xq = atof(split_line[10].c_str());
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
            data.xl = atof(split_line[14].c_str());
          } 

          // GENERATOR_S1
          if (nstr > 15) {
            data.s1 = atof(split_line[15].c_str());
          } 

          // GENERATOR_S12
          if (nstr > 16) {
            data.s12 = atof(split_line[16].c_str());
          } 
        } else if (sval == "WSIEG1") {
          // GOVERNOR_JBUS
          if (nstr > 3) {
            data.jbus = atoi(split_line[3].c_str());
          }

          // GOVERNOR_M
          if (nstr > 4) {
            data.m = atoi(split_line[4].c_str());
          }

          // GOVERNOR_K
          if (nstr > 5) {
            data.m = atof(split_line[5].c_str());
          }

          // GOVERNOR_T1
          if (nstr > 6) {
            data.t1 = atof(split_line[6].c_str());
          }

          // GOVERNOR_T2
          if (nstr > 7) {
            data.t2 = atof(split_line[7].c_str());
          }

          // GOVERNOR_T3
          if (nstr > 8) {
            data.t3 = atof(split_line[8].c_str());
          }

          // GOVERNOR_UO
          if (nstr > 9) {
            data.uo = atof(split_line[9].c_str());
          }

          // GOVERNOR_UC
          if (nstr > 10) {
            data.uc = atof(split_line[10].c_str());
          }

          // GOVERNOR_PMAX
          if (nstr > 11) {
            data.pmax = atof(split_line[11].c_str());
          }

          // GOVERNOR_PMIN
          if (nstr > 12) {
            data.pmin = atof(split_line[12].c_str());
          }

          // GOVERNOR_T4
          if (nstr > 13) {
            data.t4 = atof(split_line[13].c_str());
          }

          // GOVERNOR_K1
          if (nstr > 14) {
            data.k1 = atof(split_line[14].c_str());
          }

          // GOVERNOR_K2
          if (nstr > 15) {
            data.k2 = atof(split_line[15].c_str());
          }

          // GOVERNOR_T5
          if (nstr > 16) {
            data.t5 = atof(split_line[16].c_str());
          }

          // GOVERNOR_K3
          if (nstr > 17) {
            data.k3 = atof(split_line[17].c_str());
          }

          // GOVERNOR_K4
          if (nstr > 18) {
            data.k4 = atof(split_line[18].c_str());
          }

          // GOVERNOR_T6
          if (nstr > 19) {
            data.t6 = atof(split_line[19].c_str());
          }

          // GOVERNOR_K5
          if (nstr > 20) {
            data.k5 = atof(split_line[20].c_str());
          }

          // GOVERNOR_K6
          if (nstr > 21) {
            data.k6 = atof(split_line[21].c_str());
          }

          // GOVERNOR_T7
          if (nstr > 22) {
            data.t7 = atof(split_line[22].c_str());
          }

          // GOVERNOR_K7
          if (nstr > 23) {
            data.k7 = atof(split_line[23].c_str());
          }

          // GOVERNOR_K8
          if (nstr > 24) {
            data.k8 = atof(split_line[24].c_str());
          }

          // GOVERNOR_DB1
          if (nstr > 25) {
            data.db1 = atof(split_line[25].c_str());
          }

          // GOVERNOR_ERR
          if (nstr > 26) {
            data.err = atof(split_line[26].c_str());
          }

          // GOVERNOR_DB2
          if (nstr > 27) {
            data.db2 = atof(split_line[27].c_str());
          }

          // GOVERNOR_GV1
          if (nstr > 28) {
            data.gv1 = atof(split_line[28].c_str());
          }

          // GOVERNOR_PGV1
          if (nstr > 29) {
            data.pgv1 = atof(split_line[29].c_str());
          }

          // GOVERNOR_GV2
          if (nstr > 30) {
            data.gv2 = atof(split_line[30].c_str());
          }

          // GOVERNOR_PGV2
          if (nstr > 31) {
            data.pgv2 = atof(split_line[31].c_str());
          }

          // GOVERNOR_GV3
          if (nstr > 32) {
            data.gv3 = atof(split_line[32].c_str());
          }

          // GOVERNOR_PGV3
          if (nstr > 33) {
            data.pgv3 = atof(split_line[33].c_str());
          }

          // GOVERNOR_GV4
          if (nstr > 34) {
            data.gv4 = atof(split_line[34].c_str());
          }

          // GOVERNOR_PGV4
          if (nstr > 35) {
            data.pgv4 = atof(split_line[35].c_str());
          }

          // GOVERNOR_GV5
          if (nstr > 36) {
            data.gv5 = atof(split_line[36].c_str());
          }

          // GOVERNOR_PGV5
          if (nstr > 37) {
            data.pgv5 = atof(split_line[37].c_str());
          }

          // GOVERNOR_IBLOCK
          if (nstr > 38) {
            data.iblock = atof(split_line[38].c_str());
          }
        } else if (sval == "EXDC1" || sval == "EXDC2") {
          // EXCITER_TR
          if (nstr > 3) {
            data.tr = atof(split_line[3].c_str());
          }

          // EXCITER_KA
          if (nstr > 4) {
            data.ka = atof(split_line[4].c_str());
          }

          // EXCITER_TA
          if (nstr > 5) {
            data.ta = atof(split_line[5].c_str());
          }

          // EXCITER_TB
          if (nstr > 6) {
            data.tb = atof(split_line[6].c_str());
          }

          // EXCITER_TC
          if (nstr > 7) {
            data.tc = atof(split_line[7].c_str());
          }

          // EXCITER_VRMAX
          if (nstr > 8) {
            data.vrmax = atof(split_line[8].c_str());
          }

          // EXCITER_VRMIN
          if (nstr > 9) {
            data.vrmin = atof(split_line[9].c_str());
          }

          // EXCITER_KE
          if (nstr > 10) {
            data.ke = atof(split_line[10].c_str());
          }

          // EXCITER_TE
          if (nstr > 11) {
            data.te = atof(split_line[11].c_str());
          }

          // EXCITER_KF
          if (nstr > 12) {
            data.kf = atof(split_line[12].c_str());
          }

          // EXCITER_TF1
          if (nstr > 13) {
            data.tf1 = atof(split_line[13].c_str());
          }

          // EXCITER_SWITCH
          if (nstr > 14) {
            data.rswitch = atof(split_line[14].c_str());
          }

          // EXCITER_E1
          if (nstr > 15) {
            data.e1 = atof(split_line[15].c_str());
          }

          // EXCITER_SE1
          if (nstr > 16) {
            data.se1 = atof(split_line[16].c_str());
          }

          // EXCITER_E2
          if (nstr > 17) {
            data.e2 = atof(split_line[17].c_str());
          }

          // EXCITER_SE2
          if (nstr > 18) {
            data.se2 = atof(split_line[18].c_str());
          }
        }
        ds_vector->push_back(data);
      }
    }

    void find_uc_vector(std::ifstream & input, std::vector<uc_params> *uc_vector)
    {
      std::string          line;
      uc_vector->clear();
      // Ignore first line containing header information
      std::getline(input,line);
      while(std::getline(input,line)) {
        std::vector<std::string>  split_line;
        boost::split(split_line, line, boost::algorithm::is_any_of(","),
            boost::token_compress_on);

        uc_params data;

        int nstr = split_line.size();
        if (nstr > 1) {
          data.type = atoi(split_line[1].c_str());
        }
        if (nstr > 2) {
          data.init_level = atof(split_line[2].c_str());
        }
        if (nstr > 3) {
          data.min_gen = atof(split_line[3].c_str());
        }
        if (nstr > 4) {
          data.max_gen = atof(split_line[4].c_str());
        }
        if (nstr > 5) {
          data.max_oper = atof(split_line[5].c_str());
        }
        if (nstr > 6) {
          data.min_up = atoi(split_line[6].c_str());
        }
        if (nstr > 7) {
          data.min_down = atoi(split_line[7].c_str());
        }
        if (nstr > 8) {
          data.ramp_up = atof(split_line[8].c_str());
        }
        if (nstr > 9) {
          data.ramp_down = atof(split_line[9].c_str());
        }
        if (nstr > 10) {
          data.start_up = atof(split_line[10].c_str());
        }
        if (nstr > 11) {
          data.const_cost = atof(split_line[11].c_str());
        }
        if (nstr > 12) {
          data.lin_cost = atof(split_line[12].c_str());
        }
        if (nstr > 13) {
          data.co_2_cost = atof(split_line[13].c_str());
        }
        if (nstr > 14) {
          data.init_prd = atof(split_line[14].c_str());
        }
        if (nstr > 15) {
          data.start_cap = atof(split_line[15].c_str());
        }
        if (nstr > 16) {
          data.shut_cap = atof(split_line[16].c_str());
        }
        if (nstr > 17) {
          data.bus_id = atoi(split_line[17].c_str());
        }
        if (nstr > 18) {
          // Clean up 2 character tag for generator ID
          gridpack::utility::StringUtils util;
          std::string tag = util.clean2Char(split_line[18]);
          strcpy(data.gen_id, tag.c_str());
        }
      }
    }

    /** Store pointer to bus and branch maps
     * @param busMap pointer to map between PTI and local indices
     * @param branchMap pointer to map between PTI bus pairs and local branch
     *        indices
     */
    void setMaps(std::map<int,int> *busMap,
                 std::map<std::pair<int, int>, int> *branchMap)
    {
      p_busMap = busMap;
      p_branchMap = branchMap;
    }

    /* ************************************************************************
     **************************************************************************
     ***** OBJECT DATA
     **************************************************************************
     *********************************************************************** */
    boost::shared_ptr<_network> p_network;

    gridpack::utility::CoarseTimer *p_timer;

    // Map of PTI indices to index in p_busData
    std::map<int,int> *p_busMap;
    // Map of PTI index pair to index in p_branchData
    std::map<std::pair<int, int>, int> *p_branchMap;

}; /* end of PTI base parser */

} /* namespace parser */
} /* namespace gridpack */

#endif /* BASEPTIPARSER_HPP_ */
