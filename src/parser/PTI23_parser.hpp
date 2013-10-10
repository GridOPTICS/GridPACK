/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * PTI23parser.hpp
 *
 *  Created on: May 23, 2013
 *      Author: Kevin Glass, Bruce Palmer
 */

#ifndef PTI23_PARSER_HPP_
#define PTI23_PARSER_HPP_

#include "mpi.h"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp> // needed of is_any_of()
#include <vector>
#include <map>
#include <cstdio>
#include <cstdlib>
#include "gridpack/utilities/exception.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/network/base_network.hpp"

#define TERM_CHAR '0'
// SOURCE: http://www.ee.washington.edu/research/pstca/formats/pti.txt

namespace gridpack {
namespace parser {


template <class _network>
  class PTI23_parser {
    public:

      /**
       * Constructor
       * @param network network object that will be filled with contents of network configuration file
       */
      PTI23_parser(boost::shared_ptr<_network> network)
      {
        p_network = network;
      }

      /**
       * Destructor
       */
      virtual ~PTI23_parser()
      {
        p_busData.clear();
        p_branchData.clear();
      }

      /**
       * Parse network configuration file and create network
       * @param fileName name of network file
       */
      void parse(const std::string &fileName)
      {
        getCase(fileName);
        createNetwork();
      }

      /*
       * A case is the collection of all data associated with a PTI23 file.
       * Each case is a a vector of data_set objects the contain all the data
       * associated with a partition of the PTI file. For example, the bus
       * data in the file constitutes a data_set. Each data_set is a vector of
       * gridpack::component::DataCollection objects. Each of these objects
       * contain a single instance of the data associated with a data_set. For
       * example, each line of the bus partition corresponds to a single
       * DataCollection object.
       */
      void getCase(const std::string & fileName)
      {

        p_busData.clear();
        p_branchData.clear();
        p_busMap.clear();

        int me;
        int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &me);
        if (me == 0) {
          std::ifstream            input;
          input.open(fileName.c_str());
          if (!input.is_open()) {
            throw gridpack::Exception("failed to open case data file");
          }

          find_case(input);
          find_buses(input);
//          find_loads(input);
          find_generators(input);
          find_branches(input);
          find_transformer(input);
          find_area(input);
          find_2term(input);
          find_shunt(input);
#if 0
          find_imped_corr(input);
          find_multi_term(input);
          find_multi_section(input);
          find_zone(input);
          find_interarea(input);
          find_owner(input);
          //find_line(input);
#endif
#if 0
          // debug
          int i;
          printf("BUS data size: %d\n",p_busData.size());
          for (i=0; i<p_busData.size(); i++) {
          printf("Dumping bus: %d\n",i);
            p_busData[i]->dump();
          }
          printf("BRANCH data size: %d\n",p_branchData.size());
          for (i=0; i<p_branchData.size(); i++) {
            printf("Dumping branch: %d\n",i);
            p_branchData[i]->dump();
          }
#endif
        }
      }

      void createNetwork(void)
      {
        int me;
        int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &me);
        if (me == 0) {
          int i;
          int numBus = p_busData.size();
          for (i=0; i<numBus; i++) {
            int idx;
            p_busData[i]->getValue(BUS_NUMBER,&idx);
            p_network->addBus(idx);
            p_network->setGlobalBusIndex(i,i);
            *(p_network->getBusData(i)) = *(p_busData[i]);
            p_network->getBusData(i)->addValue(CASE_ID,p_case_id);
            p_network->getBusData(i)->addValue(CASE_SBASE,p_case_sbase);
          }
          int numBranch = p_branchData.size();
          for (i=0; i<numBranch; i++) {
            int idx1, idx2;
            p_branchData[i]->getValue(BRANCH_FROMBUS,&idx1);
            p_branchData[i]->getValue(BRANCH_TOBUS,&idx2);
            p_network->addBranch(idx1, idx2);
            p_network->setGlobalBranchIndex(i,i);
            int g_idx1, g_idx2;
            std::map<int, int>::iterator it;
            it = p_busMap.find(idx1);
            g_idx1 = it->second;
            p_network->setGlobalBusIndex1(i,g_idx1);
            p_network->setLocalBusIndex1(i,g_idx1);
            it = p_busMap.find(idx2);
            g_idx2 = it->second;
            p_network->setGlobalBusIndex2(i,g_idx2);
            p_network->setLocalBusIndex2(i,g_idx2);
            *(p_network->getBranchData(i)) = *(p_branchData[i]);
            p_network->getBranchData(i)->addValue(CASE_ID,p_case_id);
            p_network->getBranchData(i)->addValue(CASE_SBASE,p_case_sbase);
          }
#if 0
          // debug
          printf("Number of buses: %d\n",numBus);
          for (i=0; i<numBus; i++) {
          printf("Dumping bus: %d\n",i);
            p_network->getBusData(i)->dump();
          }
          printf("Number of branches: %d\n",numBranch);
          for (i=0; i<numBranch; i++) {
          printf("Dumping branch: %d\n",i);
            p_network->getBranchData(i)->dump();
          }
#endif
        }
        p_busData.clear();
        p_branchData.clear();
      }
    protected:

      void find_case(std::ifstream & input)
      {
  //      data_set                                           case_set;
        std::string                                        line;
  //      std::vector<gridpack::component::DataCollection>   case_instance;

  //      gridpack::component::DataCollection                data;

        std::getline(input, line);
        std::vector<std::string>  split_line;

        boost::algorithm::split(split_line, line, boost::algorithm::is_any_of(" "), boost::token_compress_on);

        // CASE_ID             "IC"                   ranged integer
        p_case_id = atoi(split_line[0].c_str());

        // CASE_SBASE          "SBASE"                float
        p_case_sbase = atof(split_line[1].c_str());

        /*  These do not appear in the dictionary
        // CASE_RECORD2        "RECORD2"              string
        std::getline(input, line);
        data.addValue(CASE_RECORD2, line.c_str());
        case_instance.push_back(data);

        // CASE_RECORD3        "RECORD3"              string
        std::getline(input, line);
        data.addValue(CASE_RECORD3, line.c_str());
        case_instance.push_back(data);
         */
 //       case_set.push_back(case_instance);
 //       case_data->push_back(case_set);

      }

      void find_buses(std::ifstream & input)
      {
        std::string          line;
        int                  index = 0;
        int                  o_idx;
        std::getline(input, line);
        std::getline(input, line);
        std::getline(input, line);

        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          boost::shared_ptr<gridpack::component::DataCollection>
            data(new gridpack::component::DataCollection);

          // BUS_I               "I"                   integer
          o_idx = atoi(split_line[0].c_str());
          data->addValue(BUS_NUMBER, o_idx);
          p_busData.push_back(data);
          p_busMap.insert(std::pair<int,int>(o_idx,index));

          // BUS_NAME             "NAME"                 string
          data->addValue(BUS_NAME, (char*)split_line[9].c_str());

          // BUS_BASEKV           "BASKV"               float
          data->addValue(BUS_BASEKV, atof(split_line[10].c_str()));

          // BUS_TYPE               "IDE"                   integer
          data->addValue(BUS_TYPE, atoi(split_line[1].c_str()));

          // BUS_SHUNT_GL              "GL"                  float
          data->addValue(BUS_SHUNT_GL, atof(split_line[4].c_str()));

          // BUS_SHUNT_BL              "BL"                  float
          data->addValue(BUS_SHUNT_BL, atof(split_line[5].c_str()));

          // BUS_ZONE            "ZONE"                integer
          data->addValue(BUS_AREA, atoi(split_line[11].c_str()));

          // BUS_AREA            "IA"                integer
          data->addValue(BUS_ZONE, atoi(split_line[6].c_str()));

          // BUS_VOLTAGE_MAG              "VM"                  float
          data->addValue(BUS_VOLTAGE_MAG, atof(split_line[7].c_str()));

          // BUS_VOLTAGE_ANG              "VA"                  float
          data->addValue(BUS_VOLTAGE_ANG, atof(split_line[8].c_str()));

          // BUS_OWNER              "IA"                  integer
          data->addValue(BUS_OWNER, atoi(split_line[6].c_str()));

          // LOAD_PL                "PL"                  float
          data->addValue(LOAD_PL, atof(split_line[2].c_str()));

          // LOAD_QL                "QL"                  float
          data->addValue(LOAD_QL, atof(split_line[3].c_str()));

          index++;
          std::getline(input, line);
        }
      }

      void find_loads(std::ifstream & input)
      {
        std::string          line;
        std::getline(input, line); //this should be the first line of the block

        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);

          // LOAD_BUSNUMBER               "I"                   integer
          int l_idx, o_idx;
          o_idx = atoi(split_line[0].c_str());
          std::map<int, int>::iterator it;
          it = p_busMap.find(o_idx);
          if (it != p_busMap.end()) {
            l_idx = it->second;
          } else {
            std::getline(input, line);
            continue;
          }
          p_busData[l_idx]->addValue(LOAD_BUSNUMBER, atoi(split_line[0].c_str()));

          // LOAD_ID              "ID"                  integer
          p_busData[l_idx]->addValue(LOAD_ID, atoi(split_line[1].c_str()));

          // LOAD_STATUS              "ID"                  integer
          p_busData[l_idx]->addValue(LOAD_STATUS, atoi(split_line[1].c_str()));

          // LOAD_AREA            "ZONE"                integer
          p_busData[l_idx]->addValue(LOAD_AREA, atoi(split_line[11].c_str()));

          // LOAD_ZONE            "ZONE"                integer
          p_busData[l_idx]->addValue(LOAD_ZONE, atoi(split_line[11].c_str()));

          // LOAD_PL              "PG"                  float
          p_busData[l_idx]->addValue(LOAD_PL, atof(split_line[2].c_str()));

          // LOAD_QL              "QG"                  float
          p_busData[l_idx]->addValue(LOAD_QL, atof(split_line[3].c_str()));

          // LOAD_IP              "QT"                  float
          p_busData[l_idx]->addValue(LOAD_IP, atof(split_line[4].c_str()));

          // LOAD_IQ              "QB"                  float
          p_busData[l_idx]->addValue(LOAD_IQ, atof(split_line[5].c_str()));

          // LOAD_YP              "VS"                  float
          p_busData[l_idx]->addValue(LOAD_YP, atof(split_line[6].c_str()));

          // LOAD_YQ            "IREG"                integer
          p_busData[l_idx]->addValue(LOAD_YQ, atoi(split_line[7].c_str()));

          // LOAD_OWNER              "IA"                  integer
          p_busData[l_idx]->addValue(LOAD_OWNER, atoi(split_line[6].c_str()));

          std::getline(input, line);
        }
      }

      void find_generators(std::ifstream & input)
      {
        std::string          line;
        std::getline(input, line); //this should be the first line of the block
        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);

          // GENERATOR_BUSNUMBER               "I"                   integer
          int l_idx, o_idx;
          o_idx = atoi(split_line[0].c_str());
          std::map<int, int>::iterator it;
          it = p_busMap.find(o_idx);
          if (it != p_busMap.end()) {
            l_idx = it->second;
          } else {
            std::getline(input, line);
            continue;
          }

          // Find out how many generators are already on bus
          int ngen;
          if (!p_busData[l_idx]->getValue(GENERATOR_NUMBER, &ngen)) ngen = 0;


          p_busData[l_idx]->addValue(GENERATOR_BUSNUMBER, atoi(split_line[0].c_str()), ngen);

          // GENERATOR_ID              "ID"                  integer
          p_busData[l_idx]->addValue(GENERATOR_ID, atoi(split_line[1].c_str()), ngen);

          // GENERATOR_PG              "PG"                  float
          p_busData[l_idx]->addValue(GENERATOR_PG, atof(split_line[2].c_str()),
              ngen);

          // GENERATOR_QG              "QG"                  float
          p_busData[l_idx]->addValue(GENERATOR_QG, atof(split_line[3].c_str()),
              ngen);

          // GENERATOR_QMAX              "QT"                  float
          p_busData[l_idx]->addValue(GENERATOR_QMAX,
              atof(split_line[4].c_str()), ngen);

          // GENERATOR_QMIN              "QB"                  float
          p_busData[l_idx]->addValue(GENERATOR_QMIN,
              atof(split_line[5].c_str()), ngen);

          // GENERATOR_VS              "VS"                  float
          p_busData[l_idx]->addValue(GENERATOR_VS, atof(split_line[6].c_str()),
              ngen);

          // GENERATOR_IREG            "IREG"                integer
          p_busData[l_idx]->addValue(GENERATOR_IREG,
              atoi(split_line[7].c_str()), ngen);

          // GENERATOR_MBASE           "MBASE"               float
          p_busData[l_idx]->addValue(GENERATOR_MBASE,
              atof(split_line[8].c_str()), ngen);

          // GENERATOR_ZSORCE              "ZR"                  float
          p_busData[l_idx]->addValue(GENERATOR_ZSORCE,
              atof(split_line[9].c_str()), ngen);

          // GENERATOR_XTRAN              "ZX"                  float
          p_busData[l_idx]->addValue(GENERATOR_XTRAN,
              atof(split_line[10].c_str()), ngen);

          // GENERATOR_XT              "XT"                  float
          p_busData[l_idx]->addValue(GENERATOR_XT, atof(split_line[11].c_str()),
              ngen);

          // GENERATOR_RT              "RT"                  float
          p_busData[l_idx]->addValue(GENERATOR_RT, atof(split_line[12].c_str()),
              ngen);

          // GENERATOR_GTAP              "GTAP"                  float
          p_busData[l_idx]->addValue(GENERATOR_GTAP,
              atof(split_line[13].c_str()), ngen);

          // GENERATOR_STAT              "STAT"                  float
          p_busData[l_idx]->addValue(GENERATOR_STAT,
              atoi(split_line[14].c_str()), ngen);

          // GENERATOR_RMPCT           "RMPCT"               float
          p_busData[l_idx]->addValue(GENERATOR_RMPCT,
              atof(split_line[15].c_str()), ngen);

          // GENERATOR_PMAX              "PT"                  float
          p_busData[l_idx]->addValue(GENERATOR_PMAX,
              atof(split_line[16].c_str()), ngen);

          // GENERATOR_PMIN              "PB"                  float
          p_busData[l_idx]->addValue(GENERATOR_PMIN,
              atof(split_line[17].c_str()), ngen);

          // Pick up some non-standard values for Dynamic Simulation
          if (split_line.size() >= 22) {
            // GENERATOR_REACTANCE                             float
            p_busData[l_idx]->addValue(GENERATOR_REACTANCE,
                atof(split_line[18].c_str()), ngen);

            // GENERATOR_RESISTANCE                             float
            p_busData[l_idx]->addValue(GENERATOR_RESISTANCE,
                atof(split_line[19].c_str()), ngen);

            // GENERATOR_TRANSIENT_REACTANCE                             float
            p_busData[l_idx]->addValue(GENERATOR_TRANSIENT_REACTANCE,
                atof(split_line[20].c_str()), ngen);

            // GENERATOR_SUBTRANSIENT_REACTANCE                             float
            p_busData[l_idx]->addValue(GENERATOR_SUBTRANSIENT_REACTANCE,
                atof(split_line[21].c_str()), ngen);
          }

          // Increment number of generators in data object
          if (ngen == 0) {
            ngen = 1;
            p_busData[l_idx]->addValue(GENERATOR_NUMBER,ngen);
          } else {
            ngen++;
            p_busData[l_idx]->setValue(GENERATOR_NUMBER,ngen);
          }

          std::getline(input, line);
        }
      }

      void find_branches(std::ifstream & input)
      {
        std::string line;
        int  index   = 0;
        int  o_idx1, o_idx2;

        std::getline(input, line); //this should be the first line of the block

        while(test_end(line)) {
          std::pair<int, int> branch_pair;
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          boost::shared_ptr<gridpack::component::DataCollection>
            data(new gridpack::component::DataCollection);
          
          // BRANCH_INDEX                                   integer
          data->addValue(BRANCH_INDEX, index);
          p_branchData.push_back(data);

          o_idx1 = atoi(split_line[0].c_str());
          o_idx2 = atoi(split_line[1].c_str());

          // Switch order if one of the indices is negative
          //if (o_idx1<0 || o_idx2<0) {
          //  int t_idx = o_idx2;
          //  o_idx2 = o_idx1;
          //  o_idx1 = t_idx;
            if (o_idx1 < 0) o_idx1 = -o_idx1;
            if (o_idx2 < 0) o_idx2 = -o_idx2;
          //}

          // BRANCH_FROMBUS            "I"                   integer
          data->addValue(BRANCH_FROMBUS, o_idx1);
          // BRANCH_TOBUS            "J"                   integer
          data->addValue(BRANCH_TOBUS, o_idx2);

          // record the bus pairs that form the branch for subsequent searching
          branch_pair = std::pair<int, int>(o_idx1, o_idx2);
          p_branchMap.insert(std::pair<std::pair<int, int>, int >(branch_pair, index));

          // BRANCH_CKT          "CKT"                 character
          data->addValue(BRANCH_CKT, (char*)split_line[2].c_str());

          // BRANCH_R            "R"                   float
          data->addValue(BRANCH_R, atof(split_line[3].c_str()));

          // BRANCH_X            "X"                   float
          data->addValue(BRANCH_X, atof(split_line[4].c_str()));

          // BRANCH_B            "B"                   float
          data->addValue(BRANCH_B, atof(split_line[5].c_str()));

          // BRANCH_RATING_A        "RATEA"               float
          data->addValue(BRANCH_RATING_A, atof(split_line[6].c_str()));

          // BBRANCH_RATING_        "RATEB"               float
          data->addValue(BRANCH_RATING_B, atof(split_line[7].c_str()));

          // BRANCH_RATING_C        "RATEC"               float
          data->addValue(BRANCH_RATING_C, atof(split_line[8].c_str()));

          // BRANCH_TAP        "RATIO"               float
          data->addValue(BRANCH_TAP, atof(split_line[9].c_str()));

          // BRANCH_SHIFT        "SHIFT"               float
          data->addValue(BRANCH_SHIFT, atof(split_line[10].c_str()));

          // BRANCH_SHUNT_ADMTTNC_G1        "GI"               float
          data->addValue(BRANCH_SHUNT_ADMTTNC_G1, atof(split_line[11].c_str()));

          // BRANCH_SHUNT_ADMTTNC_B1        "BI"               float
          data->addValue(BRANCH_SHUNT_ADMTTNC_B1, atof(split_line[12].c_str()));

          // BRANCH_SHUNT_ADMTTNC_G2        "GJ"               float
          data->addValue(BRANCH_SHUNT_ADMTTNC_G2, atof(split_line[13].c_str()));

          // BRANCH_SHUNT_ADMTTNC_B2        "BJ"               float
          data->addValue(BRANCH_SHUNT_ADMTTNC_B2, atof(split_line[14].c_str()));

          // BRANCH_STATUS        "STATUS"               integer
//          data->addValue(BRANCH_STATUS, atoi(split_line[15].c_str()));

          // BRANCH_LENGTH           "GI"                  float
          //data->addValue(BRANCH_LENGTH, atof(split_line[11].c_str()));

          // BRANCH_OWNER           "ST"                  integer
          //data->addValue(BRANCH_OWNER, atoi(split_line[15].c_str()));

          ++index;
          std::getline(input, line);
        }
      }

      void find_transformer(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block

        std::pair<int, int>   branch_pair;
        // Find out how many branches already exist (note that this only works
        // when all branches are read in on head node
        int index = p_branchData.size();

        // get the branch that has the same to and from buses that the transformer hadto

        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);

          // KG: I'm assuming the BRANCH_FROMBUS is the bus index we need to match
          int fromBus = atoi(split_line[0].c_str());
          if (fromBus < 0) fromBus = -fromBus;

          // KG: I'm assuming the BRANCH_TOBUS is the bus index we need to match
          int toBus = atoi(split_line[1].c_str());
          if (toBus < 0) toBus = -toBus;

          // find branch corresponding to this
          int l_idx = 0;
          branch_pair = std::pair<int,int>(fromBus, toBus);
          std::map<std::pair<int, int>, int>::iterator it;
          it = p_branchMap.find(branch_pair);

          if (it != p_branchMap.end()) {
            l_idx = it->second;
          } else {
            std::getline(input, line);
            continue;
          }

          /*
           * type: integer
           * #define BRANCH_INDEX "BRANCH_INDEX"
           */
          p_branchData[l_idx]->addValue(BRANCH_INDEX,index);

          /*
           * type: integer
           * #define TRANSFORMER_BUS1 "TRANSFORMER_BUS1"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_BUS1, atoi(split_line[0].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_BUS2 "TRANSFORMER_BUS2"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_BUS2, atoi(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_BUS3 "TRANSFORMER_BUS3"
           */
//          data->addValue(TRANSFORMER_BUS3, atoi(split_line[1].c_str()));

          /*
           * type: string
           * #define TRANSFORMER_CKT "TRANSFORMER_CKT"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_CKT, split_line[2].c_str());

          /*
           * type: integer
           * #define TRANSFORMER_CW "TRANSFORMER_CW"
           X            */
          p_branchData[l_idx]->addValue(TRANSFORMER_CW, atoi(split_line[3].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_CZ "TRANSFORMER_CZ"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_CZ, atoi(split_line[5].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_CM "TRANSFORMER_CM"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_CM, atoi(split_line[5].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_MAG1 "TRANSFORMER_MAG1"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_MAG1, atof(split_line[5].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_MAG2 "TRANSFORMER_MAG2"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_MAG2, atof(split_line[5].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_NMETR "TRANSFORMER_NMETR"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_NMETR, atoi(split_line[1].c_str()));

          /*
           * type: string
           * #define TRANSFORMER_NAME "TRANSFORMER_NAME"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_NAME, split_line[2].c_str());

          /*
           * type: integer
           * #define TRANSFORMER_STATUS "TRANSFORMER_STATUS"
           *
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_STATUS, atoi(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_OWNER "TRANSFORMER_OWNER"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_OWNER, atoi(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_R1_2 "TRANSFORMER_R1_2"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_R1_2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_X1_2 "TRANSFORMER_X1_2"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_X1_2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_SBASE1_2 "TRANSFORMER_SBASE1_2"
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_SBASE1_2, atof(split_line[1].c_str()));

          std::getline(input, line);
        }
      }

      void find_area(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block

        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);

          // AREAINTG_ISW           "ISW"                  integer
          int l_idx, o_idx;
          o_idx = atoi(split_line[1].c_str());
          std::map<int, int>::iterator it;
          it = p_busMap.find(o_idx);
          if (it != p_busMap.end()) {
            l_idx = it->second;
          } else {
            std::getline(input, line);
            continue;
          }
          p_busData[l_idx]->addValue(AREAINTG_ISW, atoi(split_line[1].c_str()));

          // AREAINTG_NUMBER             "I"                    integer
          p_busData[l_idx]->addValue(AREAINTG_NUMBER, atoi(split_line[0].c_str()));

          // AREAINTG_PDES          "PDES"                 float
          p_busData[l_idx]->addValue(AREAINTG_PDES, atof(split_line[2].c_str()));

          // AREAINTG_PTOL          "PTOL"                 float
          p_busData[l_idx]->addValue(AREAINTG_PTOL, atof(split_line[3].c_str()));

          // AREAINTG_NAME         "ARNAM"                string
          p_busData[l_idx]->addValue(AREAINTG_NAME, split_line[4].c_str());

          std::getline(input, line);
        }
      }

      void find_2term(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block

        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::getline(input, line);
        }
      }

      void find_line(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block

        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::getline(input, line);
        }
      }


      /*

       */
      void find_shunt(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block
        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);

          /*
           * type: integer
           * #define SHUNT_BUSNUMBER "SHUNT_BUSNUMBER"
           */
          int l_idx, o_idx;
          l_idx = atoi(split_line[0].c_str());
          std::map<int, int>::iterator it;
          it = p_busMap.find(l_idx);
          if (it != p_busMap.end()) {
            o_idx = it->second;
          } else {
            std::getline(input, line);
            continue;
          }
          p_busData[o_idx]->addValue(SHUNT_BUSNUMBER, atoi(split_line[0].c_str()));

          /*
           * type: integer
           * #define SHUNT_MODSW "SHUNT_MODSW"
           */
          p_busData[o_idx]->addValue(SHUNT_MODSW, atoi(split_line[1].c_str()));

          /*
           * type: real float
           * #define SHUNT_VSWHI "SHUNT_VSWHI"
           */
          p_busData[o_idx]->addValue(SHUNT_VSWHI, atof(split_line[2].c_str()));

          /*
           * type: real float
           * #define SHUNT_VSWLO "SHUNT_VSWLO"
           */
          p_busData[o_idx]->addValue(SHUNT_VSWLO, atof(split_line[3].c_str()));

          /*
           * type: integer
           * #define SHUNT_SWREM "SHUNT_SWREM"
           */
          p_busData[o_idx]->addValue(SHUNT_SWREM, atoi(split_line[4].c_str()));

          /*
           * type: real float
           * #define SHUNT_RMPCT "SHUNT_RMPCT"
           */
//          p_busData[o_idx]->addValue(SHUNT_RMPCT, atof(split_line[4].c_str()));

          /*
           * type: string
           * #define SHUNT_RMIDNT "SHUNT_RMIDNT"
           */
//          p_busData[o_idx]->addValue(SHUNT_RMIDNT, (char*)split_line[5].c_str());

          /*
           * type: real float
           * #define SHUNT_BINIT "SHUNT_BINIT"
           */
          p_busData[o_idx]->addValue(SHUNT_BINIT, atof(split_line[5].c_str()));

          /*
           * type: integer
           * #define SHUNT_N1 "SHUNT_N1"
           */
          p_busData[o_idx]->addValue(SHUNT_N1, atoi(split_line[6].c_str()));

          /*
           * type: integer
           * #define SHUNT_N2 "SHUNT_N2"
           */
          p_busData[o_idx]->addValue(SHUNT_N2, atoi(split_line[8].c_str()));

          /*
           * type: integer
           * #define SHUNT_N3 "SHUNT_N3"
           */
          p_busData[o_idx]->addValue(SHUNT_N3, atoi(split_line[10].c_str()));

          /*
           * type: integer
           * #define SHUNT_N4 "SHUNT_N4"
           */
          p_busData[o_idx]->addValue(SHUNT_N4, atoi(split_line[12].c_str()));

          /*
           * type: integer
           * #define SHUNT_N5 "SHUNT_N5"
           */
          p_busData[o_idx]->addValue(SHUNT_N5, atoi(split_line[14].c_str()));

          /*
           * type: integer
           * #define SHUNT_N6 "SHUNT_N6"
           */
          p_busData[o_idx]->addValue(SHUNT_N6, atoi(split_line[16].c_str()));

          /*
           * type: integer
           * #define SHUNT_N7 "SHUNT_N7"
           */
          p_busData[o_idx]->addValue(SHUNT_N7, atoi(split_line[18].c_str()));

          /*
           * type: integer
           * #define SHUNT_N8 "SHUNT_N8"
           */
          p_busData[o_idx]->addValue(SHUNT_N8, atoi(split_line[20].c_str()));

          /*
           * type: real float
           * #define SHUNT_B1 "SHUNT_B1"
           */
          p_busData[o_idx]->addValue(SHUNT_B1, atof(split_line[7].c_str()));

          /*
           * type: real float
           * #define SHUNT_B2 "SHUNT_B2"
           */
          p_busData[o_idx]->addValue(SHUNT_B2, atof(split_line[9].c_str()));

          /*
           * type: real float
           * #define SHUNT_B3 "SHUNT_B3"
           */
          p_busData[o_idx]->addValue(SHUNT_B3, atof(split_line[11].c_str()));

          /*
           * type: real float
           * #define SHUNT_B4 "SHUNT_B4"
           */
          p_busData[o_idx]->addValue(SHUNT_B4, atof(split_line[13].c_str()));

          /*
           * type: real float
           * #define SHUNT_B5 "SHUNT_B5"
           */
          p_busData[o_idx]->addValue(SHUNT_B5, atof(split_line[15].c_str()));

          /*
           * type: real float
           * #define SHUNT_B6 "SHUNT_B6"
           */
          p_busData[o_idx]->addValue(SHUNT_B6, atof(split_line[17].c_str()));

          /*
           * type: real float
           * #define SHUNT_B7 "SHUNT_B7"
           */
          p_busData[o_idx]->addValue(SHUNT_B7, atof(split_line[19].c_str()));

          /*
           * type: real float
           * #define SHUNT_B8 "SHUNT_B8"
           */
          p_busData[o_idx]->addValue(SHUNT_B8, atof(split_line[21].c_str()));

          std::getline(input, line);
        }
      }

      void find_imped_corr(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block

        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
#if 0
          std::vector<gridpack::component::DataCollection>   imped_corr_instance;
          gridpack::component::DataCollection          data;

          /*
           * type: integer
           * #define XFMR_CORR_TABLE_NUMBER "XFMR_CORR_TABLE_NUMBER"
           */
          data.addValue(XFMR_CORR_TABLE_NUMBER, atoi(split_line[0].c_str()));
          imped_corr_instance.push_back(data);

          /*
           * type: real float
           * #define XFMR_CORR_TABLE_Ti "XFMR_CORR_TABLE_Ti"
           */
          data.addValue(XFMR_CORR_TABLE_Ti, atoi(split_line[0].c_str()));
          imped_corr_instance.push_back(data);

          /*
           * type: real float
           * #define XFMR_CORR_TABLE_Fi "XFMR_CORR_TABLE_Fi"
           */
          data.addValue(XFMR_CORR_TABLE_Fi, atoi(split_line[0].c_str()));
          imped_corr_instance.push_back(data);

          imped_corr_set.push_back(imped_corr_instance);
#endif
          std::getline(input, line);
        }
      }

      void find_multi_section(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block

        while(test_end(line)) {
#if 0
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   multi_section_instance;
          gridpack::component::DataCollection          data;

          /*
           * type: integer
           * #define MULTI_SEC_LINE_FROMBUS "MULTI_SEC_LINE_FROMBUS"

           */
          data.addValue(MULTI_SEC_LINE_FROMBUS, atoi(split_line[0].c_str()));
          multi_section_instance.push_back(data);

          /*
           * type: integer
           * #define MULTI_SEC_LINE_TOBUS "MULTI_SEC_LINE_TOBUS"

           */
          data.addValue(MULTI_SEC_LINE_TOBUS, atoi(split_line[0].c_str()));
          multi_section_instance.push_back(data);

          /*
           * type: string
           * #define MULTI_SEC_LINE_ID "MULTI_SEC_LINE_ID"

           */
          data.addValue(MULTI_SEC_LINE_ID, (char*)split_line[0].c_str());
          multi_section_instance.push_back(data);

          /*
           * type: integer
           * #define MULTI_SEC_LINE_DUMi "MULTI_SEC_LINE_DUMi"
           */
          data.addValue(MULTI_SEC_LINE_DUMi, atoi(split_line[0].c_str()));
          multi_section_instance.push_back(data);

          multi_section.push_back(multi_section_instance);
#endif
          std::getline(input, line);
        }
      }

      void find_multi_term(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block

        while(test_end(line)) {
          // TODO: parse something here
          std::getline(input, line);
        }
      }
      /*
       * ZONE_I          "I"                       integer
       * ZONE_NAME       "NAME"                    string
       */
      void find_zone(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block
        while(test_end(line)) {
          // TODO: parse something here
          std::getline(input, line);
        }
      }

      void find_interarea(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block

        while(test_end(line)) {
#if 0
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   inter_area_instance;
          gridpack::component::DataCollection          data;

          /*
           * type: integer
           * #define INTERAREA_TRANSFER_FROM "INTERAREA_TRANSFER_FROM"
           */
          data.addValue(INTERAREA_TRANSFER_FROM, atoi(split_line[0].c_str()));
          inter_area_instance.push_back(data);

          /*
           * type: integer
           * #define INTERAREA_TRANSFER_TO "INTERAREA_TRANSFER_TO"
           */
          data.addValue(INTERAREA_TRANSFER_TO, atoi(split_line[0].c_str()));
          inter_area_instance.push_back(data);

          /*
           * type: character
           * #define INTERAREA_TRANSFER_TRID "INTERAREA_TRANSFER_TRID"
           */
          data.addValue(INTERAREA_TRANSFER_TRID, split_line[0].c_str()[0]);
          inter_area_instance.push_back(data);

          /*
           * type: real float
           * #define INTERAREA_TRANSFER_PTRAN "INTERAREA_TRANSFER_PTRAN"
           */
          data.addValue(INTERAREA_TRANSFER_PTRAN, atof(split_line[0].c_str()));
          inter_area_instance.push_back(data);

          inter_area.push_back(inter_area_instance);
#endif
          std::getline(input, line);
        }
      }

      /*
       * type: integer
       * #define OWNER_NUMBER "OWNER_NUMBER"

       * type: integer
       * #define OWNER_NAME "OWNER_NAME"
       */
      void find_owner(std::ifstream & input)
      {
        std::string          line;
        std::getline(input, line); //this should be the first line of the block

        while(test_end(line)) {
#if 0
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   owner_instance;
          gridpack::component::DataCollection          data;

          data.addValue(OWNER_NUMBER, atoi(split_line[0].c_str()));
          owner_instance.push_back(data);

          data.addValue(OWNER_NAME, (char*)split_line[1].c_str());
          owner_instance.push_back(data);

          owner.push_back(owner_instance);
#endif
          std::getline(input, line);
        }
      }

    private:
      /*
       * Test to see if string terminates a section
       * @return: false if first non-blank character is TERM_CHAR
       */
      bool test_end(std::string &str) const
      {
#if 0
        int len = str.length();
        int i=0;
        while (i<len && str[i] == ' ') {
          i++;
        }
        if (i<len && str[i] != TERM_CHAR) {
          return true;
        }
        i++;
        while (i<len && str[i] == ' '){
          i++;
        }
        if (i<len && str[i] == '/') {
          return false;
        }
        return true;
#else
        if (str[0] == '0') {
          return false;
        } else {
          return true;
        }
#endif
      }
      /*
       * The case_data is the collection of all data points in the case file.
       * Each collection in the case data contains the data associated with a given
       * type. For example, the case is the collection of data describing the
       * current case and the bus data is the collection of data associated with
       * each bus. The type data may consist of zero or more instances of the
       * given type. For example, the bus data may contain several instances of
       * a bus. These type instances are composed of a set of key value pairs.
       * Each column as an associated key and each row is an instance of a given
       * type. When the parser is reading data for a type, the value found in each
       * column associated with the key for that column in a field_data structure.
       *
       * Within the PTI file there are the following group of data sets in order:
       *     case
       *     bus
       *     generator
       *     branch
       *     transformer
       *     dc_line
       *     shunt
       *     impedence corr
       *     multi-terminal
       *     multi-section
       *     zone
       *     inter-area
       *     owner
       *     device driver
       *
       * These data sets are stored in the case data as a collection of
       * data set and each data set is a
       */
      boost::shared_ptr<_network> p_network;

      // Vector of bus data objects
      std::vector<boost::shared_ptr<gridpack::component::DataCollection> > p_busData;
      // Vector of branch data objects
      std::vector<boost::shared_ptr<gridpack::component::DataCollection> > p_branchData;
      // Map of PTI indices to index in p_busData
      std::map<int,int> p_busMap;
      // Map of PTI index pair to index in p_branchData
      std::map<std::pair<int, int>, int> p_branchMap;

      // Global variables that apply to whole network
      int p_case_id;
      double p_case_sbase;
  };

} /* namespace parser */
} /* namespace gridpack */
#endif /* PTI23PARSER_HPP_ */
