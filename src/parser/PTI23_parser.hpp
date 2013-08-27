/*
 * PTI23parser.hpp
 *
 *  Created on: May 23, 2013
 *      Author: kglass
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

      PTI23_parser(boost::shared_ptr<_network> network)
      {
        p_network = network;
      }
      virtual ~PTI23_parser()
      {
        p_busData.clear();
        p_branchData.clear();
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
#if 1
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
          }
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

        gridpack::component::DataCollection                data;

        std::getline(input, line);
        std::vector<std::string>  split_line;

        boost::algorithm::split(split_line, line, boost::algorithm::is_any_of(" "), boost::token_compress_on);

        // CASE_IC             "IC"                   ranged integer
  //      data.addValue(CASE_ID, atoi(split_line[0].c_str()));
  //      case_instance.push_back(data);

        // CASE_SBASE          "SBASE"                float
  //      data.addValue(CASE_SBASE, atof(split_line[1].c_str()));
  //      case_instance.push_back(data);

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
 //       data_set                        bus_set;
        std::string          line;
        int                  index = 0;
        int                  o_idx;
        std::getline(input, line);
        std::getline(input, line);
        std::getline(input, line);
        printf("(find_buses) Got to 1\n");

#if 0
        while(line[0] != TERM_CHAR) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   bus_instance;
          gridpack::component::DataCollection     data;

          // BUS_I               "I"                   integer
          o_idx = atoi(split_line[0].c_str());
          data.addValue(BUS_NUMBER, o_idx);
          bus_instance.push_back(data);

          // BUS_NAME             "IDE"                 ranged integer
          data.addValue(BUS_NAME, line.c_str());
          bus_instance.push_back(data);

          // BUS_BASEKV           "BASKV"               float
          data.addValue(BUS_BASEKV, atof(split_line[10].c_str()));
          bus_instance.push_back(data);

          // BUS_TYPE               "I"                   integer
          data.addValue(BUS_TYPE, atoi(split_line[0].c_str()));
          bus_instance.push_back(data);

          // BUS_SHUNT_GL              "GL"                  float
          data.addValue(BUS_SHUNT_GL, atof(split_line[4].c_str()));
          bus_instance.push_back(data);

          // BUS_SHUNT_BL              "BL"                  float
          data.addValue(BUS_SHUNT_BL, atof(split_line[5].c_str()));
          bus_instance.push_back(data);

          // BUS_AREA            "ZONE"                integer
          data.addValue(BUS_AREA, atoi(split_line[11].c_str()));
          bus_instance.push_back(data);

          // BUS_ZONE            "ZONE"                integer
          data.addValue(BUS_ZONE, atoi(split_line[11].c_str()));
          bus_instance.push_back(data);

          // BUS_VOLTAGE_MAG              "PL"                  float
          data.addValue(BUS_VOLTAGE_MAG, atof(split_line[2].c_str()));
          bus_instance.push_back(data);

          // BUS_VOLTAGE_ANG              "QL"                  float
          data.addValue(BUS_VOLTAGE_ANG, atof(split_line[3].c_str()));
          bus_instance.push_back(data);

          // BUS_OWNER              "IA"                  integer
          data.addValue(BUS_OWNER, atoi(split_line[6].c_str()));
          bus_instance.push_back(data);

          network->addBus(o_idx);
          boost::shared_ptr<gridpack::component::DataCollection> netData =
            network->getBusData(index);
          // This only works because process 0 is only process adding buses
          network->setGlobalBusIndex(index,index);

          *netData = data;
          ++index;

          bus_set.push_back(bus_instance);
          std::getline(input, line);
        }
        case_data->push_back(bus_set);
#else
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
          data->addValue(BUS_NAME, split_line[9].c_str());

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

#if 0
          network->addBus(o_idx);
          boost::shared_ptr<gridpack::component::DataCollection> netData =
            network->getBusData(index);
          // This only works because process 0 is only process adding buses
          network->setGlobalBusIndex(index,index);
#endif
          index++;
          std::getline(input, line);
        }
        printf("(find_buses) Got to 2 line: %s\n",line.c_str());
#endif
      }

      void find_loads(std::ifstream & input)
      {
//        data_set                        load_set;
        std::string          line;
        std::getline(input, line); //this should be the first line of the block
        printf("(loads) Got to 1 line: %s\n",line.c_str());
//        std::getline(input, line);
//        printf("(loads) Got to 2 line: %s\n",line.c_str());

#if 0
        while(line[0] != TERM_CHAR) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   load_instance;
          gridpack::component::DataCollection          data;

          // LOAD_BUSNUMBER               "I"                   integer
          data.addValue(LOAD_BUSNUMBER, atoi(split_line[0].c_str()));
          load_instance.push_back(data);

          // LOAD_ID              "ID"                  integer
          data.addValue(LOAD_ID, atoi(split_line[1].c_str()));
          load_instance.push_back(data);

          // LOAD_STATUS              "ID"                  integer
          data.addValue(LOAD_STATUS, atoi(split_line[1].c_str()));
          load_instance.push_back(data);

          // LOAD_AREA            "ZONE"                integer
          data.addValue(LOAD_AREA, atoi(split_line[11].c_str()));
          load_instance.push_back(data);

          // LOAD_ZONE            "ZONE"                integer
          data.addValue(LOAD_ZONE, atoi(split_line[11].c_str()));
          load_instance.push_back(data);

          // LOAD_PL              "PG"                  float
          data.addValue(LOAD_PL, atof(split_line[2].c_str()));
          load_instance.push_back(data);

          // LOAD_QL              "QG"                  float
          data.addValue(LOAD_QL, atof(split_line[3].c_str()));
          load_instance.push_back(data);

          // LOAD_IP              "QT"                  float
          data.addValue(LOAD_IP, atof(split_line[4].c_str()));
          load_instance.push_back(data);

          // LOAD_IQ              "QB"                  float
          data.addValue(LOAD_IQ, atof(split_line[5].c_str()));
          load_instance.push_back(data);

          // LOAD_YP              "VS"                  float
          data.addValue(LOAD_YP, atof(split_line[6].c_str()));
          load_instance.push_back(data);

          // LOAD_YQ            "IREG"                integer
          data.addValue(LOAD_YQ, atoi(split_line[7].c_str()));
          load_instance.push_back(data);

          // LOAD_OWNER              "IA"                  integer
          data.addValue(LOAD_OWNER, atoi(split_line[6].c_str()));
          load_instance.push_back(data);

          load_set.push_back(load_instance);
          std::getline(input, line);
        }
        case_data->push_back(load_set);
#else
        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);

          // LOAD_BUSNUMBER               "I"                   integer
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
          p_busData[o_idx]->addValue(LOAD_BUSNUMBER, atoi(split_line[0].c_str()));

          // LOAD_ID              "ID"                  integer
          p_busData[o_idx]->addValue(LOAD_ID, atoi(split_line[1].c_str()));

          // LOAD_STATUS              "ID"                  integer
          p_busData[o_idx]->addValue(LOAD_STATUS, atoi(split_line[1].c_str()));

          // LOAD_AREA            "ZONE"                integer
          p_busData[o_idx]->addValue(LOAD_AREA, atoi(split_line[11].c_str()));

          // LOAD_ZONE            "ZONE"                integer
          p_busData[o_idx]->addValue(LOAD_ZONE, atoi(split_line[11].c_str()));

          // LOAD_PL              "PG"                  float
          p_busData[o_idx]->addValue(LOAD_PL, atof(split_line[2].c_str()));

          // LOAD_QL              "QG"                  float
          p_busData[o_idx]->addValue(LOAD_QL, atof(split_line[3].c_str()));

          // LOAD_IP              "QT"                  float
          p_busData[o_idx]->addValue(LOAD_IP, atof(split_line[4].c_str()));

          // LOAD_IQ              "QB"                  float
          p_busData[o_idx]->addValue(LOAD_IQ, atof(split_line[5].c_str()));

          // LOAD_YP              "VS"                  float
          p_busData[o_idx]->addValue(LOAD_YP, atof(split_line[6].c_str()));

          // LOAD_YQ            "IREG"                integer
          p_busData[o_idx]->addValue(LOAD_YQ, atoi(split_line[7].c_str()));

          // LOAD_OWNER              "IA"                  integer
          p_busData[o_idx]->addValue(LOAD_OWNER, atoi(split_line[6].c_str()));

          std::getline(input, line);
        }
        printf("(loads) Got to 3 line: %s\n",line.c_str());
#endif
      }

      void find_generators(std::ifstream & input)
      {
//        data_set                        generator_set;
        std::string          line;
        std::getline(input, line); //this should be the first line of the block
        printf("(generators) Got to 1 line: %s\n",line.c_str());
//        std::getline(input, line);
//        printf("(generators) Got to 2 line: %s\n",line.c_str());
#if 0
        while(line[0] != TERM_CHAR) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   gen_instance;
          gridpack::component::DataCollection          data;

          // GENERATOR_BUSNUMBER               "I"                   integer
          data.addValue(GENERATOR_BUSNUMBER, atoi(split_line[0].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_ID              "ID"                  integer
          data.addValue(GENERATOR_ID, atoi(split_line[1].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_PG              "PG"                  float
          data.addValue(GENERATOR_PG, atof(split_line[2].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_QG              "QG"                  float
          data.addValue(GENERATOR_QG, atof(split_line[3].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_QMAX              "QT"                  float
          data.addValue(GENERATOR_QMAX, atof(split_line[4].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_QMIN              "QB"                  float
          data.addValue(GENERATOR_QMIN, atof(split_line[5].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_VS              "VS"                  float
          data.addValue(GENERATOR_VS, atof(split_line[6].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_IREG            "IREG"                integer
          data.addValue(GENERATOR_IREG, atoi(split_line[7].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_MBASE           "MBASE"               float
          data.addValue(GENERATOR_MBASE, atof(split_line[8].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_ZSORCE              "ZR"                  float
          data.addValue(GENERATOR_ZSORCE, atof(split_line[9].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_XTRAN              "ZX"                  float
          data.addValue(GENERATOR_XTRAN, atof(split_line[10].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_GTAP              "RT"                  float
          data.addValue(GENERATOR_GTAP, atof(split_line[11].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_XT              "XT"                  float
          data.addValue(GENERATOR_STAT, atof(split_line[12].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_RMPCT           "RMPCT"               float
          data.addValue(GENERATOR_RMPCT, atof(split_line[15].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_PMAX              "PT"                  float
          data.addValue(GENERATOR_PMAX, atof(split_line[16].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_PMIN              "PB"                  float
          data.addValue(GENERATOR_PMIN, atof(split_line[17].c_str()));
          gen_instance.push_back(data);

          // GENERATOR_OWNER              "IA"                  integer
          data.addValue(GENERATOR_OWNER, atoi(split_line[6].c_str()));
          gen_instance.push_back(data);

          generator_set.push_back(gen_instance);
          std::getline(input, line);
        }
        case_data->push_back(generator_set);
#else
        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);

          // GENERATOR_BUSNUMBER               "I"                   integer
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
          p_busData[o_idx]->addValue(GENERATOR_BUSNUMBER, atoi(split_line[0].c_str()));

          // GENERATOR_ID              "ID"                  integer
          p_busData[o_idx]->addValue(GENERATOR_ID, atoi(split_line[1].c_str()));

          // GENERATOR_PG              "PG"                  float
          p_busData[o_idx]->addValue(GENERATOR_PG, atof(split_line[2].c_str()));

          // GENERATOR_QG              "QG"                  float
          p_busData[o_idx]->addValue(GENERATOR_QG, atof(split_line[3].c_str()));

          // GENERATOR_QMAX              "QT"                  float
          p_busData[o_idx]->addValue(GENERATOR_QMAX, atof(split_line[4].c_str()));

          // GENERATOR_QMIN              "QB"                  float
          p_busData[o_idx]->addValue(GENERATOR_QMIN, atof(split_line[5].c_str()));

          // GENERATOR_VS              "VS"                  float
          p_busData[o_idx]->addValue(GENERATOR_VS, atof(split_line[6].c_str()));

          // GENERATOR_IREG            "IREG"                integer
          p_busData[o_idx]->addValue(GENERATOR_IREG, atoi(split_line[7].c_str()));

          // GENERATOR_MBASE           "MBASE"               float
          p_busData[o_idx]->addValue(GENERATOR_MBASE, atof(split_line[8].c_str()));

          // GENERATOR_ZSORCE              "ZR"                  float
          p_busData[o_idx]->addValue(GENERATOR_ZSORCE, atof(split_line[9].c_str()));

          // GENERATOR_XTRAN              "ZX"                  float
          p_busData[o_idx]->addValue(GENERATOR_XTRAN, atof(split_line[10].c_str()));

          // GENERATOR_XT              "XT"                  float
 //         p_busData[o_idx]->addValue(GENERATOR_XT, atof(split_line[11].c_str()));

          // GENERATOR_RT              "RT"                  float
 //         p_busData[o_idx]->addValue(GENERATOR_XT, atof(split_line[12].c_str()));

          // GENERATOR_GTAP              "GTAP"                  float
          p_busData[o_idx]->addValue(GENERATOR_GTAP, atof(split_line[13].c_str()));

          // GENERATOR_STAT              "STAT"                  float
          p_busData[o_idx]->addValue(GENERATOR_STAT, atof(split_line[14].c_str()));

          // GENERATOR_RMPCT           "RMPCT"               float
          p_busData[o_idx]->addValue(GENERATOR_RMPCT, atof(split_line[15].c_str()));

          // GENERATOR_PMAX              "PT"                  float
          p_busData[o_idx]->addValue(GENERATOR_PMAX, atof(split_line[16].c_str()));

          // GENERATOR_PMIN              "PB"                  float
          p_busData[o_idx]->addValue(GENERATOR_PMIN, atof(split_line[17].c_str()));

          std::getline(input, line);
        }
#endif
      printf("(generators) Got to 4 line: %s\n",line.c_str());
      }

      void find_branches(std::ifstream & input)
      {
//        data_set                        branch_set;
        std::string line;
        int  index   = 0;
        int  o_idx1, o_idx2;

        std::getline(input, line); //this should be the first line of the block
//        std::getline(input, line);

#if 0
        while(line[0] != TERM_CHAR) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   branch_instance;
          gridpack::component::DataCollection          data;

          // BRANCH_FROMBUS            "I"                   integer
          o_idx1 = atoi(split_line[0].c_str());
          data.addValue(BRANCH_FROMBUS, o_idx1);
          branch_instance.push_back(data);

          // BRANCH_TOBUS            "J"                   integer
          o_idx2 = atoi(split_line[1].c_str());
          data.addValue(BRANCH_TOBUS, o_idx2);
          branch_instance.push_back(data);

          // BRANCH_CKT          "CKT"                 character
          data.addValue(BRANCH_CKT, (split_line[2].c_str()));
          branch_instance.push_back(data);

          // BRANCH_R            "R"                   float
          data.addValue(BRANCH_R, atof(split_line[3].c_str()));
          branch_instance.push_back(data);

          // BRANCH_X            "X"                   float
          data.addValue(BRANCH_X, atof(split_line[4].c_str()));
          branch_instance.push_back(data);

          // BRANCH_B            "B"                   float
          data.addValue(BRANCH_B, atof(split_line[5].c_str()));
          branch_instance.push_back(data);

          // BRANCH_RATING_A        "RATEA"               float
          data.addValue(BRANCH_RATING_A, atof(split_line[6].c_str()));
          branch_instance.push_back(data);

          // BBRANCH_RATING_        "RATEB"               float
          data.addValue(BRANCH_RATING_B, atof(split_line[7].c_str()));
          branch_instance.push_back(data);

          // BRANCH_RATING_C        "RATEC"               float
          data.addValue(BRANCH_RATING_C, atof(split_line[8].c_str()));
          branch_instance.push_back(data);

          // BRANCH_SHUNT_ADMTTNC_G1        "RATIO"               float
          data.addValue(BRANCH_SHUNT_ADMTTNC_G1, atof(split_line[9].c_str()));
          branch_instance.push_back(data);

          // BRANCH_SHUNT_ADMTTNC_B1        "RATIO"               float
          data.addValue(BRANCH_SHUNT_ADMTTNC_B1, atof(split_line[9].c_str()));
          branch_instance.push_back(data);

          // BRANCH_SHUNT_ADMTTNC_G2        "RATIO"               float
          data.addValue(BRANCH_SHUNT_ADMTTNC_G2, atof(split_line[9].c_str()));
          branch_instance.push_back(data);

          // BRANCH_SHUNT_ADMTTNC_B2        "RATIO"               float
          data.addValue(BRANCH_SHUNT_ADMTTNC_B2, atof(split_line[9].c_str()));
          branch_instance.push_back(data);

          // BRANCH_STATUS        "ANGLE"               float
          data.addValue(BRANCH_STATUS, atoi(split_line[10].c_str()));
          branch_instance.push_back(data);

          // BRANCH_LENGTH           "GI"                  float
          data.addValue(BRANCH_LENGTH, atof(split_line[11].c_str()));
          branch_instance.push_back(data);

          // BRANCH_OWNER           "ST"                  integer
          data.addValue(BRANCH_OWNER, atoi(split_line[15].c_str()));
          branch_instance.push_back(data);

          network->addBranch(o_idx1, o_idx2);
          boost::shared_ptr<gridpack::component::DataCollection> netData =
            network->getBranchData(index);
          // This only works because process 0 is only process adding branches
          network->setGlobalBranchIndex(index,index);
          *netData = data;
          ++index;

          branch_set.push_back(branch_instance);
          std::getline(input, line);
        }
        case_data->push_back(branch_set);
#else
        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          boost::shared_ptr<gridpack::component::DataCollection>
            data(new gridpack::component::DataCollection);
          
          // BRANCH_INDEX                                   integer
          data->addValue(BRANCH_INDEX, index);
          p_branchData.push_back(data);

          // BRANCH_FROMBUS            "I"                   integer
          o_idx1 = atoi(split_line[0].c_str());
          data->addValue(BRANCH_FROMBUS, o_idx1);

          // BRANCH_TOBUS            "J"                   integer
          o_idx2 = atoi(split_line[1].c_str());
          data->addValue(BRANCH_TOBUS, o_idx2);

          // BRANCH_CKT          "CKT"                 character
          data->addValue(BRANCH_CKT, split_line[2].c_str());

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

          // BRANCH_SHUNT_ADMTTNC_G1        "RATIO"               float
//          data->addValue(BRANCH_SHUNT_ADMTTNC_G1, atof(split_line[9].c_str()));

          // BRANCH_SHUNT_ADMTTNC_B1        "ANGLE"               float
//          data->addValue(BRANCH_SHUNT_ADMTTNC_B1, atof(split_line[10].c_str()));

          // BRANCH_SHUNT_ADMTTNC_G1        "GI"               float
          data->addValue(BRANCH_SHUNT_ADMTTNC_G1, atof(split_line[11].c_str()));

          // BRANCH_SHUNT_ADMTTNC_B1        "BI"               float
          data->addValue(BRANCH_SHUNT_ADMTTNC_B1, atof(split_line[12].c_str()));

          // BRANCH_SHUNT_ADMTTNC_G2        "GJ"               float
          data->addValue(BRANCH_SHUNT_ADMTTNC_G1, atof(split_line[13].c_str()));

          // BRANCH_SHUNT_ADMTTNC_B2        "BJ"               float
          data->addValue(BRANCH_SHUNT_ADMTTNC_B1, atof(split_line[14].c_str()));

          // BRANCH_STATUS        "ANGLE"               integer
          data->addValue(BRANCH_STATUS, atoi(split_line[15].c_str()));

          // BRANCH_LENGTH           "GI"                  float
          data->addValue(BRANCH_LENGTH, atof(split_line[11].c_str()));

          // BRANCH_OWNER           "ST"                  integer
          data->addValue(BRANCH_OWNER, atoi(split_line[15].c_str()));

#if 0
          network->addBranch(o_idx1, o_idx2);
          boost::shared_ptr<gridpack::component::DataCollection> netData =
            network->getBranchData(index);
          // This only works because process 0 is only process adding branches
          network->setGlobalBranchIndex(index,index);
#endif
          ++index;
          std::getline(input, line);
        }
        printf("(branches) Got to 1 line: %s\n",line.c_str());
#endif
      }

      void find_transformer(std::ifstream & input)
      {
//        data_set                        transformer_set;
        std::string          line;

        std::getline(input, line); //this should be the first line of the block

        std::cout << "transformer block " << line << std::endl;
//        std::getline(input, line);
#if 0
        while(line[0] != TERM_CHAR) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   transformer_instance;
          gridpack::component::DataCollection          data;

          /*
           * type: integer
           * #define TRANSFORMER_BUS1 "TRANSFORMER_BUS1"
           */
          data.addValue(TRANSFORMER_BUS1, atoi(split_line[0].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_BUS2 "TRANSFORMER_BUS2"
           */
          data.addValue(TRANSFORMER_BUS2, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_BUS3 "TRANSFORMER_BUS3"
           */
          data.addValue(TRANSFORMER_BUS3, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: string
           * #define TRANSFORMER_CKT "TRANSFORMER_CKT"
           */
          data.addValue(TRANSFORMER_CKT, split_line[2].c_str());
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_CW "TRANSFORMER_CW"
           X            */
          data.addValue(TRANSFORMER_CW, atoi(split_line[3].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_CZ "TRANSFORMER_CZ"
           */
          data.addValue(TRANSFORMER_CZ, atoi(split_line[5].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_CM "TRANSFORMER_CM"
           */
          data.addValue(TRANSFORMER_CM, atoi(split_line[5].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_MAG1 "TRANSFORMER_MAG1"
           */
          data.addValue(TRANSFORMER_MAG1, atof(split_line[5].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_MAG2 "TRANSFORMER_MAG2"
           */
          data.addValue(TRANSFORMER_MAG2, atof(split_line[5].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_NMETR "TRANSFORMER_NMETR"
           */
          data.addValue(TRANSFORMER_NMETR, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: string
           * #define TRANSFORMER_NAME "TRANSFORMER_NAME"
           */
          data.addValue(TRANSFORMER_NAME, split_line[2].c_str());
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_STATUS "TRANSFORMER_STATUS"
           *
           */
          data.addValue(TRANSFORMER_STATUS, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_OWNER "TRANSFORMER_OWNER"
           */
          data.addValue(TRANSFORMER_OWNER, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_R1_2 "TRANSFORMER_R1_2"
           */
          data.addValue(TRANSFORMER_R1_2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_X1_2 "TRANSFORMER_X1_2"
           */
          data.addValue(TRANSFORMER_X1_2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_SBASE1_2 "TRANSFORMER_SBASE1_2"
           */
          data.addValue(TRANSFORMER_SBASE1_2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_R2_3 "TRANSFORMER_R2_3"
           */
          data.addValue(TRANSFORMER_R2_3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_X2_3 "TRANSFORMER_X2_3"
           */
          data.addValue(TRANSFORMER_X2_3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_SBASE2_3 "TRANSFORMER_SBASE2_3"
           */
          data.addValue(TRANSFORMER_SBASE2_3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_R3_1 "TRANSFORMER_R3_1"
           */
          data.addValue(TRANSFORMER_R3_1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_X3_1 "TRANSFORMER_X3_1"
           */
          data.addValue(TRANSFORMER_X3_1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_SBASE3_1 "TRANSFORMER_SBASE3_1"
           */
          data.addValue(TRANSFORMER_SBASE3_1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_VMSTAR "TRANSFORMER_VMSTAR"
           */
          data.addValue(TRANSFORMER_VMSTAR, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_ANSTAR "TRANSFORMER_ANSTAR"
           */
          data.addValue(TRANSFORMER_ANSTAR, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_WINDV1 "TRANSFORMER_WINDV1"
           */
          data.addValue(TRANSFORMER_WINDV1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_NOMV1 "TRANSFORMER_NOMV1"
           */
          data.addValue(TRANSFORMER_NOMV1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_ANG1 "TRANSFORMER_ANG1"
           */
          data.addValue(TRANSFORMER_ANG1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RATA1 "TRANSFORMER_RATA1"
           */
          data.addValue(TRANSFORMER_RATA1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RATB1 "TRANSFORMER_RATB1"
           */
          data.addValue(TRANSFORMER_RATB1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RATC1 "TRANSFORMER_RATC1"
           */
          data.addValue(TRANSFORMER_RATC1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_COD1 "TRANSFORMER_COD1"
           */
          data.addValue(TRANSFORMER_COD1, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_CONT1 "TRANSFORMER_CONT1"
           */
          data.addValue(TRANSFORMER_CONT1, atoi(split_line[0].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RMA1 "TRANSFORMER_RMA1"
           */
          data.addValue(TRANSFORMER_RMA1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RMI1 "TRANSFORMER_RMI1"
           */
          data.addValue(TRANSFORMER_RMI1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_VMA1 "TRANSFORMER_VMA1"
           */
          data.addValue(TRANSFORMER_VMA1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_VMI1 "TRANSFORMER_VMI1"
           */
          data.addValue(TRANSFORMER_VMI1, atof(split_line[3].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_NTP1 "TRANSFORMER_NTP1"
           */
          data.addValue(TRANSFORMER_NTP1, atoi(split_line[5].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_TAB1 "TRANSFORMER_TAB1"
           */
          data.addValue(TRANSFORMER_TAB1, atoi(split_line[5].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_CR1 "TRANSFORMER_CR1"
           */
          data.addValue(TRANSFORMER_CR1, atof(split_line[5].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_CX1 "TRANSFORMER_CX1"
           */
          data.addValue(TRANSFORMER_CX1, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_WINDV2 "TRANSFORMER_WINDV2"
           */
          data.addValue(TRANSFORMER_WINDV2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_NOMV2 "TRANSFORMER_NOMV2"
           */
          data.addValue(TRANSFORMER_NOMV2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_ANG2 "TRANSFORMER_ANG2"
           */
          data.addValue(TRANSFORMER_ANG2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RATA2 "TRANSFORMER_RATA2"
           */
          data.addValue(TRANSFORMER_RATA2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RATA2 "TRANSFORMER_RATB2"
           */
          data.addValue(TRANSFORMER_RATB2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RATC2 "TRANSFORMER_RATC2"
           */
          data.addValue(TRANSFORMER_RATC2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_COD2 "TRANSFORMER_COD2"
           */
          data.addValue(TRANSFORMER_COD2, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_CONT2 "TRANSFORMER_CONT2"
           */
          data.addValue(TRANSFORMER_CONT2, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RMA2 "TRANSFORMER_RMA2"
           */
          data.addValue(TRANSFORMER_RMA2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RMI2 "TRANSFORMER_RMI2"
           */
          data.addValue(TRANSFORMER_RMI2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_VMA2 "TRANSFORMER_VMA2"
           */
          data.addValue(TRANSFORMER_VMA2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_VMI2 "TRANSFORMER_VMI2"
           */
          data.addValue(TRANSFORMER_VMI2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_NTP2 "TRANSFORMER_NTP2"
           */
          data.addValue(TRANSFORMER_NTP2, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_TAB2 "TRANSFORMER_TAB2"
           */
          data.addValue(TRANSFORMER_TAB2, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_CR2 "TRANSFORMER_CR2"
           */
          data.addValue(TRANSFORMER_CR2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_CX2 "TRANSFORMER_CX2"
           */
          data.addValue(TRANSFORMER_CX2, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_WINDV3 "TRANSFORMER_WINDV3"
           */
          data.addValue(TRANSFORMER_WINDV3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_NOMV3 "TRANSFORMER_NOMV3"
           */
          data.addValue(TRANSFORMER_NOMV3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_ANG3 "TRANSFORMER_ANG3"
           */
          data.addValue(TRANSFORMER_ANG3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RATA3 "TRANSFORMER_RATA3"
           */
          data.addValue(TRANSFORMER_RATA3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RATB3 "TRANSFORMER_RATB3"
           */
          data.addValue(TRANSFORMER_RATB3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RATC3 "TRANSFORMER_RATC3"
           */
          data.addValue(TRANSFORMER_RATC3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_COD3 "TRANSFORMER_COD3"
           */
          data.addValue(TRANSFORMER_COD3, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_CONT3 "TRANSFORMER_CONT3"
           */
          data.addValue(TRANSFORMER_CONT3, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RMA3 "TRANSFORMER_RMA3"
           */
          data.addValue(TRANSFORMER_RMA3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_RMI3 "TRANSFORMER_RMI3"
           */
          data.addValue(TRANSFORMER_RMI3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_VMA3 "TRANSFORMER_VMA3"
           */
          data.addValue(TRANSFORMER_VMA3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_VMI3 "TRANSFORMER_VMI3"
           */
          data.addValue(TRANSFORMER_VMI3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_NTP3 "TRANSFORMER_NTP3"
           */
          data.addValue(TRANSFORMER_NTP3, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: integer
           * #define TRANSFORMER_TAB3 "TRANSFORMER_TAB3"
           */
          data.addValue(TRANSFORMER_TAB3, atoi(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_CR3 "TRANSFORMER_CR3"
           */
          data.addValue(TRANSFORMER_CR3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          /*
           * type: real float
           * #define TRANSFORMER_CX3 "TRANSFORMER_CX3"
           */
          data.addValue(TRANSFORMER_CX3, atof(split_line[1].c_str()));
          transformer_instance.push_back(data);

          transformer_set.push_back(transformer_instance);
          std::getline(input, line);
        }
        case_data->push_back(transformer_set);
#else
        // Find out how many branches already exist (note that this only works
        // when all branches are read in on head node
        int index = p_branchData.size();
        printf("(transformer) Got to 1 index: %d\n",index);
        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          boost::shared_ptr<gridpack::component::DataCollection>
            data(new gridpack::component::DataCollection);

          data->addValue(BRANCH_INDEX,index);
          p_branchData.push_back(data);
          /*
           * type: integer
           * #define TRANSFORMER_BUS1 "TRANSFORMER_BUS1"
           */
          data->addValue(TRANSFORMER_BUS1, atoi(split_line[0].c_str()));
          data->addValue(BRANCH_FROMBUS, atoi(split_line[0].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_BUS2 "TRANSFORMER_BUS2"
           */
          data->addValue(TRANSFORMER_BUS2, atoi(split_line[1].c_str()));
          data->addValue(BRANCH_TOBUS, atoi(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_BUS3 "TRANSFORMER_BUS3"
           */
//          data->addValue(TRANSFORMER_BUS3, atoi(split_line[1].c_str()));

          /*
           * type: string
           * #define TRANSFORMER_CKT "TRANSFORMER_CKT"
           */
          data->addValue(TRANSFORMER_CKT, split_line[2].c_str());

          /*
           * type: integer
           * #define TRANSFORMER_CW "TRANSFORMER_CW"
           X            */
          data->addValue(TRANSFORMER_CW, atoi(split_line[3].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_CZ "TRANSFORMER_CZ"
           */
          data->addValue(TRANSFORMER_CZ, atoi(split_line[5].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_CM "TRANSFORMER_CM"
           */
          data->addValue(TRANSFORMER_CM, atoi(split_line[5].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_MAG1 "TRANSFORMER_MAG1"
           */
          data->addValue(TRANSFORMER_MAG1, atof(split_line[5].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_MAG2 "TRANSFORMER_MAG2"
           */
          data->addValue(TRANSFORMER_MAG2, atof(split_line[5].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_NMETR "TRANSFORMER_NMETR"
           */
          data->addValue(TRANSFORMER_NMETR, atoi(split_line[1].c_str()));

          /*
           * type: string
           * #define TRANSFORMER_NAME "TRANSFORMER_NAME"
           */
          data->addValue(TRANSFORMER_NAME, split_line[2].c_str());

          /*
           * type: integer
           * #define TRANSFORMER_STATUS "TRANSFORMER_STATUS"
           *
           */
          data->addValue(TRANSFORMER_STATUS, atoi(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_OWNER "TRANSFORMER_OWNER"
           */
          data->addValue(TRANSFORMER_OWNER, atoi(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_R1_2 "TRANSFORMER_R1_2"
           */
          data->addValue(TRANSFORMER_R1_2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_X1_2 "TRANSFORMER_X1_2"
           */
          data->addValue(TRANSFORMER_X1_2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_SBASE1_2 "TRANSFORMER_SBASE1_2"
           */
          data->addValue(TRANSFORMER_SBASE1_2, atof(split_line[1].c_str()));

#if 0
          /*
           * type: real float
           * #define TRANSFORMER_R2_3 "TRANSFORMER_R2_3"
           */
          data->addValue(TRANSFORMER_R2_3, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_X2_3 "TRANSFORMER_X2_3"
           */
          data->addValue(TRANSFORMER_X2_3, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_SBASE2_3 "TRANSFORMER_SBASE2_3"
           */
          data->addValue(TRANSFORMER_SBASE2_3, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_R3_1 "TRANSFORMER_R3_1"
           */
          data->addValue(TRANSFORMER_R3_1, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_X3_1 "TRANSFORMER_X3_1"
           */
          data->addValue(TRANSFORMER_X3_1, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_SBASE3_1 "TRANSFORMER_SBASE3_1"
           */
          data->addValue(TRANSFORMER_SBASE3_1, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_VMSTAR "TRANSFORMER_VMSTAR"
           */
          data->addValue(TRANSFORMER_VMSTAR, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_ANSTAR "TRANSFORMER_ANSTAR"
           */
          data->addValue(TRANSFORMER_ANSTAR, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_WINDV1 "TRANSFORMER_WINDV1"
           */
          data->addValue(TRANSFORMER_WINDV1, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_NOMV1 "TRANSFORMER_NOMV1"
           */
          data->addValue(TRANSFORMER_NOMV1, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_ANG1 "TRANSFORMER_ANG1"
           */
          data->addValue(TRANSFORMER_ANG1, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RATA1 "TRANSFORMER_RATA1"
           */
          data->addValue(TRANSFORMER_RATA1, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RATB1 "TRANSFORMER_RATB1"
           */
          data->addValue(TRANSFORMER_RATB1, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RATC1 "TRANSFORMER_RATC1"
           */
          data->addValue(TRANSFORMER_RATC1, atof(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_COD1 "TRANSFORMER_COD1"
           */
          data->addValue(TRANSFORMER_COD1, atoi(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_CONT1 "TRANSFORMER_CONT1"
           */
          data->addValue(TRANSFORMER_CONT1, atoi(split_line[0].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RMA1 "TRANSFORMER_RMA1"
           */
          data->addValue(TRANSFORMER_RMA1, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RMI1 "TRANSFORMER_RMI1"
           */
          data->addValue(TRANSFORMER_RMI1, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_VMA1 "TRANSFORMER_VMA1"
           */
          data->addValue(TRANSFORMER_VMA1, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_VMI1 "TRANSFORMER_VMI1"
           */
          data->addValue(TRANSFORMER_VMI1, atof(split_line[3].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_NTP1 "TRANSFORMER_NTP1"
           */
          data->addValue(TRANSFORMER_NTP1, atoi(split_line[5].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_TAB1 "TRANSFORMER_TAB1"
           */
          data->addValue(TRANSFORMER_TAB1, atoi(split_line[5].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_CR1 "TRANSFORMER_CR1"
           */
          data->addValue(TRANSFORMER_CR1, atof(split_line[5].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_CX1 "TRANSFORMER_CX1"
           */
          data->addValue(TRANSFORMER_CX1, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_WINDV2 "TRANSFORMER_WINDV2"
           */
          data->addValue(TRANSFORMER_WINDV2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_NOMV2 "TRANSFORMER_NOMV2"
           */
          data->addValue(TRANSFORMER_NOMV2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_ANG2 "TRANSFORMER_ANG2"
           */
          data->addValue(TRANSFORMER_ANG2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RATA2 "TRANSFORMER_RATA2"
           */
          data->addValue(TRANSFORMER_RATA2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RATA2 "TRANSFORMER_RATB2"
           */
          data->addValue(TRANSFORMER_RATB2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RATC2 "TRANSFORMER_RATC2"
           */
          data->addValue(TRANSFORMER_RATC2, atof(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_COD2 "TRANSFORMER_COD2"
           */
          data->addValue(TRANSFORMER_COD2, atoi(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_CONT2 "TRANSFORMER_CONT2"
           */
          data->addValue(TRANSFORMER_CONT2, atoi(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RMA2 "TRANSFORMER_RMA2"
           */
          data->addValue(TRANSFORMER_RMA2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RMI2 "TRANSFORMER_RMI2"
           */
          data->addValue(TRANSFORMER_RMI2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_VMA2 "TRANSFORMER_VMA2"
           */
          data->addValue(TRANSFORMER_VMA2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_VMI2 "TRANSFORMER_VMI2"
           */
          data->addValue(TRANSFORMER_VMI2, atof(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_NTP2 "TRANSFORMER_NTP2"
           */
          data->addValue(TRANSFORMER_NTP2, atoi(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_TAB2 "TRANSFORMER_TAB2"
           */
          data->addValue(TRANSFORMER_TAB2, atoi(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_CR2 "TRANSFORMER_CR2"
           */
          data->addValue(TRANSFORMER_CR2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_CX2 "TRANSFORMER_CX2"
           */
          data->addValue(TRANSFORMER_CX2, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_WINDV3 "TRANSFORMER_WINDV3"
           */
          data->addValue(TRANSFORMER_WINDV3, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_NOMV3 "TRANSFORMER_NOMV3"
           */
          data->addValue(TRANSFORMER_NOMV3, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_ANG3 "TRANSFORMER_ANG3"
           */
          data->addValue(TRANSFORMER_ANG3, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RATA3 "TRANSFORMER_RATA3"
           */
          data->addValue(TRANSFORMER_RATA3, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RATB3 "TRANSFORMER_RATB3"
           */
          data->addValue(TRANSFORMER_RATB3, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RATC3 "TRANSFORMER_RATC3"
           */
          data->addValue(TRANSFORMER_RATC3, atof(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_COD3 "TRANSFORMER_COD3"
           */
          data->addValue(TRANSFORMER_COD3, atoi(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_CONT3 "TRANSFORMER_CONT3"
           */
          data->addValue(TRANSFORMER_CONT3, atoi(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RMA3 "TRANSFORMER_RMA3"
           */
          data->addValue(TRANSFORMER_RMA3, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_RMI3 "TRANSFORMER_RMI3"
           */
          data->addValue(TRANSFORMER_RMI3, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_VMA3 "TRANSFORMER_VMA3"
           */
          data->addValue(TRANSFORMER_VMA3, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_VMI3 "TRANSFORMER_VMI3"
           */
          data->addValue(TRANSFORMER_VMI3, atof(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_NTP3 "TRANSFORMER_NTP3"
           */
          data->addValue(TRANSFORMER_NTP3, atoi(split_line[1].c_str()));

          /*
           * type: integer
           * #define TRANSFORMER_TAB3 "TRANSFORMER_TAB3"
           */
          data->addValue(TRANSFORMER_TAB3, atoi(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_CR3 "TRANSFORMER_CR3"
           */
          data->addValue(TRANSFORMER_CR3, atof(split_line[1].c_str()));

          /*
           * type: real float
           * #define TRANSFORMER_CX3 "TRANSFORMER_CX3"
           */
          data->addValue(TRANSFORMER_CX3, atof(split_line[1].c_str()));
#endif

          std::getline(input, line);
        }
        printf("(transformer) Got to 3 line: %s\n",line.c_str());
#endif
      }

      void find_area(std::ifstream & input)
      {
//        data_set                        area_set;
        std::string          line;

        std::getline(input, line); //this should be the first line of the block
        printf("(find_area) Got to 1 line: %s\n",line.c_str());
//        std::getline(input, line);
//        printf("(find_area) Got to 2 line: %s\n",line.c_str());

#if 0
        while(line[0] != TERM_CHAR) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   area_instance;
          gridpack::component::DataCollection          data;

          // AREAINTG_NUMBER             "I"                    integer
          data.addValue(AREAINTG_NUMBER, atoi(split_line[0].c_str()));
          area_instance.push_back(data);

          // AREAINTG_ISW           "ISW"                  integer
          data.addValue(AREAINTG_ISW, atoi(split_line[1].c_str()));
          area_instance.push_back(data);

          // AREAINTG_PDES          "PDES"                 float
          data.addValue(AREAINTG_PDES, atof(split_line[2].c_str()));
          area_instance.push_back(data);

          // AREAINTG_PTOL          "PTOL"                 float
          data.addValue(AREAINTG_PTOL, atof(split_line[3].c_str()));
          area_instance.push_back(data);

          // AREAINTG_NAME         "ARNAM"                string
          data.addValue(AREAINTG_NAME, split_line[4].c_str());
          area_instance.push_back(data);

          area_set.push_back(area_instance);
          std::getline(input, line);
        }
        case_data->push_back(area_set);
#else
        while(test_end(line)) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);

          // AREAINTG_ISW           "ISW"                  integer
          int l_idx, o_idx;
          l_idx = atoi(split_line[1].c_str());
          std::map<int, int>::iterator it;
          it = p_busMap.find(l_idx);
          if (it != p_busMap.end()) {
            o_idx = it->second;
          } else {
            std::getline(input, line);
            continue;
          }
          p_busData[o_idx]->addValue(AREAINTG_ISW, atoi(split_line[1].c_str()));

          // AREAINTG_NUMBER             "I"                    integer
          p_busData[o_idx]->addValue(AREAINTG_NUMBER, atoi(split_line[0].c_str()));

          // AREAINTG_PDES          "PDES"                 float
          p_busData[o_idx]->addValue(AREAINTG_PDES, atof(split_line[2].c_str()));

          // AREAINTG_PTOL          "PTOL"                 float
          p_busData[o_idx]->addValue(AREAINTG_PTOL, atof(split_line[3].c_str()));

          // AREAINTG_NAME         "ARNAM"                string
          p_busData[o_idx]->addValue(AREAINTG_NAME, split_line[4].c_str());

          std::getline(input, line);
        }
        printf("(find_area) Got to 3 line: %s\n",line.c_str());
#endif
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
        printf("(find_2term) Got to 3 line: %s\n",line.c_str());
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
        printf("(find_line) Got to 3 line: %s\n",line.c_str());
      }


      /*

       */
      void find_shunt(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block
//        std::getline(input, line);

#if 0
        while(line[0] != TERM_CHAR) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   shunt_instance;
          gridpack::component::DataCollection          data;

          /*
           * type: integer
           * #define SHUNT_BUSNUMBER "SHUNT_BUSNUMBER"
           */
          data.addValue(SHUNT_BUSNUMBER, atoi(split_line[0].c_str()));
          shunt_instance.push_back(data);

          /*

type: integer
#define SHUNT_MODSW "SHUNT_MODSW"
           */
          data.addValue(SHUNT_MODSW, atoi(split_line[0].c_str()));
          shunt_instance.push_back(data);

          /*

type: real float
#define SHUNT_VSWHI "SHUNT_VSWHI"
           *
           */
          data.addValue(SHUNT_VSWHI, atof(split_line[1].c_str()));
          shunt_instance.push_back(data);

          /*

type: real float
#define SHUNT_VSWLO "SHUNT_VSWLO"
           *
           */
          data.addValue(SHUNT_VSWLO, atof(split_line[2].c_str()));
          shunt_instance.push_back(data);

          /*

type: integer
#define SHUNT_SWREM "SHUNT_SWREM"
           */
          data.addValue(SHUNT_SWREM, atoi(split_line[3].c_str()));
          shunt_instance.push_back(data);

          /*

type: real float
#define SHUNT_RMPCT "SHUNT_RMPCT"
           *
           */
          data.addValue(SHUNT_RMPCT, atof(split_line[4].c_str()));
          shunt_instance.push_back(data);

          /*
type: string
#define SHUNT_RMIDNT "SHUNT_RMIDNT"
           *
           */
          data.addValue(SHUNT_RMIDNT, split_line[5].c_str());
          shunt_instance.push_back(data);

          /*
           * type: real float
           * #define SHUNT_BINIT "SHUNT_BINIT"
           */
          data.addValue(SHUNT_BINIT, atof(split_line[5].c_str()));
          shunt_instance.push_back(data);

          /*
type: integer
#define SHUNT_N1 "SHUNT_N1"
           *
           */
          data.addValue(SHUNT_N1, atoi(split_line[6].c_str()));
          shunt_instance.push_back(data);

          /*
type: integer
#define SHUNT_N2 "SHUNT_N2"
           */
          data.addValue(SHUNT_N2, atoi(split_line[8].c_str()));
          shunt_instance.push_back(data);

          /*
type: integer
#define SHUNT_N3 "SHUNT_N3"
           */
          data.addValue(SHUNT_N3, atoi(split_line[10].c_str()));
          shunt_instance.push_back(data);

          /*
type: integer
#define SHUNT_N4 "SHUNT_N4"
           */
          data.addValue(SHUNT_N4, atoi(split_line[12].c_str()));
          shunt_instance.push_back(data);

          /*
type: integer
#define SHUNT_N5 "SHUNT_N5"
           */
          data.addValue(SHUNT_N5, atoi(split_line[14].c_str()));
          shunt_instance.push_back(data);

          /*
type: integer
#define SHUNT_N6 "SHUNT_N6"
           */
          data.addValue(SHUNT_N6, atoi(split_line[16].c_str()));
          shunt_instance.push_back(data);

          /*
type: integer
#define SHUNT_N7 "SHUNT_N7"
           */
          data.addValue(SHUNT_N7, atoi(split_line[18].c_str()));
          shunt_instance.push_back(data);

          /*
type: integer
#define SHUNT_N8 "SHUNT_N8"

           */
          data.addValue(SHUNT_N8, atoi(split_line[20].c_str()));
          shunt_instance.push_back(data);

          /*
type: real float
#define SHUNT_B1 "SHUNT_B1"
           */
          data.addValue(SHUNT_B1, atof(split_line[7].c_str()));
          shunt_instance.push_back(data);

          /*
type: real float
#define SHUNT_B2 "SHUNT_B2"
           */
          data.addValue(SHUNT_B2, atof(split_line[9].c_str()));
          shunt_instance.push_back(data);

          /*
type: real float
#define SHUNT_B3 "SHUNT_B3"
           */
          data.addValue(SHUNT_B3, atof(split_line[11].c_str()));
          shunt_instance.push_back(data);

          /*
type: real float
#define SHUNT_B4 "SHUNT_B4"
           */
          data.addValue(SHUNT_B4, atof(split_line[13].c_str()));
          shunt_instance.push_back(data);

          /*
type: real float
#define SHUNT_B5 "SHUNT_B5"
           */
          data.addValue(SHUNT_B5, atof(split_line[15].c_str()));
          shunt_instance.push_back(data);

          /*
type: real float
#define SHUNT_B6 "SHUNT_B6"
           */
          data.addValue(SHUNT_B6, atof(split_line[17].c_str()));
          shunt_instance.push_back(data);

          /*
type: real float
#define SHUNT_B7 "SHUNT_B7"
           */
          data.addValue(SHUNT_B7, atof(split_line[19].c_str()));
          shunt_instance.push_back(data);

          /*
type: real float
#define SHUNT_B8 "SHUNT_B8"
           */
          data.addValue(SHUNT_B8, atof(split_line[21].c_str()));
          shunt_instance.push_back(data);

          shunt_set.push_back(shunt_instance);
          std::getline(input, line);
        }
        case_data->push_back(shunt_set);
#else
        printf("(shunt) Got to 1 line: %s\n",line.c_str());
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
//          p_busData[o_idx]->addValue(SHUNT_RMIDNT, split_line[5].c_str());

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
        printf("(shunt) Got to 3 line: %s\n",line.c_str());
        }
#endif
      }

      void find_imped_corr(std::ifstream & input)
      {
//        data_set                        imped_corr_set;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
//        std::getline(input, line);

#if 0
        while(line[0] != TERM_CHAR) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
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
          std::getline(input, line);
        }
        case_data->push_back(imped_corr_set);
#else
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
#endif
      }

      void find_multi_section(std::ifstream & input)
      {
//        data_set                        multi_section;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
        std::cout << "multi section block " << line << std::endl;
//        std::getline(input, line);

#if 0
        while(line[0] != TERM_CHAR) {
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
          data.addValue(MULTI_SEC_LINE_ID, split_line[0].c_str());
          multi_section_instance.push_back(data);

          /*
           * type: integer
           * #define MULTI_SEC_LINE_DUMi "MULTI_SEC_LINE_DUMi"
           */
          data.addValue(MULTI_SEC_LINE_DUMi, atoi(split_line[0].c_str()));
          multi_section_instance.push_back(data);

          multi_section.push_back(multi_section_instance);
          std::getline(input, line);
        }
        case_data->push_back(multi_section);
#else
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
          data.addValue(MULTI_SEC_LINE_ID, split_line[0].c_str());
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
#endif
      }

      void find_multi_term(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block
        std::cout << "multi section block " << line << std::endl;

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
//        data_set                        zone;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
//        std::getline(input, line);

#if 0
        while(line[0] != TERM_CHAR) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   zone_instance;
          gridpack::component::DataCollection          data;

          /*
           *
           */
          data.addValue(ZONE_NUMBER, atoi(split_line[0].c_str()));
          zone_instance.push_back(data);

          // ZONE_NAME       "NAME"                    string
          data.addValue(ZONE_NAME, split_line[1].c_str());
          zone_instance.push_back(data);

          zone.push_back(zone_instance);
          std::getline(input, line);
        }
        case_data->push_back(zone);
#else
        while(test_end(line)) {
#if 0
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   zone_instance;
          gridpack::component::DataCollection          data;

          /*
           *
           */
          data.addValue(ZONE_NUMBER, atoi(split_line[0].c_str()));

          // ZONE_NAME       "NAME"                    string
          data.addValue(ZONE_NAME, split_line[1].c_str());

#endif
          std::getline(input, line);
        }
#endif
      }

      void find_interarea(std::ifstream & input)
      {
//        data_set                        inter_area;

        std::string          line;

        std::getline(input, line); //this should be the first line of the block
 //       std::getline(input, line);

#if 0
        while(line[0] != TERM_CHAR) {
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
          std::getline(input, line);
        }
        case_data->push_back(inter_area);
#else
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
#endif
      }

      /*
       * type: integer
       * #define OWNER_NUMBER "OWNER_NUMBER"

       * type: integer
       * #define OWNER_NAME "OWNER_NAME"
       */
      void find_owner(std::ifstream & input)
      {
//        data_set                        owner;

        std::string          line;
        std::getline(input, line); //this should be the first line of the block

 //       std::getline(input, line);

#if 0
        while(line[0] != TERM_CHAR) {
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   owner_instance;
          gridpack::component::DataCollection          data;

          data.addValue(OWNER_NUMBER, atoi(split_line[0].c_str()));
          owner_instance.push_back(data);

          data.addValue(OWNER_NAME, split_line[1].c_str());
          owner_instance.push_back(data);

          owner.push_back(owner_instance);
          std::getline(input, line);
        }
        case_data->push_back(owner);
#else
        while(test_end(line)) {
#if 0
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          std::vector<gridpack::component::DataCollection>   owner_instance;
          gridpack::component::DataCollection          data;

          data.addValue(OWNER_NUMBER, atoi(split_line[0].c_str()));
          owner_instance.push_back(data);

          data.addValue(OWNER_NAME, split_line[1].c_str());
          owner_instance.push_back(data);

          owner.push_back(owner_instance);
#endif
          std::getline(input, line);
        }
#endif
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
      // Map of PTI indices to global bus indices
      std::map<int,int> p_busMap;
      // Map of PTI indices to global branch indices
      std::map<int,int> p_branchMap;
  };

} /* namespace parser */
} /* namespace gridpack */
#endif /* PTI23PARSER_HPP_ */
