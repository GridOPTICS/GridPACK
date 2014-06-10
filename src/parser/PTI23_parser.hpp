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

#define OLD_MAP

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp> // needed of is_any_of()
#ifndef OLD_MAP
#include <boost/unordered_map.hpp>
#endif
#include <vector>
#include <map>
#include <cstdio>
#include <cstdlib>

#include "gridpack/utilities/exception.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/network/base_network.hpp"

#define TERM_CHAR '0'
// SOURCE: http://www.ee.washington.edu/research/pstca/formats/pti.txt

namespace gridpack {
namespace parser {


template <class _network>
class PTI23_parser
{
    public:
  /// Constructor 
  /**
   * 
   * @param network network object that will be filled with contents
   * of network configuration file (must be child of network::BaseNetwork<>)
   */
  PTI23_parser(boost::shared_ptr<_network> network)
    : p_network(network)
  { }

      /**
       * Destructor
       */
      virtual ~PTI23_parser()
      {
        p_busData.clear();      // unnecessary
        p_branchData.clear();
      }

      /**
       * Parse network configuration file and create network
       * @param fileName name of network file
       */
      void parse(const std::string &fileName)
      {
        p_timer = gridpack::utility::CoarseTimer::instance();
        p_timer->configTimer(false);
        int t_total = p_timer->createCategory("Parser:Total Elapsed Time");
        p_timer->start(t_total);
        getCase(fileName);
        //brdcst_data();
        createNetwork();
        p_timer->stop(t_total);
        p_timer->configTimer(true);
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

        int t_case = p_timer->createCategory("Parser:getCase");
        p_timer->start(t_case);
        p_busData.clear();
        p_branchData.clear();
        p_busMap.clear();

        int me(p_network->communicator().rank());

        if (me == 0) {
          std::ifstream            input;
          input.open(fileName.c_str());
          if (!input.is_open()) {
            char buf[512];
            sprintf(buf,"Failed to open network configuration file: %s\n\n",
                fileName.c_str());
            throw gridpack::Exception(buf);
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
          input.close();
        }
        p_timer->stop(t_case);
      }

      void createNetwork(void)
      {
        int t_create = p_timer->createCategory("Parser:createNetwork");
        p_timer->start(t_create);
        int me(p_network->communicator().rank());
        int nprocs(p_network->communicator().size());
        int i;
        // Exchange information on number of buses and branches on each
        // processor
        int sbus[nprocs], sbranch[nprocs];
        int nbus[nprocs], nbranch[nprocs];
        for (i=0; i<nprocs; i++) {
          sbus[i] = 0;
          sbranch[i] = 0;
        }
        sbus[me] = p_busData.size();
        sbranch[me] = p_branchData.size();
        MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
        int ierr;
        ierr = MPI_Allreduce(sbus,nbus,nprocs,MPI_INT,MPI_SUM,comm);
        ierr = MPI_Allreduce(sbranch,nbranch,nprocs,MPI_INT,MPI_SUM,comm);
        // Transmit CASE_ID and CASE_SBASE to all processors
        int isval, irval;
        if (me == 0) {
          isval = p_case_id;
        } else {
          isval = 0;
        }
        ierr = MPI_Allreduce(&isval,&irval,1,MPI_INT,MPI_SUM,comm);
        p_case_id = irval;
        double sval, rval;
        if (me == 0) {
          sval = p_case_sbase;
        } else {
          sval = 0.0;
        }
        ierr = MPI_Allreduce(&sval,&rval,1,MPI_DOUBLE,MPI_SUM,comm);
        p_case_sbase = rval;
        // evaluate offsets for buses and branches
        int offset_bus[nprocs], offset_branch[nprocs];
        offset_bus[0] = 0;
        offset_branch[0] = 0;
        for (i=1; i<nprocs; i++) {
          offset_bus[i] = offset_bus[i-1]+nbus[i-1];
          offset_branch[i] = offset_branch[i-1]+nbranch[i-1];
        }

        int numBus = p_busData.size();
        for (i=0; i<numBus; i++) {
          int idx;
          p_busData[i]->getValue(BUS_NUMBER,&idx);
          p_network->addBus(idx);
          p_network->setGlobalBusIndex(i,i+offset_bus[me]);
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
          p_network->setGlobalBranchIndex(i,i+offset_branch[me]);
          int g_idx1, g_idx2;
#ifdef OLD_MAP
          std::map<int, int>::iterator it;
#else
          boost::unordered_map<int, int>::iterator it;
#endif
          it = p_busMap.find(idx1);
          g_idx1 = it->second;
          it = p_busMap.find(idx2);
          g_idx2 = it->second;
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
        p_busData.clear();
        p_branchData.clear();
        p_timer->stop(t_create);
      }
    protected:

      // Clean up 2 character tags so that single quotes are removed and single
      // character tags are right-justified
      std::string clean2Char(std::string string)
      {
        std::string tag = string;
        // Find and remove single quotes
        int ntok1 = tag.find_first_not_of('\'',0);
        int ntok2 = tag.find('\'',ntok1);
        if (ntok2 == std::string::npos) ntok2 = tag.length();
        std::string clean_tag = tag.substr(ntok1,ntok2-ntok1);
        //get rid of white space
        ntok1 = clean_tag.find_first_not_of(' ',0);
        ntok2 = clean_tag.find(' ',ntok1);
        if (ntok2 == std::string::npos) ntok2 = clean_tag.length();
        tag = clean_tag.substr(ntok1,ntok2-ntok1);
        if (tag.length() == 1) {
          clean_tag = " ";
          clean_tag.append(tag);
        } else {
          clean_tag = tag;
        }
        return clean_tag;
      }

      void find_case(std::ifstream & input)
      {
  //      data_set                                           case_set;
        std::string                                        line;
  //      std::vector<gridpack::component::DataCollection>   case_instance;

  //      gridpack::component::DataCollection                data;

        std::getline(input, line);
        while (check_comment(line)) {
          std::getline(input, line);
        }
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
          data->addValue(BUS_ZONE, atoi(split_line[11].c_str()));

          // BUS_AREA            "IA"                integer
          data->addValue(BUS_AREA, atoi(split_line[6].c_str()));

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
#ifdef OLD_MAP
          std::map<int, int>::iterator it;
#else
          boost::unordered_map<int, int>::iterator it;
#endif
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
#ifdef OLD_MAP
          std::map<int, int>::iterator it;
#else
          boost::unordered_map<int, int>::iterator it;
#endif
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

          // Clean up 2 character tag
          std::string tag = clean2Char(split_line[1]);
          // GENERATOR_ID              "ID"                  integer
          p_busData[l_idx]->addValue(GENERATOR_ID, (char*)tag.c_str(), ngen);

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

          // GENERATOR_ZSOURCE                                complex
          p_busData[l_idx]->addValue(GENERATOR_ZSOURCE,
              gridpack::ComplexType(atof(split_line[9].c_str()),
                atof(split_line[10].c_str())), ngen);

          // GENERATOR_XTRAN                              complex
          p_busData[l_idx]->addValue(GENERATOR_XTRAN,
              gridpack::ComplexType(atof(split_line[11].c_str()),
                atof(split_line[12].c_str())), ngen);

          // GENERATOR_RT              "RT"                  float
          p_busData[l_idx]->addValue(GENERATOR_RT, atof(split_line[11].c_str()),
              ngen);

          // GENERATOR_XT              "XT"                  float
          p_busData[l_idx]->addValue(GENERATOR_XT, atof(split_line[12].c_str()),
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

          // Pick up some more non-standard values for Dynamic Simulation
          if (split_line.size() >= 24) {
            // GENERATOR_INERTIA_CONSTANT_H                           float
            p_busData[l_idx]->addValue(GENERATOR_INERTIA_CONSTANT_H,
                atof(split_line[22].c_str()), ngen);

            // GENERATOR_DAMPING_COEFFICIENT_0                           float
            p_busData[l_idx]->addValue(GENERATOR_DAMPING_COEFFICIENT_0,
                atof(split_line[23].c_str()), ngen);
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
        int  o_idx1, o_idx2;
        int index = 0;

        std::getline(input, line); //this should be the first line of the block

        int nelems;
        while(test_end(line)) {
          std::pair<int, int> branch_pair;
          std::vector<std::string>  split_line;
          boost::split(split_line, line, boost::algorithm::is_any_of(","), boost::token_compress_on);
          
          o_idx1 = atoi(split_line[0].c_str());
          o_idx2 = atoi(split_line[1].c_str());

          // Switch sign if indices are negative
          if (o_idx1 < 0) o_idx1 = -o_idx1;
          if (o_idx2 < 0) o_idx2 = -o_idx2;

          // Check to see if pair has already been created
          int l_idx = 0;
          branch_pair = std::pair<int,int>(o_idx1, o_idx2);
#ifdef OLD_MAP
          std::map<std::pair<int, int>, int>::iterator it;
#else
          boost::unordered_map<std::pair<int, int>, int>::iterator it;
#endif
          it = p_branchMap.find(branch_pair);

          bool switched = false;
          if (it != p_branchMap.end()) {
            l_idx = it->second;
            p_branchData[l_idx]->getValue(BRANCH_NUM_ELEMENTS,&nelems);
          } else {
            // Check to see if from and to buses have been switched
            std::pair<int, int> new_branch_pair;
            new_branch_pair = std::pair<int,int>(o_idx2, o_idx1);
            it = p_branchMap.find(new_branch_pair);
            if (it != p_branchMap.end()) {
              printf("Found multiple lines with switched buses 1: %d 2: %d\n",
                  o_idx1,o_idx2);
              l_idx = it->second;
              p_branchData[l_idx]->getValue(BRANCH_NUM_ELEMENTS,&nelems);
              switched = true;
            } else {
              boost::shared_ptr<gridpack::component::DataCollection>
                data(new gridpack::component::DataCollection);
              l_idx = p_branchData.size();
              p_branchData.push_back(data);
              nelems = 0;
              p_branchData[l_idx]->addValue(BRANCH_NUM_ELEMENTS,nelems);
            }
          }

          if (nelems == 0) {
            // BRANCH_INDEX                                   integer
            p_branchData[l_idx]->addValue(BRANCH_INDEX, index);

            // BRANCH_FROMBUS            "I"                   integer
            p_branchData[l_idx]->addValue(BRANCH_FROMBUS, o_idx1);

            // BRANCH_TOBUS            "J"                   integer
            p_branchData[l_idx]->addValue(BRANCH_TOBUS, o_idx2);

            // add pair to branch map
            p_branchMap.insert(std::pair<std::pair<int, int>, int >(branch_pair,
                  index));
            index++;
          }
          
          // BRANCH_SWITCHED
          p_branchData[l_idx]->addValue(BRANCH_SWITCHED, switched, nelems);

          // Clean up 2 character tag
          std::string tag = clean2Char(split_line[2]);
          // BRANCH_CKT          "CKT"                 character
          p_branchData[l_idx]->addValue(BRANCH_CKT, (char*)tag.c_str(),
              nelems);

          // BRANCH_R            "R"                   float
          p_branchData[l_idx]->addValue(BRANCH_R, atof(split_line[3].c_str()),
              nelems);

          // BRANCH_X            "X"                   float
          p_branchData[l_idx]->addValue(BRANCH_X, atof(split_line[4].c_str()),
              nelems);

          // BRANCH_B            "B"                   float
          p_branchData[l_idx]->addValue(BRANCH_B, atof(split_line[5].c_str()),
              nelems);

          // BRANCH_RATING_A        "RATEA"               float
          p_branchData[l_idx]->addValue(BRANCH_RATING_A,
              atof(split_line[6].c_str()), nelems);

          // BBRANCH_RATING_        "RATEB"               float
          p_branchData[l_idx]->addValue(BRANCH_RATING_B,
              atof(split_line[7].c_str()), nelems);

          // BRANCH_RATING_C        "RATEC"               float
          p_branchData[l_idx]->addValue(BRANCH_RATING_C,
              atof(split_line[8].c_str()), nelems);

          // BRANCH_TAP        "RATIO"               float
          p_branchData[l_idx]->addValue(BRANCH_TAP, atof(split_line[9].c_str()), nelems);

          // BRANCH_SHIFT        "SHIFT"               float
          p_branchData[l_idx]->addValue(BRANCH_SHIFT,
              atof(split_line[10].c_str()), nelems);

          // BRANCH_SHUNT_ADMTTNC_G1        "GI"               float
          p_branchData[l_idx]->addValue(BRANCH_SHUNT_ADMTTNC_G1,
              atof(split_line[11].c_str()), nelems);

          // BRANCH_SHUNT_ADMTTNC_B1        "BI"               float
          p_branchData[l_idx]->addValue(BRANCH_SHUNT_ADMTTNC_B1,
              atof(split_line[12].c_str()), nelems);

          // BRANCH_SHUNT_ADMTTNC_G2        "GJ"               float
          p_branchData[l_idx]->addValue(BRANCH_SHUNT_ADMTTNC_G2,
              atof(split_line[13].c_str()), nelems);

          // BRANCH_SHUNT_ADMTTNC_B2        "BJ"               float
          p_branchData[l_idx]->addValue(BRANCH_SHUNT_ADMTTNC_B2,
              atof(split_line[14].c_str()), nelems);

          // BRANCH_STATUS        "STATUS"               integer
          p_branchData[l_idx]->addValue(BRANCH_STATUS,
              atoi(split_line[15].c_str()), nelems);

          nelems++;
          p_branchData[l_idx]->setValue(BRANCH_NUM_ELEMENTS,nelems);
          std::getline(input, line);
        }
      }

      // TODO: This code is NOT handling these elements correctly. Need to bring
      // it in line with find_branch routine and the definitions in the
      // ex_pti_file
      void find_transformer(std::ifstream & input)
      {
        std::string          line;

        std::getline(input, line); //this should be the first line of the block

        std::pair<int, int>   branch_pair;

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

          // find branch corresponding to this transformer line
          int l_idx = 0;
          branch_pair = std::pair<int,int>(fromBus, toBus);
#ifdef OLD_MAP
          std::map<std::pair<int, int>, int>::iterator it;
#else
          boost::unordered_map<std::pair<int, int>, int>::iterator it;
#endif
          it = p_branchMap.find(branch_pair);

          if (it != p_branchMap.end()) {
            l_idx = it->second;
          } else {
            std::getline(input, line);
            continue;
          }
          
          // Find number of transmission elements on this branch and then
          // compare the circuit number for this transformer adjustment with
          // BRANCH_CKT values
          int nelems = 0;
          p_branchData[l_idx]->getValue(BRANCH_NUM_ELEMENTS,&nelems);
          std::string b_ckt(clean2Char(split_line[2]));
          int i;
          int idx = -1;
          for (i=0; i<nelems; i++) {
            std::string t_ckt;
            p_branchData[l_idx]->getValue(BRANCH_CKT,&t_ckt,i);
            if (b_ckt == t_ckt) {
              idx = i;
              break;
            }
          }
          if (idx == -1) {
            printf("No match for transformer from %d to %d\n",
                fromBus,toBus);
            continue;
          }

          /*
           * type: integer
           * TRANSFORMER_CONTROL
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_CONTROL,
              atoi(split_line[3].c_str()),idx);

          /*
           * type: float
           * TRANSFORMER_RMA
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_RMA,
              atof(split_line[4].c_str()),idx);

          /*
           * type: float
           * TRANSFORMER_RMI
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_RMI,
              atof(split_line[5].c_str()),idx);

          /*
           * type: float
           * TRANSFORMER_VMA
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_VMA,
              atof(split_line[6].c_str()),idx);

          /*
           * type: float
           * TRANSFORMER_VMI
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_VMI,
              atof(split_line[7].c_str()),idx);

          /*
           * type: float
           * TRANSFORMER_STEP
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_STEP,
              atof(split_line[8].c_str()),idx);

          /*
           * type: float
           * TRANSFORMER_TABLE
           */
          p_branchData[l_idx]->addValue(TRANSFORMER_TABLE,
              atof(split_line[9].c_str()),idx);

          // This stuff is probably all wrong
#if 0

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
#endif

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
#ifdef OLD_MAP
          std::map<int, int>::iterator it;
#else
          boost::unordered_map<int, int>::iterator it;
#endif
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
#ifdef OLD_MAP
          std::map<int, int>::iterator it;
#else
          boost::unordered_map<int, int>::iterator it;
#endif
          it = p_busMap.find(l_idx);
          if (it != p_busMap.end()) {
            o_idx = it->second;
          } else {
            std::getline(input, line);
            continue;
          }
          int nval = split_line.size();

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
          if (8<nval) 
            p_busData[o_idx]->addValue(SHUNT_N2, atoi(split_line[8].c_str()));

          /*
           * type: integer
           * #define SHUNT_N3 "SHUNT_N3"
           */
          if (10<nval) 
            p_busData[o_idx]->addValue(SHUNT_N3, atoi(split_line[10].c_str()));

          /*
           * type: integer
           * #define SHUNT_N4 "SHUNT_N4"
           */
          if (12<nval) 
            p_busData[o_idx]->addValue(SHUNT_N4, atoi(split_line[12].c_str()));

          /*
           * type: integer
           * #define SHUNT_N5 "SHUNT_N5"
           */
          if (14<nval) 
            p_busData[o_idx]->addValue(SHUNT_N5, atoi(split_line[14].c_str()));

          /*
           * type: integer
           * #define SHUNT_N6 "SHUNT_N6"
           */
          if (16<nval) 
            p_busData[o_idx]->addValue(SHUNT_N6, atoi(split_line[16].c_str()));

          /*
           * type: integer
           * #define SHUNT_N7 "SHUNT_N7"
           */
          if (18<nval) 
            p_busData[o_idx]->addValue(SHUNT_N7, atoi(split_line[18].c_str()));

          /*
           * type: integer
           * #define SHUNT_N8 "SHUNT_N8"
           */
          if (20<nval) 
            p_busData[o_idx]->addValue(SHUNT_N8, atoi(split_line[20].c_str()));

          /*
           * type: real float
           * #define SHUNT_B1 "SHUNT_B1"
           */
          if (7<nval) 
            p_busData[o_idx]->addValue(SHUNT_B1, atof(split_line[7].c_str()));

          /*
           * type: real float
           * #define SHUNT_B2 "SHUNT_B2"
           */
          if (9<nval) 
            p_busData[o_idx]->addValue(SHUNT_B2, atof(split_line[9].c_str()));

          /*
           * type: real float
           * #define SHUNT_B3 "SHUNT_B3"
           */
          if (11<nval) 
            p_busData[o_idx]->addValue(SHUNT_B3, atof(split_line[11].c_str()));

          /*
           * type: real float
           * #define SHUNT_B4 "SHUNT_B4"
           */
          if (13<nval) 
            p_busData[o_idx]->addValue(SHUNT_B4, atof(split_line[13].c_str()));

          /*
           * type: real float
           * #define SHUNT_B5 "SHUNT_B5"
           */
          if (15<nval) 
            p_busData[o_idx]->addValue(SHUNT_B5, atof(split_line[15].c_str()));

          /*
           * type: real float
           * #define SHUNT_B6 "SHUNT_B6"
           */
          if (17<nval) 
            p_busData[o_idx]->addValue(SHUNT_B6, atof(split_line[17].c_str()));

          /*
           * type: real float
           * #define SHUNT_B7 "SHUNT_B7"
           */
          if (19<nval) 
            p_busData[o_idx]->addValue(SHUNT_B7, atof(split_line[19].c_str()));

          /*
           * type: real float
           * #define SHUNT_B8 "SHUNT_B8"
           */
          if (21<nval) 
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

      // Distribute data uniformly on processors
      void brdcst_data(void)
      {
        int t_brdcst = p_timer->createCategory("Parser:brdcst_data");
        int t_serial = p_timer->createCategory("Parser:data packing and unpacking");
        MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
        int me(p_network->communicator().rank());
        int nprocs(p_network->communicator().size());
        if (nprocs == 1) return;
        p_timer->start(t_brdcst);

        // find number of buses and branches and broadcast this information to
        // all processors
        int sbus, sbranch;
        if (me == 0) {
          sbus = p_busData.size();
          sbranch = p_branchData.size();
        } else {
          sbus = 0;
          sbranch = 0;
        }
        int ierr, nbus, nbranch;
        ierr = MPI_Allreduce(&sbus, &nbus, 1, MPI_INT, MPI_SUM, comm);
        ierr = MPI_Allreduce(&sbranch, &nbranch, 1, MPI_INT, MPI_SUM, comm);
        double rprocs = static_cast<double>(nprocs);
        double rme = static_cast<double>(me);
        int n, i;
        std::vector<gridpack::component::DataCollection>recvV;
        // distribute buses
        if (me == 0) {
          for (n=0; n<nprocs; n++) {
            double rn = static_cast<double>(n);
            int istart = static_cast<int>(static_cast<double>(nbus)*rn/rprocs);
            int iend = static_cast<int>(static_cast<double>(nbus)*(rn+1.0)/rprocs);
            if (n != 0) {
              p_timer->start(t_serial);
              std::vector<gridpack::component::DataCollection> sendV;
              for (i=istart; i<iend; i++) {
                sendV.push_back(*(p_busData[i]));
              }
              p_timer->stop(t_serial);
              static_cast<boost::mpi::communicator>(p_network->communicator()).send(n,n,sendV);
            } else {
              p_timer->start(t_serial);
              for (i=istart; i<iend; i++) {
                recvV.push_back(*(p_busData[i]));
              }
              p_timer->stop(t_serial);
            }
          }
        } else {
          int istart = static_cast<int>(static_cast<double>(nbus)*rme/rprocs);
          int iend = static_cast<int>(static_cast<double>(nbus)*(rme+1.0)/rprocs)-1;
          static_cast<boost::mpi::communicator>(p_network->communicator()).recv(0,me,recvV);
        }
        int nsize = recvV.size();
        p_busData.clear();
        p_timer->start(t_serial);
        for (i=0; i<nsize; i++) {
          boost::shared_ptr<gridpack::component::DataCollection> data(new
              gridpack::component::DataCollection);
          *data = recvV[i];
          p_busData.push_back(data);
        }
        p_timer->stop(t_serial);
        recvV.clear();
        // distribute branches
        if (me == 0) {
          for (n=0; n<nprocs; n++) {
            double rn = static_cast<double>(n);
            int istart = static_cast<int>(static_cast<double>(nbranch)*rn/rprocs);
            int iend = static_cast<int>(static_cast<double>(nbranch)*(rn+1.0)/rprocs);
            if (n != 0) {
              p_timer->start(t_serial);
              std::vector<gridpack::component::DataCollection> sendV;
              for (i=istart; i<iend; i++) {
                sendV.push_back(*(p_branchData[i]));
              }
              p_timer->stop(t_serial);
              static_cast<boost::mpi::communicator>(p_network->communicator()).send(n,n,sendV);
            } else {
              p_timer->start(t_serial);
              for (i=istart; i<iend; i++) {
                recvV.push_back(*(p_branchData[i]));
              }
              p_timer->stop(t_serial);
            }
          }
        } else {
          int istart = static_cast<int>(static_cast<double>(nbranch)*rme/rprocs);
          int iend = static_cast<int>(static_cast<double>(nbranch)*(rme+1.0)/rprocs)-1;
          static_cast<boost::mpi::communicator>(p_network->communicator()).recv(0,me,recvV);
        }
        nsize = recvV.size();
        p_branchData.clear();
        p_timer->start(t_serial);
        for (i=0; i<nsize; i++) {
          boost::shared_ptr<gridpack::component::DataCollection> data(new
              gridpack::component::DataCollection);
          *data = recvV[i];
          p_branchData.push_back(data);
        }
        p_timer->stop(t_serial);
#if 0
        // debug
        printf("p[%d] BUS data size: %d\n",me,p_busData.size());
        for (i=0; i<p_busData.size(); i++) {
          printf("p[%d] Dumping bus: %d\n",me,i);
          p_busData[i]->dump();
        }
        printf("p[%d] BRANCH data size: %d\n",me,p_branchData.size());
        for (i=0; i<p_branchData.size(); i++) {
          printf("p[%d] Dumping branch: %d\n",me,i);
          p_branchData[i]->dump();
        }
#endif
        p_timer->stop(t_brdcst);
      }

    private:
      /**
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

      /**
       * Test to see if string is a comment line. Check to see if first
       * non-blank characters are "//"
       */
      bool check_comment(std::string &str) const
      {
        int ntok = str.find_first_not_of(' ',0);
        if (ntok != std::string::npos && ntok+1 != std::string::npos &&
            str[ntok] == '/' && str[ntok+1] == '/') {
          return true;
        } else {
          return false;
        }
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
#ifdef OLD_MAP
      std::map<int,int> p_busMap;
#else
      boost::unordered_map<int, int> p_busMap;
#endif
      // Map of PTI index pair to index in p_branchData
#ifdef OLD_MAP
      std::map<std::pair<int, int>, int> p_branchMap;
#else
      boost::unordered_map<std::pair<int, int>, int> p_branchMap;
#endif

      // Global variables that apply to whole network
      int p_case_id;
      double p_case_sbase;
      gridpack::utility::CoarseTimer *p_timer;
  };

} /* namespace parser */
} /* namespace gridpack */
#endif /* PTI23PARSER_HPP_ */
