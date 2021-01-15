/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/*
 * MAT_parser.hpp
 *
 *  Created on: May 20, 2020
 *      Author: Bruce Palmer
 */

#ifndef MAT_PARSER_HPP_
#define MAT_PARSER_HPP_


#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp> // needed of is_any_of()
#include <vector>
#include <map>
#include <cstdio>
#include <cstdlib>

#include "gridpack/utilities/exception.hpp"
#include "gridpack/timer/coarse_timer.hpp"
#include "gridpack/parser/dictionary.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/parser/base_parser.hpp"
#include "gridpack/parser/base_pti_parser.hpp"
#include "gridpack/parser/hash_distr.hpp"

#define TERM_CHAR '0'
// SOURCE: http://www.ee.washington.edu/research/pstca/formats/pti.txt

namespace gridpack {
namespace parser {


template <class _network>
class MAT_parser : public BaseParser<_network>
{
  public:
    /// Constructor 
    /**
     * 
     * @param network network object that will be filled with contents
     * of network configuration file (must be child of network::BaseNetwork<>)
     */
    MAT_parser(boost::shared_ptr<_network> network)
      : p_network(network)
    { 
      p_network_data = network->getNetworkData();
      this->setNetwork(network);
    }

    /**
     * Destructor
     */
    virtual ~MAT_parser()
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

      p_busData.clear();
      p_branchData.clear();
      p_busMap.clear();

      p_case_id = 0;
      p_base_case = 0.0;
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
        std::string line;
        std::string function, fbase, fbus, fgen, fbranch, farea, fcost;
        while (std::getline(input, line)) {
          if (!check_comment(line)) {
            // Find out what block is being parsed
            std::vector<std::string>  split_line = p_util.blankTokenizer(line);
            if (split_line.size() >= 4 && split_line[0] == "function") {
              function = split_line[1];
              fbase = function;
              fbase.append(".baseMVA");
              fbus = function;
              fbus.append(".bus");
              fgen = function;
              fgen.append(".gen");
              fbranch = function;
              fbranch.append(".branch");
              farea = function;
              farea.append(".areas");
              fcost = function;
              fcost.append(".gencost");
            } else if (split_line[0] == fbase) {
              printf("Parsing base case\n");
              parse_case(line);
            } else if (split_line[0] == fbus) {
              printf("Parsing buses\n");
              parse_buses(input);
            } else if (split_line[0] == fgen) {
              printf("Parsing generators\n");
              parse_generators(input);
            } else if (split_line[0] == fbranch) {
              printf("Parsing branches\n");
              parse_branches(input);
            } else if (split_line[0] == farea) {
              printf("Parsing areas\n");
              parse_areas(input);
            } else if (split_line[0] == fcost) {
              printf("Parsing costs\n");
              parse_gen_costs(input);
            }
          }
        }
      }
      MPI_Comm comm = static_cast<MPI_Comm>(p_network->communicator());
      double sval, rval;
      sval = p_base_case;
      int ierr = MPI_Allreduce(&sval,&rval,1,MPI_DOUBLE,MPI_SUM,comm);
      p_base_case = rval;
      int isval, irval;
      isval = p_case_id;
      ierr = MPI_Allreduce(&isval,&irval,1,MPI_INT,MPI_SUM,comm);
      p_case_id = irval;
      this->setCaseSBase(p_base_case);
      this->setCaseID(p_case_id);
      p_network->broadcastNetworkData(0);
      p_network_data = p_network->getNetworkData();
      this->createNetwork(p_busData,p_branchData);
      printf("p[%d] Total buses: %d Total branches: %d\n",p_network->communicator().rank(),
          p_network->numBuses(),p_network->numBranches());
      p_timer->stop(t_total);
      p_timer->configTimer(true);
    }

  protected:

    void parse_case(std::string &line)
    {
      removeSemiColon(line);
      std::vector<std::string>  split_line;
      boost::algorithm::split(split_line, line, boost::algorithm::is_any_of(" "), boost::token_compress_on);

      // Just default to 0
      // CASE_ID             "IC"                   ranged integer
      p_case_id = 0;

      // CASE_SBASE          "SBASE"                float
      p_base_case = atof(split_line[2].c_str());

      p_network_data->addValue(CASE_SBASE, p_base_case);
      p_network_data->addValue(CASE_ID, p_case_id);
    }

    void parse_buses(std::ifstream & input)
    {
      std::string         line;
      int                 index = 0;
      int                 o_idx;
      double pl,ql,bl,gl;

      std::getline(input, line);
      while(!test_end(line)) {
        removeSemiColon(line);
        std::vector<std::string>  split_line = p_util.blankTokenizer(line);
        boost::shared_ptr<gridpack::component::DataCollection>
          data(new gridpack::component::DataCollection);
        int nstr = split_line.size();

        // BUS_I               "I"                   integer
        o_idx = atoi(split_line[0].c_str());
        data->addValue(BUS_NUMBER, o_idx);
        p_busData.push_back(data);
        p_busMap.insert(std::pair<int,int>(o_idx,index));

        // Add case parameters to data set
        data->addValue(CASE_SBASE, p_base_case);
        data->addValue(CASE_ID, p_case_id);

        // Construct a name
        char name[32];
        sprintf(name,"BUS_%d",o_idx);
        // BUS_NAME             "NAME"                 string
        data->addValue(BUS_NAME, name);

        // BUS_TYPE               "IDE"                   integer
        if (nstr > 1) data->addValue(BUS_TYPE, atoi(split_line[1].c_str()));

        // BUS_BASEKV           "BASKV"               float
        if (nstr > 9) data->addValue(BUS_BASEKV, atof(split_line[9].c_str()));

        // BUS_SHUNT_GL              "GL"                  float
        gl = 0.0;
        if (nstr > 4) {
          gl = atof(split_line[4].c_str());
        }

        // BUS_SHUNT_BL              "BL"                  float
        bl = 0.0;
        if (nstr > 5) {
          bl = atof(split_line[5].c_str());
        }
        if (gl != 0.0 || bl != 0.0) {
          data->addValue(BUS_SHUNT_GL, atof(split_line[4].c_str()));
          data->addValue(BUS_SHUNT_GL, atof(split_line[4].c_str()),0);
          data->addValue(BUS_SHUNT_BL, atof(split_line[5].c_str()));
          data->addValue(BUS_SHUNT_BL, atof(split_line[5].c_str()),0);
          data->addValue(SHUNT_BUSNUMBER,o_idx);
          int ival = 1;
          data->addValue(SHUNT_NUMBER,ival);
          data->addValue(SHUNT_ID, "1 ", 0);
        }

        // BUS_ZONE            "ZONE"                integer
        if (nstr > 10) data->addValue(BUS_ZONE, atoi(split_line[10].c_str()));

        // BUS_AREA            "IA"                integer
        if (nstr > 6) data->addValue(BUS_AREA, atoi(split_line[6].c_str()));

        // BUS_VOLTAGE_MAG              "VM"                  float
        if (nstr > 7) data->addValue(BUS_VOLTAGE_MAG, atof(split_line[7].c_str()));

        // BUS_VOLTAGE_ANG              "VA"                  float
        if (nstr > 8) data->addValue(BUS_VOLTAGE_ANG, atof(split_line[8].c_str()));

        // LOAD_PL                "PL"                  float
        pl = 0.0;
        if (nstr > 2) {
          pl = atof(split_line[2].c_str());
        }

        // LOAD_QL                "QL"                  float
        ql = 0.0;
        if (nstr > 3) {
          ql = atof(split_line[3].c_str());
        }
        if (pl != 0.0 || ql != 0.0) {
          data->addValue(LOAD_PL, atof(split_line[2].c_str()));
          data->addValue(LOAD_PL, atof(split_line[2].c_str()),0);
          std::string tmp("1 ");
          data->addValue(LOAD_ID,tmp.c_str(),0);
          data->addValue(LOAD_QL, atof(split_line[3].c_str()));
          data->addValue(LOAD_QL, atof(split_line[3].c_str()),0);
          int ival = 1;
          data->addValue(LOAD_NUMBER,ival);
          data->addValue(LOAD_STATUS,ival,0);
          data->addValue(LOAD_BUSNUMBER,o_idx);
        }

        // BUS_VOLTAGE_MAX "VOLTAGE_MAX" float
        if (nstr > 11) data->addValue(BUS_VOLTAGE_MAX, atof(split_line[11].c_str()));

        // BUS_VOLTAGE_MIN "VOLTAGE_MIN" float
        if (nstr > 12) data->addValue(BUS_VOLTAGE_MIN, atof(split_line[12].c_str()));

        index++;
        std::getline(input, line);
      }
    }

    void parse_generators(std::ifstream & input)
    {
      std::string          line;
      std::getline(input, line);
      std::multimap<int, int> genmap;
      std::multimap<int, int>::iterator gen_it;
      int index = 0;
      while(!test_end(line)) {
        this->removeSemiColon(line);
        std::vector<std::string>  split_line = p_util.blankTokenizer(line);

        // GENERATOR_BUSNUMBER               "I"                   integer
        int l_idx, o_idx;
        o_idx = atoi(split_line[0].c_str());
        std::map<int, int>::iterator it;
        int nstr = split_line.size();
        it = p_busMap.find(o_idx);
        if (it != p_busMap.end()) {
          l_idx = it->second;
        } else {
          std::getline(input, line);
          continue;
        }

        gen_it = genmap.find(o_idx);
        int ngen = 0;
        while (gen_it != genmap.end()) {
          ngen++;
          gen_it++;
        }
        genmap.insert(std::pair<int, int>(o_idx,ngen));

        p_busData[l_idx]->addValue(GENERATOR_BUSNUMBER, atoi(split_line[0].c_str()), ngen);

        p_busData[l_idx]->addValue(GENERATOR_GLOBAL_IDX,index,ngen);
        p_gen_bus.push_back(l_idx);
        p_gen_dev.push_back(ngen);

        // Generate 2 character tag
        std::string tag = get_char_id(ngen+1);
        // GENERATOR_ID              "ID"                  integer
        p_busData[l_idx]->addValue(GENERATOR_ID, tag.c_str(), ngen);

        // GENERATOR_PG              "PG"                  float
        if (nstr > 1) p_busData[l_idx]->addValue(GENERATOR_PG, atof(split_line[1].c_str()),
            ngen);

        // GENERATOR_QG              "QG"                  float
        if (nstr > 2) p_busData[l_idx]->addValue(GENERATOR_QG, atof(split_line[2].c_str()),
            ngen);

        // GENERATOR_QMAX              "QT"                  float
        if (nstr > 3) p_busData[l_idx]->addValue(GENERATOR_QMAX,
            atof(split_line[3].c_str()), ngen);

        // GENERATOR_QMIN              "QB"                  float
        if (nstr > 4) p_busData[l_idx]->addValue(GENERATOR_QMIN,
            atof(split_line[4].c_str()), ngen);

        // GENERATOR_VS              "VS"                  float
        if (nstr > 5) p_busData[l_idx]->addValue(GENERATOR_VS, atof(split_line[5].c_str()),
            ngen);

        // Just set this to 0?
        p_busData[l_idx]->addValue(GENERATOR_IREG, 0, ngen);

        // GENERATOR_MBASE           "MBASE"               float
        if (nstr > 6) p_busData[l_idx]->addValue(GENERATOR_MBASE,
            atof(split_line[6].c_str()), ngen);

        // GENERATOR_STAT              "STAT"                  float
        if (nstr > 7)  p_busData[l_idx]->addValue(GENERATOR_STAT,
            atoi(split_line[7].c_str()), ngen);

        // GENERATOR_PMAX              "PT"                  float
        if (nstr > 8) p_busData[l_idx]->addValue(GENERATOR_PMAX,
            atof(split_line[8].c_str()), ngen);

        // GENERATOR_PMIN              "PB"                  float
        if (nstr > 9) p_busData[l_idx]->addValue(GENERATOR_PMIN,
            atof(split_line[9].c_str()), ngen);

        // set generator number on bus
        if (ngen == 0) {
          p_busData[l_idx]->addValue(GENERATOR_NUMBER,ngen+1);
        } else {
          p_busData[l_idx]->setValue(GENERATOR_NUMBER,ngen+1);
        }

        index++;
        std::getline(input, line);
      }
      p_ngen = index;
    }

    void parse_branches(std::ifstream & input)
    {
      std::string          line;
      int  o_idx1, o_idx2;
      std::getline(input, line);
      std::multimap<int, int> genmap;
      std::multimap<int, int>::iterator gen_it;
      int nelems;
      int index = 0;
      while(!test_end(line)) {
        std::pair<int, int> branch_pair;
        this->removeSemiColon(line);
        std::vector<std::string>  split_line = p_util.blankTokenizer(line);

        o_idx1 = atoi(split_line[0].c_str());
        o_idx2 = atoi(split_line[1].c_str());

        // Check to see if pair has already been created
        int l_idx = 0;
        branch_pair = std::pair<int,int>(o_idx1, o_idx2);
        std::map<std::pair<int, int>, int>::iterator it;
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

        // Generate 2 character ID for line
        std::string tag = get_char_id(nelems+1);
        // BRANCH_CKT          "CKT"                 character
        p_branchData[l_idx]->addValue(BRANCH_CKT, tag.c_str(), nelems);

        // BRANCH_R            "R"                   float
        p_branchData[l_idx]->addValue(BRANCH_R, atof(split_line[2].c_str()),
            nelems);

        // BRANCH_X            "X"                   float
        p_branchData[l_idx]->addValue(BRANCH_X, atof(split_line[3].c_str()),
            nelems);

        // BRANCH_B            "B"                   float
        p_branchData[l_idx]->addValue(BRANCH_B, atof(split_line[4].c_str()),
            nelems);

        // BRANCH_RATING_A        "RATEA"               float
        p_branchData[l_idx]->addValue(BRANCH_RATING_A,
            atof(split_line[5].c_str()), nelems);

        // BBRANCH_RATING_        "RATEB"               float
        p_branchData[l_idx]->addValue(BRANCH_RATING_B,
            atof(split_line[6].c_str()), nelems);

        // BRANCH_RATING_C        "RATEC"               float
        p_branchData[l_idx]->addValue(BRANCH_RATING_C,
            atof(split_line[7].c_str()), nelems);

        // BRANCH_TAP        "RATIO"               float
        p_branchData[l_idx]->addValue(BRANCH_TAP, atof(split_line[8].c_str()), nelems);

        // BRANCH_SHIFT        "SHIFT"               float
        p_branchData[l_idx]->addValue(BRANCH_SHIFT,
            atof(split_line[9].c_str()), nelems);

        // BRANCH_STATUS        "STATUS"               integer
        p_branchData[l_idx]->addValue(BRANCH_STATUS,
            atoi(split_line[10].c_str()), nelems);

        p_branchData[l_idx]->addValue(BRANCH_SHUNT_ADMTTNC_G1, 0.0, nelems);
        p_branchData[l_idx]->addValue(BRANCH_SHUNT_ADMTTNC_B1, 0.0, nelems);
        p_branchData[l_idx]->addValue(BRANCH_SHUNT_ADMTTNC_G2, 0.0, nelems);
        p_branchData[l_idx]->addValue(BRANCH_SHUNT_ADMTTNC_B2, 0.0, nelems);

        nelems++;
        p_branchData[l_idx]->setValue(BRANCH_NUM_ELEMENTS,nelems);
        std::getline(input, line);
      }
    }


    void parse_areas(std::ifstream & input)
    {
      std::string          line;
      std::getline(input, line);
      int ncnt = 0;
      while(!test_end(line)) {
        this->removeSemiColon(line);
        std::vector<std::string>  split_line = p_util.blankTokenizer(line);

        // AREAINTG_ISW           "ISW"                  integer
        p_network_data->addValue(AREAINTG_ISW, atoi(split_line[1].c_str()),ncnt);

        // AREAINTG_NUMBER             "I"                    integer
        p_network_data->addValue(AREAINTG_NUMBER, atoi(split_line[0].c_str()),ncnt);

        ncnt++;
        std::getline(input, line);
      }
    }


    void parse_gen_costs(std::ifstream & input)
    {
      std::string          line;
      std::getline(input, line);
      int gcnt = 0;
      char buf[32];
      int i;
      while(!test_end(line)) {
        this->removeSemiColon(line);
        std::vector<std::string>  split_line = p_util.blankTokenizer(line);
        if (gcnt < p_ngen) {
          int ibus = p_gen_bus[gcnt];
          int igen = p_gen_dev[gcnt];
          int model = atoi(split_line[0].c_str()); 
          int npar = atoi(split_line[3].c_str());
          int tpar = split_line.size()-4;
          if (model == 1) {
            if (2*npar > tpar) continue;
            for (i=0; i<npar; i++) {
              sprintf(buf,"GENCOST_PARAM_P%d",i+1);
              p_busData[ibus]->addValue(buf,atof(split_line[2*i+4].c_str()),igen);
              sprintf(buf,"GENCOST_PARAM_F%d",i+1);
              p_busData[ibus]->addValue(buf,atof(split_line[2*i+1+4].c_str()),igen);
            }
          } else if (model == 2) {
            if (npar > tpar) continue;
            for (i=0; i<npar; i++) {
              sprintf(buf,"GENCOST_PARAM_C%d",npar-i);
              p_busData[ibus]->addValue(buf,atof(split_line[i+4].c_str()),igen);
            }
          } else {
            printf("Parsing unknown cost model: %d\n",model);
            continue;
          }
        } else {
          int ibus = p_gen_bus[gcnt-p_ngen];
          int igen = p_gen_dev[gcnt-p_ngen];
          int model = atoi(split_line[0].c_str()); 
          int npar = atoi(split_line[3].c_str());
          int tpar = split_line.size()-4;
          if (model == 1) {
            if (2*npar > tpar) continue;
            for (i=0; i<npar; i++) {
              sprintf(buf,"GENCOST_PARAM_R_P%d",i+1);
              p_busData[ibus]->addValue(buf,atof(split_line[2*i+4].c_str()),igen);
              sprintf(buf,"GENCOST_PARAM_R_F%d",i+1);
              p_busData[ibus]->addValue(buf,atof(split_line[2*i+1+4].c_str()),igen);
            }
          } else if (model == 2) {
            if (npar > tpar) continue;
            for (i=0; i<npar; i++) {
              sprintf(buf,"GENCOST_PARAM_R_C%d",npar-i);
              p_busData[ibus]->addValue(buf,atof(split_line[i+4].c_str()),igen);
            }
          } else {
            printf("Parsing unknown cost model: %d\n",model);
            continue;
          }
        }
        gcnt++;
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
      p_timer->stop(t_brdcst);
    }

  private:
    /**
     * Test to see if string terminates a section
     * @return: false if first non-blank character is TERM_CHAR
     */
    bool test_end(std::string &line) const
    {
      int ntok1 = line.find_first_of(']',0);
      if (ntok1 != std::string::npos) {
        int ntok2 = line.find_first_of(';',ntok1);
        if (ntok2 != std::string::npos) return true;
      }
      return false;
    }

    /**
     * Test to see if string is a comment line. Check to see if first
     * non-blank characters are "//"
     */
    bool check_comment(std::string &str) const
    {
      int ntok = str.find_first_not_of(' ',0);
      if (ntok != std::string::npos && str[ntok] == '%') {
        return true;
      } else if (ntok == std::string::npos) {
        return true;
      } else {
        return false;
      }
      return false;
    }

    /**
     * Remove semicolon at end of line
     * @param line text string for semicolon is being removed
     * @param returns true if semicolon found, false otherwise
     */
    bool removeSemiColon(std::string &line)
    {
      int ntok = line.find_first_of(';',0);
      if (ntok != std::string::npos) {
        line.erase(ntok);
        p_util.trim(line);
        return true;
      }
    }

    /**
     * Routine for generating character string identifiers for generators and
     * lines. These are not provided by the Mat Power files
     * @param idx index of object (generator or line)
     * @return two character string for object
     */
    std::string get_char_id(int idx)
    {
      std::string ret;
      if (idx < 10) {
        char c = '0'+idx;
        ret = c;
      } else if (idx < 36) {
        char c = 'A'+idx-10;
        ret = c;
      } else if (idx < 100) {
        char c[3];
        sprintf(c,"%d",idx);
        ret = c;
      }
      return p_util.clean2Char(ret);
    }

    boost::shared_ptr<_network> p_network;

    // Vector of bus data objects
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > p_busData;
    // Vector of branch data objects
    std::vector<boost::shared_ptr<gridpack::component::DataCollection> > p_branchData;
    // Map of PTI indices to index in p_busData
    std::map<int,int> p_busMap;
    // Map of PTI index pair to index in p_branchData
    std::map<std::pair<int, int>, int> p_branchMap;
    // Vectors that map back to generators
    std::vector<int> p_gen_bus;
    std::vector<int> p_gen_dev;
    int p_ngen;

    // Global variables that apply to whole network
    int p_case_id;
    double p_base_case;
    gridpack::utility::CoarseTimer *p_timer;

    gridpack::utility::StringUtils p_util;
    /**
     * Data collection object associated with network as a whole
     */
    boost::shared_ptr<gridpack::component::DataCollection> p_network_data;
};

} /* namespace parser */
} /* namespace gridpack */
#endif /* MAT_PARSER_HPP_ */
